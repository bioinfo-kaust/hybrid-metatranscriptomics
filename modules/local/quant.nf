/*
 * Stage 5 — quantification. Decoy-aware salmon index (the genome is the decoy), selective
 * alignment with GC and sequence bias correction.
 */

process SALMON_INDEX {
    label 'process_high'
    tag "$label"

    input:
    tuple val(label), path(fasta)
    path genome

    output:
    path 'salmon_index', emit: index
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: '-k 31'
    // no genome = no decoys: the index is the transcripts alone
    def decoy_arg = genome.name == 'NO_FILE_GENOME' ? '' : '-d decoys.txt'
    """
    grep '^>' ${genome} | sed 's/^>//; s/ .*//' > decoys.txt || true
    cat ${fasta} ${genome} > gentrome.fa
    salmon index -t gentrome.fa ${decoy_arg} -i salmon_index ${args} -p ${task.cpus}
    rm -f gentrome.fa
    echo "[${label}] salmon index over \$(grep -c '^>' ${fasta}) transcripts + \$(wc -l < decoys.txt) decoys"
    emit_versions.sh "${task.process}" salmon
    """

    stub:
    """
    mkdir -p salmon_index && touch salmon_index/info.json
    emit_versions.sh "${task.process}" salmon
    """
}

process SALMON_QUANT {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}"), emit: quant
    path "${meta.id}/aux_info/meta_info.json", emit: meta_info
    path "versions.yml", emit: versions, topic: versions

    script:
    // -l A auto-detects the library type, and the detection is CHECKED against params.strandedness
    // rather than printed: a per-sample flip (a swapped R1/R2, a mislabelled library) otherwise
    // passes silently and every fold change for that sample is wrong.
    def args = task.ext.args ?: '--gcBias --seqBias'
    def expect = params.strandedness == 'reverse'  ? 'ISR'
               : params.strandedness == 'forward'  ? 'ISF'
               : 'IU'
    """
    salmon quant -i ${index} -l A \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -p ${task.cpus} ${args} -o ${meta.id}

    python3 - <<'PY_CHECK'
import json, sys
d = json.load(open('${meta.id}/aux_info/meta_info.json'))
lib, mapped = d['library_types'][0], d['percent_mapped']
print('[${meta.id}] library %s, mapped %.2f%%' % (lib, mapped))
if '${params.strandedness_check}' == 'true' and lib != '${expect}':
    sys.exit("[ERROR] ${meta.id}: salmon detected library type %s but --strandedness "
             "'${params.strandedness}' implies ${expect}. Fix the samplesheet or set "
             "--strandedness_check false to override." % lib)
PY_CHECK
    emit_versions.sh "${task.process}" salmon python3
    """

    stub:
    """
    mkdir -p ${meta.id}/aux_info && touch ${meta.id}/quant.sf && echo "{}" > ${meta.id}/aux_info/meta_info.json
    emit_versions.sh "${task.process}" salmon python3
    """
}

process QUANT_SUMMARY {
    label 'process_single'

    input:
    path quants, stageAs: 'quants/*'

    output:
    path 'mapping_rates.tsv', emit: table
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    quant_summary.py --dir quants --out mapping_rates.tsv --min-mapping-rate ${params.min_mapping_rate}
    cat mapping_rates.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch mapping_rates.tsv
    emit_versions.sh "${task.process}" python3
    """
}

/*
 * Expression-filtered reference. The ladder's cross-source union carries every near-duplicate and
 * every unexpressed fragment into the index — which is what 31% BUSCO duplication looks like —
 * so the first quantification pass is used to drop what the experiment never saw, and the
 * survivors are re-indexed and re-quantified.
 */
process FILTER_BY_EXPRESSION {
    label 'process_medium'
    tag "$label"

    input:
    path quants, stageAs: 'quants/*'
    tuple val(label), path(fasta)
    tuple val(label2), path(tx2gene)
    tuple val(label3), path(annotation)

    output:
    tuple val("${label}_expressed"), path("${label}_expressed.fasta")           , emit: fasta
    tuple val("${label}_expressed"), path("tx2gene.${label}_expressed.tsv")     , emit: tx2gene
    tuple val("${label}_expressed"), path("annotation.${label}_expressed.tsv")  , emit: annotation
    path "${label}_expressed.filter.txt"                                   , emit: log
    path "${label}_expressed.de_eligible_tx.txt"                           , emit: expressed_tx
    path "versions.yml", emit: versions, topic: versions

    script:
    def ann = annotation.name != 'NO_FILE_ANN'
        ? "--annotation ${annotation} --out-annotation annotation.${label}_expressed.tsv" : ''
    // the rescue needs the annotation, so it silently degrades to expression-only without one
    def rescue = (params.expr_filter_rescue && annotation.name != 'NO_FILE_ANN')
        ? '--rescue-annotated' : ''
    """
    filter_by_expression.py \\
        --quant-dir quants --fasta ${fasta} --tx2gene ${tx2gene} ${ann} ${rescue} \\
        --min-tpm ${params.expr_filter_tpm} --min-samples ${params.expr_filter_samples} \\
        --out-fasta ${label}_expressed.fasta --out-tx2gene tx2gene.${label}_expressed.tsv \\
        --out-expressed-tx ${label}_expressed.de_eligible_tx.txt \\
        --log ${label}_expressed.filter.txt
    # a reference without annotation still needs the file to exist for the channel
    [ -s annotation.${label}_expressed.tsv ] || touch annotation.${label}_expressed.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch ${label}_expressed.fasta tx2gene.${label}_expressed.tsv annotation.${label}_expressed.tsv ${label}_expressed.filter.txt ${label}_expressed.de_eligible_tx.txt
    emit_versions.sh "${task.process}" python3
    """
}

/*
 * A third gene map, built from read sharing rather than sequence identity.
 *
 * String-stripped ids apply a per-source rule and reconcile nothing across sources; the locus
 * rebuild groups by genomic position and over-merges on a fragmented draft genome. Corset uses
 * salmon's fragment equivalence classes — which contigs the same reads are compatible with, and
 * how those counts move across samples — so contigs that are one locus to the data end up in one
 * cluster whatever assembler produced them. Its own gene-level counts are kept beside the map.
 */
process CORSET {
    label 'process_medium'
    label 'error_ignore'
    tag "$label"

    input:
    path eq_classes, stageAs: 'eq/*'
    tuple val(label), path(tx2gene)
    path samplesheet

    output:
    tuple val("${label}_corset"), path("tx2gene.${label}_corset.tsv"), emit: tx2gene
    path "corset_${label}"                                           , emit: results
    path "${label}_corset.summary.txt"                               , emit: summary
    path "versions.yml", emit: versions, topic: versions

    script:
    def bin = params.corset_bin ? "export PATH=${params.corset_bin}:\$PATH" : ''
    def args = task.ext.args ?: '-D 99999999'
    """
    ${bin}
    mkdir -p corset_${label}

    # one -g group and -n name per library, in the order the equivalence-class files are staged
    SRC=\$(ls -d eq/*/aux_info 2>/dev/null || true)
    [ -n "\$SRC" ] || { echo "[WARN] no salmon aux_info staged — is --dumpEq on?" \\
                        | tee ${label}_corset.summary.txt; \\
                        cut -f1,2 ${tx2gene} > tx2gene.${label}_corset.tsv; exit 0; }

    # salmon writes eq_classes.txt.gz. Corset opens the gzip without complaining and reads
    # "0 transcripts in 0 equivalence classes" from it, so decompress first rather than trust it.
    mkdir -p eqs
    for d in \$SRC; do
        smp=\$(basename \$(dirname \$d))
        mkdir -p eqs/\$smp
        if   [ -s "\$d/eq_classes.txt.gz" ]; then gunzip -c "\$d/eq_classes.txt.gz" > eqs/\$smp/eq_classes.txt
        elif [ -s "\$d/eq_classes.txt" ];    then cp       "\$d/eq_classes.txt"      eqs/\$smp/eq_classes.txt
        fi
    done
    FILES=\$(ls eqs/*/eq_classes.txt 2>/dev/null || true)
    [ -n "\$FILES" ] || { echo "[WARN] no equivalence classes found under aux_info" \\
                          | tee ${label}_corset.summary.txt; \\
                          cut -f1,2 ${tx2gene} > tx2gene.${label}_corset.tsv; exit 0; }

    # NB: not GROUPS. That is a bash special array; assigning to it returns non-zero and
    # `set -e` kills the task on that line, before corset is ever reached.
    SMP=\$(for f in \$FILES; do basename \$(dirname \$f); done | paste -sd, -)
    # Label-free by default: one group for every library. With -g <condition> corset's LRT splits
    # clusters on the very labels DESeq2 then tests; on a label-shuffle null that gave 89 false
    # calls per shuffle against 56 for the locus map (60 label-free), and its genes matched BUSCO /
    # Gacesa gene models far worse (pair recall 0.35 vs 0.74 label-free). scripts/50_genemap_eval.sh
    if [ "${params.corset_use_groups}" = "true" ]; then
        GRP=\$(for f in \$FILES; do
                   s=\$(basename \$(dirname \$f))
                   awk -F',' -v s="\$s" 'NR>1 && \$1==s {print \$2; exit}' ${samplesheet}
                done | paste -sd, -)
    else
        GRP=\$(for f in \$FILES; do echo all; done | paste -sd, -)
    fi
    echo "[corset] samples: \$SMP"
    echo "[corset] groups : \$GRP"

    corset -i salmon_eq_classes ${args} -f true -g "\$GRP" -n "\$SMP" \\
           -p corset_${label}/${label} \$FILES > corset.log 2>&1

    # corset writes <prefix>-clusters.txt as transcript<TAB>cluster — already a tx2gene map
    cp corset_${label}/${label}-clusters.txt tx2gene.${label}_corset.tsv


    {
      echo "Corset gene map — ${label}"
      echo "  transcripts clustered : \$(wc -l < tx2gene.${label}_corset.tsv)"
      echo "  clusters (genes)      : \$(cut -f2 tx2gene.${label}_corset.tsv | sort -u | wc -l)"
      echo "  string-stripped genes : \$(cut -f2 ${tx2gene} | sort -u | wc -l)"
    } | tee ${label}_corset.summary.txt
    emit_versions.sh "${task.process}" corset
    """

    stub:
    """
    mkdir -p corset_${label}
    touch tx2gene.${label}_corset.tsv ${label}_corset.summary.txt
    emit_versions.sh "${task.process}" corset
    """
}

/*
 * Re-key an annotation table onto a different transcript->gene map.
 *
 * Its own process on purpose: folding this into CORSET made corset's 37-minute clustering
 * re-run every time the annotation or the remap changed, for a result corset does not even read.
 * Anything downstream of a gene map can now change without touching the expensive step.
 */
process REMAP_ANNOTATION {
    label 'process_single'
    tag "$label"

    input:
    tuple val(label), path(tx2gene)
    tuple val(label2), path(annotation)

    output:
    tuple val(label), path("annotation.${label}.tsv"), emit: annotation
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    # The enrichment keys on the annotation's gene column while the DE table is keyed on the new
    # map's ids: leave them unmatched and the universe comes out empty.
    awk -F'\t' -v OFS='\t' '
        NR==FNR { g[\$1] = \$2; next }
        FNR==1  { print; next }
        { if (\$1 in g) \$2 = g[\$1]; print }
    ' ${tx2gene} ${annotation} > annotation.${label}.tsv

    echo "[remap] ${label}: \$(( \$(wc -l < annotation.${label}.tsv) - 1 )) rows re-keyed onto \$(cut -f2 ${tx2gene} | sort -u | wc -l) genes"
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch annotation.${label}.tsv
    emit_versions.sh "${task.process}" python3
    """
}
