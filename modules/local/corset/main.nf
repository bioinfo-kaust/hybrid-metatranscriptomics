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
    tag "${meta.id}"
    label 'process_medium'
    label 'error_ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/corset_gawk_gzip:72b219419508ba4f'
        : 'community.wave.seqera.io/library/corset_gawk_gzip:0f5abd91c07a51dc' }"

    input:
    path eq_classes, stageAs: 'eq/*'
    tuple val(meta), path(tx2gene)
    path samplesheet

    output:
    tuple val(meta_out), path("tx2gene.${meta.id}_corset.tsv"), emit: tx2gene
    path "corset_${meta.id}"                                  , emit: results
    path "${meta.id}_corset.summary.txt"                      , emit: summary
    path "versions.yml"                                       , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    meta_out  = meta + [id: "${meta.id}_corset"]
    def label = meta.id
    def bin   = params.corset_bin ? "export PATH=${params.corset_bin}:\$PATH" : ''
    def args  = task.ext.args ?: ''
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
    meta_out  = meta + [id: "${meta.id}_corset"]
    """
    mkdir -p corset_${meta.id}
    touch tx2gene.${meta.id}_corset.tsv ${meta.id}_corset.summary.txt
    emit_versions.sh "${task.process}" corset
    """
}
