/*
 * SQANTI3 structural QC of the long-read models against the short-read annotation.
 *
 * The categories (FSM / ISM / NIC / NNC / antisense / intergenic) and the RT-switching and
 * intra-priming flags say which collapse products are supported by an independent assembly and
 * which look like artefacts. This reports; it does not filter the reference.
 */
process SQANTI3_QC {
    tag "${meta.id}"
    label 'process_high'
    label 'error_ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/sqanti3:6.0.2--7cea524dfbe448e2'
        : 'community.wave.seqera.io/library/sqanti3:6.0.2--82505825316f6ee9' }"

    input:
    tuple val(meta), path(isoforms_gtf)
    tuple val(meta2), path(reference_gtf)
    path genome

    output:
    path "sqanti3_${meta.id}"            , emit: results, optional: true
    path "sqanti3_${meta.id}.summary.txt", emit: summary
    path "versions.yml"                  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: '--skipORF --report skip'
    def label = meta.id
    """
    # A site install can be shadowed by ~/.local site-packages; ignore user site-packages.
    export PYTHONNOUSERSITE=1

    # SQANTI3 >= 5.3 takes named inputs; earlier releases are positional and need -g for a GTF
    if sqanti3_qc.py --help 2>&1 | grep -q -- '--isoforms'; then
        INPUTS="--isoforms ${isoforms_gtf} --refGTF ${reference_gtf} --refFasta ${genome}"
    else
        INPUTS="${isoforms_gtf} ${reference_gtf} ${genome} -g"
    fi
    sqanti3_qc.py \$INPUTS -t ${task.cpus} -d sqanti3_${label} -o ${label} ${args} > sqanti3.log 2>&1 \\
      || echo "[WARN] sqanti3_qc.py exited non-zero, see sqanti3.log"

    CLS=\$(find sqanti3_${label} -maxdepth 1 -name '*classification.txt' 2>/dev/null | sort | sed -n 1p)
    {
      echo "SQANTI3 structural categories — ${label} vs the short-read annotation"
      if [ -n "\$CLS" ] && [ -s "\$CLS" ]; then
          awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) h[\$i]=i; next}
                       { c[\$h["structural_category"]]++; n++ }
                       END { for (k in c) printf "  %-28s %8d  %5.1f%%\\n", k, c[k], 100*c[k]/n;
                             printf "  %-28s %8d\\n", "total", n }' "\$CLS" | sort -k2,2nr
      else
          echo "  [WARN] SQANTI3 produced no classification table — see sqanti3.log"
          tail -20 sqanti3.log
      fi
    } | tee sqanti3_${label}.summary.txt
    emit_versions.sh "${task.process}" sqanti3_qc.py
    """

    stub:
    """
    mkdir -p sqanti3_${meta.id} && touch sqanti3_${meta.id}.summary.txt
    emit_versions.sh "${task.process}" sqanti3_qc.py
    """
}
