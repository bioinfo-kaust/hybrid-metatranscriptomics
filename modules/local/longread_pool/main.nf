/* Pool the (usually few, unreplicated) long-read samples into one reference-building input. */
process LONGREAD_POOL {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/seqkit_gawk:8a03f2b8f8f6d92a'
        : 'community.wave.seqera.io/library/seqkit_gawk:787aa427251d2eb0' }"

    input:
    path(fastas, stageAs: 'in/*')

    output:
    tuple val(meta_out), path('longreads_pooled.fasta'), emit: seqs
    path "versions.yml"                                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    meta_out = [id: 'longreads']
    """
    cat in/* > pooled.raw.fasta
    # duplicate names across pooled samples would silently collide downstream
    if seqkit rmdup -n pooled.raw.fasta -o longreads_pooled.fasta 2> rmdup.log; then
        rm -f pooled.raw.fasta
    else
        mv pooled.raw.fasta longreads_pooled.fasta
    fi
    cat rmdup.log
    echo "[pooled] long reads: \$(grep -c '^>' longreads_pooled.fasta || true)"
    emit_versions.sh "${task.process}" seqkit
    """

    stub:
    meta_out = [id: 'longreads']
    """
    touch longreads_pooled.fasta
    emit_versions.sh "${task.process}" seqkit
    """
}
