/*
 * salmon runs with -l A and auto-detects the library type; this CHECKS the detection against
 * --strandedness rather than printing it. A per-sample flip (a swapped R1/R2, a mislabelled
 * library) otherwise passes silently and every fold change for that sample is wrong. The quant
 * directory passes through, so nothing downstream sees a sample that failed the check.
 */
process SALMON_LIBTYPE_CHECK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(quant)
    val strandedness
    val enforce

    output:
    tuple val(meta), path("${quant}", includeInputs: true), emit: quant
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def expect = strandedness == 'reverse' ? 'ISR'
               : strandedness == 'forward' ? 'ISF'
               : 'IU'
    """
    python3 - <<'PY_CHECK'
import json, sys
d = json.load(open('${quant}/aux_info/meta_info.json'))
lib, mapped = d['library_types'][0], d['percent_mapped']
print('[${meta.id}] library %s, mapped %.2f%%' % (lib, mapped))
if '${enforce}' == 'true' and lib != '${expect}':
    sys.exit("[ERROR] ${meta.id}: salmon detected library type %s but --strandedness "
             "'${strandedness}' implies ${expect}. Fix the samplesheet or set "
             "--strandedness_check false to override." % lib)
PY_CHECK
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    emit_versions.sh "${task.process}" python3
    """
}
