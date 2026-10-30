// ─────────────────────────────────────────────────────────────────────────────
// MODULE 4: Dereplication of the BGCs predicted by modules 1-3
// Input:  module-1-antismash, module-2-deepbgc and module-3-gecco output dirs of one sample (whichever ran)
// Output: module-4-dereplication/<sample>/dereplicated/{longest,core,extended}/ plus the audit dirs
//         {overlapped,partially_overlapped,non_overlapped}/
//         output/<sample>: link to module-4-dereplication/<sample>/dereplicated (workflow.onComplete
//         in nextflow.config)
// ─────────────────────────────────────────────────────────────────────────────

process MODULE4_DEREPLICATION {

    container "ghcr.io/pereiramemo/sedna/module-4:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: params.publish_mode,
           enabled: params.full_output.toString().toBoolean() || (params.stop_at_module as Integer) == 4

    tag "${sample_name}"

    input:
    // Each of modules 1-3 emits a dir called <sample>; stageAs numbers them
    // (tool_1, tool_2, ...) in the order sedna.nf sorted them: antiSMASH, deepBGC, GECCO.
    // That order fixes the pairwise merge order, which can change which of two
    // equal-length BGCs is kept.
    tuple val(sample_name), path(fasta), path(tool_dirs, stageAs: "tool_*")
    path(synonyms_file)
    // opts: the params this script uses (built in sedna.nf). They come in as an
    // input so that they are part of the task hash: Nextflow does not hash
    // param values read inside the script, and -resume would reuse stale results
    val(opts)

    output:
    tuple val(sample_name), path("module-4-dereplication/${sample_name}"), emit: derep

    script:
    def input_dirs = [tool_dirs].flatten().join(' ')
    """
    sedna-module-4.py \
        --input_dirs     ${input_dirs} \
        --input_fasta    ${fasta} \
        --overlap_thres  ${opts.overlap_thres} \
        --synonyms_file  ${synonyms_file} \
        --output_dir     module-4-dereplication/${sample_name} \
        --overwrite
    """
}
