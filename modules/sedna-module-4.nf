// ─────────────────────────────────────────────────────────────────────────────
// MODULE 4: Dereplication of the BGCs predicted by modules 1-3
// Input:  module-1/2/3 output dirs of one sample (whichever ran)
// Output: module-4/<sample>/dereplicated/{longest,core,extended}/ plus the audit dirs
//         {overlapped,partially_overlapped,non_overlapped}/
// ─────────────────────────────────────────────────────────────────────────────

process MODULE4 {

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

    output:
    tuple val(sample_name), path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}"), emit: derep

    script:
    def input_dirs = [tool_dirs].flatten().join(' ')
    """
    sedna-module-4.py \
        --input_dirs     ${input_dirs} \
        --input_fasta    ${fasta} \
        --overlap_thres  ${params.overlap_thres} \
        --synonyms_file  ${synonyms_file} \
        --output_dir     ${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name} \
        --overwrite
    """
}
