// ─────────────────────────────────────────────────────────────────────────────
// MODULE 3: BGC annotation with GECCO
// Input:  assembly / MAG FASTA (per sample)
// Output: module-3/<sample>/{gbks/,annot_metadata.tsv}
// ─────────────────────────────────────────────────────────────────────────────

process MODULE3 {

    container "ghcr.io/pereiramemo/sedna/module-3:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: params.publish_mode,
           enabled: params.full_output.toString().toBoolean() || (params.stop_at_module as Integer) <= 3

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(fasta)

    output:
    tuple val(sample_name), path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}"), emit: annot

    script:
    """
    sedna-module-3.py \
        --input_fasta      ${fasta} \
        --sample_name      ${sample_name} \
        --score_thres      ${params.gecco_score_thres} \
        --cds_count_thres  ${params.gecco_cds_count_thres} \
        --nslots           ${task.cpus} \
        --output_dir       ${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name} \
        --overwrite
    """
}
