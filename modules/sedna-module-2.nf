// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2: BGC annotation with deepBGC
// Input:  assembly / MAG FASTA (per sample)
// Output: module-2/<sample>/{gbks/,gbks_removed/,annot_metadata.tsv}
// ─────────────────────────────────────────────────────────────────────────────

process MODULE2 {

    container "ghcr.io/pereiramemo/sedna/module-2:latest"
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
    sedna-module-2.py \
        --input_fasta      ${fasta} \
        --sample_name      ${sample_name} \
        --deepbgc_db       ${params.deepbgc_db ?: "${params.db_dir}/deepbgc/data"} \
        --score_thres      ${params.deepbgc_score_thres} \
        --cds_count_thres  ${params.deepbgc_cds_count_thres} \
        --minlength        ${params.minlength} \
        --nslots           ${task.cpus} \
        --output_dir       ${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name} \
        --overwrite
    """
}
