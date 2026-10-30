// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2: BGC annotation with deepBGC
// Input:  assembly / MAG FASTA (per sample)
// Output: module-2-deepbgc/<sample>/{gbks/,gbks_removed/,annot_metadata.tsv}
// ─────────────────────────────────────────────────────────────────────────────

process MODULE2_DEEPBGC {

    container "ghcr.io/pereiramemo/sedna/module-2:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: params.publish_mode,
           enabled: params.full_output.toString().toBoolean() || (params.stop_at_module as Integer) <= 3

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(fasta)
    // opts: the params this script uses (built in sedna.nf). They come in as an
    // input so that they are part of the task hash: Nextflow does not hash
    // param values read inside the script, and -resume would reuse stale results
    val(opts)

    output:
    tuple val(sample_name), path("module-2-deepbgc/${sample_name}"), emit: annot

    script:
    """
    sedna-module-2.py \
        --input_fasta      ${fasta} \
        --sample_name      ${sample_name} \
        --deepbgc_db       ${opts.deepbgc_db} \
        --score_thres      ${opts.score_thres} \
        --cds_count_thres  ${opts.cds_count_thres} \
        --minlength        ${opts.minlength} \
        --nslots           ${task.cpus} \
        --output_dir       module-2-deepbgc/${sample_name} \
        --overwrite
    """
}
