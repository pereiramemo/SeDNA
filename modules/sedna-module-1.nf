// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1: BGC annotation with antiSMASH
// Input:  assembly / MAG FASTA (per sample)
// Output: module-1/<sample>/{gbks/,annot_metadata.tsv}
// ─────────────────────────────────────────────────────────────────────────────

process MODULE1 {

    container "ghcr.io/pereiramemo/sedna/module-1:latest"
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
    sedna-module-1.py \
        --input_fasta       ${fasta} \
        --sample_name       ${sample_name} \
        --antismash_db      ${params.antismash_db ?: "${params.db_dir}/antismash"} \
        --taxon             ${params.taxon} \
        --genefinding_tool  ${params.genefinding_tool} \
        --minlength         ${params.minlength} \
        --nslots            ${task.cpus} \
        --output_dir        ${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name} \
        --overwrite
    """
}
