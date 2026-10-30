// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1: BGC annotation with antiSMASH
// Input:  assembly / MAG FASTA (per sample)
// Output: module-1-antismash/<sample>/{gbks/,annot_metadata.tsv}
// ─────────────────────────────────────────────────────────────────────────────

process MODULE1_ANTISMASH {

    container "ghcr.io/pereiramemo/sedna/module-1:latest"
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
    tuple val(sample_name), path("module-1-antismash/${sample_name}"), emit: annot

    script:
    """
    sedna-module-1.py \
        --input_fasta       ${fasta} \
        --sample_name       ${sample_name} \
        --antismash_db      ${opts.antismash_db} \
        --taxon             ${opts.taxon} \
        --genefinding_tool  ${opts.genefinding_tool} \
        --minlength         ${opts.minlength} \
        --nslots            ${task.cpus} \
        --output_dir        module-1-antismash/${sample_name} \
        --overwrite
    """
}
