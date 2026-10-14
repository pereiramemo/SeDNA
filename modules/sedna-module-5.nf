// ─────────────────────────────────────────────────────────────────────────────
// MODULE 5: Pfam domain annotation of the dereplicated BGCs
// (pyHMMER + cath-resolve-hits)
// Input:  module-4 dereplicated/<version> dir (per sample; sedna.nf picks the
//         version from --derep_version)
// Output: module-5/<sample>/<bgc_id>_annotdoms_resolved.tsv
// ─────────────────────────────────────────────────────────────────────────────

process MODULE5 {

    container "ghcr.io/pereiramemo/sedna/module-5:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: params.publish_mode,
           enabled: params.full_output.toString().toBoolean() || (params.stop_at_module as Integer) == 5

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(derep_dir)
    // opts: the params this script uses (built in sedna.nf). They come in as an
    // input so that they are part of the task hash: Nextflow does not hash
    // param values read inside the script, and -resume would reuse stale results
    val(opts)

    output:
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}"), emit: domains

    script:
    """
    sedna-module-5.py \
        --input_dir     ${derep_dir} \
        --cds_id        ${opts.cds_id} \
        --hmm_db        ${opts.hmm_db} \
        --evalue_thres  ${opts.evalue_thres} \
        ${opts.cut_ga ? "--cut_ga" : "--no-cut_ga"} \
        --nslots        ${task.cpus} \
        --output_dir    ${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name} \
        --overwrite
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// MODULE 5 (precomputed input, --input_mode from_annot): reuse an existing
// genome-wide hmmsearch domtblout + GFF instead of re-running pyHMMER. Reuses
// the same container/script as MODULE5 and produces the same output shape, so
// MODULE6 cannot tell which entry path ran.
// Input:  module-4 dereplicated/<version> dir + GFF + domtblout (per sample)
// Output: module-5/<sample>/<bgc_id>_annotdoms_resolved.tsv
// ─────────────────────────────────────────────────────────────────────────────

process MODULE5_PRECOMPUTED {

    container "ghcr.io/pereiramemo/sedna/module-5:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: params.publish_mode,
           enabled: params.full_output.toString().toBoolean() || (params.stop_at_module as Integer) == 5

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(derep_dir), path(gff), path(domtblout)

    output:
    path("module-5/${sample_name}"), emit: domains

    script:
    """
    sedna-module-5.py \
        --input_dir   ${derep_dir} \
        --gff         ${gff} \
        --domtblout   ${domtblout} \
        --nslots      ${task.cpus} \
        --output_dir  module-5/${sample_name} \
        --overwrite
    """
}
