// ─────────────────────────────────────────────────────────────────────────────
// MODULE 5: Pfam domain annotation of the dereplicated BGCs
// (pyHMMER + cath-resolve-hits)
// Input:  module-4 output dir (per sample)
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

    output:
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}"), emit: domains

    script:
    def cut_ga = params.cut_ga.toBoolean() ? "--cut_ga" : "--no-cut_ga"
    """
    sedna-module-5.py \
        --input_dir     ${derep_dir}/dereplicated \
        --cds_id        ${params.cds_id} \
        --hmm_db        ${params.pfam_hmm ?: "${params.db_dir}/pfam/Pfam-A.hmm"} \
        --evalue_thres  ${params.evalue_thres} \
        ${cut_ga} \
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
// Input:  module-4 output dir + GFF + domtblout (per sample)
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
        --input_dir   ${derep_dir}/dereplicated \
        --gff         ${gff} \
        --domtblout   ${domtblout} \
        --nslots      ${task.cpus} \
        --output_dir  module-5/${sample_name} \
        --overwrite
    """
}
