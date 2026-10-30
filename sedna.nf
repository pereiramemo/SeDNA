#!/usr/bin/env -S nextflow run

// Include modules
include { MODULE1_ANTISMASH } from './modules/sedna-module-1.nf'
include { MODULE2_DEEPBGC } from './modules/sedna-module-2.nf'
include { MODULE3_GECCO } from './modules/sedna-module-3.nf'
include { MODULE4_DEREPLICATION } from './modules/sedna-module-4.nf'

workflow {

    main:
    if (params.help) {
        log.info """
        SeDNA: Biosynthetic Gene Cluster bioprospecting of metagenomic assemblies and MAGs

        Usage: nextflow run sedna.nf [options]

        General:
          --output_dir      DIR   Output directory (default: ${params.output_dir})
          --nslots          INT   CPU threads per per-sample task (default: ${params.nslots})
          --stop_at_module  INT   Stop after module N, 1-4; 1-3 all mean "after BGC annotation" (default: ${params.stop_at_module})
          --full_output     BOOL  Publish all intermediate outputs (default: ${params.full_output})
          --publish_mode    STR   publishDir mode: copy | symlink | rellink | link | move (default: ${params.publish_mode})
          --maxForks        INT   Max per-sample tasks run in parallel (default: ${params.maxForks})
          --db_dir          DIR   Root of the external databases, bind-mounted into the containers (default: ${params.db_dir})

        Threading: modules 1-4 run one task per sample, each using --nslots
        threads; --maxForks caps each module separately, not the pipeline as a
        whole. Modules 1, 2 and 3 are independent branches and run at the same
        time.

        Input:
          --from_fasta_tsv  TSV   sample_name, fasta (default: ${params.from_fasta_tsv})

        MODULES 1-3 — BGC annotation:
          --skip_antismash  BOOL  Skip MODULE1_ANTISMASH (default: ${params.skip_antismash})
          --skip_deepbgc    BOOL  Skip MODULE2_DEEPBGC (default: ${params.skip_deepbgc})
          --skip_gecco      BOOL  Skip MODULE3_GECCO (default: ${params.skip_gecco})
          --minlength       INT   Minimum sequence length, antiSMASH and deepBGC (default: ${params.minlength})
          --taxon           STR   antiSMASH taxon: bacteria | fungi (default: ${params.taxon})
          --genefinding_tool STR  antiSMASH gene finder (default: ${params.genefinding_tool})
          --antismash_db    PATH  antiSMASH databases dir (default: <db_dir>/antismash)
          --deepbgc_db      PATH  deepBGC downloads dir (default: <db_dir>/deepbgc/data)
          --deepbgc_score_thres     NUM  Minimum deepBGC score (default: ${params.deepbgc_score_thres})
          --deepbgc_cds_count_thres INT  Minimum CDSs per deepBGC BGC (default: ${params.deepbgc_cds_count_thres})
          --gecco_score_thres       NUM  GECCO probability threshold (default: ${params.gecco_score_thres})
          --gecco_cds_count_thres   INT  Minimum CDSs per GECCO BGC (default: ${params.gecco_cds_count_thres})

        MODULE4_DEREPLICATION — Dereplication:
          --overlap_thres   NUM   Overlap fraction separating overlapped from partially overlapped (default: ${params.overlap_thres})
          --synonyms_file   PATH  BGC class synonyms table (default: ${params.synonyms_file})
        """.stripIndent()
        exit 0
    }

    stop = params.stop_at_module as Integer
    log.info "Pipeline will run up to module ${stop}"

    // Per-module options. Each process receives the params its script uses as a
    // `val` input instead of reading params directly: Nextflow's task hash covers
    // the script text and the inputs but not param values, so a param read in
    // the script would not invalidate -resume (bugs/003)
    module1_opts = [antismash_db:     (params.antismash_db ?: "${params.db_dir}/antismash").toString(),
                    taxon:            params.taxon.toString(),
                    genefinding_tool: params.genefinding_tool.toString(),
                    minlength:        params.minlength.toString()]
    module2_opts = [deepbgc_db:       (params.deepbgc_db ?: "${params.db_dir}/deepbgc/data").toString(),
                    score_thres:      params.deepbgc_score_thres.toString(),
                    cds_count_thres:  params.deepbgc_cds_count_thres.toString(),
                    minlength:        params.minlength.toString()]
    module3_opts = [score_thres:      params.gecco_score_thres.toString(),
                    cds_count_thres:  params.gecco_cds_count_thres.toString()]
    module4_opts = [overlap_thres:    params.overlap_thres.toString()]

    // Input samplesheet
    samplesheet = params.from_fasta_tsv

    fasta_ch = channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_name, file(row.fasta, checkIfExists: true)) }

    n_samples = file(samplesheet).readLines().findAll { it.trim() }.size() - 1
    log.info "Pipeline will run on ${n_samples} samples from ${samplesheet}"
    log.info "Outputs will be written to ${params.output_dir}"

    // MODULES 1-3: BGC annotation with antiSMASH, deepBGC and GECCO. The three
    // are independent branches of the DAG and run concurrently (always run,
    // unless skipped)
    annot_ch = channel.empty()
    n_tools = 0
    if (!params.skip_antismash.toBoolean()) {
        annot_ch = annot_ch.mix(MODULE1_ANTISMASH(fasta_ch, module1_opts).annot.map { s, d -> tuple(s, 1, d) })
        n_tools += 1
    }
    if (!params.skip_deepbgc.toBoolean()) {
        annot_ch = annot_ch.mix(MODULE2_DEEPBGC(fasta_ch, module2_opts).annot.map { s, d -> tuple(s, 2, d) })
        n_tools += 1
    }
    if (!params.skip_gecco.toBoolean()) {
        annot_ch = annot_ch.mix(MODULE3_GECCO(fasta_ch, module3_opts).annot.map { s, d -> tuple(s, 3, d) })
        n_tools += 1
    }
    if (n_tools == 0) {
        error "--skip_antismash, --skip_deepbgc and --skip_gecco are all set; at least one BGC annotation tool must run"
    }

    // MODULE4_DEREPLICATION: Dereplicate the per-tool predictions of each
    // sample. groupTuple with size: releases a sample as soon as all its active
    // tools finished; the dirs are then put back in tool order (antiSMASH,
    // deepBGC, GECCO), which is the pairwise merge order and must not depend on
    // which tool finished first.
    if (stop >= 4) {
        module4_in_ch = annot_ch
            .groupTuple(by: 0, size: n_tools)
            .map { sample_name, ranks, dirs ->
                tuple(sample_name, [ranks, dirs].transpose().sort { it[0] }.collect { it[1] }) }
            .join(fasta_ch, by: 0)
            .map { sample_name, dirs, fasta -> tuple(sample_name, fasta, dirs) }
            // -> tuple(sample_name, fasta, [tool dirs in tool order])
        module4_out   = MODULE4_DEREPLICATION(module4_in_ch, file(params.synonyms_file, checkIfExists: true), module4_opts)
    }

}
