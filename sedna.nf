#!/usr/bin/env -S nextflow run

// Include modules
include { MODULE1 } from './modules/sedna-module-1.nf'
include { MODULE2 } from './modules/sedna-module-2.nf'
include { MODULE3 } from './modules/sedna-module-3.nf'
include { MODULE4 } from './modules/sedna-module-4.nf'
include { MODULE5; MODULE5_PRECOMPUTED } from './modules/sedna-module-5.nf'
include { MODULE6 } from './modules/sedna-module-6.nf'

workflow {

    main:
    if (params.help) {
        log.info """
        SeDNA: Biosynthetic Gene Cluster bioprospecting of metagenomic assemblies and MAGs

        Usage: nextflow run sedna.nf [options]

        General:
          --output_dir      DIR   Output directory (default: ${params.output_dir})
          --nslots          INT   CPU threads per per-sample task (default: ${params.nslots})
          --stop_at_module  INT   Stop after module N, 1-6; 1-3 all mean "after BGC annotation" (default: ${params.stop_at_module})
          --full_output     BOOL  Publish all intermediate outputs (default: ${params.full_output})
          --publish_mode    STR   publishDir mode: copy | symlink | rellink | link | move (default: ${params.publish_mode})
          --maxForks        INT   Max per-sample tasks run in parallel (default: ${params.maxForks})
          --db_dir          DIR   Root of the external databases, bind-mounted into the containers (default: ${params.db_dir})

        Threading: modules 1-5 run one task per sample, each using --nslots
        threads; --maxForks caps each module separately, not the pipeline as a
        whole. Modules 1, 2 and 3 are independent branches and run at the same
        time. Module 6 runs a single task over all samples, taking
        --nslots_collected threads (${params.nslots_collected}).

        Input mode (selects which TSV samplesheet is read):
          --input_mode      STR   from_fasta | from_annot (default: ${params.input_mode})
          --from_fasta_tsv  TSV   [from_fasta]  sample_name, fasta                    (default: ${params.from_fasta_tsv})
          --from_annot_tsv  TSV   [from_annot]  sample_name, fasta, gff, domtblout    (default: ${params.from_annot_tsv})
                                  from_annot reuses a genome-wide hmmsearch domtblout in MODULE5
                                  instead of running pyHMMER

        MODULES 1-3 — BGC annotation:
          --skip_antismash  BOOL  Skip MODULE1 antiSMASH (default: ${params.skip_antismash})
          --skip_deepbgc    BOOL  Skip MODULE2 deepBGC (default: ${params.skip_deepbgc})
          --skip_gecco      BOOL  Skip MODULE3 GECCO (default: ${params.skip_gecco})
          --minlength       INT   Minimum sequence length, antiSMASH and deepBGC (default: ${params.minlength})
          --taxon           STR   antiSMASH taxon: bacteria | fungi (default: ${params.taxon})
          --genefinding_tool STR  antiSMASH gene finder (default: ${params.genefinding_tool})
          --antismash_db    PATH  antiSMASH databases dir (default: <db_dir>/antismash)
          --deepbgc_db      PATH  deepBGC downloads dir (default: <db_dir>/deepbgc/data)
          --deepbgc_score_thres     NUM  Minimum deepBGC score (default: ${params.deepbgc_score_thres})
          --deepbgc_cds_count_thres INT  Minimum CDSs per deepBGC BGC (default: ${params.deepbgc_cds_count_thres})
          --gecco_score_thres       NUM  GECCO probability threshold (default: ${params.gecco_score_thres})
          --gecco_cds_count_thres   INT  Minimum CDSs per GECCO BGC (default: ${params.gecco_cds_count_thres})

        MODULE4 — Dereplication:
          --overlap_thres   NUM   Overlap fraction separating overlapped from partially overlapped (default: ${params.overlap_thres})
          --synonyms_file   PATH  BGC class synonyms table (default: ${params.synonyms_file})
          --derep_version   STR   Catalog version fed to MODULE5-6 (default: ${params.derep_version}):
                                    longest  - the longest prediction of each locus
                                    core     - the region supported by the most tools
                                               (all tools' intersection when they all overlap)
                                    extended - leftmost start to rightmost end of all predictions
                                  All three are always written to module-4/<sample>/dereplicated/

        MODULE5 — Domain annotation:
          --pfam_hmm        PATH  Pfam-A HMM file (default: <db_dir>/pfam/Pfam-A.hmm)
          --cds_id          STR   GBK qualifier used as CDS ID (default: ${params.cds_id})
          --evalue_thres    NUM   hmmsearch e-value threshold (default: ${params.evalue_thres})
          --cut_ga          BOOL  Use Pfam gathering cutoffs (default: ${params.cut_ga})

        MODULE6 — Embedding + clustering:
          --embedding_method STR  sif | mean (default: ${params.embedding_method})
          --birch_threshold NUM   BIRCH threshold, i.e. GCF radius (default: ${params.birch_threshold})
          --alpha           NUM   SIF a parameter (default: ${params.alpha})
          --embeddings_file PATH  Domain embeddings pkl (default: ${params.embeddings_file})
          --counts_file     PATH  Domain counts pkl (default: ${params.counts_file})
        """.stripIndent()
        exit 0
    }

    stop = params.stop_at_module as Integer
    if (!(params.derep_version in ["longest", "core", "extended"])) {
        error "Unknown --derep_version '${params.derep_version}'; expected one of: longest, core, extended"
    }
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
    module5_opts = [hmm_db:           (params.pfam_hmm ?: "${params.db_dir}/pfam/Pfam-A.hmm").toString(),
                    cds_id:           params.cds_id.toString(),
                    evalue_thres:     params.evalue_thres.toString(),
                    cut_ga:           params.cut_ga.toString().toBoolean()]
    module6_opts = [embedding_method: params.embedding_method.toString(),
                    threshold:        params.birch_threshold.toString(),
                    alpha:            params.alpha.toString()]

    // Input samplesheet, selected by --input_mode
    if (params.input_mode == "from_fasta") {
        samplesheet = params.from_fasta_tsv
    } else if (params.input_mode == "from_annot") {
        if (!params.from_annot_tsv) {
            error "--from_annot_tsv is required when --input_mode is 'from_annot' " +
                  "(TSV columns: sample_name, fasta, gff, domtblout)"
        }
        samplesheet = params.from_annot_tsv
    } else {
        error "Unknown --input_mode '${params.input_mode}'; expected one of: from_fasta, from_annot"
    }

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
        annot_ch = annot_ch.mix(MODULE1(fasta_ch, module1_opts).annot.map { s, d -> tuple(s, 1, d) })
        n_tools += 1
    }
    if (!params.skip_deepbgc.toBoolean()) {
        annot_ch = annot_ch.mix(MODULE2(fasta_ch, module2_opts).annot.map { s, d -> tuple(s, 2, d) })
        n_tools += 1
    }
    if (!params.skip_gecco.toBoolean()) {
        annot_ch = annot_ch.mix(MODULE3(fasta_ch, module3_opts).annot.map { s, d -> tuple(s, 3, d) })
        n_tools += 1
    }
    if (n_tools == 0) {
        error "--skip_antismash, --skip_deepbgc and --skip_gecco are all set; at least one BGC annotation tool must run"
    }

    // MODULE4: Dereplicate the per-tool predictions of each sample. groupTuple
    // with size: releases a sample as soon as all its active tools finished; the
    // dirs are then put back in tool order (antiSMASH, deepBGC, GECCO), which is
    // the pairwise merge order and must not depend on which tool finished first.
    if (stop >= 4) {
        module4_in_ch = annot_ch
            .groupTuple(by: 0, size: n_tools)
            .map { sample_name, ranks, dirs ->
                tuple(sample_name, [ranks, dirs].transpose().sort { it[0] }.collect { it[1] }) }
            .join(fasta_ch, by: 0)
            .map { sample_name, dirs, fasta -> tuple(sample_name, fasta, dirs) }
            // -> tuple(sample_name, fasta, [tool dirs in tool order])
        module4_out   = MODULE4(module4_in_ch, file(params.synonyms_file, checkIfExists: true), module4_opts)
    }

    // MODULE5: Pfam domain annotation of each dereplicated catalog, or -- when
    // --input_mode from_annot supplied gff/domtblout -- import of those
    if (stop >= 5) {
        // The catalog version is picked here, as an input path, not inside the
        // process script: Nextflow's task hash ignores param values used only in
        // the script, so -resume would otherwise reuse another version's results
        derep_version_ch = module4_out.derep
            .map { sample_name, derep_dir ->
                tuple(sample_name, derep_dir.resolve("dereplicated/${params.derep_version}")) }

        if (params.input_mode == "from_annot") {

            precomputed_annot_ch = channel.fromPath(params.from_annot_tsv, checkIfExists: true)
                .splitCsv(header: true, sep: '\t')
                .map { row -> tuple(row.sample_name,
                                    file(row.gff, checkIfExists: true),
                                    file(row.domtblout, checkIfExists: true)) }

            module5_in_ch = derep_version_ch
                .join(precomputed_annot_ch, by: 0) // -> tuple(sample_name, derep_dir, gff, domtblout)
            module5_out = MODULE5_PRECOMPUTED(module5_in_ch)

        } else {
            module5_out = MODULE5(derep_version_ch, module5_opts)
        }
    }

    // MODULE6: Embed and cluster all BGCs of all samples into GCFs
    if (stop >= 6) {
        module6_out = MODULE6(module5_out.domains.collect(),
                              file(params.embeddings_file, checkIfExists: true),
                              file(params.counts_file, checkIfExists: true),
                              module6_opts)
    }

}
