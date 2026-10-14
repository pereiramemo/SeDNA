// ─────────────────────────────────────────────────────────────────────────────
// MODULE 6: BGC embedding (SIF or mean over word2vec domain vectors) and BIRCH
// clustering into Gene Cluster Families
// Input:  module-5 output dirs of all samples (collected) + embedding resources
// Output: module-6/{bgc_embeddings-<method>.tsv, bgc_clust-<method>-thres<T>.tsv}
// ─────────────────────────────────────────────────────────────────────────────

process MODULE6 {

    container "ghcr.io/pereiramemo/sedna/module-6:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: params.publish_mode,
           enabled: true

    input:
    // Every module 5 task emits a dir called <sample>; stageAs numbers them
    // (domains_1, domains_2, ...) so they can sit side by side.
    path(domain_dirs, stageAs: "domains_*")
    path(embeddings_file)
    path(counts_file)

    output:
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}")

    script:
    """
    sedna-module-6.py \
        --input_dirs        ${domain_dirs} \
        --embedding_method  ${params.embedding_method} \
        --threshold         ${params.birch_threshold} \
        --alpha             ${params.alpha} \
        --embeddings_file   ${embeddings_file} \
        --counts_file       ${counts_file} \
        --output_dir        ${task.process.toLowerCase().replaceFirst('module', 'module-')} \
        --overwrite
    """
}
