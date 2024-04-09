/*
========================================================================================
    PROCESSES
========================================================================================
*/

/*
------------------------------
    Integration (batch) metrics
------------------------------
*/

process METRIC_BATCHPURITY {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "batchPurity.tsv" }

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-batchPurity.tsv")

    script:
        """
        metric-batchPurity.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-batchPurity.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-batchPurity.tsv"
        """
}

process METRIC_MIXING {
    conda "envs/seurat.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "mixing.tsv" }

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(io_functions)
        path(metric_functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-mixing.tsv")

    script:
        """
        metric-mixing.R \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-mixing.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-mixing.tsv"
        """
}

process METRIC_KBET {
    conda "envs/scib-kBET.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "kBET.tsv" }

    memory { get_memory(reference.size(), "96.GB", task.attempt, "12.GB") }

    time { check_max_time(8.h * task.attempt) }

    label "error_retry"

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-kBET.tsv")

    script:
        """
        metric-kBET-install.R
        metric-kBET.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-kBET.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-kBET.tsv"
        """
}

process METRIC_ILISI {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "iLISI.tsv" }

    label "process_tiny"

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-iLISI.tsv")

    script:
        """
        metric-iLISI.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-iLISI.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-iLISI.tsv"
        """
}

process METRIC_BATCHPCR {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "batchPCR.tsv" }

    memory { get_memory(reference.size(), "6.GB", task.attempt, "12.GB") }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-batchPCR.tsv")

    script:
        """
        metric-batchPCR.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --exprs "${reference_exprs}" \\
            --out-file "${dataset}-${method}-${integration}-batchPCR.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-batchPCR.tsv"
        """
}

process METRIC_BATCHASW {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "batchASW.tsv" }

    memory { get_memory(reference.size(), "6.GB", task.attempt, "12.GB") }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-batchASW.tsv")

    script:
        """
        metric-batchASW.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-batchASW.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-batchASW.tsv"
        """
}

process METRIC_GRAPHCONNECTIVITY {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "graphConnectivity.tsv" }

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-graphConnectivity.tsv")

    script:
        """
        metric-graphConnectivity.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-graphConnectivity.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-graphConnectivity.tsv"
        """
}

process METRIC_CMS {
    conda "envs/CellMixS.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "cms.tsv" }

    label "process_medium"

    memory { get_memory(reference.size(), "12.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(io_functions)
        path(metric_functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-cms.tsv")

    script:
        """
        metric-cms.R \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-cms.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-cms.tsv"
        """
}

/*
------------------------------
    Integration (bio) metrics
------------------------------
*/

process METRIC_LOCALSTRUCTURE {
    conda "envs/seurat.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "localStructure.tsv" }

    memory { get_memory(reference_exprs.size(), "12.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(io_functions)
        path(metric_functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-localStructure.tsv")

    script:
        """
        metric-localStructure.R \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --exprs "${reference_exprs}" \\
            --out-file "${dataset}-${method}-${integration}-localStructure.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-localStructure.tsv"
        """
}

process METRIC_ARI {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "ari.tsv" }

    label "process_low"

    memory { get_memory(reference.size(), "4.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-ari.tsv")

    script:
        """
        metric-ari.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-ari.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-ari.tsv"
        """
}

process METRIC_BARI {
    conda "envs/balanced-clustering.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "bARI.tsv" }

    label "process_low"

    memory { get_memory(reference.size(), "4.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-bARI.tsv")

    script:
        """
        metric-bARI.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-bARI.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-bARI.tsv"
        """
}

process METRIC_NMI {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "nmi.tsv" }

    label "process_low"

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-nmi.tsv")

    script:
        """
        metric-nmi.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-nmi.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-nmi.tsv"
        """
}

process METRIC_BNMI {
    conda "envs/balanced-clustering.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "bNMI.tsv" }

    label "process_low"

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-bNMI.tsv")

    script:
        """
        metric-bNMI.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-bNMI.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-bNMI.tsv"
        """
}

process METRIC_LABELASW {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "labelASW.tsv" }

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-labelASW.tsv")

    script:
        """
        metric-labelASW.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-labelASW.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-labelASW.tsv"
        """
}

process METRIC_CLISI {
    conda "envs/scib.yml"

    label "process_tiny"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "cLISI.tsv" }

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-cLISI.tsv")

    script:
        """
        metric-cLISI.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-cLISI.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-cLISI.tsv"
        """
}

process METRIC_ISOLATEDLABELSF1 {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "isolatedLabelsF1.tsv" }

    label "process_low"

    memory { get_memory(reference.size(), "8.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-isolatedLabelsF1.tsv")

    script:
        """
        metric-isolatedLabels.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --cluster \\
            --out-file "${dataset}-${method}-${integration}-isolatedLabelsF1.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-isolatedLabelsF1.tsv"
        """
}

process METRIC_ISOLATEDLABELSASW {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "isolatedLabelsASW.tsv" }

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-isolatedLabelsASW.tsv")

    script:
        """
        metric-isolatedLabels.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-isolatedLabelsASW.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-isolatedLabelsASW.tsv"
        """
}

process METRIC_CELLCYCLE {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "cellCycle.tsv" }

    memory { get_memory(reference_exprs.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(tirosh_genes)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-cellCycle.tsv")

    script:
        """
        metric-cellCycle.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --exprs "${reference_exprs}" \\
            --cc-genes "${tirosh_genes}" \\
            --out-file "${dataset}-${method}-${integration}-cellCycle.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-cellCycle.tsv"
        """
}

process METRIC_LDFDIFF {
    conda "envs/CellMixS.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "ldfDiff.tsv" }

    label "process_low"

    memory { get_memory(reference_exprs.size(), "12.GB", task.attempt, "12.GB") }

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(reference_exprs)
        path(io_functions)
        path(metric_functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-ldfDiff.tsv")

    script:
        """
        metric-ldfDiff.R \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --exprs "${reference_exprs}" \\
            --out-file "${dataset}-${method}-${integration}-ldfDiff.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-ldfDiff.tsv"
        """
}

/*
------------------------------
    Mapping metrics
------------------------------
*/

process METRIC_MLISI {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "mLISI.tsv" }

    label "process_tiny"

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference, stageAs: "reference.h5ad"), path(query, stageAs: "query.h5ad"), path(query_dir), path(reference_exprs), path(query_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-mLISI.tsv")

    script:
        """
        metric-mLISI.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --reference reference.h5ad \\
            --out-file "${dataset}-${method}-${integration}-mLISI.tsv" \\
            query.h5ad
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-mLISI.tsv"
        """
}

process METRIC_QLISI {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "qLISI.tsv" }

    label "process_tiny"

    memory { get_memory(query.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference, stageAs: "reference.h5ad"), path(query, stageAs: "query.h5ad"), path(query_dir), path(reference_exprs), path(query_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-qLISI.tsv")

    script:
        """
        metric-qLISI.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-qLISI.tsv" \\
            query.h5ad
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-qLISI.tsv"
        """
}

process METRIC_CELLDIST {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "cellDist.tsv" }

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference, stageAs: "reference.h5ad"), path(query, stageAs: "query.h5ad"), path(query_dir), path(reference_exprs), path(query_exprs)
        path(functions)
        path(distance_functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-cellDist.tsv")

    script:
        """
        metric-cellDist.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --reference reference.h5ad \\
            --out-file "${dataset}-${method}-${integration}-cellDist.tsv" \\
            query.h5ad
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-cellDist.tsv"
        """
}

process METRIC_LABELDIST {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "labelDist.tsv" }

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference, stageAs: "reference.h5ad"), path(query, stageAs: "query.h5ad"), path(query_dir), path(reference_exprs), path(query_exprs)
        path(functions)
        path(distance_functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-labelDist.tsv")

    script:
        """
        metric-labelDist.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --reference reference.h5ad \\
            --out-file "${dataset}-${method}-${integration}-labelDist.tsv" \\
            query.h5ad
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-labelDist.tsv"
        """
}

process METRIC_RECONSTRUCTION {
    conda "envs/scvi-tools.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "reconstruction.tsv" }

    label "process_high"

    memory { get_memory(query_exprs.size(), "8.GB", task.attempt, "8.GB") }

    input:
        tuple val(dataset), val(method), val(integration), path("reference.h5ad"), path(query), path(query_dir), path(reference_exprs), path(query_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-reconstruction.tsv")

    script:
        """
        metric-reconstruction.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --exprs "${query_exprs}" \\
            --out-file "${dataset}-${method}-${integration}-reconstruction.tsv" \\
            "${query_dir}/adata.h5ad"
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-reconstruction.tsv"
        """
}

process METRIC_KNNCORR {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "kNNcorr.tsv" }

    memory { get_memory(query_exprs.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path("reference.h5ad"), path("query.h5ad"), path(query_dir), path(reference_exprs), path(query_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-kNNcorr.tsv")

    script:
        """
        metric-kNNcorr.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --exprs "${query_exprs}" \\
            --out-file "${dataset}-${method}-${integration}-kNNcorr.tsv" \\
            query.h5ad
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-kNNcorr.tsv"
        """
}

/*
------------------------------
    Classification metrics
------------------------------
*/

process METRIC_ACCURACY {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "accuracy.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-accuracy.tsv")

    script:
        """
        metric-accuracy.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-accuracy.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-accuracy.tsv"
        """
}

process METRIC_RAREACCURACY {
    conda "envs/tidyverse.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "rareAccuracy.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(io_functions)
        path(metric_functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-rareAccuracy.tsv")

    script:
        """
        metric-rareAccuracy.R \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-rareAccuracy.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-rareAccuracy.tsv"
        """
}

process METRIC_MCC {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "MCC.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-MCC.tsv")

    script:
        """
        metric-MCC.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-MCC.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-MCC.tsv"
        """
}

process METRIC_F1_MICRO {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "F1Micro.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-F1Micro.tsv")

    script:
        """
        metric-f1score.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --average "micro" \\
            --out-file "${dataset}-${method}-${integration}-F1Micro.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-F1Micro.tsv"
        """
}

process METRIC_F1_MACRO {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "F1Macro.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-F1Macro.tsv")

    script:
        """
        metric-f1score.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --average "macro" \\
            --out-file "${dataset}-${method}-${integration}-F1Macro.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-F1Macro.tsv"
        """
}

process METRIC_F1_RARITY {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "F1Rarity.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-F1Rarity.tsv")

    script:
        """
        metric-f1score.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --average "rarity" \\
            --out-file "${dataset}-${method}-${integration}-F1Rarity.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-F1Rarity.tsv"
        """
}

process METRIC_JACCARDINDEX_MICRO {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "JaccardIndexMicro.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-JaccardIndexMicro.tsv")

    script:
        """
        metric-jaccardIndex.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --average "micro" \\
            --out-file "${dataset}-${method}-${integration}-JaccardIndexMicro.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-JaccardIndexMicro.tsv"
        """
}

process METRIC_JACCARDINDEX_MACRO {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "JaccardIndexMacro.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-JaccardIndexMacro.tsv")

    script:
        """
        metric-jaccardIndex.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --average "macro" \\
            --out-file "${dataset}-${method}-${integration}-JaccardIndexMacro.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-JaccardIndexMacro.tsv"
        """
}

process METRIC_JACCARDINDEX_RARITY {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "JaccardIndexRarity.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-JaccardIndexRarity.tsv")

    script:
        """
        metric-jaccardIndex.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --average "rarity" \\
            --out-file "${dataset}-${method}-${integration}-JaccardIndexRarity.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-JaccardIndexRarity.tsv"
        """
}

process METRIC_AUPRC {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "AUPRC.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-AUPRC.tsv")

    script:
        """
        metric-auprc.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-AUPRC.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-AUPRC.tsv"
        """
}


/*
------------------------------
    Unseen metrics
------------------------------
*/

process METRIC_UNSEEN_CELLDIST {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "unseenCellDist.tsv" }

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference, stageAs: "reference.h5ad"), path(query, stageAs: "query.h5ad"), path(query_dir), path(reference_exprs), path(query_exprs)
        path(functions)
        path(distance_functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-unseenCellDist.tsv")

    script:
        """
        metric-unseenCellDist.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --reference reference.h5ad \\
            --out-file "${dataset}-${method}-${integration}-unseenCellDist.tsv" \\
            query.h5ad
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-unseenCellDist.tsv"
        """
}

process METRIC_UNSEEN_LABELDIST {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "unseenLabelDist.tsv" }

    memory { get_memory(reference.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(reference, stageAs: "reference.h5ad"), path(query, stageAs: "query.h5ad"), path(query_dir), path(reference_exprs), path(query_exprs)
        path(functions)
        path(distance_functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-unseenLabelDist.tsv")

    script:
        """
        metric-unseenLabelDist.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --reference reference.h5ad \\
            --out-file "${dataset}-${method}-${integration}-unseenLabelDist.tsv" \\
            query.h5ad
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-unseenLabelDist.tsv"
        """
}

process METRIC_MILO {
    conda "envs/milopy.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "MILO.tsv" }

    label "process_low"

    memory { get_memory(reference.size(), "12.GB", task.attempt, "12.GB") }

    input:
        tuple val(dataset), val(method), val(integration), path(reference, stageAs: "reference.h5ad"), path(query, stageAs: "query.h5ad"), path(query_dir), path(reference_exprs), path(query_exprs)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-MILO.tsv")

    script:
        """
        metric-milo.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --reference reference.h5ad \\
            --out-file "${dataset}-${method}-${integration}-MILO.tsv" \\
            query.h5ad
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-MILO.tsv"
        """
}

process METRIC_UNCERTAINTY {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "uncertainty.tsv" }

    memory { get_memory(labels.size(), "2.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-uncertainty.tsv")

    script:
        """
        metric-uncertainty.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-uncertainty.tsv" \\
            ${labels}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-uncertainty.tsv"
        """
}

/*
------------------------------
    Other processes
------------------------------
*/

process COMBINE_METRICS {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/", mode: "copy"

    input:
        path(metrics)

    output:
        tuple path("all-metrics.tsv"), path("missing-summary.tsv"), path("missing-metrics.txt")

    script:
        """
        combine-metrics.py \\
            --out-file "all-metrics.tsv" \\
            --missing-values "missing-summary.tsv" \\
            --missing-files "missing-metrics.txt" \\
            .
        """

    stub:
        """
        touch "all-metrics.tsv"
        touch "missing-summary.tsv"
        touch "missing-metrics.txt"
        """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow METRICS {

    take:
        reference_ch
        query_ch
        labels_ch
        tirosh_genes_ch

    main:

        metric_names = params.metrics.collect{metric -> metric.name}

        // Function file paths
        py_metrics_funcs = file(params.bindir + "/functions/metrics.py")
        py_distances_funcs = file(params.bindir + "/functions/distances.py")
        r_io_funcs = file(params.bindir + "/functions/io.R")
        r_metrics_funcs = file(params.bindir + "/functions/metrics.R")

        // Integration (batch) metrics
        batchPurity_ch = metric_names.contains("batchPurity") ?
            METRIC_BATCHPURITY(reference_ch, py_metrics_funcs) :
            Channel.empty()
        mixing_ch = metric_names.contains("mixing") ?
            METRIC_MIXING(reference_ch, r_io_funcs, r_metrics_funcs) :
            Channel.empty()
        kBET_ch = metric_names.contains("kBET") ?
            METRIC_KBET(reference_ch, py_metrics_funcs) :
            Channel.empty()
        batchPCR_ch = metric_names.contains("batchPCR") ?
            METRIC_BATCHPCR(reference_ch, py_metrics_funcs) :
            Channel.empty()
        batchASW_ch = metric_names.contains("batchASW") ?
            METRIC_BATCHASW(reference_ch, py_metrics_funcs) :
            Channel.empty()
        iLISI_ch = metric_names.contains("iLISI") ?
            METRIC_ILISI(reference_ch,py_metrics_funcs) :
            Channel.empty()
        graphConnectivity_ch = metric_names.contains("graphConnectivity") ?
            METRIC_GRAPHCONNECTIVITY(reference_ch, py_metrics_funcs) :
            Channel.empty()
        cms_ch = metric_names.contains("CMS") ?
            METRIC_CMS(reference_ch, r_io_funcs, r_metrics_funcs) :
            Channel.empty()

        // Integration (bio) metrics
        localStructure_ch = metric_names.contains("localStructure") ?
            METRIC_LOCALSTRUCTURE(reference_ch, r_io_funcs, r_metrics_funcs) :
            Channel.empty()
        cLISI_ch = metric_names.contains("cLISI") ?
            METRIC_CLISI(reference_ch, py_metrics_funcs) :
            Channel.empty()
        ari_ch = metric_names.contains("ari") ?
            METRIC_ARI(reference_ch, py_metrics_funcs) :
            Channel.empty()
        bARI_ch = metric_names.contains("bARI") ?
            METRIC_BARI(reference_ch, py_metrics_funcs) :
            Channel.empty()
        nmi_ch = metric_names.contains("nmi") ?
            METRIC_NMI(reference_ch, py_metrics_funcs) :
            Channel.empty()
        bNMI_ch = metric_names.contains("bNMI") ?
            METRIC_BNMI(reference_ch, py_metrics_funcs) :
            Channel.empty()
        labelASW_ch = metric_names.contains("labelASW") ?
            METRIC_LABELASW(reference_ch, py_metrics_funcs) :
            Channel.empty()
        isolatedLabelsF1_ch = metric_names.contains("isolatedLabelsF1") ?
            METRIC_ISOLATEDLABELSF1(reference_ch, py_metrics_funcs) :
            Channel.empty()
        isolatedLabelsASW_ch = metric_names.contains("isolatedLabelsASW") ?
            METRIC_ISOLATEDLABELSASW(reference_ch, py_metrics_funcs) :
            Channel.empty()
        cellCycle_ch = metric_names.contains("cellCycle") ?
            METRIC_CELLCYCLE(reference_ch, tirosh_genes_ch, py_metrics_funcs) :
            Channel.empty()
        ldfDiff_ch = metric_names.contains("ldfDiff") ?
            METRIC_LDFDIFF(reference_ch, r_io_funcs, r_metrics_funcs) :
            Channel.empty()

        // Mapping metrics
        mLISI_ch = metric_names.contains("mLISI") ?
            METRIC_MLISI(query_ch, py_metrics_funcs) :
            Channel.empty()
        qLISI_ch = metric_names.contains("qLISI") ?
            METRIC_QLISI(query_ch,py_metrics_funcs) :
            Channel.empty()
        cellDist_ch = metric_names.contains("cellDist") ?
            METRIC_CELLDIST(query_ch, py_metrics_funcs, py_distances_funcs) :
            Channel.empty()
        labelDist_ch = metric_names.contains("labelDist") ?
            METRIC_LABELDIST(query_ch, py_metrics_funcs, py_distances_funcs) :
            Channel.empty()
        reconstruction_ch = metric_names.contains("reconstruction") ?
            METRIC_RECONSTRUCTION(query_ch, py_metrics_funcs) :
            Channel.empty()
        kNNcorr_ch = metric_names.contains("kNNcorr") ?
            METRIC_KNNCORR(query_ch, py_metrics_funcs) :
            Channel.empty()

        // Classification metrics
        accuracy_ch = metric_names.contains("accuracy") ?
            METRIC_ACCURACY(labels_ch, py_metrics_funcs) :
            Channel.empty()
        rareAccuracy_ch = metric_names.contains("rareAccuracy") ?
            METRIC_RAREACCURACY(labels_ch, r_io_funcs, r_metrics_funcs) :
            Channel.empty()
        f1_micro_ch = metric_names.contains("f1Micro") ?
            METRIC_F1_MICRO(labels_ch, py_metrics_funcs) :
            Channel.empty()
        f1_macro_ch = metric_names.contains("f1Macro") ?
            METRIC_F1_MACRO(labels_ch, py_metrics_funcs) :
            Channel.empty()
        f1_rarity_ch = metric_names.contains("f1Rarity") ?
            METRIC_F1_RARITY(labels_ch, py_metrics_funcs) :
            Channel.empty()
		jaccard_micro_ch = metric_names.contains("jaccardIndexMicro") ?
            METRIC_JACCARDINDEX_MICRO(labels_ch, py_metrics_funcs) :
            Channel.empty()
        jaccard_macro_ch = metric_names.contains("jaccardIndexMacro") ?
            METRIC_JACCARDINDEX_MACRO(labels_ch, py_metrics_funcs) :
            Channel.empty()
        jaccard_rarity_ch = metric_names.contains("jaccardIndexRarity") ?
            METRIC_JACCARDINDEX_RARITY(labels_ch, py_metrics_funcs) :
            Channel.empty()
        mcc_ch = metric_names.contains("MCC") ?
            METRIC_MCC(labels_ch, py_metrics_funcs) :
            Channel.empty()
        auprc_ch = metric_names.contains("AUPRC") ?
            METRIC_AUPRC(labels_ch, py_metrics_funcs) :
            Channel.empty()

        // Unseen metrics
        unseen_cellDist_ch = metric_names.contains("unseenCellDist") ?
            METRIC_UNSEEN_CELLDIST(query_ch, py_metrics_funcs, py_distances_funcs) :
            Channel.empty()
        unseen_labelDist_ch = metric_names.contains("unseenLabelDist") ?
            METRIC_UNSEEN_LABELDIST(query_ch, py_metrics_funcs, py_distances_funcs) :
            Channel.empty()
        milo_ch = metric_names.contains("MILO") ?
            METRIC_MILO(query_ch, py_metrics_funcs) :
            Channel.empty()
        uncertainty_ch = metric_names.contains("uncertainty") ?
            METRIC_UNCERTAINTY(labels_ch, py_metrics_funcs) :
            Channel.empty()

        metrics_ch = Channel.empty()
            .mix(
                // Integration (batch) metrics
                batchPurity_ch,
                mixing_ch,
                kBET_ch,
                batchPCR_ch,
                batchASW_ch,
                iLISI_ch,
                graphConnectivity_ch,
                cms_ch,
                // Integration (bio) metrics
                localStructure_ch,
                cLISI_ch,
                ari_ch,
                bARI_ch,
                nmi_ch,
                bNMI_ch,
                labelASW_ch,
                isolatedLabelsF1_ch,
                isolatedLabelsASW_ch,
                cellCycle_ch,
                ldfDiff_ch,
                // Mapping metrics
                mLISI_ch,
                qLISI_ch,
                cellDist_ch,
                labelDist_ch,
                reconstruction_ch,
                kNNcorr_ch,
                // Classification metrics
                accuracy_ch,
                rareAccuracy_ch,
                f1_micro_ch,
                f1_macro_ch,
                f1_rarity_ch,
                jaccard_micro_ch,
                jaccard_macro_ch,
                jaccard_rarity_ch,
                mcc_ch,
                auprc_ch,
                // Unseen metrics
                unseen_cellDist_ch,
                unseen_labelDist_ch,
                milo_ch,
                uncertainty_ch
            )
            .map {it -> file(it[3])}
            .toList()

        COMBINE_METRICS(metrics_ch)

    emit:
        combined_metrics_ch = COMBINE_METRICS.out
            .map {it -> file(it[0])}
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def get_memory(file_size, mem_per_gb = "1.GB", attempt = 1, overhead = "8.GB") {
    file_mem = new nextflow.util.MemoryUnit(file_size)
    mem_per_gb = new nextflow.util.MemoryUnit(mem_per_gb)
    overhead = new nextflow.util.MemoryUnit(overhead)
    max_mem = new nextflow.util.MemoryUnit(params.max_memory)

    file_gb = file_mem.toGiga()
    file_gb = file_gb > 0 ? file_gb : 1
    mem_use = (mem_per_gb * file_gb * attempt) + overhead

    if (mem_use > max_mem) {
        mem_use = max_mem
    }

    return mem_use
}

def check_max_time(obj) {
    try {
        if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
            return params.max_time as nextflow.util.Duration
        else
            return obj
    } catch (all) {
        println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
        return obj
    }
}


/*
========================================================================================
    THE END
========================================================================================
*/
