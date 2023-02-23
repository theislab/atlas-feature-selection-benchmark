/*
========================================================================================
    PROCESSES
========================================================================================
*/

/*
------------------------------
    Integration metrics
------------------------------
*/

process METRIC_BATCHPURITY {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "batchPurity.tsv" }

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
        path(functions)

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

process METRIC_LOCALSTRUCTURE {
    conda "envs/seurat.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "localStructure.tsv" }

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-localStructure.tsv")

    script:
        """
        metric-localStructure.R \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
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

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

process METRIC_KBET {
    conda "envs/scib-kBET.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "kBET.tsv" }

    label "process_low"

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

process METRIC_LABELASW {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "labelASW.tsv" }

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "cLISI.tsv" }

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

process METRIC_ILISI {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "iLISI.tsv" }

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

process METRIC_ISOLATEDLABELSF1 {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "isolatedLabelsF1.tsv" }

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

process METRIC_BATCHPCR {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "batchPCR.tsv" }

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-batchPCR.tsv")

    script:
        """
        metric-batchPCR.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-batchPCR.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-batchPCR.tsv"
        """
}

process METRIC_GRAPHCONNECTIVITY {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "graphConnectivity.tsv" }

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
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

process METRIC_CELLCYCLE {
    conda "envs/scib.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "cellCycle.tsv" }

    input:
        tuple val(dataset), val(method), val(integration), path(reference)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), path("${dataset}-${method}-${integration}-cellCycle.tsv")

    script:
        """
        metric-cellCycle.py \\
            --dataset "${dataset}" \\
            --method "${method}" \\
            --integration "${integration}" \\
            --out-file "${dataset}-${method}-${integration}-cellCycle.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "${dataset}-${method}-${integration}-cellCycle.tsv"
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

    input:
        tuple val(dataset), val(method), val(integration), path("reference.h5ad"), path("query.h5ad")
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

    input:
        tuple val(dataset), val(method), val(integration), path("reference.h5ad"), path("query.h5ad")
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

    input:
        tuple val(dataset), val(method), val(integration), path("reference.h5ad"), path("query.h5ad")
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

    input:
        tuple val(dataset), val(method), val(integration), path("reference.h5ad"), path("query.h5ad")
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

/*
------------------------------
    Classification metrics
------------------------------
*/

process METRIC_ACCURACY {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "accuracy.tsv" }

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

    input:
        tuple val(dataset), val(method), val(integration), path(query), path(labels)
        path(functions)

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

/*
------------------------------
    Unseen metrics
------------------------------
*/

process METRIC_UNSEEN_CELLDIST {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/metrics/${dataset}/${method}/${integration}",
        saveAs: { filename -> "unseenCellDist.tsv" }

    input:
        tuple val(dataset), val(method), val(integration), path("reference.h5ad"), path("query.h5ad")
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

    input:
        tuple val(dataset), val(method), val(integration), path("reference.h5ad"), path("query.h5ad")
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
        path("all-metrics.tsv")

    script:
        """
        combine-metrics.py \\
            --out-file "all-metrics.tsv" \\
            ${metrics}
        """

    stub:
        """
        touch "all-metrics.tsv"
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

    main:

        metric_names = params.metrics.collect{metric -> metric.name}

        // Integration metrics
        batchPurity_ch = metric_names.contains("batchPurity") ?
            METRIC_BATCHPURITY(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        mixing_ch = metric_names.contains("mixing") ?
            METRIC_MIXING(reference_ch, file(params.bindir + "/functions/functions.R")) :
            Channel.empty()
        localStructure_ch = metric_names.contains("localStructure") ?
            METRIC_LOCALSTRUCTURE(reference_ch, file(params.bindir + "/functions/functions.R")) :
            Channel.empty()
        kBET_ch = metric_names.contains("kBET") ?
            METRIC_KBET(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        cLISI_ch = metric_names.contains("cLISI") ?
            METRIC_CLISI(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        ari_ch = metric_names.contains("ari") ?
            METRIC_ARI(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        bARI_ch = metric_names.contains("bARI") ?
            METRIC_BARI(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        nmi_ch = metric_names.contains("nmi") ?
            METRIC_NMI(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        bNMI_ch = metric_names.contains("bNMI") ?
            METRIC_BNMI(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        labelASW_ch = metric_names.contains("labelASW") ?
            METRIC_LABELASW(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        batchPCR_ch = metric_names.contains("batchPCR") ?
            METRIC_BATCHPCR(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
		iLISI_ch = metric_names.contains("iLISI") ?
            METRIC_ILISI(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        graphConnectivity_ch = metric_names.contains("graphConnectivity") ?
            METRIC_GRAPHCONNECTIVITY(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        isolatedLabelsF1_ch = metric_names.contains("isolatedLabelsF1") ?
            METRIC_ISOLATEDLABELSF1(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        isolatedLabelsASW_ch = metric_names.contains("isolatedLabelsASW") ?
            METRIC_ISOLATEDLABELSASW(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        cellCycle_ch = metric_names.contains("cellCycle") ?
            METRIC_CELLCYCLE(reference_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()

        // Mapping metrics
        mLISI_ch = metric_names.contains("mLISI") ?
            METRIC_MLISI(query_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        qLISI_ch = metric_names.contains("qLISI") ?
            METRIC_QLISI(query_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        cellDist_ch = metric_names.contains("cellDist") ?
            METRIC_CELLDIST(
                query_ch,
                file(params.bindir + "/functions/functions.py"),
                file(params.bindir + "/functions/distances.py")
            ) :
            Channel.empty()
        labelDist_ch = metric_names.contains("labelDist") ?
            METRIC_LABELDIST(
                query_ch,
                file(params.bindir + "/functions/functions.py"),
                file(params.bindir + "/functions/distances.py")
            ) :
            Channel.empty()

        // Classification metrics
        accuracy_ch = metric_names.contains("accuracy") ?
            METRIC_ACCURACY(labels_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        rareAccuracy_ch = metric_names.contains("rareAccuracy") ?
            METRIC_RAREACCURACY(labels_ch, file(params.bindir + "/functions/functions.R")) :
            Channel.empty()
        f1_micro_ch = metric_names.contains("f1Micro") ?
            METRIC_F1_MICRO(labels_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        f1_macro_ch = metric_names.contains("f1Macro") ?
            METRIC_F1_MACRO(labels_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        f1_rarity_ch = metric_names.contains("f1Rarity") ?
            METRIC_F1_RARITY(labels_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
		jaccard_micro_ch = metric_names.contains("jaccardIndexMicro") ?
            METRIC_JACCARDINDEX_MICRO(labels_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        jaccard_macro_ch = metric_names.contains("jaccardIndexMacro") ?
            METRIC_JACCARDINDEX_MACRO(labels_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        jaccard_rarity_ch = metric_names.contains("jaccardIndexRarity") ?
            METRIC_JACCARDINDEX_RARITY(labels_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()
        mcc_ch = metric_names.contains("MCC") ?
            METRIC_MCC(labels_ch, file(params.bindir + "/functions/functions.py")) :
            Channel.empty()

        // Unseen metrics
        unseen_cellDist_ch = metric_names.contains("unseenCellDist") ?
            METRIC_UNSEEN_CELLDIST(
                query_ch,
                file(params.bindir + "/functions/functions.py"),
                file(params.bindir + "/functions/distances.py")
            ) :
            Channel.empty()
        unseen_labelDist_ch = metric_names.contains("unseenLabelDist") ?
            METRIC_UNSEEN_LABELDIST(
                query_ch,
                file(params.bindir + "/functions/functions.py"),
                file(params.bindir + "/functions/distances.py")
            ) :
            Channel.empty()

        metrics_ch = batchPurity_ch
            .mix(
                mixing_ch,
                localStructure_ch,
                kBET_ch,
                accuracy_ch,
                rareAccuracy_ch,
                f1_micro_ch,
                f1_macro_ch,
                f1_rarity_ch,
				cLISI_ch,
                ari_ch,
                bARI_ch,
                nmi_ch,
                bNMI_ch,
                labelASW_ch,
				isolatedLabelsF1_ch,
                isolatedLabelsASW_ch,
                jaccard_micro_ch,
                jaccard_macro_ch,
                jaccard_rarity_ch,
                mcc_ch,
				graphConnectivity_ch,
				batchPCR_ch,
				iLISI_ch,
                cellCycle_ch,
                mLISI_ch,
                qLISI_ch,
                cellDist_ch,
                labelDist_ch,
                unseen_cellDist_ch,
                unseen_labelDist_ch
            )
            .map {it -> file(it[3])}
            .toList()

        COMBINE_METRICS(metrics_ch)

    emit:
        combined_metrics_ch = COMBINE_METRICS.out
}

/*
========================================================================================
    THE END
========================================================================================
*/
