
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

process METRICS_REPORT {
    conda "envs/tidyverse.yml"

    publishDir "$params.outdir/metrics/", mode: "copy"

    stageInMode "copy"

    input:
        tuple path(all_metrics), path(rmd), path(functions)

    output:
        path("metrics.html")

    script:
        """
        render-rmarkdown.R \\
            --params "metrics_file=${all_metrics},functions_file=${functions}" \\
            --out-file "metrics.html" \\
            ${rmd}
        """

    stub:
        """
        touch "metrics.html"
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

    main:

        // Integration metrics
        METRIC_BATCHPURITY(reference_ch)
        METRIC_MIXING(reference_ch)

        // Classifcation metrics
        METRIC_ACCURACY(query_ch)
        METRIC_RAREACCURACY(query_ch)

        metrics_ch = METRIC_BATCHPURITY.out
            .mix(
                METRIC_MIXING.out,
                METRIC_ACCURACY.out,
                METRIC_RAREACCURACY.out
            )
            .map {it -> file(it[3])}
            .toList()

        COMBINE_METRICS(metrics_ch)

        report_ch = COMBINE_METRICS.out
            .map {it ->
                tuple(
                    it,
                    file(params.reportsdir + "/metrics.Rmd"),
                    file(params.reportsdir + "/functions.R")
                )
            }
        METRICS_REPORT(report_ch)
}

/*
========================================================================================
    THE END
========================================================================================
*/
