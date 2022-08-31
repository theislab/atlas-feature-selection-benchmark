
/*
========================================================================================
    PROCESSES
========================================================================================
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
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow METRICS {

    take:
        reference_ch
        query_ch

    main:

        METRIC_BATCHPURITY(reference_ch)
        METRIC_MIXING(reference_ch)
}

/*
========================================================================================
    THE END
========================================================================================
*/
