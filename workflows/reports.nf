/*
========================================================================================
    PROCESSES
========================================================================================
*/

process METRICS_REPORT {
    conda "envs/tidyverse.yml"

    publishDir "$params.outdir/reports/", mode: "copy"

    stageInMode "copy"

    input:
        path(all_metrics)
        path(rmd)
        path(functions)

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

workflow REPORTS {

    take:
        combined_metrics_ch

    main:

        METRICS_REPORT(
            combined_metrics_ch,
            file(params.reportsdir + "/metrics.Rmd"),
            file(params.reportsdir + "/functions.R")
        )
}

/*
========================================================================================
    THE END
========================================================================================
*/
