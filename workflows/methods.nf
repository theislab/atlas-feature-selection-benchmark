
/*
========================================================================================
    PROCESSES
========================================================================================
*/

process METHOD_ALL {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${name}", mode: "copy"

    input:
        tuple val(name), path(reference), path(query)

    output:
        tuple val(name), val("all"), path("all.tsv")

    script:
        """
        method-all.py \\
            --out-file "all.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "all.tsv"
        """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow METHODS {

    take:
        prepared_datasets_ch

    main:

        METHOD_ALL(prepared_datasets_ch)

        selected_features_ch = METHOD_ALL.out

    emit:
        selected_features_ch
}

/*
========================================================================================
    THE END
========================================================================================
*/
