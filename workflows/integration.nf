
/*
========================================================================================
    PROCESSES
========================================================================================
*/

process INTEGRATE_SCVI {
    conda "envs/scvi-tools.yml"

    publishDir "$params.outdir/integration-models/${dataset}/scVI", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query), val(method), path(features)

    output:
        tuple val(dataset), val(method), val("scVI"), path("${method}")

    script:
        """
        integrate-scvi.py \\
            --features "${features}" \\
            --out-dir "${method}" \\
            ${reference}
        """

    stub:
        """
        mkdir ${method}
        touch ${method}/adata.h5ad
        touch ${method}/model.pt
        """
}
/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow INTEGRATION {

    take:
        datasets_features_ch

    main:

        INTEGRATE_SCVI(datasets_features_ch)

    // emit:
    //     datasets_features_ch = prepared_datasets_ch.combine(selected_features_ch, by: 0)
}

/*
========================================================================================
    THE END
========================================================================================
*/
