
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
        tuple val(dataset), val(method), val("scVI"), path("${method}"), path(query)

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

process INTEGRATE_SCANVI {
    conda "envs/scvi-tools.yml"

    publishDir "$params.outdir/integration-models/${dataset}/scANVI", mode: "copy"

    input:
        tuple val(dataset), val(method), val(integration), path(model), path(query)

    output:
        tuple val(dataset), val(method), val("scANVI"), path("${method}"), path(query)

    script:
        """
        integrate-scanvi.py \\
            --out-dir "${method}" \\
            ${model}
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
        INTEGRATE_SCANVI(INTEGRATE_SCVI.out)

    // emit:
    //     datasets_features_ch = prepared_datasets_ch.combine(selected_features_ch, by: 0)
}

/*
========================================================================================
    THE END
========================================================================================
*/
