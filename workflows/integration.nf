
/*
========================================================================================
    PROCESSES
========================================================================================
*/

process INTEGRATE_SCVI {
    conda "envs/scvi-tools.yml"

    publishDir "$params.outdir/integration-models/${dataset}/${method}",
        mode: "copy",
        pattern: "scVI-reference"

    input:
        tuple val(dataset), path(reference), path(query), val(method), path(features)

    output:
        tuple val(dataset), val(method), val("scVI"), path("scVI-reference"), path(query)

    script:
        """
        integrate-scvi.py \\
            --features "${features}" \\
            --out-dir scVI-reference \\
            ${reference}
        """

    stub:
        """
        mkdir scVI-reference
        touch scVI-reference/adata.h5ad
        touch scVI-reference/model.pt
        """
}

process INTEGRATE_SCANVI {
    conda "envs/scvi-tools.yml"

    publishDir "$params.outdir/integration-models/${dataset}/${method}",
        mode: "copy",
        pattern: "scANVI-reference"

    input:
        tuple val(dataset), val(method), val(integration), path(scVI), path(query)

    output:
        tuple val(dataset), val(method), val("scANVI"), path("scANVI-reference"), path(query)

    script:
        """
        integrate-scanvi.py \\
            --out-dir "scANVI-reference" \\
            ${scVI}
        """

    stub:
        """
        mkdir scANVI-reference
        touch scANVI-reference/adata.h5ad
        touch scANVI-reference/model.pt
        """
}

process MAP_SCVI {
    conda "envs/scvi-tools.yml"

    publishDir "$params.outdir/integration-models/${dataset}/${method}",
        mode: "copy",
        pattern: "scVI-mapped"

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(query)

    output:
        tuple val(dataset), val(method), val(integration), path(reference), path("scVI-mapped")

    script:
        """
        map-scvi.py \\
            --reference "${reference}" \\
            --out-dir scVI-mapped \\
            ${query}
        """

    stub:
        """
        mkdir scVI-mapped
        touch scVI-mapped/adata.h5ad
        touch scVI-mapped/model.pt
        """
}

process MAP_SCANVI {
    conda "envs/scvi-tools.yml"

    publishDir "$params.outdir/integration-models/${dataset}/${method}",
        mode: "copy",
        pattern: "scANVI-mapped"

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(query)

    output:
        tuple val(dataset), val(method), val(integration), path(reference), path("scANVI-mapped")

    script:
        """
        map-scanvi.py \\
            --reference "${reference}" \\
            --out-dir scANVI-mapped \\
            ${query}
        """

    stub:
        """
        mkdir scANVI-mapped
        touch scANVI-mapped/adata.h5ad
        touch scANVI-mapped/model.pt
        """
}

process PREDICT_LABELS {
    conda "envs/lightgbm.yml"

    publishDir "$params.outdir/integration-models/${dataset}/${method}",
        mode: "copy"

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(query)

    output:
        tuple val(dataset), val(method), val(integration), val("${reference}/adata.h5ad"), val("${query}/adata.h5ad"), path("${integration}-labels.tsv")

    script:
        """
        predict-labels.py \\
            --reference "${reference}/adata.h5ad" \\
            --out-file "${integration}-labels.tsv" \\
            ${query}/adata.h5ad
        """

    stub:
        """
        mkdir "${integration}-labels.tsv"
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

        MAP_SCVI(INTEGRATE_SCVI.out)
        MAP_SCANVI(INTEGRATE_SCANVI.out)

        mapped_ch = MAP_SCVI.out.mix(MAP_SCANVI.out)

        PREDICT_LABELS(mapped_ch)

    // emit:
    //     datasets_features_ch = prepared_datasets_ch.combine(selected_features_ch, by: 0)
}

/*
========================================================================================
    THE END
========================================================================================
*/
