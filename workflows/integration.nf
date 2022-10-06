
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

    label "process_medium"

    input:
        tuple val(dataset), path(reference), path(query), val(method), path(features)
        path(functions)

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
        path(functions)

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
        path(functions)

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
        path(functions)

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
        mode: "copy",
        pattern: "*-labels.tsv"

    label "process_medium"

    input:
        tuple val(dataset), val(method), val(integration), path(reference), path(query)

    output:
        tuple val(dataset), val(method), val(integration), path("${query}/adata.h5ad"), path("${integration}-labels.tsv")

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

        INTEGRATE_SCVI(datasets_features_ch, file(params.bindir + "/_functions.py"))
        INTEGRATE_SCANVI(INTEGRATE_SCVI.out, file(params.bindir + "/_functions.py"))

        MAP_SCVI(INTEGRATE_SCVI.out, file(params.bindir + "/_functions.py"))
        MAP_SCANVI(INTEGRATE_SCANVI.out, file(params.bindir + "/_functions.py"))

        mapped_ch = MAP_SCVI.out.mix(MAP_SCANVI.out)

        PREDICT_LABELS(mapped_ch)

    emit:
        reference_ch = INTEGRATE_SCVI.out
            .mix(INTEGRATE_SCANVI.out)
            .map { it ->
                tuple(
                    it[0],                       // Dataset name
                    it[1],                       // Method name
                    it[2],                       // Integration name
                    file(it[3] + "/adata.h5ad"), // Path to reference H5AD
                )
            }
        query_ch = PREDICT_LABELS.out
}

/*
========================================================================================
    THE END
========================================================================================
*/
