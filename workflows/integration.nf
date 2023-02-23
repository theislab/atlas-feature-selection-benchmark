
/*
========================================================================================
    PROCESSES
========================================================================================
*/

process INTEGRATE_SCVI {
    conda "envs/scvi-tools.yml"

    publishDir "$params.outdir/integration-models/${dataset}/${method}",
        mode: "copy",
        pattern: "scVI-reference",
        saveAs: { pathname -> pathname + "-${seed}" }

    label "process_medium"

    input:
        tuple val(dataset), path(reference), path(query), val(method), path(features), val(seed)
        path(functions)

    output:
        tuple val(dataset), val(method), val("scVI"), val(seed), path("scVI-reference"), path(query)

    script:
        """
        integrate-scvi.py \\
            --features "${features}" \\
            --out-dir scVI-reference \\
            --seed "${seed}" \\
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
        pattern: "scANVI-reference",
        saveAs: { pathname -> pathname + "-${seed}" }

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(scVI), path(query)
        path(functions)

    output:
        tuple val(dataset), val(method), val("scANVI"), val(seed), path("scANVI-reference"), path(query)

    script:
        """
        integrate-scanvi.py \\
            --out-dir "scANVI-reference" \\
            --seed "${seed}" \\
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
        pattern: "scVI-mapped",
        saveAs: { pathname -> pathname + "-${seed}" }

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(query)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path("scVI-mapped")

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
        pattern: "scANVI-mapped",
        saveAs: { pathname -> pathname + "-${seed}" }

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(query)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path("scANVI-mapped")

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

process OPTIMISE_CLASSIFIER {
    conda "envs/lightgbm.yml"

    publishDir "$params.outdir/predicted-labels/${dataset}",
        mode: "copy"

    label "process_medium"

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(query)

    output:
        tuple val(dataset), path("parameters.tsv")

    script:
        """
        optimise-classifier.py \\
            --out-file "parameters.tsv" \\
            ${reference}/adata.h5ad
        """

    stub:
        """
        mkdir "parameters.tsv"
        """
}

process PREDICT_LABELS {
    conda "envs/lightgbm.yml"

    publishDir "$params.outdir/predicted-labels/${dataset}/${method}",
        mode: "copy",
        pattern: "*.tsv"

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(query), path(parameters)

    output:
        tuple val(dataset), val(method), val("${integration}-${seed}"), path("${query}/adata.h5ad"), path("${method}-${integration}-${seed}.tsv")

    script:
        """
        predict-labels.py \\
            --reference "${reference}/adata.h5ad" \\
            --params "${parameters}" \\
            --out-file "${method}-${integration}-${seed}.tsv" \\
            ${query}/adata.h5ad
        """

    stub:
        """
        mkdir "${method}-${integration}-${seed}.tsv"
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

        scvi_ch = datasets_features_ch
            .combine(Channel.fromList(params.integration.seeds))

        INTEGRATE_SCVI(scvi_ch, file(params.bindir + "/functions/integration.py"))
        INTEGRATE_SCANVI(INTEGRATE_SCVI.out, file(params.bindir + "/functions/integration.py"))

        MAP_SCVI(INTEGRATE_SCVI.out, file(params.bindir + "/functions/integration.py"))
        MAP_SCANVI(INTEGRATE_SCANVI.out, file(params.bindir + "/functions/integration.py"))

        // Use one scVI integration with all features for each dataset to
        // optimise classifier hyperparameters
        scvi_all_ch = INTEGRATE_SCVI.out
            .filter { it[1] == "all" }
            .filter { it[3] == params.integration.seeds.min()}

        OPTIMISE_CLASSIFIER(scvi_all_ch)

        mapped_ch = MAP_SCVI.out.mix(MAP_SCANVI.out)
            .combine(OPTIMISE_CLASSIFIER.out, by: 0)

        PREDICT_LABELS(mapped_ch)

    emit:
        reference_ch = INTEGRATE_SCVI.out
            .mix(INTEGRATE_SCANVI.out)
            .map { it ->
                tuple(
                    it[0],                       // Dataset name
                    it[1],                       // Method name
                    it[2] + "-" + it[3],         // Integration name
                    file(it[4] + "/adata.h5ad"), // Path to reference H5AD
                )
            }

        query_ch = MAP_SCVI.out
            .mix(MAP_SCANVI.out)
            .map { it ->
                tuple(
                    it[0],                       // Dataset name
                    it[1],                       // Method name
                    it[2] + "-" + it[3],         // Integration name
                    file(it[4] + "/adata.h5ad"), // Path to reference H5AD
                    file(it[5] + "/adata.h5ad"), // Path to query H5AD
                )
            }

        labels_ch = PREDICT_LABELS.out
}

/*
========================================================================================
    THE END
========================================================================================
*/
