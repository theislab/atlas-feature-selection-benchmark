
/*
========================================================================================
    PROCESSES
========================================================================================
*/

process INTEGRATE_SCVI {
    conda "envs/scvi-tools.yml"

    publishDir "$params.outdir/integration-models/${dataset}/${method}",
        pattern: "scVI-reference",
        saveAs: { pathname -> pathname + "-${seed}" }

    label "process_medium"

    memory { get_memory(reference.size(), "6.GB", task.attempt) }

    input:
        tuple val(dataset), path(reference), path(query), val(method), path(features), val(seed)
        path(functions)

    output:
        tuple val(dataset), val(method), val("scVI"), val(seed), path(reference), path("scVI-reference"), path(query)

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
        pattern: "scANVI-reference",
        saveAs: { pathname -> pathname + "-${seed}" }

    label "process_medium"

    memory { get_memory(reference.size(), "12.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(scVI), path(query)
        path(functions)

    output:
        tuple val(dataset), val(method), val("scANVI"), val(seed), path(reference), path("scANVI-reference"), path(query)

    script:
        """
        integrate-scanvi.py \\
            --scvi "${scVI}" \\
            --out-dir "scANVI-reference" \\
            --seed "${seed}" \\
            ${reference}
        """

    stub:
        """
        mkdir scANVI-reference
        touch scANVI-reference/adata.h5ad
        touch scANVI-reference/model.pt
        """
}

process INTEGRATE_SYMPHONY {
    conda "envs/symphonypy.yml"

    publishDir "$params.outdir/integration-models/${dataset}/${method}",
        pattern: "symphony-reference",
        saveAs: { pathname -> pathname + "-${seed}" }

    label "process_medium"

    memory { get_memory(reference.size(), "16.GB", task.attempt, "16.GB") }

    input:
        tuple val(dataset), path(reference), path(query), val(method), path(features), val(seed)
        path(functions)

    output:
        tuple val(dataset), val(method), val("Symphony"), val(seed), path(reference), path("symphony-reference"), path(query)

    script:
        """
        integrate-symphony.py \\
            --features "${features}" \\
            --out-dir symphony-reference \\
            --seed "${seed}" \\
            ${reference}
        """

    stub:
        """
        mkdir symphony-reference
        touch symphony-reference/adata.h5ad
        """
}

process MAP_SCVI {
    conda "envs/scvi-tools.yml"

    publishDir "$params.outdir/integration-models/${dataset}/${method}",
        pattern: "scVI-mapped",
        saveAs: { pathname -> pathname + "-${seed}" }

    label "process_high"

    memory { get_memory(reference.size(), "16.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(reference_model), path(query)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(reference_model), path(query), path("scVI-mapped")

    script:
        """
        map-scvi.py \\
            --reference "${reference}" \\
            --reference-model "${reference_model}" \\
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
        pattern: "scANVI-mapped",
        saveAs: { pathname -> pathname + "-${seed}" }

    label "process_high"

    memory { get_memory(reference.size(), "16.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(reference_model), path(query)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(reference_model), path(query), path("scANVI-mapped")

    script:
        """
        map-scanvi.py \\
            --reference "${reference}" \\
            --reference-model "${reference_model}" \\
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

process MAP_SYMPHONY {
    conda "envs/symphonypy.yml"

    publishDir "$params.outdir/integration-models/${dataset}/${method}",
        pattern: "symphony-mapped",
        saveAs: { pathname -> pathname + "-${seed}" }

    label "process_medium"

    memory { get_memory(reference.size(), "16.GB", task.attempt) }

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(reference_model), path(query)
        path(functions)

    output:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(reference_model), path(query), path("symphony-mapped")

    script:
        """
        map-symphony.py \\
            --reference "${reference}" \\
            --reference-model "${reference_model}" \\
            --out-dir symphony-mapped \\
            ${query}
        """

    stub:
        """
        mkdir symphony-mapped
        touch symphony-mapped/adata.h5ad
        """
}

process OPTIMISE_CLASSIFIER {
    conda "envs/lightgbm.yml"

    publishDir "$params.outdir/predicted-labels/${dataset}",
        mode: "copy"

    label "process_high"

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(reference_model), path(query)

    output:
        tuple val(dataset), path("parameters.tsv")

    script:
        """
        optimise-classifier.py \\
            --out-file "parameters.tsv" \\
            ${reference_model}/adata.h5ad
        """

    stub:
        """
        mkdir "parameters.tsv"
        """
}

process PREDICT_LABELS {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/predicted-labels/${dataset}/${method}",
        mode: "copy",
        pattern: "*.tsv"

    input:
        tuple val(dataset), val(method), val(integration), val(seed), path(reference), path(reference_model), path(query), path(query_model)

    output:
        tuple val(dataset), val(method), val("${integration}-${seed}"), path("${query_model}/adata.h5ad"), path("${method}-${integration}-${seed}.tsv")

    script:
        """
        predict-labels.py \\
            --reference "${reference_model}/adata.h5ad" \\
            --out-file "${method}-${integration}-${seed}.tsv" \\
            ${query_model}/adata.h5ad
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

        integration_ch = datasets_features_ch
            .combine(Channel.fromList(params.integration.seeds))

        INTEGRATE_SCVI(integration_ch, file(params.bindir + "/functions/integration.py"))
        INTEGRATE_SCANVI(INTEGRATE_SCVI.out, file(params.bindir + "/functions/integration.py"))
        INTEGRATE_SYMPHONY(integration_ch, file(params.bindir + "/functions/integration.py"))

        MAP_SCVI(INTEGRATE_SCVI.out, file(params.bindir + "/functions/integration.py"))
        MAP_SCANVI(INTEGRATE_SCANVI.out, file(params.bindir + "/functions/integration.py"))
        MAP_SYMPHONY(INTEGRATE_SYMPHONY.out, file(params.bindir + "/functions/integration.py"))

        // Use one scVI integration with all features for each dataset to
        // optimise classifier hyperparameters
        // scvi_all_ch = INTEGRATE_SCVI.out
        //     .filter { it[1] == "all" }
        //     .filter { it[3] == params.integration.seeds.min()}

        // OPTIMISE_CLASSIFIER(scvi_all_ch)

        mapped_ch = MAP_SCVI.out
            .mix(MAP_SCANVI.out)
            .mix(MAP_SYMPHONY.out)
            // .combine(OPTIMISE_CLASSIFIER.out, by: 0)

        PREDICT_LABELS(mapped_ch)

    emit:
        reference_ch = INTEGRATE_SCVI.out
            .mix(INTEGRATE_SCANVI.out)
            .mix(INTEGRATE_SYMPHONY.out)
            .map { it ->
                tuple(
                    it[0],                       // Dataset name
                    it[1],                       // Method name
                    it[2] + "-" + it[3],         // Integration name
                    file(it[5] + "/adata.h5ad"), // Path to reference H5AD
                    it[4]                        // Path to reference expression H5AD
                )
            }

        query_ch = MAP_SCVI.out
            .mix(MAP_SCANVI.out)
            .mix(MAP_SYMPHONY.out)
            .map { it ->
                tuple(
                    it[0],                       // Dataset name
                    it[1],                       // Method name
                    it[2] + "-" + it[3],         // Integration name
                    file(it[5] + "/adata.h5ad"), // Path to reference H5AD
                    file(it[7] + "/adata.h5ad"), // Path to query H5AD
                    it[7],                       // Path to query model directory
                    it[4],                       // Path to reference expression H5AD
                    it[6]                        // Path to query expression H5AD
                )
            }

        labels_ch = PREDICT_LABELS.out
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def get_memory(file_size, mem_per_gb = "1.GB", attempt = 1, overhead = "8.GB") {
    file_mem = new nextflow.util.MemoryUnit(file_size)
    mem_per_gb = new nextflow.util.MemoryUnit(mem_per_gb)
    overhead = new nextflow.util.MemoryUnit(overhead)
    max_mem = new nextflow.util.MemoryUnit(params.max_memory)

    file_gb = file_mem.toGiga()
    file_gb = file_gb > 0 ? file_gb : 1
    mem_use = (mem_per_gb * file_gb * attempt) + overhead

    if (mem_use > max_mem) {
        mem_use = max_mem
    }

    return mem_use
}

/*
========================================================================================
    THE END
========================================================================================
*/
