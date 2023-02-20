
/*
========================================================================================
    PROCESSES
========================================================================================
*/


process DATASET_TINYSIM {
    conda "envs/splatter.yml"

    publishDir "$params.outdir/datasets-raw/", mode: "copy"

    input:
        path(functions)

    output:
        tuple val("tinySim"), path("tinySim.h5ad")

    script:
        """
        dataset-tinySim.R --out-file "tinySim.h5ad"
        """

    stub:
        """
        touch tinySim.h5ad
        """
}

process DATASET_TINYSIM2 {
    conda "envs/splatter.yml"

    publishDir "$params.outdir/datasets-raw/", mode: "copy"

    input:
        path(functions)

    output:
        tuple val("tinySim2"), path("tinySim2.h5ad")

    script:
        """
        dataset-tinySim2.R --out-file "tinySim2.h5ad"
        """

    stub:
        """
        touch tinySim2.h5ad
        """
}

process DATASET_SCIBPANCREAS {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/datasets-raw/", mode: "copy"

    output:
        tuple val("scIBPancreas"), path("scIBPancreas.h5ad")

    script:
        """
        dataset-scIBPancreas.py --out-file "scIBPancreas.h5ad"
        """

    stub:
        """
        touch scIBPancreas.h5ad
        """
}

process DATASET_NEURIPS {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/datasets-raw/", mode: "copy"

    input:
        path(functions)

    output:
        tuple val("neurips"), path("neurips.h5ad")

    script:
        """
        dataset-neurips.py --out-file "neurips.h5ad"
        """

    stub:
        """
        touch neurips.h5ad
        """
}

process PREPARE_DATASET {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/datasets-prepped/"

    input:
        tuple val(name), val(batch_col), val(query_batches), val(label_col), val(unseen_labels), val(species), path(file)

    output:
        tuple val(name), path("${name}-reference.h5ad"), path("${name}-query.h5ad")

    script:
        """
        prepare-dataset.py \\
            --name "${name}" \\
            --batch-col "${batch_col}" \\
            --query-batches "${query_batches}" \\
            --label-col "${label_col}" \\
            --unseen-labels "${unseen_labels}" \\
            --species "${species}" \\
            --reference-out "${name}-reference.h5ad" \\
            --query-out "${name}-query.h5ad" \\
            ${file}
        """

    stub:
        """
        touch "${name}-reference.h5ad"
        touch "${name}-query.h5ad"
        """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow DATASETS {
    main:

        dataset_names = params.datasets.collect{dataset -> dataset.name}

        tinySim_ch  = dataset_names.contains("tinySim")  ?
            DATASET_TINYSIM(file(params.bindir + "/_functions.R"))  :
            Channel.empty()
        tinySim2_ch = dataset_names.contains("tinySim2") ?
            DATASET_TINYSIM2(file(params.bindir + "/_functions.R")) :
            Channel.empty()
        scIBPancreas_ch = dataset_names.contains("scIBPancreas") ?
            DATASET_SCIBPANCREAS() :
            Channel.empty()
		neurips_ch = dataset_names.contains("neurips") ?
            DATASET_NEURIPS(file(params.bindir + "/_functions.R")) :
            Channel.empty()

        raw_datasets_ch = tinySim_ch
            .mix(
                tinySim2_ch,
                scIBPancreas_ch,
				neurips_ch
            )

        datasets_ch = Channel
            .fromList(params.datasets)
            .map { dataset ->
                tuple(
                    dataset.name,
                    dataset.batch_col,
                    dataset.query_batches,
                    dataset.label_col,
                    dataset.unseen_labels,
                    dataset.species
                )
            }
            .join(raw_datasets_ch)

        PREPARE_DATASET(datasets_ch)

    emit:
        prepared_datasets_ch = PREPARE_DATASET.out
}

/*
========================================================================================
    THE END
========================================================================================
*/
