
/*
========================================================================================
    PROCESSES
========================================================================================
*/


process DATASET_TINYSIM {
    conda "envs/splatter.yml"

    publishDir "$params.outdir/datasets-raw/", mode: "copy"

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

process PREPARE_DATASET {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/datasets-prepped/"

    input:
        tuple val(name), val(batch_col), val(label_col), val(query_batches), path(file)

    output:
        tuple val(name), path("${name}-reference.h5ad"), path("${name}-query.h5ad")

    script:
        """
        prepare-dataset.py \\
            --name "${name}" \\
            --batch-col "${batch_col}" \\
            --label-col "${label_col}" \\
            --query-batches "${query_batches}" \\
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

        tinySim_ch  = dataset_names.contains("tinySim")  ? DATASET_TINYSIM()  : Channel.empty()
        tinySim2_ch = dataset_names.contains("tinySim2") ? DATASET_TINYSIM2() : Channel.empty()

        raw_datasets_ch = tinySim_ch
            .mix(tinySim2_ch)

        datasets_ch = Channel
            .fromList(params.datasets)
            .map { dataset ->
                tuple(
                    dataset.name,
                    dataset.batch_col,
                    dataset.label_col,
                    dataset.query_batches
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