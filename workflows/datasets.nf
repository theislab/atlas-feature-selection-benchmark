
/*
========================================================================================
    PROCESSES
========================================================================================
*/


process DATASET_TINYSIM {
    conda "envs/splatter.yml"

    publishDir "$params.outdir/datasets-raw/", mode: "copy"

    output:
        path("tinySim.h5ad")

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
        path("tinySim2.h5ad")

    script:
        """
        dataset-tinySim2.R --out-file "tinySim2.h5ad"
        """

    stub:
        """
        touch tinySim2.h5ad
        """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow DATASETS {
    main:

        tinySim_ch  = params.datasets.contains("tinySim")  ? DATASET_TINYSIM()  : Channel.empty()
        tinySim2_ch = params.datasets.contains("tinySim2") ? DATASET_TINYSIM2() : Channel.empty()

        raw_datasets_ch = tinySim_ch
            .mix(tinySim2_ch)

        raw_datasets_ch.view()

    // emit:
    //     prep_ch = PREP_SAMPLE.out
}

/*
========================================================================================
    THE END
========================================================================================
*/
