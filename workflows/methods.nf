
/*
========================================================================================
    PROCESSES
========================================================================================
*/


/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow METHODS {

    take:
        prepared_datasets_ch

    main:

        prepared_datasets_ch.view()

    // emit:
    //     prepared_datasets_ch = PREPARE_DATASET.out
}

/*
========================================================================================
    THE END
========================================================================================
*/