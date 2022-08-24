#!/usr/bin/env nextflow
/*
========================================================================================
    Pig pancreas analysis
========================================================================================
    GitHub : https://github.com/theislab/atlas-feature-selection-benchmark
    Website:
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    WORKFLOWS
========================================================================================
*/

include { DATASETS } from './workflows/datasets'
include { METHODS } from './workflows/methods'

//
// WORKFLOW: Run main analysis pipeline, prints a message and ends
//
workflow WF_MAIN {
    println '\n'
    println "==== ATLAS FEATURE SELECTION BENCHMARKING PIPELINE ===="
    println '\n'
    println "Specify the workflow to run using the '-entry' option."
    println '\n'
    println 'Current workflows are:'
    println '\n'
    println '* WF_DATASETS - Download and prepare datasets'
    println '* WF_RUN - Run feature selection methods'
    println '* WF_ALL - Run all analysis steps in order'
    println '\n'
    println 'Stopping.'
    println '\n'
    println "======================================================="
    println '\n'
}

//
// WORKFLOW: Download and prepare datasets
//
workflow WF_DATASETS {
    DATASETS()
}

//
// WORKFLOW: Run all analysis steps in order
//
workflow WF_METHODS {

    prepared_datasets_ch = Channel
        .fromList(params.datasets)
        .map { dataset ->
            tuple(
                dataset.name,
                file(params.outdir + "/datasets-prepped/" + dataset.name + "-reference.h5ad"),
                file(params.outdir + "/datasets-prepped/" + dataset.name + "-query.h5ad")
            )
        }

    METHODS(prepared_datasets_ch)
}

//
// WORKFLOW: Run all analysis steps in order
//
workflow WF_ALL {
    DATASETS()
    METHODS(DATASETS.out.prepared_datasets_ch)
}

/*
========================================================================================
    DEFAULT WORKFLOW
========================================================================================
*/

//
// Implicit workflow executed when no entry is specified
//
workflow {
    WF_MAIN ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
