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
    println '* WF_ALL - Run all analysis steps in order'
    println '\n'
    println 'Stopping.'
    println '\n'
    println "======================================================="
    println '\n'
}

//
// WORKFLOW: Run all analysis steps in order
//
workflow WF_ALL {

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
