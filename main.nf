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
include { INTEGRATION } from './workflows/integration'
include { METRICS } from './workflows/metrics'
include { REPORTS } from './workflows/reports'

prepared_datasets_ch = Channel
    .fromList(params.datasets)
    .map { dataset ->
        tuple(
            dataset.name,
            file(params.outdir + "/datasets-prepped/" + dataset.name + "-reference.h5ad"),
            file(params.outdir + "/datasets-prepped/" + dataset.name + "-query.h5ad")
        )
    }

features_ch = Channel
    .fromList(params.methods)
    .map { method ->
        tuple(
            method.name
        )
    }

datasets_features_ch = prepared_datasets_ch
    .combine(features_ch)
    .map { input ->
        tuple(
            input[0], // Dataset name
            input[1], // Path to reference file
            input[2], // Path to query file
            input[3], // Method name
            file(params.outdir + "/selected-features/" + input[0] + "/" + input[3] + ".tsv")
        )
    }

integrations_ch = Channel.fromList(["scVI", "scANVI"])

reference_ch = datasets_features_ch
    .combine(integrations_ch)
    .map { input ->
        tuple(
            input[0], // Dataset name
            input[3], // Method name
            input[5], // Integration name
            file(params.outdir + "/" + input[0] + "/" + input[3] + "/" + input[5] + "-reference/adata.h5ad")
        )
    }

query_ch = reference_ch
    .map { input ->
        tuple(
            input[0], // Dataset name
            input[1], // Method name
            input[2], // Integration name
            file(params.outdir + "/" + input[0] + "/" + input[1] + "/" + input[2] + "-mapped/adata.h5ad"),
            file(params.outdir + "/" + input[0] + "/" + input[1] + "/" + input[2] + "-labels.tsv")
        )
    }

combined_metrics_ch = Channel.fromList([file(params.outdir + "/metrics/all-metrics.tsv")])

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
    println '* WF_METHODS - Run feature selection methods'
    println '* WF_INTEGRATION - Run integration steps'
    println '* WF_METRICS - Run evaluation metrics'
    println '* WF_REPORTS - Run reports'
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
// WORKFLOW: Run feature selection methods
//
workflow WF_METHODS {
    METHODS(prepared_datasets_ch)
}

//
// WORKFLOW: Run integration
//
workflow WF_INTEGRATION {
    INTEGRATION(datasets_features_ch)
}

//
// WORKFLOW: Run metrics
//
workflow WF_METRICS {
    METRICS(reference_ch, query_ch)
}

//
// WORKFLOW: Run reports
//
workflow WF_REPORTS {
    REPORTS(combined_metrics_ch)
}

//
// WORKFLOW: Run all analysis steps in order
//
workflow WF_ALL {
    DATASETS()
    METHODS(DATASETS.out.prepared_datasets_ch)
    INTEGRATION(METHODS.out.datasets_features_ch)
    METRICS(INTEGRATION.out.reference_ch, INTEGRATION.out.query_ch)
    REPORTS(METRICS.out)
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
