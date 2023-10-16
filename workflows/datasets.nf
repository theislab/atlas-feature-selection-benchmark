
/*
========================================================================================
    PROCESSES
========================================================================================
*/


process DATASET_TIROSH_GENES {
    conda "envs/biomaRt.yml"

    publishDir "$params.outdir/datasets-raw/"

    output:
        path("tirosh-genes.tsv")

    script:
        """
        dataset-tirosh-genes.R --out-file "tirosh-genes.tsv"
        """

    stub:
        """
        touch tirosh-genes.tsv
        """
}

process DATASET_HUMAN_TFS {
    conda "envs/biomaRt.yml"

    publishDir "$params.outdir/datasets-raw/"

    output:
        path("human-tfs.tsv")

    script:
        """
        dataset-human-tfs.R --out-file "human-tfs.tsv"
        """

    stub:
        """
        touch human-tfs.tsv
        """
}

process DATASET_TINYSIM {
    conda "envs/splatter.yml"

    publishDir "$params.outdir/datasets-raw/"

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

    publishDir "$params.outdir/datasets-raw/"

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

    publishDir "$params.outdir/datasets-raw/"

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

    publishDir "$params.outdir/datasets-raw/"

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

process DATASET_FETALLIVER {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/datasets-raw/"

    output:
        tuple val("fetalLiver"), path("fetalLiver.h5ad")

    script:
        """
        dataset-fetalLiver.py --out-file "fetalLiver.h5ad"
        """

    stub:
        """
        touch fetalLiver.h5ad
        """
}

process DATASET_REEDBREAST {
    conda "envs/cellxgene-census.yml"

    publishDir "$params.outdir/datasets-raw/"

    output:
        tuple val("reedBreast"), path("reedBreast.h5ad")

    script:
        """
        dataset-reedBreast.py --out-file "reedBreast.h5ad"
        """

    stub:
        """
        touch reedBreast.h5ad
        """
}

process DATASET_SCEIAD {
    conda "envs/seurat.yml"

    publishDir "$params.outdir/datasets-raw/"

    label "process_medium"

    input:
        path(functions)

    output:
        tuple val("scEiaD"), path("scEiaD.h5ad")

    script:
        """
        dataset-scEiaD.R --out-file "scEiaD.h5ad"
        """

    stub:
        """
        touch scEiaD.h5ad
        """
}

process DATASET_HUMANENDODERM {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/datasets-raw/"

    output:
        tuple val("humanEndoderm"), path("humanEndoderm.h5ad")

    script:
        """
        dataset-humanEndoderm.py --out-file "humanEndoderm.h5ad"
        """

    stub:
        """
        touch humanEndoderm.h5ad
        """
}

process DATASET_HLCA {
    conda "envs/cellxgene-census.yml"

    publishDir "$params.outdir/datasets-raw/"

    output:
        tuple val("HLCA"), path("HLCA.h5ad")

    script:
        """
        dataset-HLCA.py --out-file "HLCA.h5ad"
        """

    stub:
        """
        touch HLCA.h5ad
        """
}

process DATASET_HLCAIMMUNE {
    conda "envs/cellxgene-census.yml"

    publishDir "$params.outdir/datasets-raw/"

    output:
        tuple val("HLCAImmune"), path("HLCAImmune.h5ad")

    script:
        """
        dataset-HLCAImmune.py --out-file "HLCAImmune.h5ad"
        """

    stub:
        """
        touch HLCAImmune.h5ad
        """
}

process PREPARE_DATASET {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/datasets-prepped/"

    memory { get_memory(file.size(), "2.GB", task.attempt) }

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

process DATASET_SPLAT {
    conda "envs/splatter.yml"

    publishDir "$params.outdir/datasets-raw/"

    input:
        path(functions)

    output:
        tuple val("splat"), path("splat.h5ad")

    script:
        """
        dataset-splat.R --out-file "splat.h5ad"
        """

    stub:
        """
        touch splat.h5ad
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
            DATASET_TINYSIM(file(params.bindir + "/functions/io.R")) :
            Channel.empty()
        tinySim2_ch = dataset_names.contains("tinySim2") ?
            DATASET_TINYSIM2(file(params.bindir + "/functions/io.R")) :
            Channel.empty()
        scIBPancreas_ch = dataset_names.contains("scIBPancreas") ?
            DATASET_SCIBPANCREAS() :
            Channel.empty()
		neurips_ch = dataset_names.contains("neurips") ?
            DATASET_NEURIPS() :
            Channel.empty()
        fetalLiver_ch = dataset_names.contains("fetalLiver") ?
            DATASET_FETALLIVER() :
            Channel.empty()
        reedBreast_ch = dataset_names.contains("reedBreast") ?
            DATASET_REEDBREAST() :
            Channel.empty()
        scEiaD_ch = dataset_names.contains("scEiaD") ?
            DATASET_SCEIAD(file(params.bindir + "/functions/io.R")) :
            Channel.empty()
        humanEndoderm_ch = dataset_names.contains("humanEndoderm") ?
            DATASET_HUMANENDODERM() :
            Channel.empty()
        hlca_ch = dataset_names.contains("HLCA") ?
            DATASET_HLCA() :
            Channel.empty()
        hlcaImmune_ch = dataset_names.contains("HLCAImmune") ?
            DATASET_HLCAIMMUNE() :
            Channel.empty()
        splat_ch = dataset_names.contains("splat") ?
            DATASET_SPLAT() :
            Channel.empty()

        raw_datasets_ch = tinySim_ch
            .mix(
                tinySim2_ch,
                scIBPancreas_ch,
				neurips_ch,
                fetalLiver_ch,
                reedBreast_ch,
                scEiaD_ch,
                humanEndoderm_ch,
                hlca_ch,
                hlcaImmune_ch,
                splat_ch
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

        DATASET_TIROSH_GENES()

        DATASET_HUMAN_TFS()

    emit:
        prepared_datasets_ch = PREPARE_DATASET.out
        tirosh_genes_ch = DATASET_TIROSH_GENES.out
        human_tfs_ch = DATASET_HUMAN_TFS.out
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
