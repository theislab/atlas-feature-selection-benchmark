
/*
========================================================================================
    PROCESSES
========================================================================================
*/

process METHOD_ALL {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)

    output:
        tuple val(dataset), val("all"), path("all.tsv")

    script:
        """
        method-all.py \\
            --out-file "all.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "all.tsv"
        """
}

process METHOD_RANDOM {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query), val(n_features), val(seed)

    output:
        tuple val(dataset), val("random-N${n_features}-S${seed}"), path("random_N${n_features}_S${seed}.tsv")

    script:
        """
        method-random.py \\
            --n-features ${n_features} \\
            --seed ${seed} \\
            --out-file "random_N${n_features}_S${seed}.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "random_N${n_features}_S${seed}.tsv"
        """
}

process METHOD_SCANPY {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query), val(flavor), val(n_features), val(batch)

    output:
        tuple val(dataset), val("scanpy-${flavor}-N${n_features}-Batch${batch}"), path("scanpy_${flavor}_N${n_features}_Batch${batch}.tsv")

    script:
        """
        method-scanpy.py \\
            --flavor ${flavor} \\
            --n-features ${n_features} \\
            --batch ${batch} \\
            --out-file "scanpy_${flavor}_N${n_features}_Batch${batch}.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "scanpy_${flavor}_N${n_features}_Batch${batch}.tsv"
        """
}

process METHOD_TRIKU {
    conda "envs/triku.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    label "process_high"
    label "error_ignore"

    input:
        tuple val(dataset), path(reference), path(query)

    output:
        tuple val(dataset), val("triku"), path("triku.tsv")

    script:
        """
        method-triku.py \\
            --out-file "triku.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "triku.tsv"
        """
}

process METHOD_HOTSPOT {
    conda "envs/hotspot.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)

    output:
        tuple val(dataset), val("hotspot"), path("hotspot.tsv")

    script:
        """
        method-hotspot.py \\
            --n-features 500 \\
            --out-file "hotspot.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "hotspot.tsv"
        """
}

process METHOD_SEURAT {
    conda "envs/seurat.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query), val(method), val(n_features)
        path(functions)

    output:
        tuple val(dataset), val("seurat-${method}-N${n_features}"), path("seurat_${method}_N${n_features}.tsv")

    script:
        """
        method-seurat.R \\
            --method ${method} \\
            --n-features ${n_features} \\
            --out-file "seurat_${method}_N${n_features}.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "seurat_${method}_N${n_features}.tsv"
        """
}

process METHOD_SCSEGINDEX {
    conda "envs/scmerge.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)
        path(functions)

    output:
        tuple val(dataset), val("scsegindex"), path("scsegindex.tsv")

    script:
        """
        method-scSEGIndex.R \\
            --out-file "scsegindex.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "scsegindex.tsv"
        """
}

process METHOD_NBUMI {
    conda "envs/m3drop.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    label "process_medium"

    input:
        tuple val(dataset), path(reference), path(query)
        path(functions)

    output:
        tuple val(dataset), val("nbumi"), path("nbumi.tsv")

    script:
        """
        method-NBumi.R \\
            --out-file "nbumi.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "nbumi.tsv"
        """
}

process METHOD_OSCA {
    conda "envs/osca.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)
        path(functions)

    output:
        tuple val(dataset), val("osca"), path("osca.tsv")

    script:
        """
        method-osca.R \\
            --out-file "osca.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "osca.tsv"
        """
}

process METHOD_DUBSTEPR {
    conda "envs/dubstepr.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)
        path(functions)

    output:
        tuple val(dataset), val("dubstepr"), path("dubstepr.tsv")

    script:
        """
        method-DUBStepR.R \\
            --out-file "dubstepr.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "dubstepr.tsv"
        """
}

process METHOD_SCRY {
    conda "envs/scry.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    label "process_low"

    input:
        tuple val(dataset), path(reference), path(query)
        path(functions)

    output:
        tuple val(dataset), val("scry"), path("scry.tsv")

    script:
        """
        method-scry.R \\
            --out-file "scry.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "scry.tsv"
        """
}

process METHOD_SINGLECELLHAYSTACK {
    conda "envs/singleCellHaystack.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)
        path(functions)

    output:
        tuple val(dataset), val("singleCellHaystack"), path("singleCellHaystack.tsv")

    script:
        """
        method-singleCellHaystack.R \\
            --out-file "singleCellHaystack.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "singleCellHaystack.tsv"
        """
}

process METHOD_BRENNECKE {
    conda "envs/osca.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)
        path(functions)

    output:
        tuple val(dataset), val("Brennecke"), path("Brennecke.tsv")

    script:
        """
        method-Brennecke.R \\
            --out-file "Brennecke.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "Brennecke.tsv"
        """
}

process METHOD_WILCOXON {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)

    output:
        tuple val(dataset), val("wilcoxon"), path("wilcoxon.tsv")

    script:
        """
        method-wilcoxon.py \\
            --out-file "wilcoxon.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "wilcoxon.tsv"
        """
}

process METHOD_STATISTIC {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query), val(statistic), val(n_features)

    output:
        tuple val(dataset), val("statistic-${statistic}-N${n_features}"), path("statistic_${statistic}_N${n_features}.tsv")

    script:
        """
        method-statistic.py \\
            --statistic ${statistic} \\
            --n-features ${n_features} \\
            --out-file "statistic_${statistic}_N${n_features}.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "statistic_${statistic}_N${n_features}.tsv"
        """
}

process METHOD_SCPNMF {
    conda "envs/scPNMF.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    label "process_high"
    label "error_ignore"

    input:
        tuple val(dataset), path(reference), path(query)
        path(functions)

    output:
        tuple val(dataset), val("scPNMF"), path("scPNMF.tsv")

    script:
        """
        method-scPNMF.R \\
            --out-file "scPNMF.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "scPNMF.tsv"
        """
}

process METHOD_ANTICOR {
    conda "envs/anticor.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)

    output:
        tuple val(dataset), val("anticor"), path("anticor.tsv")

    script:
        """
        method-anticor.py \\
            --out-file "anticor.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "anticor.tsv"
        """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow METHODS {

    take:
        prepared_datasets_ch

    main:

        method_names = params.methods.collect{method -> method.name}

        // Function file paths
        r_io_funcs = file(params.bindir + "/functions/io.R")

        all_ch = METHOD_ALL(prepared_datasets_ch)
        triku_ch = method_names.contains("triku") ?
            METHOD_TRIKU(prepared_datasets_ch) :
            Channel.empty()
        hotspot_ch = method_names.contains("hotspot") ?
            METHOD_HOTSPOT(prepared_datasets_ch) :
            Channel.empty()
        scsegindex_ch = method_names.contains("scsegindex") ?
            METHOD_SCSEGINDEX(prepared_datasets_ch, r_io_funcs) :
            Channel.empty()
        dubstepr_ch = method_names.contains("dubstepr") ?
            METHOD_DUBSTEPR(prepared_datasets_ch, r_io_funcs) :
            Channel.empty()
        nbumi_ch = method_names.contains("nbumi") ?
            METHOD_NBUMI(prepared_datasets_ch, r_io_funcs) :
            Channel.empty()
        osca_ch = method_names.contains("osca") ?
            METHOD_OSCA(prepared_datasets_ch, r_io_funcs) :
            Channel.empty()
        scry_ch = method_names.contains("scry") ?
            METHOD_SCRY(prepared_datasets_ch, r_io_funcs) :
            Channel.empty()
        singleCellHaystack_ch = method_names.contains("singleCellHaystack") ?
            METHOD_SINGLECELLHAYSTACK(prepared_datasets_ch, r_io_funcs) :
            Channel.empty()
        brennecke_ch = method_names.contains("Brennecke") ?
            METHOD_BRENNECKE(prepared_datasets_ch, r_io_funcs) :
            Channel.empty()
        wilcoxon_ch = method_names.contains("wilcoxon") ?
            METHOD_WILCOXON(prepared_datasets_ch) :
            Channel.empty()
        scpnmf_ch = method_names.contains("scPNMF") ?
            METHOD_SCPNMF(prepared_datasets_ch, r_io_funcs) :
            Channel.empty()
        anticor_ch = method_names.contains("anticor") ?
            METHOD_ANTICOR(prepared_datasets_ch) :
            Channel.empty()

        if (method_names.contains("random")) {
            random_params_ch = Channel
                .fromList(params.methods[method_names.indexOf("random")].settings)
                .map { settings ->
                    tuple(
                        settings.n_features,
                        settings.seed
                    )
                }
            random_ch = METHOD_RANDOM(prepared_datasets_ch.combine(random_params_ch))
        } else {
            random_ch = Channel.empty()
        }

        if (method_names.contains("scanpy")) {
            scanpy_params_ch = Channel
                .fromList(params.methods[method_names.indexOf("scanpy")].settings)
                .map { settings ->
                    tuple(
                        settings.flavor,
                        settings.n_features,
                        settings.batch
                    )
                }
            scanpy_ch = METHOD_SCANPY(prepared_datasets_ch.combine(scanpy_params_ch))
        } else {
            scanpy_ch = Channel.empty()
        }

        if (method_names.contains("seurat")) {
            seurat_params_ch = Channel
                .fromList(params.methods[method_names.indexOf("seurat")].settings)
                .map { settings ->
                    tuple(
                        settings.method,
                        settings.n_features
                    )
                }
            seurat_ch = METHOD_SEURAT(
                prepared_datasets_ch.combine(seurat_params_ch),
                r_io_funcs
            )
        } else {
            seurat_ch = Channel.empty()
        }

        if (method_names.contains("statistic")) {
            statistic_params_ch = Channel
                .fromList(params.methods[method_names.indexOf("statistic")].settings)
                .map { settings ->
                    tuple(
                        settings.statistic,
                        settings.n_features
                    )
                }
            statistic_ch = METHOD_STATISTIC(prepared_datasets_ch.combine(statistic_params_ch))
        } else {
            statistic_ch = Channel.empty()
        }

        selected_features_ch = all_ch
            .mix(
                random_ch,
                scanpy_ch,
                triku_ch,
                hotspot_ch,
                nbumi_ch,
                scsegindex_ch,
                dubstepr_ch,
                osca_ch,
                seurat_ch,
                scry_ch,
                singleCellHaystack_ch,
                brennecke_ch,
                wilcoxon_ch,
                statistic_ch,
                scpnmf_ch,
                anticor_ch
            )

    emit:
        datasets_features_ch = prepared_datasets_ch.combine(selected_features_ch, by: 0)
}

/*
========================================================================================
    THE END
========================================================================================
*/
