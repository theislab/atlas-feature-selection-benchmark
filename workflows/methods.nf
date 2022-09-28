
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

    label "process_low"

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

        all_ch            = method_names.contains("all")            ? METHOD_ALL(prepared_datasets_ch)            : Channel.empty()
        triku_ch          = method_names.contains("triku")          ? METHOD_TRIKU(prepared_datasets_ch)          : Channel.empty()
        hotspot_ch        = method_names.contains("hotspot")        ? METHOD_HOTSPOT(prepared_datasets_ch)        : Channel.empty()
        scsegindex_ch     = method_names.contains("scsegindex")     ? METHOD_SCSEGINDEX(prepared_datasets_ch, file(params.bindir + "/_functions.R"))     : Channel.empty()
        nbumi_ch          = method_names.contains("nbumi")          ? METHOD_NBUMI(prepared_datasets_ch, file(params.bindir + "/_functions.R"))          : Channel.empty()

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
            random_ch = Channel.empty()
        }

        selected_features_ch = all_ch
            .mix(
                random_ch,
                scanpy_ch,
                triku_ch,
                hotspot_ch,
                nbumi_ch,
                scsegindex_ch,
            )

    emit:
        datasets_features_ch = prepared_datasets_ch.combine(selected_features_ch, by: 0)
}

/*
========================================================================================
    THE END
========================================================================================
*/
