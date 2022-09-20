
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

process METHOD_RANDOM_N500 {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)

    output:
        tuple val(dataset), val("random_N500"), path("random_N500.tsv")

    script:
        """
        method-random.py \\
            --n-features 500 \\
            --out-file "random_N500.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "random_N500.tsv"
        """
}

process METHOD_RANDOM_N1000 {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)

    output:
        tuple val(dataset), val("random_N1000"), path("random_N1000.tsv")

    script:
        """
        method-random.py \\
            --n-features 1000 \\
            --out-file "random_N1000.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "random_N1000.tsv"
        """
}

process METHOD_RANDOM_N2000 {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)

    output:
        tuple val(dataset), val("random_N2000"), path("random_N2000.tsv")

    script:
        """
        method-random.py \\
            --n-features 2000 \\
            --out-file "random_N2000.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "random_N2000.tsv"
        """
}

process METHOD_RANDOM_N5000 {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(reference), path(query)

    output:
        tuple val(dataset), val("random_N5000"), path("random_N5000.tsv")

    script:
        """
        method-random.py \\
            --n-features 5000 \\
            --out-file "random_N5000.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "random_N5000.tsv"
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

        all_ch          = method_names.contains("all")          ? METHOD_ALL(prepared_datasets_ch)          : Channel.empty()
        random_n500_ch  = method_names.contains("random-N500")  ? METHOD_RANDOM_N500(prepared_datasets_ch)  : Channel.empty()
        random_n1000_ch = method_names.contains("random-N1000") ? METHOD_RANDOM_N1000(prepared_datasets_ch) : Channel.empty()
        random_n2000_ch = method_names.contains("random-N2000") ? METHOD_RANDOM_N2000(prepared_datasets_ch) : Channel.empty()
        random_n5000_ch = method_names.contains("random-N5000") ? METHOD_RANDOM_N5000(prepared_datasets_ch) : Channel.empty()
        hotspot_ch      = method_names.contains("hotspot")      ? METHOD_HOTSPOT(prepared_datasets_ch)      : Channel.empty()

        selected_features_ch = all_ch
            .mix(
                random_n500_ch,
                random_n1000_ch,
                random_n2000_ch,
                random_n5000_ch,
                hotspot_ch
            )

    emit:
        datasets_features_ch = prepared_datasets_ch.combine(selected_features_ch, by: 0)
}

/*
========================================================================================
    THE END
========================================================================================
*/
