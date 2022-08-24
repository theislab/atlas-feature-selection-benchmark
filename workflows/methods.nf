
/*
========================================================================================
    PROCESSES
========================================================================================
*/

process METHOD_ALL {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${name}", mode: "copy"

    input:
        tuple val(name), path(reference), path(query)

    output:
        tuple val(name), val("all"), path("all.tsv")

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

    publishDir "$params.outdir/selected-features/${name}", mode: "copy"

    input:
        tuple val(name), path(reference), path(query)

    output:
        tuple val(name), val("random-N500"), path("random-N500.tsv")

    script:
        """
        method-random.py \\
            --n-features 500 \\
            --out-file "random-N500.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "random-N500.tsv"
        """
}

process METHOD_RANDOM_N1000 {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${name}", mode: "copy"

    input:
        tuple val(name), path(reference), path(query)

    output:
        tuple val(name), val("random-N1000"), path("random-N1000.tsv")

    script:
        """
        method-random.py \\
            --n-features 1000 \\
            --out-file "random-N1000.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "random-N1000.tsv"
        """
}

process METHOD_RANDOM_N2000 {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${name}", mode: "copy"

    input:
        tuple val(name), path(reference), path(query)

    output:
        tuple val(name), val("random-N2000"), path("random-N2000.tsv")

    script:
        """
        method-random.py \\
            --n-features 2000 \\
            --out-file "random-N2000.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "random-N2000.tsv"
        """
}

process METHOD_RANDOM_N5000 {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/selected-features/${name}", mode: "copy"

    input:
        tuple val(name), path(reference), path(query)

    output:
        tuple val(name), val("random-N5000"), path("random-N5000.tsv")

    script:
        """
        method-random.py \\
            --n-features 5000 \\
            --out-file "random-N5000.tsv" \\
            ${reference}
        """

    stub:
        """
        touch "random-N5000.tsv"
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

        METHOD_ALL(prepared_datasets_ch)
        METHOD_RANDOM_N500(prepared_datasets_ch)
        METHOD_RANDOM_N1000(prepared_datasets_ch)
        METHOD_RANDOM_N2000(prepared_datasets_ch)
        METHOD_RANDOM_N5000(prepared_datasets_ch)

        selected_features_ch = METHOD_ALL.out
            .mix(
                METHOD_RANDOM_N500.out,
                METHOD_RANDOM_N1000.out,
                METHOD_RANDOM_N2000.out,
                METHOD_RANDOM_N5000.out
            )

    emit:
        datasets_features_ch = prepared_datasets_ch.combine(selected_features_ch, by: 0)
}

/*
========================================================================================
    THE END
========================================================================================
*/
