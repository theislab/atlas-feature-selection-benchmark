#!/usr/bin/env Rscript

"
Create a simulated dataset using the splat simulation

Usage:
    dataset-splat.R --out-file=<path> [options]
Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

# Load libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

# Source functions
suppressMessages({
    source("io.R")
})

#' Install splatter from GitHub using remotes
install_splatter <- function() {
    if (!requireNamespace("splatter", quietly = TRUE)) {
        message("Installing splatter...")
        remotes::install_github(
            "Oshlack/splatter@v1.25.1",
            dependencies = FALSE
        )
    } else {
        message("splatter already installed")
    }
}

#' Simulate a tiny dataset
#'
#' @returns SingleCellExperiment with simulated data
simulate_dataset <- function() {

    set.seed(1)

    final_cells <- 100000

    message("Defining batches...")
    batch_params <- data.frame(
        Batch = c(
            "Batch1", "Batch2", "Batch3", "Batch4", "Batch5", "Batch6"
        ),
        RelProp      = c(0.90, 1.00, 0.22, 0.20, 2.60, 2.80),
        FacLoc       = c(0.08, 0.05, 0.12, 0.18, 0.15, 0.25),
        FacScale     = c(0.06, 0.08, 0.12, 0.10, 0.16, 0.18),
        RelLibSize   = c(0.22, 0.18, 0.90, 1.05, 0.05, 0.04)
    )
    batch_params$Prop <- batch_params$RelProp / sum(batch_params$RelProp)
    batch_params$Cells <- round(batch_params$Prop * final_cells)
    rownames(batch_params) <- batch_params$Batch

    message("Defining groups...")
    group_params = data.frame(
        Group = c(
            "Path1", "Path2", "Path3", "Path4", "Path5",
            "Path6", "Path7", "Path8", "Path9", "Path10"
        ),
        Name = c(
            "Progenitor", "Intermediate", "Trajectory1", "Trajectory2",
            "Discrete1", "Common", "Discrete2", "Perturbed", "Discrete3",
            "Rare"
        ),
        #              1     2     3     4     5     6     7     8     9    10
        DEProb  = c(0.00, 0.05, 0.10, 0.15, 0.05, 0.20, 0.10, 0.05, 0.08, 0.02),
        DELoc   = c(0.00, 0.40, 1.00, 0.80, 0.60, 0.20, 1.00, 0.60, 1.20, 2.50),
        DEScale = c(0.00, 0.80, 0.20, 0.40, 0.20, 0.60, 0.40, 0.20, 0.20, 0.20),
        From    = c(   0,    1,    2,    0,    4,    0,    0,    7,    6,    0),
        Steps   = c(   1,   20,   20,   60,    2,    2,    2,    2,    2,    2),
        Skew    = c( 0.5,  0.5,  0.3,  0.7,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0),
        Dropout = c(  -1, -0.1, -0.5, -0.2,   -2, -0.6,   -3, -0.8, -0.4, -0.9),
        Noise   = c(0.04, 0.03, 0.04, 0.05, 0.08, 0.15, 0.07, 0.10, 0.11, 0.06)
    )

    message("Defining batch group proportions...")
    batch_groups <- data.frame(
        Batch1 = c(0.02, 0.00, 0.15, 0.13, 0.08, 0.35, 0.14, 0.00, 0.13, 0.00),
        Batch2 = c(0.01, 0.00, 0.12, 0.16, 0.06, 0.38, 0.12, 0.00, 0.15, 0.00),
        Batch3 = c(0.04, 0.00, 0.15, 0.10, 0.17, 0.24, 0.14, 0.00, 0.16, 0.00),
        Batch4 = c(0.05, 0.00, 0.18, 0.10, 0.12, 0.26, 0.16, 0.00, 0.13, 0.00),
        Batch5 = c(0.04, 0.05, 0.11, 0.16, 0.06, 0.21, 0.05, 0.15, 0.14, 0.03),
        Batch6 = c(0.03, 0.07, 0.09, 0.15, 0.09, 0.20, 0.04, 0.13, 0.17, 0.03)
    )
    group_maxs <- apply(as.matrix(batch_groups), 1, max)
    initial_group_probs <- group_maxs / sum(group_maxs)

    n_groups <- nrow(group_params)

    initial_batch_cells <- sapply(seq_len(nrow(batch_params)), function(batch) {
        batch_props <- batch_groups[[batch]]
        less_than_initial <- batch_props < initial_group_probs
        batch_props[less_than_initial] <- initial_group_probs[less_than_initial]
        batch_factor <- sum(batch_props)
        ceiling(batch_params$Cells[batch] * batch_factor * 1.5)
    })

    message("Setting initial parameters...")
    initial_params <- splatter::newSplatParams(
        lib.loc        = 12,
        lib.scale      = 0.3,
        # Initial batches of equal size
        batchCells     = initial_batch_cells,
        batch.facLoc   = batch_params$FacLoc,
        batch.facScale = batch_params$FacScale,
        # Initial groups with estimated max probabilities
        group.prob     = initial_group_probs,
        path.from      = group_params$From,
        path.nSteps    = group_params$Steps,
        path.skew      = group_params$Skew,
        # Differential expression by group
        de.prob        = group_params$DEProb,
        de.facLoc      = group_params$DELoc,
        de.facScale    = group_params$DEScale,
        # Dropout by group
        dropout.mid    = rep(1, n_groups),
        dropout.shape  = group_params$Dropout,
        dropout.type   = "group",
        # Seed
        seed           = 1
    )

    message(
        "Simulating initial dataset with ",
        sum(initial_batch_cells), " cells..."
    )
    message("Getting parameters...")
    params <- splatter:::expandParams(initial_params)

    # Set random seed
    seed <- splatter::getParam(params, "seed")
    withr::with_seed(seed, {

        # Get the parameters we are going to use
        nCells <- splatter::getParam(params, "nCells")
        nGenes <- splatter::getParam(params, "nGenes")
        nBatches <- splatter::getParam(params, "nBatches")
        batch.cells <- splatter::getParam(params, "batchCells")
        nGroups <- splatter::getParam(params, "nGroups")
        group.prob <- splatter::getParam(params, "group.prob")

        # Set up name vectors
        message("Creating simulation object...")
        cell.names <- paste0("Cell", seq_len(nCells))
        gene.names <- paste0("Gene", seq_len(nGenes))
        batch.names <- paste0("Batch", seq_len(nBatches))
        group.names <- paste0("Path", seq_len(nGroups))

        # Create SingleCellExperiment to store simulation
        cells <-  data.frame(Cell = cell.names)
        rownames(cells) <- cell.names
        features <- data.frame(Gene = gene.names)
        rownames(features) <- gene.names
        sim <- SingleCellExperiment(rowData = features, colData = cells)

        # Make batches vector which is the index of param$batchCells repeated
        # params$batchCells[index] times
        batches <- lapply(seq_len(nBatches), function(i, b) {rep(i, b[i])},
                          b = batch.cells)
        batches <- unlist(batches)
        colData(sim)$Batch <- batch.names[batches]

        groups <- sample(seq_len(nGroups), nCells, prob = group.prob,
                         replace = TRUE)
        colData(sim)$Group <- factor(group.names[groups], levels = group.names)

        gc()
        message("Simulating library sizes...")
        sim <- splatter:::splatSimLibSizes(sim, params)
        gc()
        message("Simulating gene means...")
        sim <- splatter:::splatSimGeneMeans(sim, params)
        rowData(sim)$BaseGeneMean <- NULL
        rowData(sim)$OutlierFactor <- NULL
        gc()
        message("Simulating batch effects...")
        sim <- splatter:::splatSimBatchEffects(sim, params)
        gc()
        sim <- splatter:::splatSimBatchCellMeans(sim, params)
        rowData(sim)$GeneMean <- NULL
        rowData(sim) <- rowData(sim)[
            , !grepl("BatchFacBatch", colnames(rowData(sim))), drop = FALSE
        ]
        gc()
        message("Simulating path endpoints...")
        sim <- splatter:::splatSimPathDE(sim, params)
        gc()
        message("Simulating path steps...")
        sim <- splatter:::splatSimPathCellMeans(sim, params)
        colData(sim)$ExpLibSize <- NULL
        assay(sim, "BatchCellMeans") <- NULL
        rowData(sim) <- rowData(sim)[
            , !grepl("DEFacPath", colnames(rowData(sim))), drop = FALSE
        ]
        rowData(sim) <- rowData(sim)[
            , !grepl("LocalDEFacPath", colnames(rowData(sim))), drop = FALSE
        ]
        rowData(sim) <- rowData(sim)[
            , !grepl("SigmaFacPath", colnames(rowData(sim))), drop = FALSE
        ]
        colData(sim)$Step <- NULL
        gc()
        message("Simulating BCV...")
        sim <- splatter:::splatSimBCVMeans(sim, params)
        assay(sim, "BaseCellMeans") <- NULL
        assay(sim, "BCV") <- NULL
        gc()
        message("Adding additional noise...")
        groups_numeric <- as.numeric(sim$Group)
        noise_factors <- matrix(
            rnorm(
                ncol(sim) * nrow(sim),
                mean = 1,
                sd = group_params$Noise[groups_numeric]
            ),
            ncol = ncol(sim), nrow = nrow(sim),
            byrow = TRUE
        )
        noise_factors[noise_factors < 0] <- 0.01
        assay(sim, "CellMeans") <- assay(sim, "CellMeans") * noise_factors
        rm(noise_factors)
        gc()
        message("Simulating counts...")
        sim <- splatter:::splatSimTrueCounts(sim, params)
        gc()
        message("Simulating dropout (if needed)...")
        sim <- splatter:::splatSimDropout(sim, params)
        assay(sim, "CellMeans") <- NULL
        assay(sim, "TrueCounts") <- NULL
        assay(sim, "DropProb") <- NULL
        assay(sim, "Dropout") <- NULL
        gc()

        message("Sparsifying assays...")
        assays(sim) <- splatter:::sparsifyMatrices(
            assays(sim),
            auto = TRUE, verbose = TRUE
        )
        gc()
    })
    sim$Label <- factor(sim$Group, levels = group_params$Group,
                        labels = group_params$Name)

    message("Adjusting group proportions for each batch...")
    # Get the current counts of each group in each batch
    initial_batch_group_count <- table(sim$Batch, sim$Group)

    # Get the final number of each group in each batch
    batch_group_ncells <- lapply(seq_along(batch_groups), function(idx) {
        ncells <- batch_params$Cells[idx] * batch_groups[[idx]]
        if (any(ncells > initial_batch_group_count[idx, ])) {
            stop("Not enough cells for these proportions in batch ", idx)
        }
        ncells
    })

    message("Selecting ", final_cells, " cells...")
    # Select a subset of cells in each batch to adjust group proportions
    selected <- lapply(seq_along(batch_group_ncells), function(batch) {
        group_cells <- batch_group_ncells[[batch]]
        is_batch <- sim$Batch == paste0("Batch", batch)
        selected_batch <- sapply(seq_along(group_cells), function(group) {
            is_group <- sim$Group == paste0("Path", group)
            sample(sim$Cell[is_batch & is_group], group_cells[group])
        })
        unlist(selected_batch)
    })
    selected <- as.character(unlist(selected))

    # Subset SingleCellExperiment
    sim <- sim[, selected]
    gc()

    message("Downsampling counts to adjust sequencing depth...")
    cell_count_props <- batch_params[sim$Batch, "RelLibSize"]

    counts(sim) <- scuttle::downsampleMatrix(
        counts(sim),
        cell_count_props,
        bycol = TRUE
    )

    message("Performing quick cell filtering...")
    sim <- scuttle::quickPerCellQC(sim)

    return(sim)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    out_file <- args[["--out-file"]]

    install_splatter()

    output <- simulate_dataset()
    print(output)
    write_h5ad(output, out_file, X_name = "counts")
    message("Done!")
}

if (sys.nframe() == 0) {
    main()
}
