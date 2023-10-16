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
    # Initial number needs to be larger as this is sampled to get the final
    # number with adjusted group proportions
    initial_cells <- final_cells * 4

    message("Defining batches...")
    batches <- data.frame(
        Batch = c(
            "Batch1", "Batch2", "Batch3", "Batch4", "Batch5", "Batch6"
        ),
        RelProp    = c(0.90, 1.00, 0.22, 0.20, 2.60, 2.80),
        FacLoc     = c(0.08, 0.05, 0.12, 0.18, 0.15, 0.25),
        FacScale   = c(0.06, 0.08, 0.12, 0.10, 0.16, 0.18),
        RelLibSize = c(0.22, 0.18, 0.90, 1.05, 0.05, 0.04)
    )
    batches$Prop <- batches$RelProp / sum(batches$RelProp)
    batches$Cells <- round(batches$Prop * final_cells)
    rownames(batches) <- batches$Batch

    message("Defining groups...")
    groups = data.frame(
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
        Skew    = c( 0.5,  0.5,  0.3,  0.7,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0)
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

    n_batches <- nrow(batches)
    n_groups <- nrow(groups)

    message("Setting initial parameters...")
    initial_params <- splatter::newSplatParams(
        lib.loc        = 11,
        # Initial batches of equal size
        batchCells     = round(rep(initial_cells, n_batches) / n_batches),
        batch.facLoc   = batches$FacLoc,
        batch.facScale = batches$FacScale,
        # Initial groups with estimated max probabilities
        group.prob     = initial_group_probs,
        path.from      = groups$From,
        path.nSteps    = groups$Steps,
        path.skew      = groups$Skew,
        # Differential expression by group
        de.prob        = groups$DEProb,
        de.facLoc      = groups$DELoc,
        de.facScale    = groups$DEScale,
        # Seed
        seed           = 1
    )

    message("Simulating initial dataset...")
    initial_sim <- splatter::splatSimulatePaths(initial_params)
    initial_sim$Label <- factor(initial_sim$Group, levels = groups$Group,
                                labels = groups$Name)

    message("Adjusting group proportions for each batch...")
    # Get the current counts of each group in each batch
    initial_batch_group_count <- table(initial_sim$Batch, initial_sim$Group)

    # Get the final number of each group in each batch
    batch_group_ncells <- lapply(seq_along(batch_groups), function(idx) {
        ncells <- batches$Cells[idx] * batch_groups[[idx]]
        if (any(ncells > initial_batch_group_count[idx, ])) {
            stop("Not enough cells for these proportions in batch ", idx)
        }
        ncells
    })

    # Select a subset of cells in each batch to adjust group proportions
    selected <- lapply(seq_along(batch_group_ncells), function(batch) {
        group_cells <- batch_group_ncells[[batch]]
        is_batch <- initial_sim$Batch == paste0("Batch", batch)
        selected_batch <- sapply(seq_along(group_cells), function(group) {
            is_group <- initial_sim$Group == paste0("Path", group)
            sample(initial_sim$Cell[is_batch & is_group], group_cells[group])
        })
        unlist(selected_batch)
    })
    selected <- as.character(unlist(selected))

    # Subset SingleCellExperiment
    sim <- initial_sim[, selected]

    message("Downsampling counts to adjust sequencing depth...")
    cell_count_props <- batches[sim$Batch, "RelLibSize"]

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
