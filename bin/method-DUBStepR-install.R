#!/usr/bin/env R
install_DUBStepR <- function() {
  if (!requireNamespace("DUBStepR", quietly = TRUE)) {
      message("Installing DUBStepR...")
      remotes::install_github(
          "prabhakarlab/DUBStepR@76aa39485742d5f5bcfb86346a8a1784ee08f6b9",
          dependencies = FALSE
      )
  } else {
      message("DUBStepR already installed")
  }
}
install_DUBStepR()
