#!/usr/bin/env R
install_kBET <- function() {
  if (!requireNamespace("kBET", quietly = TRUE)) {
      message("Installing kBET...")
      remotes::install_github(
          "theislab/kBET@f35171dfb04c7951b8a09ac778faf7424c4b6bc0",
          dependencies = FALSE
      )
  } else {
      message("kBET already installed")
  }
}
install_kBET()
