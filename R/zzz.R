#' R/zzz.R 
#' 
#' install required python dependencies for scvelo_py.py
.onLoad <- function(libname, pkgname) {
  library(reticulate)

  # required Python packages for scVelo
  packages <- c("scvelo", "anndata", "loompy", "matplotlib", "numba", "numpy", "pandas", "scanpy", "scikit-learn", "scipy", "scvi-tools", "umap-learn")
 
  if (reticulate::conda_available()) {
    # Check if Conda is available
    message("Conda is available, proceeding with Conda installation.")
    # Attempt to install necessary dependencies via conda
    tryCatch({
      reticulate::conda_install(envname = "scvelo", packages = packages, channel = "conda-forge")
      message("Required Conda packages installed successfully.")
    }, error = function(e) {
      message("Failed to install Conda packages. Please check your Conda setup.")
      stop("Unable to install Conda packages.")
    }) 
  } else if (reticulate::py_available()) {
    # Check if Python is available (without Conda)
    message("Python is available, proceeding with Python package installation.")
    # Attempt to install necessary dependencies via pip
    tryCatch({
      reticulate::py_install(packages)
      message("Required Python packages installed successfully.")
    }, error = function(e) {
      message("Failed to install Python packages. Please check your Python setup.")
      stop("Unable to install Python packages.")
    })
  } else {
    stop("Neither Conda nor Python is installed. Please install either Conda or Python to continue.")
  }


}
