#' R/zzz.R 
#' 
#' install required python dependencies for scvelo_workflow.py
.onLoad <- function(libname, pkgname) {

  library(reticulate)

  # Check if Python is available without initializing
  reticulate::py_available(initialize = FALSE)

  # Set the path to the Miniconda installation if Python is not available
  miniconda_path <- "~/.local/share/r-miniconda"
  # Check if Miniconda is installed at the specified path
  if (!reticulate::py_available(initialize = FALSE)) {
    if (!file.exists(miniconda_path)) {
      # If Miniconda is not installed, install it at specific path
      reticulate::install_miniconda()
    } 
    # use the minoconda at specific path
    reticulate::use_miniconda(miniconda_path, required = TRUE)
  }

  env_name <- "scvelo"
  env_exists <- env_name %in% reticulate::conda_list()$name 
  # # Check if the conda environment exists
  if (!env_exists) {
    reticulate::conda_create(env_name)
  }
  # Use the conda environment
  reticulate::use_condaenv(env_name, required = TRUE)

  packages=c("scvelo", "anndata", "loompy", "matplotlib", "numba", "numpy", "pandas", "scanpy", "scikit-learn", "scipy", "scvi-tools", "umap-learn")
  # Install the required packages 
  for (pkg in packages) {
    installed_packages <- py_module_available(pkg)
    if (!installed_packages) {
      cat("Installing Conda package:", pkg, "\n")
      reticulate::conda_install(envname = env_name, packages = pkg, channel = "conda-forge")
    }
  }


}
