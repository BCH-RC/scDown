#' R/zzz.R 
#' 
#' install required python dependencies for scvelo_workflow.py
.onLoad <- function(env_path=NULL, packages=c("scvelo", "anndata", "loompy", "matplotlib", "numba", "numpy", "pandas", "scanpy", "scikit-learn", "scipy", "scvi-tools", "umap-learn")) {
  library(reticulate)

  env_name <- "scvelo"
  env_exists <- reticulate::conda_list()$name %in% env_name

  reticulate::py_available(initialize = FALSE)
  if (!is.null(env_path) & dir.exists(env_path)) {
    cat("The specified Conda environment exists. Activating it...\n")
    # Use the existing Conda environment
    reticulate::use_condaenv(env_path, required = TRUE)
  } 

  if (!reticulate::py_available(initialize = FALSE)) {
    reticulate::install_miniconda()
  }

  if (env_exists)
    reticulate::use_condaenv(env_name)
  } else {  
    reticulate::conda_create(env_name)
    reticulate::use_condaenv(env_name)
  }


  # Install the required packages 
  for (pkg in required_packages) {
    cat("Installing Conda package:", pkg, "\n")
    reticulate::conda_install(envname = env_name, packages = pkg, channel = "conda-forge")
  }

    reticulate::use_condaenv(env_name)

}
