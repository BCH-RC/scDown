#' R/zzz.R 
#' 
#' install required python dependencies for scvelo_workflow.py for run_scvelo_full function
#'
#' @noRd

.onLoad <- function(libname, pkgname) {

  # Check if Python is available without initializing
  reticulate::py_available(initialize = FALSE)

  # Detect if running inside Docker or Singularity (check if /opt/micromamba exists)
  docker_singularity_scvelo_path <- "/opt/micromamba/envs/scvelo/bin/python"
  if (file.exists(docker_singularity_scvelo_path)) {
    # Force Singularity users to use the correct environment
    Sys.setenv(RETICULATE_PYTHON = "/opt/micromamba/envs/scvelo/bin/python")
    Sys.setenv(R_MINICONDA_PATH = "/opt/micromamba")
    Sys.setenv(RETICULATE_MINICONDA_PATH = "/opt/micromamba")
    reticulate::py_module_available("scvelo")
    message("Using scvelo from /opt/micromamba/envs/scvelo/.")
  } else {
    # Set the path to the Miniforge installation if Python is not available
    conda_path <- reticulate::conda_binary()
    # Check if Miniconda is installed 
    if (!reticulate::py_available(initialize = FALSE)) {
      if (is.null(conda_path)) {
        # If Miniconda is not installed, install it
        miniconda_path <- "~/.local/share/r-miniconda"
        reticulate::install_miniconda()
      }
    }
    conda_path <- reticulate::conda_binary()
    conda_path_pre=gsub("/conda$","",conda_path)

    env_name <- "scvelo"
    env_file <- "inst/env/environment.yml"  # conda env for Linux, Windows or Intel-based Mac

    env_exists <- env_name %in% reticulate::conda_list()$name 
    # Check if the conda environment exists
    if (!env_exists) {
      system(paste0(conda_path_pre, "/conda env create -f ", env_file))
    }
    # Use the conda environment
    print(conda_path_pre)
    envs <- system(paste(conda_path, "info --envs"), intern = TRUE)
    scvelo_path <- envs[grepl("^scvelo", envs) | grepl("/scvelo", envs)]
    print(scvelo_path)
    reticulate::use_condaenv(env_name, conda = conda_path, required = TRUE)

  }

}
