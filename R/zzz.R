#' R/zzz.R
#'
#' install required python dependencies
#'
#' @noRd

.onLoad <- function(libname, pkgname) {

  # Check if Python is available without initializing
  reticulate::py_available(initialize = FALSE)

  # Detect if running inside Docker or Singularity (check if /opt/micromamba exists)
  docker_singularity_scdownv2_path <- "/opt/micromamba/envs/scdownv2/bin/python"
  if (file.exists(docker_singularity_scdownv2_path)) {
    # Force Singularity users to use the correct environment
    Sys.setenv(RETICULATE_PYTHON = "/opt/micromamba/envs/scdownv2/bin/python")
    Sys.setenv(R_MINICONDA_PATH = "/opt/micromamba")
    Sys.setenv(RETICULATE_MINICONDA_PATH = "/opt/micromamba")
    reticulate::py_module_available("scdownv2")
    message("Using scdownv2 from /opt/micromamba/envs/scdownv2/.")
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

    env_name <- "scdownv2"

    # Use the conda environment
    print(conda_path_pre)
    envs <- system(paste(conda_path, "info --envs"), intern = TRUE)
    scdownv2_path <- envs[grepl("^scdownv2", envs) | grepl("/scdownv2", envs)]
    print(scdownv2_path)
    reticulate::use_condaenv(env_name, conda = conda_path, required = TRUE)

  }

}
