#' R/zzz.R 
#' 
#' install required python dependencies for scvelo_workflow.py for run_scvelo_full function
#'
#' @noRd

.onLoad <- function(libname, pkgname) {

  # Check if Python is available without initializing
  reticulate::py_available(initialize = FALSE)

  # Detect if running inside Docker or Singularity (check if /opt/miniconda exists)
  docker_singularity_scvelo_path <- "/opt/miniconda/envs/scvelo/bin/python"
  if (file.exists(docker_singularity_scvelo_path)) {
    # Force Singularity users to use the correct environment
    Sys.setenv(RETICULATE_PYTHON = "/opt/miniconda/envs/scvelo/bin/python")
    Sys.setenv(R_MINICONDA_PATH = "/opt/miniconda")
    Sys.setenv(RETICULATE_MINICONDA_PATH = "/opt/miniconda")
    message("Using scvelo from /opt/miniconda.")
  } else {
    # Set the path to the Miniconda installation if Python is not available
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
    conda_path2=gsub("/conda$","",conda_path)

    env_name <- "scvelo"
    
    if(Sys.info()["machine"] != "aarch64") {
      env_file <- "inst/env/environment.yml"        # conda env for Linux, Windows or Intel-based Mac
    } else{
      env_file <- "inst/arm64_env/environment.yml"  # conda env for Apple-silicon Mac M1/M2/M3
    }
    
    env_exists <- env_name %in% reticulate::conda_list()$name 
    # # Check if the conda environment exists
    if (!env_exists) {  
      system(paste0(conda_path2, "/conda env create -f ", env_file))
      #reticulate::conda_create(env_name, envfile = env_file, dependencies = c("python=3.10.9"))
    }
    # Use the conda environment
    print(conda_path2)
    envs <- system(paste(conda_path, "info --envs"), intern = TRUE)
    scvelo_path <- envs[grepl("^scvelo", envs) | grepl("/scvelo", envs)]
    print(scvelo_path)
    reticulate::use_condaenv(env_name, conda = conda_path, required = TRUE)

  }

}
