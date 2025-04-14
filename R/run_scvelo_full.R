#' Function to run the scVelo pipeline using python script
#'
#' This function performs RNA velocity calculations from .h5ad file using scVelo python script.
#' Workflow of this function: 
#' 1. calculate RNA velocity using scVelo workflow
#' 2. cluster-specific differential velocity genes
#' 3. trajectory inference using PAGA
#' This function outputs basic figures such as spliced/unspliced count proportion, projected RNA velocity 
#' on umap, the phase portrait (ratio of spliced/unspliced RNA abundance) for top differential genes, and 
#' directed graphs of predicted lineages from PAGA trajectory inference 
#'
#' @param h5ad_file input h5ad file path and name, if running after `run_scvelo()`, this object has a fixed 
#' name and does not need to be changed
#' @param output_dir A character vector specifying the output directory
#' @param annotation_column A character variable specifying which metadata column of the h5ad object contains 
#' cell type annotations
#' @param mode Mode to conduct scvelo velocity calculation, either 'stochastic (default)', 'deterministic', 
#' or 'dynamical (slowest)'
#' @param top_gene The number of top differential velocity genes to plot phase portrait for
#' @param groups A list of character vectors representing groups of conditions or time points used to 
#' calculate RNA velocity separately, default: NULL
#' @param group_column A string specifying the name of the metadata column in the h5ad object that should be 
#' used for subsetting each group in `groups`
#'
#' @return A list of scVelo data objects
#'
#' @export
#'
#' Estimate RNA velocity for spliced and unspliced counts of scRNA-seq data


run_scvelo_full <- function(h5ad_file="scvelo/rds/obj_spliced_unspliced.h5ad", 
                        output_dir=".", 
                        annotation_column = 'ID', 
                        mode = 'stochastic', 
                        top_gene = 5,
                        groups=NULL, 
                        group_column=NULL){

# create subdirectories in the output directory
setwd(output_dir)
subdirectories <- c("scvelo",
                    "scvelo/csv",
                    "scvelo/rds",
                    "scvelo/images")

for(i in subdirectories){
    dir.create(file.path(output_dir,i), showWarnings = F, recursive = T)
}


checkmate::assert_string(h5ad_file, null.ok = FALSE)
checkmate::assert_file_exists(h5ad_file)
checkmate::assert_string(annotation_column, null.ok = FALSE)
checkmate::assert_string(mode, null.ok = FALSE)
checkmate::assert_numeric(top_gene, null.ok = FALSE, any.missing=FALSE)

checkmate::assert_list(groups, types = c("numeric","integer","character"), null.ok = TRUE)
if(!is.null(groups)){
  groups <- lapply(groups, function(element) {
    if (is.integer(element) || is.numeric(element)) {
      return(as.character(element))
    }
    return(element)
  })
}
if(!is.null(groups)){
  checkmate::assert_string(group_column, null.ok = FALSE)
}

# Call the main python function from scvelo_workflow.py with parameters
library(reticulate)

py_script <- system.file("python", "scvelo_workflow.py", package = "scDown")
if (py_script == "") {
  stop("Python script not found in the installed scDown package")
}
reticulate::source_python(py_script)
#reticulate::source_python("inst/python/scvelo_workflow.py")

# RNA velocity for the entire object
run_scvelo_workflow(h5ad_file,annotation_column,mode,top_gene,group_label="ALL")
system("stty echo")

# RNA velocity for specified conditions or time points, if any
if (length(groups) != 0){
  library(anndata)
  file_base <- gsub(".h5ad", "", basename(h5ad_file))
  input_dir <- dirname(h5ad_file)
  for (group in groups){        
      group_label=paste(group, collapse="_")
      file_name=paste0(file_base,"_",group_label,".h5ad")
      h5ad_group_file_input=file.path(input_dir, file_name)
      h5ad_group_file_output=file.path(paste0(output_dir,"/scvelo/rds"), file_name)
      if(file.exists(h5ad_group_file_input)){
        h5ad_group_file <- h5ad_group_file_input
      } else {
        h5ad_group_file <- h5ad_group_file_output
        if (!file.exists(h5ad_group_file_output)){
          adata <- anndata::read_h5ad(h5ad_file)
          subset_adata <- adata[adata$obs[[group_column]] %in% group, ]
          subset_adata$write_h5ad(h5ad_group_file)
        }
      }
      run_scvelo_workflow(h5ad_group_file,annotation_column,mode,top_gene,group_label)
  }
}
system("stty sane")
#system("stty echo")


sessionInfo()

}
