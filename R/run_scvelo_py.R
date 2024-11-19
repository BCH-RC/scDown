#' Function to run the scVelo pipeline using python script
#'
#' This function performs RNA velocity calculations from .h5ad file using scVelo python script.
#' Workflow of this function: 
#' 1. calculate RNA velocity using scVelo workflow
#' 2. cluster-specific differential velocity genes
#' 3. trajectory inference using PAGA
#' This function outputs basic figures such as spliced/unspliced count proportion, projected RNA velocity on umap,
#' the phase portrait (ratio of spliced/unspliced RNA abundance) for top differential genes, 
#' and directed graphs of predicted lineages from PAGA trajectory inference 
#'
#' @param h5ad_file input h5ad file, if running after run_scvelo(), this object will have a fixed name and does not need to be changed
#' @param output_dir A character vector specifying the output directory
#' @param annotation_column A character variable specifying which metadata column of the h5ad object contains cell type annotations
#' @param mode Mode to conduct scvelo velocity calculation, either 'stochastic (default)', 'deterministic', or 'dynamical (slowest)'
#' @param top_gene The number of top differential velocity genes to plot phase portrait for
#'
#' @return A list of scVelo data objects
#'
#' @export
#'
#' Estimate RNA velocity for spliced and unspliced counts of scRNA-seq data


run_scvelo_py <- function(h5ad_file="scvelo/rds/obj_spliced_unspliced.h5ad", 
                        output_dir=".", 
                        annotation_column = 'ID', 
                        mode = 'stochastic', 
                        top_gene = 5){

# create subdirectories in the output directory
setwd(output_dir)
subdirectories <- c("rds",
                    "csv",
                    "images",
                    "scvelo/csv",
                    "scvelo/rds",
                    "scvelo/images")

for(i in subdirectories){
    dir.create(file.path(output_dir,i), showWarnings = F, recursive = T)
}

# Call the main python function from scvelo.py with parameters
library(reticulate)
reticulate::source_python("inst/python/scvelo.py")

run_scvelo_workflow(h5ad_file,annotation_column,mode,top_gene)



sessionInfo()

}
