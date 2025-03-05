#' Function to run scProportionTest pipeline 
#'
#' This function runs scProportionTest for all the pairwise comparision conditions
#' in parallel and generate figures and table of results.
#' The user can aslo chose specific conditions for comparsion rather than all the comparision.
#' 
#' @param output_dir Directory path where output figures will be saved.
#' @param seurat_obj A Seurat object containing count data and metadata.
#' @param annotation_column The name of the metadata column in `seurat_obj` containing cluster labels or cell type names.
#' @param group_column The name of the metadata column in `seurat_obj` that contains sample identifiers.
#' @param comparision1 Optional: name of the first group for comparison (default NULL, which compares all pairs).
#' @param comparision2 Optional: name of the second group for comparison (default NULL, which compares all pairs).
#' @param output_format Format of the output figure
#' @param verbose Print the processing stept
#' @param cores  Numenr of the cores for parallel computation
#'
#' @return NULL Save comparison figures and results in the specified directory.
#' 
#' @export
#' 
#'
run_scproportion <- function(seurat_obj,annotation_column,group_column,comparision1=NULL,comparision2=NULL,
                             output_dir=".",output_format = "png",verbose = TRUE,cores = 8){

  # check the input data format 
  check_required_variables(seurat_obj,species=NULL,output_dir,annotation_column,group_column)
  checkmate::expect_choice(comparision1,unique(seurat_obj@meta.data[,group_column]),label = "comparision1",null.ok = TRUE)
  checkmate::expect_choice(comparision2,unique(seurat_obj@meta.data[,group_column]),label = "comparision2",null.ok = TRUE)
  checkmate::expect_choice(output_format,c("png","pdf","jpeg"),label = "output_format")
  meta <- seurat_obj@meta.data
  if(!group_column %in% colnames(meta)){
    stop("Sample name does not exits in the seurat object meta data!")
  }
  if(!annotation_column %in% colnames(meta)){
    stop("Cluster column does not exits in the seurat object meta data!")
  }
  
  conditions <- as.character(unique(meta[,group_column]))
  if (!is.null(comparision1) && !is.null(comparision2)) {
    if (!(comparision1 %in% conditions) || !(comparision2 %in% conditions)) {
      stop("One or both specified comparison conditions do not exist in the meta data!")
    }
    comparisons_condition <- matrix(c(comparision1, comparision2), ncol = 2, byrow = TRUE)
  } else {
    comparisons_condition <- gtools::permutations(length(conditions), 2, conditions)
  }
  
  create_dir(output_dir)
  subdirectories <- c(file.path("scproportion", "images"),file.path("scproportion","results"))
  
  for(dir.i in subdirectories){
    dir.create(dir.i, showWarnings = F, recursive = T)
  }
  #output_dir <- file.path(output_dir,"scproportion")
  prop_test <- scProportionTest::sc_utils(seurat_obj)
  
  # Function to process each comparison
  library(scProportionTest)
  library(ggplot2)
  library(dplyr)
  process_comparison <- function(i) {
    if (verbose) message("Running comparison: ", comparisons_condition[i, 2], " vs ", comparisons_condition[i, 1])
    
    prop_test_i <- scProportionTest::permutation_test(prop_test,
                                    cluster_identity = annotation_column,
                                    sample_1 = comparisons_condition[i, 1],
                                    sample_2 = comparisons_condition[i, 2],
                                    sample_identity = group_column)
    
    # save the figure
    generate_figure(prop_test_i, output_format,comparisons_condition, output_dir, annotation_column, i)
    
    # save the results
    stat_res(prop_test_i,comparisons_condition, output_dir, i)
    
    if (verbose) message("Completed comparison: ", comparisons_condition[i, 2], " vs ", comparisons_condition[i, 1])
  }
  
  parallel::mclapply(1:nrow(comparisons_condition), process_comparison, mc.cores = cores)
  if (verbose) message("scProportion test completed.")
}
