#' Function to run scProportionTest pipeline 
#'
#' This function runs scProportionTest for all the pairwise comparision conditions
#' in parallel and generate figures and table of results.
#' The user can aslo chose specific conditions for comparsion rather than all the comparision.
#' 
#' @param dir_scproportion Directory path where output figures will be saved.
#' @param seurat_obj A Seurat object containing count data and metadata.
#' @param cluster_col The name of the metadata column in `seurat_obj` containing cluster labels or cell type names.
#' @param sample_col The name of the metadata column in `seurat_obj` that contains sample identifiers.
#' @param comparision1 Optional: name of the first group for comparison (default NULL, which compares all pairs).
#' @param comparision2 Optional: name of the second group for comparison (default NULL, which compares all pairs).
#' @param output.format Format of the output figure
#' @param verbose Print the processing stept
#' @param cores  Numenr of the cores for parallel computation
#'
#' @return NULL Save comparison figures and results in the specified directory.
#' 
#' @export
#' 
#'
run_scproportion <- function(dir_scproportion=".",seurat_obj,cluster_col,sample_col,comparision1=NULL,
                             comparision2=NULL,output.format = "png",verbose = TRUE,cores = detectCores() - 1){

  # check the input data format 
  #checkmate::expect_class(seurat_obj,"Seurat",label="seurat_obj")
  #checkmate::expect_numeric(cores, min.len = 1, max.len = detectCores() - 1, any.missing = FALSE,label="cores")
  # check the sample_col and cluster_col are in the meta data
  meta <- seurat_obj@meta.data
  if(!sample_col %in% colnames(meta)){
    stop("Sample name does not exits in the seurat object meta data!")
  }
  if(!cluster_col %in% colnames(meta)){
    stop("Cluster column does not exits in the seurat object meta data!")
  }
  
  conditions <- unique(meta[,sample_col])
  if (!is.null(comparision1) && !is.null(comparision2)) {
    if (!(comparision1 %in% conditions) || !(comparision2 %in% conditions)) {
      stop("One or both specified comparison conditions do not exist in the meta data!")
    }
    comparisons_condition <- matrix(c(comparision1, comparision2), ncol = 2, byrow = TRUE)
  } else {
    comparisons_condition <- permutations(length(conditions), 2, conditions)
  }
  
  create_dir(dir_scproportion)
  #subdirectories <- c(file.path("scproportion", "images"),file.path("scproportion","results"))
  
  #for(dir.i in subdirectories){
  #  dir.create(dir.i, showWarnings = F, recursive = T)
  #}
  #dir_scproportion <- file.path(dir_scproportion,"scproportion")
  prop_test <- scProportionTest::sc_utils(seurat_obj)
  
  # Function to process each comparison
  process_comparison <- function(i) {
    if (verbose) message("Running comparison: ", comparisons_condition[i, 1], " vs ", comparisons_condition[i, 2])
    
    prop_test_i <- scProportionTest::permutation_test(prop_test,
                                    cluster_identity = cluster_col,
                                    sample_1 = comparisons_condition[i, 1],
                                    sample_2 = comparisons_condition[i, 2],
                                    sample_identity = sample_col)
    
    # save the figure
    generate_figure(prop_test_i, output.format,comparisons_condition, dir_scproportion, cluster_col, i)
    
    # save the results
    stat_res(prop_test_i,comparisons_condition, dir_scproportion, i)
    
    if (verbose) message("Completed comparison: ", comparisons_condition[i, 1], " vs ", comparisons_condition[i, 2])
  }
  
  parallel::mclapply(1:nrow(comparisons_condition), process_comparison, mc.cores = cores)
  if (verbose) message("scProportion test completed.")
}
