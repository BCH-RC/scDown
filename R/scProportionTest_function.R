library(scProportionTest)
library(Seurat)
library(SeuratObject)
library(gtools)
library(parallel)
library(ggplot2)

#' Create result directories for scProportionTest figures
#'
#' @param dir_scproportion Folder path for scProportionTest figures
#' @return NULL
#'
create_dir_scproportion <- function(dir_scproportion) {
  subdirectories <- c("scproportion","scproportion/images","scproportion/results")
  for(dir.i in subdirectories){
    dir.create(paste0(dir_scproportion, dir.i), showWarnings = F, recursive = T)
  }
}

#' Generate plot for each comparision
#' @param prop_test.i scProportion object
#' @param output.format  The format of output figure
#' @param comparisons_condition
#' @param dir_scproportion
#' @param cluster_col
#' @param i The index of the comparision condition
#'
#' @return NULL Saves comparison figures in the specified directory.
#'
generate_figure <- function(prop_test.i, output.format,comparisons_condition, dir_scproportion, cluster_col, i) {
  p <- permutation_plot(prop_test.i) +
    theme_bw(base_size = 12) +
    labs(title = paste0(comparisons_condition[i, 1], " vs ", comparisons_condition[i, 2]),
         x = "log2(FD)",
         y = cluster_col)

  output_format <- match.arg(output.format, choices = c("png", "pdf", "jpeg"))
  file_extension <- switch(output_format, png = "png", pdf = "pdf", jpeg = "jpg")

  output_path <- paste0(dir_scproportion, "scproportion/images/scProportiontest_",
                        comparisons_condition[i, 1], "vs", comparisons_condition[i, 2], ".", file_extension)
  if (output_format == "png") {
    png(output_path, width = 2000, height = 1250, res = 350)
  } else if (output_format == "pdf") {
    pdf(output_path, width = 10, height = 6.25)
  } else if (output_format == "jpeg") {
    jpeg(output_path, width = 2000, height = 1250, res = 350)
  }
  print(p)
  dev.off()
}

#' Generate table of stat results for each comparision
#' @param prop_test.i scProportion object
#' @param comparisons_condition table of comparisions
#' @param dir_scproportion directory to save the results table
#' @param i The index of the comparison condition
#'
#' @return NULL Saves comparison figures in the specified directory.
#'
stat_res <- function(prop_test.i,comparisons_condition, dir_scproportion, i){
  res_tab <- prop_test.i@results %>% as.data.frame()
  output_path <- paste0(dir_scproportion, "scproportion/results/scProportiontest_",
                        comparisons_condition[i, 1], "vs", comparisons_condition[i, 2], ".csv")
  write.csv(res_tab, output_path, row.names = FALSE)
}

#' Run scProportionTest
#'
#' @param dir_scproportion Directory path where output figures will be saved.
#' @param seurat_obj A Seurat object containing count data and metadata.
#' @param cluster_col The name of the metadata column in `seurat_obj` containing cluster labels or cell type names.
#' @param sample_col The name of the metadata column in `seurat_obj` that contains sample identifiers.
#' @param comparision1 Optional: name of the first group for comparison (default NULL, which compares all pairs).
#' @param comparision2 Optional: name of the second group for comparison (default NULL, which compares all pairs).
#' @param output.format Format of the output figure
#' @param verbose Print the processing stept
#' @param num_cores  Number of the cores for parallel computation
#'
#' @return Saves comparison figures in the specified directory.
#'
run_scproportion <- function(dir_scproportion,seurat_obj,cluster_col,sample_col,comparision1=NULL,
                             comparision2=NULL,output.format = "png",verbose = TRUE,
                             num_cores = detectCores() - 1){
  # Capture the start time
  start_time <- Sys.time()
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

  create_dir_scproportion(dir_scproportion)
  prop_test <- sc_utils(seurat_obj)

  # Function to process each comparison in parallel
  process_comparison <- function(i) {
    if (verbose) message("Running comparison: ", comparisons_condition[i, 1], " vs ", comparisons_condition[i, 2])

    # Run the proportion test for the current pair
    prop_test_i <- permutation_test(prop_test,
                                    cluster_identity = cluster_col,
                                    sample_1 = comparisons_condition[i, 1],
                                    sample_2 = comparisons_condition[i, 2],
                                    sample_identity = sample_col)

    # Generate and save the figure
    generate_figure(prop_test_i, output.format,comparisons_condition, dir_scproportion, cluster_col, i)

    # Save the results to a CSV file
    stat_res(prop_test_i,comparisons_condition, dir_scproportion, i)

    if (verbose) message("Completed comparison: ", comparisons_condition[i, 1], " vs ", comparisons_condition[i, 2])
  }

  # Run comparisons in parallel
  mclapply(1:nrow(comparisons_condition), process_comparison, mc.cores = num_cores)
  # Capture the end time and calculate duration
  end_time <- Sys.time()
  duration <- end_time - start_time
  # Print execution time
  print(paste("scProportion test completed in", duration))
  if (verbose) message("scProportion test completed.")
}

