library(scProportionTest)
library(Seurat)
library(SeuratObject)
library(gtools)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyverse)

#' Create result directories for scProportionTest figures
#' @param output_dir Folder path for scProportionTest figures
#' @return NULL
#' @noRd

create_dir <- function(output_dir) {

  subdirectories <- c(file.path(output_dir, "scproportion", "images"),
                      file.path(output_dir, "scproportion", "results"))
  
  # Create each directory if it doesn't exist
  for (dir.i in subdirectories) {
    dir.create(dir.i, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Return the path to the main directory for further use if needed
  return(file.path(output_dir, "scproportion"))
}

#' Generate plot for each comparison
#' @param prop_test.i scProportion object
#' @param output_format  The format of output figure
#' @param comparisons_condition table of all the pairwise comparison
#' @param output_dir path to save the figures 
#' @param annotation_column the column name in the meta, cluster or celltype
#' @param i index of the comparison condition
#' @return NULL saves comparison figures in the specified directory.
#' @noRd

generate_figure <- function(prop_test.i, output_format="png",comparisons_condition, output_dir=".", annotation_column, i) {
  p <- permutation_plot(prop_test.i) +
    theme_bw(base_size = 12) +
    labs(title = paste0(comparisons_condition[i, 2], " vs ", comparisons_condition[i, 1]),
         x = annotation_column,
         y = "log2(FD)") +
    theme(legend.text = element_text(size = 8)) +
    scale_shape_manual(name = "significance",labels = c("FDR < 0.05 &\nabs(Log2FD) > 0.58", "n.s."),values = c(16, 1)) +
    scale_color_manual(name = "significance",labels = c("FDR < 0.05 &\nabs(Log2FD) > 0.58", "n.s."),values = c("red", "grey") )
  
  output_format <- match.arg(output_format, choices = c("png", "pdf", "jpeg"))
  file_extension <- switch(output_format, png = "png", pdf = "pdf", jpeg = "jpg")
  
  output_path <- paste0(output_dir, "/scproportion/images/scProportiontest_",
                        comparisons_condition[i, 2], "vs", comparisons_condition[i, 1], ".", file_extension)
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
#' @param output_dir directory to save the results table
#' @param i index of the comparison condition
#' @return NULL saves comparison figures in the specified directory.
#' @noRd

stat_res <- function(prop_test.i,comparisons_condition, output_dir, i){
  res_tab <- prop_test.i@results %>% as.data.frame()
  output_path <- paste0(output_dir, "/scproportion/results/scProportiontest_",
                        comparisons_condition[i, 2], "vs", comparisons_condition[i, 1], ".csv")
  write.csv(res_tab, output_path, row.names = FALSE)
}
