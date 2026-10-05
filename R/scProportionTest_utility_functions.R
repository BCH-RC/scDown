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
#'
create_dir <- function(output_dir) {

  subdirectories <- c(file.path(output_dir, "scproportion", "plots"),
                      file.path(output_dir, "scproportion", "tables"))

  # Create each directory if it doesn't exist
  for (dir.i in subdirectories) {
    dir.create(dir.i, showWarnings = FALSE, recursive = TRUE)
  }

  # Return the path to the main directory for further use if needed
  return(file.path(output_dir, "scproportion"))
}

#' Generate plot for each comparison
#' @param prop_test.i scProportion object
#' @param output_format The format of output figure
#' @param comparisons_condition table of all the pairwise comparison
#' @param output_dir path to save the figures
#' @param annotation_column the column name in the meta, cluster or celltype
#' @param i index of the comparison condition
#' @param meta metadata information
#' @param fdr_threshold Threshold of significance for plotting purpose.
#' @param log2fd_threshold Threshold of log2 fold difference for plotting purpose.
#' @return NULL saves comparison figures in the specified directory.
#' @noRd
#'
generate_figure <- function(prop_test.i,
                            output_format="png",
                            comparisons_condition,
                            output_dir=".",
                            annotation_column,
                            i,
                            meta,
                            significance_threshold,
                            log2fd_threshold) {
  p <- permutation_plot(prop_test.i) +
    theme_bw(base_size = 12) +
    labs(title = paste0(comparisons_condition[i, 2], " vs ", comparisons_condition[i, 1]),
         x = annotation_column,
         y = "log2(FD)") +
    theme(legend.text = element_text(size = 8)) +
    scale_shape_manual(name = "significance",labels = c(paste0("FDR < ",significance_threshold," &\nabs(Log2FD) > ",log2fd_threshold), "n.s."),values = c(16, 1)) +
    scale_color_manual(name = "significance",labels = c(paste0("FDR < ",significance_threshold," &\nabs(Log2FD) > ",log2fd_threshold), "n.s."),values = c("red", "grey") )

  output_format <- match.arg(output_format, choices = c("png", "pdf", "jpg", "jpeg"))
  if (output_format == "jpeg") {
    output_format <- "jpg"
  }
  file_extension <- switch(output_format, png = "png", pdf = "pdf", jpg = "jpg")

  # set figure height and width
  label_vec <- as.character(meta[[annotation_column]])
  n_anno <- length(unique(label_vec))
  max_label_nchar <- max(nchar(unique(label_vec)))

  # adaptive height
  img_height <- max(1200, 250 + 50 * n_anno)

  # adaptive width
  img_width <- max(2100, 900 + 18 * max_label_nchar)

  output_path <- paste0(output_dir, "/scproportion/plots/scProportiontest_",
                        comparisons_condition[i, 2], "vs", comparisons_condition[i, 1], ".", file_extension)
  if (output_format == "png") {
    png(output_path, width = img_width, height = img_height, res = 350)
  } else if (output_format == "pdf") {
    pdf(output_path, width = 10, height = 6.25)
  } else if (output_format == "jpg") {
    jpeg(output_path, width = img_width, height = img_height, res = 350)
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
#'
stat_res <- function(prop_test.i,comparisons_condition, output_dir, i){
  res_tab <- prop_test.i@results %>% as.data.frame()
  output_path <- paste0(output_dir, "/scproportion/tables/scProportiontest_",
                        comparisons_condition[i, 2], "vs", comparisons_condition[i, 1], ".txt")
  write.csv(res_tab, output_path, sep='\t', row.names = FALSE)
}
