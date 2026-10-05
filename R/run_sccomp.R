#' Function to run scComp pipeline
#'
#' This function runs scComp to perform differential composition and variability analysis
#' and generate figures and table of results.
#' The user can optionally choose a control condition to compare all others against it.
#'
#' @param seurat_obj A Seurat object containing count data and metadata.
#' @param annotation_column The name of the metadata column in `seurat_obj` containing cluster labels or cell type names.
#' @param sample_column The name of the metadata column in `seurat_obj` that contains sample identifiers.
#' @param group_column The name of the metadata column in `seurat_obj` that contains conditions for comparisons.
#' @param output_dir Directory path where output figures will be saved.
#' @param output_format Format of the output figure.
#' @param control_cond Control condition group label.
#' @param significance_threshold Threshold of significance for plotting purpose.
#' @param additional_effects Additional effects to be added to formula.
#' @param verbose Print the processing steps.
#' @param n_jobs A numeric variable to set the number of cores to be used.
#'
#' @return NULL Save comparison figures and results in the specified directory.
#'
#' @export
#'
run_sccomp <- function(seurat_obj,
                       annotation_column,
                       sample_column,
                       group_column,
                       output_dir=".",
                       output_format="png",
                       control_cond=NULL,
                       significance_threshold=0.05,
                       additional_effects=NULL,
                       verbose=TRUE,
                       n_jobs=1){

  # Check the input data format
  checkmate::expect_class(seurat_obj,"Seurat",label="seurat_obj")
  checkmate::expect_choice(annotation_column, colnames(seurat_obj@meta.data),label="annotation_column",null.ok = TRUE)
  checkmate::expect_choice(sample_column, colnames(seurat_obj@meta.data),label="sample_column",null.ok = TRUE)
  checkmate::expect_choice(group_column, colnames(seurat_obj@meta.data),label="group_column",null.ok = TRUE)
  checkmate::expect_directory(output_dir,access="rw",label = "output_dir")
  checkmate::expect_choice(output_format,c("png","pdf","jpg","jpeg"),label = "output_format")
  meta <- seurat_obj@meta.data

  if(!sample_column %in% colnames(meta)){
    stop("Sample name does not exits in the seurat object meta data!")
  }
  if(!annotation_column %in% colnames(meta)){
    stop("Cluster column does not exits in the seurat object meta data!")
  }
  if(!group_column %in% colnames(meta)){
    stop("Group (condition) column does not exits in the seurat object meta data!")
  }

  raw_groups <- as.character(seurat_obj@meta.data[[group_column]])
  clean_groups <- gsub("[^[:alnum:]]", "_", raw_groups)
  seurat_obj@meta.data$sccomp <- factor(clean_groups)

  all_conditions <- levels(seurat_obj@meta.data$sccomp)
  if (length(all_conditions) < 2) stop("Need at least two conditions to run comparisons.")

  # Determine which groups will act as the control baseline
  controls_to_run <- if (is.null(control_cond)) all_conditions else gsub("[^[:alnum:]]", "_", control_cond)
  formula_string <- if (is.null(additional_effects)) "~ sccomp" else paste("~ sccomp", additional_effects)

  create_dir_sccomp(output_dir)

  # Function to run scComp
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(ggrepel)
  library(cmdstanr)
  library(sccomp)

  # Dynamically locate and set the CmdStan path inside the container
  cmdstan_dir <- list.files('/opt/cmdstan', full.names=TRUE, pattern='cmdstan-')[1]
  if (length(cmdstan_dir) > 0) set_cmdstan_path(cmdstan_dir) else stop("CmdStan not found")
  # Force sccomp to use the shared container cache
  utils::assignInNamespace("sccomp_stan_models_cache_dir", "/opt/sccomp_models", ns = "sccomp")

  # Cycle through each control condition
  for (current_control in controls_to_run) {
    if (verbose) message(sprintf("\nTesting all groups against control: %s", current_control))

    current_output_dir <- file.path(output_dir, "sccomp", paste0("vs_", current_control))
    dir.create(file.path(current_output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(current_output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

    # Relevel the factor so the current_control is the Intercept baseline
    seurat_obj@meta.data$sccomp <- fct_relevel(seurat_obj@meta.data$sccomp, current_control)

    set.seed(42)
    # Re-estimate for each control because Bayesian prior changes with each Intercept
    sccomp_estimate_obj = seurat_obj |>
      sccomp_estimate(
        formula_composition = as.formula(formula_string),
        formula_variability = as.formula(formula_string),
        sample = sample_column,
        cell_group = annotation_column,
        cores = n_jobs,
        bimodal_mean_variability_association = ifelse(sample_column==group_column, FALSE, TRUE),
        output_directory = paste0(output_dir, "/sccomp/sccomp_draws_files"),
        mcmc_seed = 42,
        verbose = TRUE
      ) |>
      sccomp_remove_outliers(cores = n_jobs, verbose = FALSE, max_sampling_iterations = 2000)

    res_main <- sccomp_estimate_obj |> sccomp_test()
    current_meta <- seurat_obj@meta.data
    current_meta$sccomp <- fct_relevel(current_meta$sccomp, current_control)

    main_suffix <- paste0("all_vs_", current_control, ".")
    # Save the main figure
    generate_figure_sccomp(
      res_main,
      output_format,
      current_output_dir,
      main_suffix,
      annotation_column,
      group_column,
      current_meta,
      significance_threshold,
      plot_title = sprintf("All Individual Conditions vs %s", current_control),
      plot_subtitle = sprintf("Global Baseline: %s", current_control),
      base_name = current_control
    )

    # Save the main results
    stat_res_sccomp(res_main, current_output_dir, main_suffix)

    # Pairwise comparisons for non-control conditions
    treatments <- setdiff(all_conditions, current_control)

    if (length(treatments) >= 2) {
      if (verbose) message("Generating pairwise contrasts...")
      pairs <- combn(treatments, 2, simplify = FALSE)

      for (p in pairs) {
        cond1 <- p[1]
        cond2 <- p[2]

        contrast_str <- sprintf("`sccomp%s` - `sccomp%s`", cond1, cond2)
        res_pair <- sccomp_estimate_obj |> sccomp_test(contrasts = contrast_str)

        # Rename the math string and relevel inside the results table
        res_pair$parameter <- ifelse(res_pair$parameter == contrast_str, paste0("sccomp", cond1), res_pair$parameter)
        res_pair$factor <- ifelse(res_pair$parameter == paste0("sccomp", cond1), "sccomp", res_pair$factor)

        # Clone and relevel the meta data so the plotter calculates sizes and baselines correctly
        pair_meta <- current_meta
        pair_meta$sccomp <- fct_relevel(pair_meta$sccomp, cond2)

        pair_suffix <- paste0(cond1, "_vs_", cond2, ".")
        # Save the pair figure
        generate_figure_sccomp(
          res_pair,
          output_format,
          current_output_dir,
          pair_suffix,
          annotation_column,
          group_column,
          pair_meta,
          significance_threshold,
          plot_title = sprintf("Pairwise Contrast: %s vs %s", cond1, cond2),
          plot_subtitle = sprintf("Global Baseline: %s", current_control),
          trt_name = cond1,
          ctrl_name = cond2,
          base_name = current_control
        )

        # Save the pair results
        stat_res_sccomp(res_pair, current_output_dir, pair_suffix)
      }
    }
  }

  # Clean up tmp folders
  unlink(paste0(output_dir, "/sccomp/sccomp_draws_files"), recursive = TRUE)

  if (verbose) message("scComp test completed.")
}
