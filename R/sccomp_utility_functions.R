#' Create result directories for scComp
#' @param output_dir Folder path for scComp figures
#' @return NULL
#' @noRd
#'
create_dir_sccomp <- function(output_dir) {

  subdirectories <- c(file.path(output_dir, "sccomp"))

  # Create each directory if it doesn't exist
  for (dir.i in subdirectories) {
    dir.create(dir.i, showWarnings = FALSE, recursive = TRUE)
  }

  # Return the path to the main directory for further use if needed
  return(file.path(output_dir, "sccomp"))
}

#' Generate and save plots for scComp
#' @param sccomp_res scComp test result object
#' @param output_format The format of output figure
#' @param output_dir Path to save the figures
#' @param annotation_column The column name in the meta, cluster or celltype
#' @param meta Metadata information
#' @param significance_threshold Threshold of significance for plotting purpose.
#' @param plot_title Title of plots
#' @param plot_subtitle Subtitle of plots
#' @param trt_name Treatment condition name
#' @param ctrl_name Control condition name
#' @param base_name Baseline condition name
#' @return NULL saves comparison figures in the specified directory.
#' @noRd
#'
generate_figure_sccomp <- function(sccomp_res,
                                   output_format="png",
                                   output_dir=".",
                                   output_suffix="",
                                   annotation_column,
                                   group_column,
                                   meta,
                                   significance_threshold,
                                   plot_title = NULL,
                                   plot_subtitle = NULL,
                                   trt_name = NULL,
                                   ctrl_name = NULL,
                                   base_name = NULL) {
  output_format <- match.arg(output_format, choices = c("png", "pdf", "jpg", "jpeg"))
  if (output_format == "jpeg") {
    output_format <- "jpg"
  }
  file_extension <- switch(output_format, png = "png", pdf = "pdf", jpg = "jpg")

  # Set figure height and width
  label_vec <- as.character(meta[[annotation_column]])
  n_anno <- length(unique(label_vec))
  max_label_nchar <- max(nchar(unique(label_vec)))
  num_conditions <- length(unique(meta[[group_column]]))
  dynamic_height <- 2200 + (num_conditions * 50)

  # Make box plot with customization
  set.seed(42)
  pbox = sccomp_boxplot(sccomp_res, factor = "sccomp", remove_unwanted_effects = TRUE, significance_threshold=significance_threshold)
  # Renames specific groups
  if (!is.null(trt_name) && !is.null(ctrl_name) && !is.null(base_name)) {
    pbox <- pbox & scale_x_discrete(labels = function(x) {
      custom_labels <- setNames(
        c(paste0(trt_name, "\n(Treatment)"), paste0(ctrl_name, "\n(Control)"), paste0(base_name, "\n(Baseline)")),
        c(trt_name, ctrl_name, base_name)
      )
      ifelse(x %in% names(custom_labels), custom_labels[x], x)
    })
  }

  fixed_boxplot <- pbox &
    theme(
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5)
    )&
    guides(fill = guide_legend(ncol = 3, byrow = TRUE))

  anno_args <- list(theme = theme(plot.margin = margin(t = 5, r = 10, b = 10, l = 10, unit = "pt")))
  if (!is.null(plot_title)) anno_args$title <- plot_title
  if (!is.null(plot_subtitle)) anno_args$subtitle <- plot_subtitle

  fixed_boxplot <- fixed_boxplot + do.call(patchwork::plot_annotation, anno_args)

  output_path <- paste0(output_dir, "/plots/scComp_boxplot.", output_suffix, file_extension)
  if (output_format %in% c("png", "jpg")) {
    ggsave(
      filename = output_path,
      plot = fixed_boxplot,
      width = 2600,
      height = dynamic_height,
      units = "px",
      dpi = 350
    )
  } else if (output_format == "pdf") {
    ggsave(
      filename = output_path,
      plot = fixed_boxplot,
      width = 10,
      height = 8,
      units = "in"
    )
  }

  # Save 1D interval plot
  num_params <- sccomp_res %>% 
    filter(str_starts(parameter, "sccomp")) %>% 
    pull(parameter) %>% 
    unique() %>% 
    length()

  base_height_px <- max(1200, 300 + (60 * n_anno))
  base_height_in <- max(6, 1.5 + (0.3 * n_anno))

  img_height_px <- base_height_px * num_params + 200
  img_height_in <- base_height_in * num_params
  img_width_px  <- max(2500, 900 + 20 * max_label_nchar)

  output_path <- paste0(output_dir, "/plots/scComp_1D_interval.", output_suffix, file_extension)
  p1d = plot_1D_intervals(sccomp_res, show_fdr_message=FALSE, significance_threshold=significance_threshold) &
    scale_color_manual(
      values = c("FALSE" = "grey50", "TRUE" = "red"),
      limits = c("FALSE", "TRUE"),
      na.value = "grey50"
    ) &
    theme(
      plot.margin = margin(t = 5, r = 25, b = 10, l = 5, unit = "pt"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8, face = "plain", hjust = 0)
    )

  anno_args_1d <- list()
  if (!is.null(plot_title)) anno_args_1d$title <- plot_title
  if (!is.null(plot_subtitle)) anno_args_1d$subtitle <- plot_subtitle

  if (length(anno_args_1d) > 0) {
    p1d <- p1d + do.call(patchwork::plot_annotation, anno_args_1d)
  }
  if (output_format %in% c("png", "jpg")) {
    ggsave(
      filename = output_path,
      plot = p1d,
      width = img_width_px,
      height = img_height_px,
      units = "px",
      dpi = 350
    )
  } else if (output_format == "pdf") {
    ggsave(
      filename = output_path,
      plot = p1d,
      width = 10,
      height = img_height_in,
      units = "in"
    )
  }

  # Save regression and interval plots separately. Fix label misalignment issue.
  # Save scatter plot with regressions only for main tests (requires Intercept)
  if ("(Intercept)" %in% unique(sccomp_res$parameter)) {
    p2d <- plot_2D_intervals(sccomp_res, show_fdr_message=FALSE, significance_threshold=significance_threshold)

    # Keep only the Intercept panel
    keep_params <- "(Intercept)" 

    p2d$data <- p2d$data %>%
      filter(parameter %in% keep_params) %>%
      mutate(parameter = factor(parameter, levels = keep_params))

    for (i in seq_along(p2d$layers)) {
      # Check if the layer actually has data
      if (!is.null(p2d$layers[[i]]$data)) {
        # Check if that data contains the 'parameter' column
        if ("parameter" %in% colnames(p2d$layers[[i]]$data)) {
          # Filter it out and drop the levels
          p2d$layers[[i]]$data <- p2d$layers[[i]]$data %>%
            filter(parameter %in% keep_params) %>%
            mutate(parameter = factor(parameter, levels = keep_params))
        }
      }
    }

    p2d <- p2d + facet_wrap(~ parameter, scales = "free") + 
      labs(
        title = "Mean-Variability Association",
        subtitle = if (!is.null(base_name)) {
          sprintf("The (Intercept) represents the absolute abundance of the baseline: %s", base_name)
        } else {
          "The (Intercept) represents the absolute abundance of the baseline"
        }
      ) +
      theme(
        plot.title = element_text(size = 9, face = "bold", hjust=0),
        plot.subtitle = element_text(size = 6, color = "grey30", hjust=0)
      )

    output_path <- paste0(output_dir, "/plots/scComp_2D_scatter.", output_suffix, file_extension)

    if (output_format %in% c("png", "jpg")) {
      ggsave(
        filename = output_path,
        plot = p2d,
        width = 1250,
        height = 1250,
        units = "px",
        dpi = 350
      )
    } else if (output_format == "pdf") {
      ggsave(
        filename = output_path,
        plot = p2d,
        width = 6,
        height = 6,
        units = "in"
      )
    }
  }

  # Prepare and plot the 2D credible intervals.
  plot_data <- sccomp_res %>%
    filter(str_starts(parameter, "sccomp")) %>%
    # Define what significant means
    mutate(c_sig = c_FDR < significance_threshold,
           v_sig = v_FDR < significance_threshold,
           significant = c_FDR < significance_threshold | v_FDR < significance_threshold)

  num_panels <- length(unique(plot_data$parameter))
  img_width_px <- 500 + (1250 * num_panels)  # e.g., 1 panel = 1750px, 2 panels = 3000px
  img_width_in <- 2 + (5 * num_panels)

  p2di = ggplot(plot_data, aes(x = c_effect, y = v_effect)) +
    geom_errorbar(aes(ymin = v_lower, ymax = v_upper, color = v_sig),
                  linewidth = 0.4,
                  alpha = 0.7,
                  orientation = 'x') +
    geom_errorbar(aes(xmin = c_lower, xmax = c_upper, color = c_sig),
                  linewidth = 0.4,
                  alpha = 0.7,
                  orientation = 'y') +
    geom_point(color = "black", size = 0.8, alpha = 1) +
    # Add the dashed boundaries
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey70") +
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey70") +
    # Add the labels using ggrepel
    geom_text_repel(
      data = subset(plot_data, significant==TRUE),
      aes(label = .data[[annotation_column]]),
      size = 2.5,
      max.overlaps = Inf,
      box.padding = 0.2,
      point.padding = 0.2,
      segment.color = "grey50"
    ) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
    theme_bw() +
    facet_wrap(~ parameter, scales = "free", nrow = 1) +
    labs(
      title = if(!is.null(plot_title)) plot_title else "scComp 2D credible intervals",
      subtitle = if(!is.null(plot_subtitle)) plot_subtitle else NULL,
      x = "c_effect (Abundance effect)",
      y = "v_effect (Variability effect)",
      color = paste0("FDR < ", significance_threshold)
    ) +
    theme(
      axis.title = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 6)
    )
  # Save the plot
  output_path <- paste0(output_dir, "/plots/scComp_2D_interval.", output_suffix, file_extension)
  if (output_format %in% c("png", "jpg")) {
    ggsave(
      filename = output_path,
      plot = p2di,
      width = img_width_px,
      height = 1500,
      units = "px",
      dpi = 350
    )
  } else if (output_format == "pdf") {
    ggsave(
      filename = output_path,
      plot = p2di,
      width = img_width_in,
      height = 7.5,
      units = "in"
    )
  }
}

#' Generate table of stat results scComp
#' @param sccomp_res scComp test result object
#' @param output_dir directory to save the results table
#' @return NULL saves comparison figures in the specified directory.
#' @noRd
#'
stat_res_sccomp <- function(sccomp_res, output_dir=".", output_suffix=""){
  res_tab <- sccomp_res %>% as.data.frame()
  output_path <- paste0(output_dir, "/tables/scComp_result.", output_suffix, "txt")
  write.table(res_tab, output_path, sep='\t', row.names = FALSE)
}
