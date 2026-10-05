#' Run clusterProfiler
#'
#' This function performs differential gene expression (DEG) and pathway enrichment analysis 
#' (KEGG, GO_BP, GO_CC, GO_MF) on a Seurat object. It generates volcano plots, pathway barplots, 
#' and comprehensive CSV tables as output. The function allows for the specification of various 
#' parameters, including the directory for storing results, the metadata column for cell annotations, 
#' and specific cell types of interest. The function also supports adding covariates for MAST-based DEG testing, 
#' specifying the number of top pathways to visualize, and utilizing parallel processing.
#'
#' @param seurat_obj Seurat object containing the single-cell data and metadata.
#' @param species The species of the data, either 'human' or 'mouse', required.
#' @param annotation_column Name of the metadata column in the Seurat object that contains cell annotations or cluster IDs.
#' @param group_column Name of the metadata column in the Seurat object that defines conditions or groups (e.g., Treatment vs Control).
#' @param control_cond String representing the baseline/control condition. If provided, all other groups are compared against this. If NULL, all unique pairwise comparisons are performed.
#' @param selected_cell_types A character vector of specific cell types to analyze. If not provided, all cell types in the annotation column will be considered.
#' @param output_dir Path to the folder where clusterProfiler results will be saved (figures and tables).
#' @param covariates Optional character vector of metadata columns to use as latent variables for Seurat FindMarkers (switches to MAST test).
#' @param n_jobs A numeric variable to set the number of cores to be used for parallel processing.
#' @param n_top_terms The number of top signaling pathways/terms to analyze and visualize in the barplots.
#' @param significance_threshold The adjusted p-value threshold used for defining significant DEGs and enriched pathways (default: 0.05).
#' @param output_format Format of output figures: "png", "pdf", "jpg", or "jpeg" (default: "png").
#'
#' @return NULL
#'
#' @export
#'
run_clusterProfiler <- function(seurat_obj, 
                                species, 
                                annotation_column=NULL, 
                                group_column=NULL, 
                                control_cond=NULL, 
                                selected_cell_types=NULL, 
                                output_dir=".", 
                                covariates=NULL, 
                                n_jobs=2, 
                                n_top_terms=20,
                                significance_threshold=0.05, 
                                output_format="png") {
  library(Seurat)
  library(clusterProfiler)
  library(enrichplot)
  library(EnhancedVolcano)
  library(ggplot2)
  library(stringr)
  library(foreach)
  library(doParallel)
  library(checkmate)
  
  # Set global seed for reproducibility
  set.seed(42)
  
  # Check the input data format
  checkmate::expect_class(seurat_obj, "Seurat", label="seurat_obj")
  checkmate::expect_choice(annotation_column, colnames(seurat_obj@meta.data), label="annotation_column", null.ok = TRUE)
  checkmate::expect_choice(group_column, colnames(seurat_obj@meta.data), label="group_column", null.ok = TRUE)
  checkmate::expect_directory(output_dir, access="rw", label = "output_dir")
  checkmate::expect_choice(output_format, c("png","pdf","jpg","jpeg"), label = "output_format")
  
  meta <- seurat_obj@meta.data
  
  if(!annotation_column %in% colnames(meta)){
    stop("Cluster column does not exist in the seurat object meta data!")
  }
  if(!group_column %in% colnames(meta)){
    stop("Group (condition) column does not exist in the seurat object meta data!")
  }
  
  # Pre-create all analysis subdirectories for tables and plots
  analysis_subdirs <- c("DEGs", "KEGG", "GO_BP", "GO_CC", "GO_MF")
  for (subdir in analysis_subdirs) {
    dir.create(file.path(output_dir, "clusterprofiler", "tables", subdir), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, "clusterprofiler", "plots", subdir), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Retrieve configurations based on species
  config <- get_species_config(species)
  OrgDb <- config$OrgDb
  organism <- config$organism
  keggReplaceString <- config$keggReplaceString
  
  # Ensure the necessary OrgDb is loaded dynamically
  if (!requireNamespace(OrgDb, quietly = TRUE)) {
    stop(paste0("The package '", OrgDb, "' is required but not installed. Please install it via BiocManager."))
  }
  library(OrgDb, character.only = TRUE)
  
  # Extract unique clusters and conditions
  all_clusters <- as.character(unique(seurat_obj[[annotation_column]][,1]))
  conditions_vector <- as.character(unique(seurat_obj[[group_column]][,1]))
  
  # Filter by selected cell types if provided
  if (!is.null(selected_cell_types)) {
    clusters <- intersect(all_clusters, selected_cell_types)
    if (length(clusters) == 0) {
      stop("None of the selected_cell_types were found in the annotation_column.")
    }
  } else {
    clusters <- all_clusters
  }
  
  if (length(conditions_vector) < 2) {
    stop("There are not enough conditions to perform a differential gene expression analysis.")
  }
  
  # Create a combined identity in the Seurat object (e.g., "Cluster1_ConditionA")
  seurat_obj$celltype.group <- paste(seurat_obj[[annotation_column]][,1], seurat_obj[[group_column]][,1], sep = "_")
  Idents(seurat_obj) <- "celltype.group"
  
  # Determine comparisons (Control-based vs Pairwise)
  if (!is.null(control_cond)) {
    if (!control_cond %in% conditions_vector) {
      stop(paste("control_cond", control_cond, "not found in the group_column."))
    }
    
    comparisons_list <- list()
    for (cond in conditions_vector) {
      if (cond != control_cond) {
        pair_id <- paste(cond, control_cond, sep="_")
        comparisons_list[[pair_id]] <- c(cond, control_cond)
      }
    }
    pairwise_comparisons <- comparisons_list
  } else {
    pairwise_comparisons <- pairwise(conditions_vector)
  }
  
  # Pre-cache the KEGG database
  message("Checking KEGG database...")
  prepare_KEGG(organism, "KEGG", "kegg", file.path(output_dir, "clusterprofiler"))
  
  # Register cores for parallel computing
  cl <- makeCluster(n_jobs, outfile="")
  registerDoParallel(cl)
  
  on.exit(stopCluster(cl), add = TRUE)
  
  # Loop through each cluster and perform DEG and pathway analysis
  foreach(cluster_name = clusters, 
    .packages = c("Seurat", "clusterProfiler", "enrichplot", "EnhancedVolcano", "stringr", "ggplot2", OrgDb), 
    .export = c("run_functional_analysis", "generate_pathway_plots", "enrichKEGG_custom", "prepare_KEGG", "pairwise", "open_plot_device")) %dopar% {
    
    for (pairwise_comparison in pairwise_comparisons) {
      cond1 <- pairwise_comparison[1]
      cond2 <- pairwise_comparison[2]
      
      ident1 <- paste(cluster_name, cond1, sep = "_")
      ident2 <- paste(cluster_name, cond2, sep = "_")
      
      # Check if both idents exist in the seurat object to prevent errors
      if(ident1 %in% unique(Idents(seurat_obj)) && ident2 %in% unique(Idents(seurat_obj))) {
        
        # Check cell counts (need at least 3 cells per condition)
        count1 <- sum(Idents(seurat_obj) == ident1)
        count2 <- sum(Idents(seurat_obj) == ident2)
        
        if (count1 >= 3 && count2 >= 3) {
          
          print(paste("======== PROCESSING:", cluster_name, "|", cond1, "vs", cond2, "========"))
          
          # Create a non-redundant comparison name (e.g., "Astrocytes_Treatment_vs_Control")
          comparison_name <- paste0(cluster_name, "_", cond1, "_vs_", cond2)
          
          # Compute DEGs
          if (is.null(covariates)) {
            DEG <- Seurat::FindMarkers(seurat_obj, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0.1, min.pct = 0.01, random.seed = 42)
          } else {
            DEG <- Seurat::FindMarkers(seurat_obj, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0.1, min.pct = 0.01, test.use = "MAST", latent.vars = covariates, random.seed = 42)
          }
          
          # Write DEGs to tables
          write.csv(DEG, file = file.path(output_dir, "clusterprofiler", "tables", "DEGs", paste0("DEG_ALL_", comparison_name, ".csv")))
          DEG_sig <- DEG[which(DEG$p_val_adj < significance_threshold),]
          write.csv(DEG_sig, file = file.path(output_dir, "clusterprofiler", "tables", "DEGs", paste0("DEG_SIG_", comparison_name, ".csv")))
          
          # Run functional analysis 
          run_functional_analysis(DEG_sig, comparison_name, organism, OrgDb, keggReplaceString, n_top_terms, significance_threshold, output_format, output_dir)
          
          # Volcano plot
          volcano_file <- file.path(output_dir, "clusterprofiler", "plots", "DEGs", paste0("Volcano_", comparison_name, ".", output_format))
          open_plot_device(volcano_file, 2600, 2200, output_format)
          volcano_plot <- suppressWarnings(EnhancedVolcano(DEG,
                                          lab = rownames(DEG),
                                          x = 'avg_log2FC',
                                          y = 'p_val',
                                          title = comparison_name,
                                          pCutoff = significance_threshold,
                                          FCcutoff = log2(1.5)))
          print(volcano_plot)
          dev.off()
          
        } else {
          message(paste0("Skipping ", ident1, " vs ", ident2, ". Not enough cells (requires >= 3 in both)."))
        }
      }
    }
  }
  print(paste("clusterProfiler pathway analysis completed. Results saved in:", file.path(output_dir, "clusterprofiler")))
}
