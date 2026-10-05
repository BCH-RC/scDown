# Utility functions for Seurat DEG and Custom clusterProfiler Pathway Analysis

#' Get Species Configuration
#'
#' Retrieves default configurations (OrgDb, organism ID, and regex patterns) 
#' based on the specified species.
#'
#' @param species Character string specifying the species ("human" or "mouse").
#' @return A list containing configuration variables specific to the species.
#' @noRd
get_species_config <- function(species) {
  species_lower <- tolower(species)
  
  if (species_lower == "mouse") {
    return(list(
      OrgDb = "org.Mm.eg.db",
      organism = "mmu",
      mt_pattern = "^mt-",
      mt_threshold = 0.05,
      keggReplaceString = " - Mus musculus \\(house mouse\\)"
    ))
  } else if (species_lower == "human") {
    return(list(
      OrgDb = "org.Hs.eg.db",
      organism = "hsa",
      mt_pattern = "^MT-",
      mt_threshold = 0.10,
      keggReplaceString = " - Homo sapiens \\(human\\)"
    ))
  } else {
    stop("Species '", species, "' is not supported. Please configure it in get_species_config().")
  }
}

#' Generate Pairwise Comparisons
#'
#' Creates all unique pairwise permutations of conditions for differential expression analysis.
#'
#' @param conditions_vector A character vector of unique conditions/groups.
#' @return A list of character vectors, where each vector contains two conditions to compare.
#' @noRd
pairwise <- function(conditions_vector) {
  comparisons <- list()
  # Ensure unique pairs without duplicates (e.g., gets A_B but skips B_A)
  for (i in 1:(length(conditions_vector) - 1)) {
    for (j in (i + 1):length(conditions_vector)) {
      cond1 <- conditions_vector[i]
      cond2 <- conditions_vector[j]
      pair_id <- paste(cond1, cond2, sep="_")
      comparisons[[pair_id]] <- c(cond1, cond2)
    }
  }
  return(comparisons)
}

#' Prepare and Cache KEGG Database
#'
#' Downloads the latest KEGG pathway mappings from the API or loads them from 
#' a local RDS cache to prevent redundant API calls during parallel processing.
#'
#' @param species Character string specifying the organism code.
#' @param KEGG_Type Character string specifying the database type (default: "KEGG").
#' @param keyType Character string specifying the gene ID key type (default: "kegg").
#' @param workdir Path to the directory where the cache should be stored.
#' @return KEGG data.
#' @noRd
prepare_KEGG <- function(species, KEGG_Type="KEGG", keyType="kegg", workdir=".") {
  kegg_dir <- file.path(workdir, "kegg_database")
  if(!dir.exists(kegg_dir)) dir.create(kegg_dir, recursive = TRUE)
  
  kegg_db_file <- file.path(kegg_dir, paste0("kegg_",species,".RDS")) 
  if(file.exists(kegg_db_file)) {
    kegg <- readRDS(kegg_db_file)
  } else {
    kegg <- download_KEGG(species, KEGG_Type, keyType)
    saveRDS(kegg, kegg_db_file)
  }
  return(kegg)
}

#' Custom KEGG Enrichment Analysis
#'
#' Runs hypergeometric testing for KEGG pathway enrichment using the locally cached database.
#'
#' @param gene A vector of Entrez gene IDs.
#' @param organism Character string for the KEGG organism code (e.g., "hsa").
#' @param keyType Target key type.
#' @param pvalueCutoff Adjusted p-value cutoff.
#' @param pAdjustMethod Method for multiple testing correction.
#' @param universe Background gene universe.
#' @param minGSSize Minimum gene set size.
#' @param maxGSSize Maximum gene set size.
#' @param qvalueCutoff Q-value cutoff.
#' @param workdir Directory path for loading the database cache.
#' @return An `enrichResult` instance.
#' @noRd
enrichKEGG_custom <- function(gene, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = NULL, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, workdir=".") {
  # load the pre-cached KEGG data
  kegg <- prepare_KEGG(organism, "KEGG", keyType, workdir)
  term2gene <- as.data.frame(kegg$KEGGPATHID2EXTID)
  term2name <- as.data.frame(kegg$KEGGPATHID2NAME)
  
  # run enricher function directly
  res <- clusterProfiler::enricher(
    gene = gene, 
    pvalueCutoff = pvalueCutoff, 
    pAdjustMethod = pAdjustMethod, 
    universe = universe, 
    minGSSize = minGSSize, 
    maxGSSize = maxGSSize, 
    qvalueCutoff = qvalueCutoff, 
    TERM2GENE = term2gene,
    TERM2NAME = term2name
  )
  
  if (is.null(res)) return(res)
  
  # assign KEGG metadata
  res@ontology <- "KEGG"
  res@organism <- organism
  res@keytype <- keyType
  return(res)
}

#' Open Graphic Device
#'
#' Opens a graphics device (PNG, JPEG, or PDF) with calculated dimensions.
#'
#' @param filename Full path to the output file.
#' @param w_px Target width in pixels (or equivalent for PDF).
#' @param h_px Target height in pixels (or equivalent for PDF).
#' @param fmt Output format ("png", "jpeg", "jpg", or "pdf").
#' @return None. Opens a graphic device.
#' @noRd
open_plot_device <- function(filename, w_px, h_px, fmt) {
  if (fmt == "pdf") {
    pdf(filename, width = w_px/300, height = h_px/300)
  } else if (fmt %in% c("jpg", "jpeg")) {
    jpeg(filename, width = w_px, height = h_px, res = 300)
  } else {
    png(filename, width = w_px, height = h_px, res = 300)
  }
}

#' Generate Standard Pathway Visualizations
#'
#' Takes an enrichResult object and generates a barplot, dotplot, cnetplot, 
#' and treeplot, saving them directly to disk.
#'
#' @param enrich_obj The enrichResult object.
#' @param db_name String for the database (e.g., "KEGG", "GO_BP").
#' @param direction String for regulation direction ("ALL", "UP", "DOWN").
#' @param comparison_name String for the comparison.
#' @param n_top_terms Number of terms to show.
#' @param output_format Image format ("png", "pdf", etc.).
#' @param output_dir Base directory for output.
#' @param foldChange_vector Named numeric vector of log2 fold changes.
#' @param keggReplaceString Optional regex string to clean KEGG descriptions.
#' @noRd
generate_pathway_plots <- function(enrich_obj, db_name, direction, comparison_name, n_top_terms, output_format, output_dir, foldChange_vector, keggReplaceString = NULL) {
  
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) return(NULL)
  
  if (db_name == "KEGG" && !is.null(keggReplaceString)) {
    enrich_obj@result$Description <- gsub(keggReplaceString, "", enrich_obj@result$Description)
  }
  
  # use the same top N terms across all plots
  n_terms_plot <- min(n_top_terms, nrow(enrich_obj@result))
  enrich_obj@result <- enrich_obj@result[1:n_terms_plot, ]
  
  # truncate descriptions to 50 characters
  enrich_obj@result$Description <- stringr::str_trunc(enrich_obj@result$Description, 50)
  
  # use information from the longest text name result
  max_desc <- max(sapply(enrich_obj$Description, nchar))
  height_calc <- (60 + (n_terms_plot) * 10)
  width_calc <- ifelse(max_desc <= 50, (300 + (max_desc) * 5), 550)
  width_calc <- max(width_calc, height_calc)
  
  plot_base_dir <- file.path(output_dir, "clusterprofiler", "plots", db_name)
  if(!dir.exists(plot_base_dir)) dir.create(plot_base_dir, recursive = TRUE)
  
  base_title <- paste(db_name, direction, "-", comparison_name)
  base_filename <- paste0(db_name, "_", direction, "_top_", n_top_terms, "_", comparison_name)
  
  # define basic theme
  custom_theme <- theme(
    plot.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm") 
  )
  
  # Barplot (Ranked by p.adjust)
  open_plot_device(file.path(plot_base_dir, paste0(base_filename, "_barplot.", output_format)), width_calc * 4, height_calc * 4, output_format)
  p1 <- suppressWarnings(suppressMessages(
    barplot(enrich_obj, 
            x = "GeneRatio", 
            showCategory = n_top_terms, 
            title = paste(base_title, "(Barplot - Ranked by p.adjust)"), 
            font.size = 8,
            label_format = 100) +
      scale_x_continuous(expand = expansion(mult = c(0, .1))) + 
      custom_theme
  ))
  print(p1)
  dev.off()
  
  # Dotplot (Ranked by GeneRatio)
  open_plot_device(file.path(plot_base_dir, paste0(base_filename, "_dotplot.", output_format)), width_calc * 4, height_calc * 4, output_format)
  p2 <- suppressWarnings(suppressMessages(
    dotplot(enrich_obj, 
            showCategory = n_top_terms, 
            title = paste(base_title, "(Dotplot - Ranked by GeneRatio)"), 
            font.size = 8, 
            label_format = 100) +
      custom_theme
  ))
  print(p2)
  dev.off()
}

#' Run Functional Enrichment Analysis and Plotting
#'
#' Executes KEGG and GO (BP, CC, MF) pathway enrichment on ALL, UP, and DOWN subsets of DEGs.
#'
#' @param DEG Data frame containing differential expression results from Seurat.
#' @param comparison_name String containing the active comparison identifier.
#' @param organism KEGG organism code.
#' @param OrgDb Bioconductor OrgDb package name.
#' @param keggReplaceString String (regex) to strip redundant species tags from KEGG descriptions.
#' @param n_top_terms Number of top pathways to visualize in the plots.
#' @param significance_threshold P-value cutoff for filtering DEGs and enriched pathways.
#' @param output_format Format for the output plots ("png", "pdf", etc.).
#' @param output_dir Base directory path to save tables and plots.
#' @return NULL. Outputs are saved directly to disk.
#' @noRd
run_functional_analysis <- function(DEG, comparison_name, organism, OrgDb, keggReplaceString, n_top_terms, significance_threshold, output_format, output_dir) {
  
  # fold-change vector named by gene symbol for cnetplot
  foldChange_vector <- DEG$avg_log2FC
  names(foldChange_vector) <- rownames(DEG)
  
  # define the three sets of genes
  gene_sets <- list(
    "ALL"  = rownames(subset(DEG, p_val_adj < significance_threshold)),
    "UP"   = rownames(subset(DEG, p_val_adj < significance_threshold & avg_log2FC > 0)),
    "DOWN" = rownames(subset(DEG, p_val_adj < significance_threshold & avg_log2FC < 0))
  )
  
  for (direction in names(gene_sets)) {
    symbol <- gene_sets[[direction]]
    
    if(length(symbol) <= 4) {
      message(paste("Skipping", direction, "for", comparison_name, "- not enough significant genes."))
      next
    }
    
    entrezID = suppressWarnings(suppressMessages(clusterProfiler::bitr(symbol, fromType="SYMBOL", toType ="ENTREZID", OrgDb=OrgDb)$ENTREZID))
    
    # run analyses (enrichGO output gene symbols by default)
    kegg_result <- enrichKEGG_custom(gene = entrezID, organism = organism, pvalueCutoff = significance_threshold, workdir = file.path(output_dir, "clusterprofiler"))
    go_cc_all <- enrichGO(gene = entrezID, OrgDb = OrgDb, ont = "CC", pvalueCutoff = significance_threshold, readable = TRUE)
    go_bp_all <- enrichGO(gene = entrezID, OrgDb = OrgDb, ont = "BP", pvalueCutoff = significance_threshold, readable = TRUE)
    go_mf_all <- enrichGO(gene = entrezID, OrgDb = OrgDb, ont = "MF", pvalueCutoff = significance_threshold, readable = TRUE)
    
    # translate Entrez IDs to gene symbols
    if(!is.null(kegg_result)) {
      kegg_result <- suppressMessages(clusterProfiler::setReadable(kegg_result, OrgDb = OrgDb, keyType = "ENTREZID"))
    }
    
    # write tables to disk
    keggFile <- file.path(output_dir, "clusterprofiler", "tables", "KEGG", paste0("KEGG_", direction, "_", comparison_name, ".csv"))
    ccFile <- file.path(output_dir, "clusterprofiler", "tables", "GO_CC", paste0("GO_CC_", direction, "_", comparison_name, ".csv"))
    bpFile <- file.path(output_dir, "clusterprofiler", "tables", "GO_BP", paste0("GO_BP_", direction, "_", comparison_name, ".csv"))
    mfFile <- file.path(output_dir, "clusterprofiler", "tables", "GO_MF", paste0("GO_MF_", direction, "_", comparison_name, ".csv"))
    
    if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
      write.csv(as.data.frame(kegg_result), file=keggFile, row.names = FALSE)
    }
    if (!is.null(go_cc_all) && nrow(as.data.frame(go_cc_all)) > 0) {
      write.csv(as.data.frame(go_cc_all), file=ccFile, row.names = FALSE) 
    }
    if (!is.null(go_bp_all) && nrow(as.data.frame(go_bp_all)) > 0) {
      write.csv(as.data.frame(go_bp_all), file=bpFile, row.names = FALSE) 
    }
    if (!is.null(go_mf_all) && nrow(as.data.frame(go_mf_all)) > 0) {
      write.csv(as.data.frame(go_mf_all), file=mfFile, row.names = FALSE) 
    }
    
    # generate and save all plot types
    generate_pathway_plots(kegg_result, "KEGG", direction, comparison_name, n_top_terms, output_format, output_dir, foldChange_vector, keggReplaceString)
    generate_pathway_plots(go_cc_all, "GO_CC", direction, comparison_name, n_top_terms, output_format, output_dir, foldChange_vector)
    generate_pathway_plots(go_bp_all, "GO_BP", direction, comparison_name, n_top_terms, output_format, output_dir, foldChange_vector)
    generate_pathway_plots(go_mf_all, "GO_MF", direction, comparison_name, n_top_terms, output_format, output_dir, foldChange_vector)
  }
  
  return(NULL)
}
