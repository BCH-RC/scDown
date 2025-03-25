#' Run CellChat
#'
#' run CellChat analysis
#' This function performs CellChat analysis on a Seurat object to investigate cell-to-cell communication. 
#' It generates figures, tables, and RDS files as output. The function allows for the specification of 
#' various parameters, including the directory for storing results, the metadata column for cell annotations, 
#' and the cell annotations of interest. Additionally, users can provide species information, group metadata 
#' for group comparisons, and define pairwise group for differential analysis. The function also supports 
#' specifying the number of top pathways to analyze and visualize.
#'
#' @param output_dir Path to the folder where CellChat results will be saved (figures, tables, and RDS files).
#' @param seurat_obj Seurat object containing the UMI counts and metadata for the single-cell data.
#' @param sample_column Name of the metadata column in the Seurat object that contains sample information.
#' @param annotation_column Name of the metadata column in the Seurat object that contains cell annotations.
#' @param annotation_selected A vector of cell annotations of interest for running the CellChat analysis. If not provided, all cell types will be considered.
#' @param species The species of the data, either 'human' or 'mouse'.
#' @param group_column Name of the metadata column in the Seurat object that defines conditions or groups.
#' @param group_cmp A list of pairwise condition comparisons for differential CellChat analysis (e.g., comparing different groups or conditions).
#' @param top_n The number of top signaling pathways to analyze and visualize.
#' 
#' @return NULL

run_cellchatV2 <- function(output_dir, seurat_obj, sample_column = NULL, annotation_column, annotation_selected = NULL, species, group_column = NULL, group_cmp = NULL, top_n = 10, num_cores = 4) {

  library(Seurat)
  library(CellChat)
  library(ComplexHeatmap)
  library(checkmate)
  library(foreach)
  library(doParallel)

  # Checking the cellchat inputs
  check_required_variables(seurat_obj, species, output_dir, annotation_column, group_column)
  # Check sample_column
  checkmate::expect_choice(sample_column, colnames(seurat_obj@meta.data), label = "sample_column", null.ok = TRUE)
  # Check annotation_selected
  checkmate::expect_subset(annotation_selected, as.vector(unique(seurat_obj@meta.data[, annotation_column])), label = "annotation_selected")
  # Check group_cmp
  checkmate::expect_subset(unique(unlist(group_cmp)), as.vector(unique(seurat_obj@meta.data[, group_column])), label = "group_cmp")

  # Create cell chat folders
  create_dir_cellchat(dir_cellchat = output_dir)

  # Set the cell annotations as the seurat object Idents
  Idents(seurat_obj) <- seurat_obj@meta.data[, annotation_column]

  # Create a sample column
  if (!is.null(sample_column)) {
    seurat_obj$samples <- as.factor(seurat_obj@meta.data[, sample_column])
  } else {
    seurat_obj$samples <- as.factor("sample")
  }
  
  # Run CellChat analysis by conditions
  if (is.null(group_column)) {
    # No condition or group information are specified
    cellchat_object <- doCellCom(X = seurat_obj, species = species)
    
    # Subset the cell types if there's any
    if (is.null(annotation_selected)) {
      # No cell type of interest specified
      seurat_obj_condition <- seurat_obj
    } else {
      cellchat_object <- subsetCellChatMod(cellchat_object, idents.use = annotation_selected)
      cellchat_object <- netAnalysis_computeCentrality(cellchat_object)
      seurat_obj_condition <- seurat_obj[, seurat_obj@meta.data[, annotation_column] %in% annotation_selected]
    }
    
    saveRDS(cellchat_object, file = paste0(output_dir, "/cellchat/rds/cellchat_obj_ALL.rds"))
    saveRDS(seurat_obj_condition, file = paste0(output_dir, "/cellchat/rds/seurat_obj_ALL.rds"))
    
  } else {
    conditions <- unique(as.vector(seurat_obj@meta.data[, group_column]))
    metadata_cond <- FetchData(object = seurat_obj, vars = group_column)

    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    # for (condition in conditions) {
    foreach(condition = conditions, .packages = c("Seurat", "CellChat"), .export = c("doCellCom", "subsetCellChatMod", "netAnalysis_computeCentrality")) %dopar% {
      seurat_obj_condition <- seurat_obj[, which(x = (metadata_cond == condition))]
      cellchat_object_condition <- doCellCom(seurat_obj_condition, species)
      seurat_obj_condition <- seurat_obj[, which(metadata_cond == condition)]
      
      # Subset the cell types if there's any
      if (!is.null(annotation_selected)) {
        cellchat_object_condition <- subsetCellChatMod(cellchat_object_condition, idents.use = annotation_selected)
        cellchat_object_condition <- netAnalysis_computeCentrality(cellchat_object_condition)
        seurat_obj_condition <- seurat_obj_condition[, seurat_obj_condition@meta.data[, annotation_column] %in% annotation_selected]
      }
      
      saveRDS(cellchat_object_condition, file = paste0(output_dir, "/cellchat/rds/cellchat_obj_", condition, ".rds"))
      saveRDS(seurat_obj_condition, file = paste0(output_dir, "/cellchat/rds/seurat_obj_", condition, ".rds"))
    }
    # Stop the parallel cluster
    stopCluster(cl)

  }
  
  # CellChat visualization at the aggregated level
  if (is.null(group_column)) {
    condition <- "ALL"
    cellchat_object <- readRDS(paste0(output_dir, "/cellchat/rds/cellchat_obj_", condition, ".rds"))
    aggregate_visu(X = cellchat_object, condition = condition, dir_cellchat = output_dir)
  } else {
    for (condition in conditions) {
      cellchat_object <- readRDS(paste0(output_dir, "/cellchat/rds/cellchat_obj_", condition, ".rds"))
      aggregate_visu(X = cellchat_object, condition = condition, dir_cellchat = output_dir)
    }
  }
  
  # CellChat visualization at signaling pathway level
  if (is.null(group_column)) {
    condition <- "ALL"
    # Get the top n pathways
    cellchat_object <- readRDS(paste0(output_dir, "/cellchat/rds/cellchat_obj_", condition, ".rds"))
    pathways_top <- cellchat_object@netP$pathways[1:top_n]
    pathways_top <- pathways_top[!is.na(pathways_top)]
    message(paste0("Top ", top_n, " Pathways for condition: ", condition, " (Not specified)"))
    message("If fewer than ", top_n, " pathways are displayed, only the available pathways have been found.")
    cat(pathways_top, sep = ";\n")
    
    if (length(pathways_top) == 0) {
      warning("No top pathways found for condition: ", condition)
    } else {
      seurat_obj <- readRDS(paste0(output_dir, "/cellchat/rds/seurat_obj_", condition, ".rds"))
      doCellComVisu(X = cellchat_object, Y = seurat_obj, pathways_to_show = pathways_top, condition = condition, dir_cellchat = output_dir, species = species)
    }
    
    # Save all inferred cell-cell communications at the level of: 1) ligands/receptors, 2) pathways
    # cellchat_object_net <- subsetCommunication(cellchat_object)
    cellchat_object_netp <- subsetCommunication(cellchat_object, slot.name = "netP")
    # write.csv(cellchat_object_net, file=paste0(output_dir, "/cellchat/csv/", condition, "_enriched_LRpairs.csv", sep=""), quote = F, row.names = F)
    write.csv(cellchat_object_netp, file=paste0(output_dir, "/cellchat/csv/", condition, "_enriched_pathways.csv", sep=""), quote = F, row.names = F)
    
    # Save all inferred cell-cell communications by pathways
    info_flow <- calc_infoflow(X = cellchat_object, condition)
    write.csv(info_flow, file=paste0(output_dir, "/cellchat/csv/", condition, "_enriched_pathways_aggregated.csv", sep=""), quote = F, row.names = F)
    
  } else {
    for (condition in conditions) {
      # Get the top n pathways
      cellchat_object <- readRDS(paste0(output_dir, "/cellchat/rds/cellchat_obj_", condition, ".rds"))
      pathways_top <- cellchat_object@netP$pathways[1:top_n]
      pathways_top <- pathways_top[!is.na(pathways_top)]
      message(paste0("Top ", top_n, " Pathways for condition: ", condition))
      message("If fewer than ", top_n, " pathways are displayed, only the available pathways have been found.")
      cat(pathways_top, sep = ";\n")
      
      if (length(pathways_top) == 0) {
        warning("No top pathways found for condition: ", condition)
      } else {
        seurat_obj <- readRDS(paste0(output_dir, "/cellchat/rds/seurat_obj_", condition, ".rds"))
        doCellComVisu(X = cellchat_object, Y = seurat_obj, pathways_to_show = pathways_top, condition = condition, dir_cellchat = output_dir, species = species)
      }
      
      # Save all inferred cell-cell communications at the level of: 1) ligands/receptors, 2) pathways
      # cellchat_object_net <- subsetCommunication(cellchat_object)
      cellchat_object_netp <- subsetCommunication(cellchat_object, slot.name = "netP")
      # write.csv(cellchat_object_net, file=paste0(output_dir, "/cellchat/csv/", condition, "_enriched_LRpairs.csv", sep=""), quote = F, row.names = F)
      write.csv(cellchat_object_netp, file=paste0(output_dir, "/cellchat/csv/", condition, "_enriched_pathways.csv", sep=""), quote = F, row.names = F)
      
      # Save all inferred cell-cell communications by pathways
      info_flow <- calc_infoflow(X = cellchat_object, condition)
      write.csv(info_flow, file=paste0(output_dir, "/cellchat/csv/", condition, "_enriched_pathways_aggregated.csv", sep=""), quote = F, row.names = F)
    }
  }
  
  # Perform pairwise CellChat comparisons
  if (!is.null(group_cmp)) {
    for (cmp in group_cmp) {
      # Conditions for pair wise comparison
      cond_1 <- cmp[1]
      cond_2 <- cmp[2]
      if (cond_1 == cond_2) {
        warning("Skip the condition comparison since two conditions are identical: ", cond_1)
        next()
      }
      
      # Load the seurat objects and cellchat objects
      # Condition 1
      seurat_obj_cond1 <- readRDS(paste0(output_dir, "/cellchat/rds/seurat_obj_", cond_1, ".rds"))
      cellchat_obj_cond1 <- readRDS(paste0(output_dir, "/cellchat/rds/cellchat_obj_", cond_1, ".rds"))
      # Condition 2
      seurat_obj_cond2 <- readRDS(paste0(output_dir, "/cellchat/rds/seurat_obj_", cond_2, ".rds"))
      cellchat_obj_cond2 <- readRDS(paste0(output_dir, "/cellchat/rds/cellchat_obj_", cond_2, ".rds"))
      
      run_cellchatV2_cmp(dir_cellchat = output_dir, seurat_obj_cond1 = seurat_obj_cond1, cellchat_obj_cond1 = cellchat_obj_cond1, seurat_obj_cond2 = seurat_obj_cond2, cellchat_obj_cond2 = cellchat_obj_cond2, condition_col = group_column, condition_1 = cond_1, condition_2 = cond_2, top_n = top_n)
    }
  }
}

#' Run CellChatV2 Pathway Visualization
#'
#' This function performs pathway visualization for the CellChatV2 analysis. 
#' It generates various visualizations, including signaling strength plots, heatmaps, scatter plots, 
#' LR gene expression violin plots, and communication probability bubble plots for specified pathways. 
#' Users can define the number of top pathways to visualize and specify pathways of interest that may not 
#' have been included in the top_n most significant pathways already visualized.
#'
#' @param output_dir Path to the folder where the visualization results will be saved (figures and tables).
#' @param species The species of the data, either 'human' or 'mouse'.
#' @param pathway_to_show A specific pathways to visualize.
#' 
#' @return NULL

cellchatV2_path_visu <- function(output_dir, species, pathway_to_show) {
  # Identify the conditions
  rds_files <- list.files(paste0(output_dir, "/cellchat/rds/"), pattern = "\\.rds$", full.names = TRUE)
  rds_files <- rds_files[which(regexpr("seurat_obj", rds_files)>0)]

  # Check if the folder is empty
  if (length(rds_files) == 0) {
    stop("Error: The folder is empty. No .rds files found.")
  }
  
  # Get the conditions
  conditions <- unique(sub(".*seurat_obj_(.*?)\\.rds", "\\1", rds_files))
  
  # Load the corresponding cellchat object and perform pathway analysis
  for (condition in conditions) {
    seurat_obj <- readRDS(paste0(output_dir, "/cellchat/rds/seurat_obj_", condition, ".rds"))
    cellchat_object <- readRDS(paste0(output_dir, "/cellchat/rds/cellchat_obj_", condition, ".rds"))
    
    cat("Pathway visualization for pathway:", pathway_to_show, "; condition:", condition, "\n")
    
    # Significant pathways
    pathways_sig <- cellchat_object@netP$pathways[!is.na(cellchat_object@netP$pathways)]
    
    if (pathway_to_show %in% pathways_sig) {
      pathway_visu(X = cellchat_object, Y = seurat_obj, pathway_to_show, condition, dir_cellchat = output_dir, species)
    } else {
      cat(pathway_to_show, "is not a significant pathway in the current CellChat Analysis for condition:", condition, "\n")
    }
  }
  
  cat(paste0("Pathway plots are updated in the folder: ", output_dir, "/cellchat/images/pathway"))
}


