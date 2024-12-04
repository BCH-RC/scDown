#' Run CellChat
#'
#' run CellChat analysis and generate figures, tables and rds files.
#'
#' @param dir_cellchat Path to the folder for storing CellChat results.
#' @param seurat_obj Seurat object containing UMI counts and metadata.
#' @param celltype_col Name of the metadata column used for cell type annotation in the Seurat object.
#' @param celltypes Cell types of interest for running CellChat analysis.
#' @param species Species for the data, either 'human' or 'mouse'.
#' @param condition_col Name of the metadata column for conditions or groups in the Seurat object.
#' @param conditions_cmp List containing the pairwise condition comparisons for CellChat analysis.
#' @param top_n the number of pathways.
#' 
#' @return NULL

run_cellchatV2 <- function(dir_cellchat, seurat_obj, celltype_col, celltypes = "ALL", species, condition_col = NULL, conditions_cmp = NULL, top_n = 10) {
  
  # Checking the cellchat inputs
  cellchat_input_check(dir_cellchat, seurat_obj, celltype_col, celltypes, species, condition_col)
  
  # Set the cell types as the seurat object Idents
  Idents(seurat_obj) <- seurat_obj@meta.data[, celltype_col]
  
  # Run CellChat analysis by conditions
  if (is.null(condition_col)) {
    # No condition or group information are specified
    cellchat_object <- doCellCom(X = seurat_obj, species = species)
    
    # Subset the cell types if there's any
    if (celltypes == "ALL") {
      # No cell type of interest specified
      seurat_obj_condition <- seurat_obj
    } else {
      celltypes_of_interest <- intersect(celltypes, unique(seurat_obj@meta.data[, celltype_col]))
      cellchat_object <- subsetCellChatMod(cellchat_object, idents.use = celltypes_of_interest)
      cellchat_object <- netAnalysis_computeCentrality(cellchat_object)
      seurat_obj_condition <- seurat_obj[, seurat_obj@meta.data[, celltype_col] %in% celltypes_of_interest]
    }
    
    saveRDS(cellchat_object, file = paste0(dir_cellchat, "cellchat/rds/cellchat_obj_ALL.rds"))
    saveRDS(seurat_obj_condition, file = paste0(dir_cellchat, "cellchat/rds/seurat_obj_ALL.rds"))
    
  } else {
    conditions <- unique(seurat_obj@meta.data[, condition_col])
    metadata_cond <- FetchData(object = seurat_obj, vars = condition_col)
    for (condition in conditions) {
      seurat_obj_condition <- seurat_obj[, which(x = (metadata_cond == condition))]
      cellchat_object_condition <- doCellCom(seurat_obj_condition, species)
      seurat_obj_condition <- seurat_obj[, which(metadata_cond == condition)]
      
      # Subset the cell types if there's any
      if (celltypes == "ALL") {
        
      } else {
        celltypes_of_interest <- intersect(celltypes, unique(seurat_obj_condition@meta.data[, celltype_col]))
        cellchat_object_condition <- subsetCellChatMod(cellchat_object_condition, idents.use = celltypes_of_interest)
        cellchat_object_condition <- netAnalysis_computeCentrality(cellchat_object_condition)
        seurat_obj_condition <- seurat_obj_condition[, seurat_obj_condition@meta.data[, celltype_col] %in% celltypes_of_interest]
      }
      
      saveRDS(cellchat_object_condition, file = paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", condition, ".rds"))
      saveRDS(seurat_obj_condition, file = paste0(dir_cellchat, "cellchat/rds/seurat_obj_", condition, ".rds"))
    }
  }
  
  # CellChat visualization at the aggregated level
  if (is.null(condition_col)) {
    condition <- "ALL"
    cellchat_object <- readRDS(paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", condition, ".rds"))
    aggregate_visu(X = cellchat_object, condition = condition)
  } else {
    for (condition in conditions) {
      cellchat_object <- readRDS(paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", condition, ".rds"))
      aggregate_visu(X = cellchat_object, condition = condition)
    }
  }
  
  # CellChat visualization at signaling pathway level
  if (is.null(condition_col)) {
    condition <- "ALL"
    # Get the top n pathways
    cellchat_object <- readRDS(paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", condition, ".rds"))
    pathways_top <- cellchat_object@netP$pathways[1:top_n]
    pathways_top <- pathways_top[!is.na(pathways_top)]
    message(paste0("Top ", top_n, " Pathways for condition: ", condition, " (Not specified)"))
    message("If fewer than ", top_n, " pathways are displayed, only the available pathways have been found.")
    cat(pathways_top, sep = ";\n")
    
    if (length(pathways_top) == 0) {
      warning("No top pathways found for condition: ", condition)
    } else {
      seurat_obj <- readRDS(paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", condition, ".rds"))
      doCellComVisu(X = cellchat_object, Y = seurat_obj, pathways_to_show = pathways_top, condition = condition)
    }
  } else {
    for (condition in conditions) {
      # Get the top n pathways
      cellchat_object <- readRDS(paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", condition, ".rds"))
      pathways_top <- cellchat_object@netP$pathways[1:top_n]
      pathways_top <- pathways_top[!is.na(pathways_top)]
      message(paste0("Top ", top_n, " Pathways for condition: ", condition))
      message("If fewer than ", top_n, " pathways are displayed, only the available pathways have been found.")
      cat(pathways_top, sep = ";\n")
      
      if (length(pathways_top) == 0) {
        warning("No top pathways found for condition: ", condition)
      } else {
        seurat_obj <- readRDS(paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", condition, ".rds"))
        doCellComVisu(X = cellchat_object, Y = seurat_obj, pathways_to_show = pathways_top, condition = condition)
      }
    }
  }
  
  # Perform pairwise CellChat comparisons
  if (!is.null(conditions_cmp)) {
    for (cmp in conditions_cmp) {
      # Conditions for pair wise comparison
      cond_1 <- cmp[1]
      cond_2 <- cmp[2]
      if (cond_1 == cond_2) {
        warning("Skip the condition comparison since two conditions are identical: ", cond_1)
        next()
      }
      
      # Load the seurat objects and cellchat objects
      # Condition 1
      seurat_obj_cond1 <- readRDS(paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", cond_1, ".rds"))
      cellchat_obj_cond1 <- readRDS(paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", cond_1, ".rds"))
      # Condition 2
      seurat_obj_cond2 <- readRDS(paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", cond_2, ".rds"))
      cellchat_obj_cond2 <- readRDS(paste0(dir_cellchat, "cellchat/rds/cellchat_obj_", cond_2, ".rds"))
      
      run_cellchatV2_cmp(dir_cellchat, seurat_obj_cond1, cellchat_obj_cond1, seurat_obj_cond2, cellchat_obj_cond2, condition_col, condition_1 = cond_1, condition_2 = cond_2, top_n)
      
    }
  }
}
