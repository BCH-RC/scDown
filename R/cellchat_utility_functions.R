#' Create directories for storing CellChat results.
#'
#' This function creates a set of subdirectories under the specified folder path 
#' to organize CellChat results, including directories for figures, tables, and RDS files.
#'
#' @param dir_cellchat The folder path where CellChat results will be stored, including subdirectories for RDS files, figures, and tables.
#' @return NULL
#' 
#' @noRd

create_dir_cellchat <- function(dir_cellchat) {
  # Create folders for storing rds files, figures and tables
  subdirectories <- c("/cellchat",
                      "/cellchat/rds",
                      "/cellchat/csv",
                      "/cellchat/images",
                      "/cellchat/images/aggregate",
                      "/cellchat/images/pathway",
                      "/cellchat/images/pathway/LR_gene",
                      "/cellchat/images/comparison",
                      "/cellchat/images/comparison/Net",
                      "/cellchat/images/comparison/infoFlow",
                      "/cellchat/images/comparison/sidebyside")
  
  for(dir.i in subdirectories){
    dir.create(paste0(dir_cellchat, dir.i), showWarnings = F, recursive = T)
  }
}

#' Check the input data and report error if conditions are not satisfied
#'
#' Take the input for running CellChat analysis, verify if all conditions are met, 
#' and report any errors. Provide guidance on resolving the error messages if 
#' any issues are detected.
#'
#' @param dir_cellchat Path to the folder for storing CellChat results.
#' @param seurat_obj Seurat object containing UMI counts and metadata.
#' @param celltype_col Name of the metadata column used for cell type annotation in the Seurat object.
#' @param celltypes Cell types of interest for running CellChat analysis.
#' @param species Species for the data, either 'human' or 'mouse'.
#' @param condition_col Name of the metadata column for conditions or groups in the Seurat object.
#' @param conditions_cmp List containing the pairwise condition comparisons for CellChat analysis.
#' @return NULL
#' 
#' @noRd

cellchat_input_check <- function(dir_cellchat, seurat_obj, celltype_col, celltypes = "ALL", species, condition_col = NULL, conditions_cmp = NULL) {
  
  # Create cellchat directory
  create_dir_cellchat(dir_cellchat)
  
  # Check the meta data columns
  meta <- seurat_obj@meta.data
  if (!celltype_col %in% colnames(meta)) {
    stop(celltype_col, " cell type column does not exist in the Seurat object metadata. Please specify the correct metadata column name for cell type annotation.")
  }
  
  # Check cell types of interest
  if (celltypes == "ALL") {
    # No cell types of interest are specified
  } else {
    celltypes <- intersect(celltypes, unique(meta.data[, celltype_col]))
    celltypes_missing <- setdiff(celltypes, unique(meta.data[, celltype_col]))
    if (length(celltypes_missing)>0) {
      warning(celltypes_missing, " provided are not found in the cell types!")
    }
    if (length(celltypes)<2) {
      stop("There are not enough cell types! Need to provide at least two cell types of interest!")
    }
  }
  
  # Check species
  if (!species %in% c("human", "mouse")) {
    stop("Please use 'human' or 'mouse' data!")
  }
  
  # Check condition
  if (is.null(condition_col)) {
    # No condition or group information are specified!
    if (!is.null(conditions_cmp)) {
      stop("Please specify the condition column in the Seurat object metadata if you are performing condition comparisons for CellChat analysis! Otherwise, set conditions_cmp = NULL if no comparisons between conditions or groups are needed.")
    }
  } else {
    if (!condition_col %in% colnames(meta)) {
      stop(condition_col, " condition column does not exist in the Seurat object metadata! Please specify the correct metadata column name for condition/group information.")
    } else {
      missing_conditions <- setdiff(unique(unlist(conditions_cmp)), unique(meta[, condition_col]))
      if (length(missing_conditions)>0) {
        stop(missing_conditions, " do not exist in the Seurat object condition metadata column!")
      }
    }
  }
}

#' Do cell-cell communication analysis using cellchat and create a cellchat V2 object
#'
#' Takes as input a Seurat object with cell-type labels as identity of the cells.
#' This Seurat object should have identities populated at Idents(X) and counts
#' normalized in X@assays$RNA@data.
#'
#' @param X a Seurat object.
#' @param species Species for the data, either 'human' or 'mouse'.
#' @return A CellChat object with cell-cell communication analysis results.
#' 
#' @noRd

doCellCom <- function(X, species) {
  ccMetaData <- data.frame(label = Idents(X))
  ccMetaData <- cbind(ccMetaData, X@meta.data)
  ccX <- createCellChat(X@assays$RNA@data, meta = ccMetaData, group.by = 'label')
  if (species == "mouse"){
    ccDB <- CellChatDB.mouse
  } else if (species == "human"){
    ccDB <- CellChatDB.human
  } else {
    print("Other species currently not supported.")
  }
  ccX@DB <- ccDB
  ccX <- subsetData(ccX)
  ccX <- identifyOverExpressedGenes(ccX)
  ccX <- identifyOverExpressedInteractions(ccX)
  ccX <- computeCommunProb(ccX)
  ccX <- filterCommunication(ccX, min.cells = 10)
  ccX <- computeCommunProbPathway(ccX)
  ccX <- aggregateNet(ccX)
  ccX <- netAnalysis_computeCentrality(ccX)
  return(ccX)
}

#' Run CellChat visualization at the aggregated level.
#' 
#' Take as input a CellChat object and graph the aggregated
#' cell-cell communication network.
#'
#' The CellChat object needs to have centrality scores calculated.
#'
#' @param X a CellChat object.
#' @param condition the condition of the object used for naming files.
#' 
#' @noRd

aggregate_visu <- function(X, condition, dir_cellchat){
  
  groupSize <- as.numeric(table(X@idents))
  
  # Circle plot: interaction strength and total interactions for all cell types
  # According to https://github.com/sqjin/CellChat/issues/499, position of vertex labels cannot be changed?
  png(paste0(dir_cellchat, "/cellchat/images/aggregate/", condition, "_net_interaction_and_weight.png", sep=""), height = 600*2, width = 800*3, res=300)
  par(mfrow = c(1, 2), xpd=TRUE)
  netVisual_circle(X@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(X@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
  # Circle plot: interaction strength for each individual cell type
  mat <- X@net$weight
  png(paste0(dir_cellchat, "/cellchat/images/aggregate/", condition, "_net_weight_per_celltype.png", sep=""), height = 600*3*ceiling(length(groupSize)/4), width = 600*4*3, res = 300)
  par(mfrow = c(ceiling(length(groupSize)/4),4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
  
  # Signaling role analysis on the aggregated communication network from all signaling pathways
  p1 <- netAnalysis_signalingRole_scatter(X)
  ggsave(file=paste0(dir_cellchat, "/cellchat/images/aggregate/", condition, "_signaling_role.png", sep=""), plot=p1, height = 6, width = 6)
  
  # Signals contributing most to outgoing or incoming signaling of cell types, need to load ComplexHeatmap library 
  tryCatch(
    {
      pathway_num <- length(X@netP$pathways) # number of pathways that will be shown in heatmap, use this to tune figure height
      # "/50" is used here because even with really large datasets, number of pathways normally won't exceed 100.
      # "/30" is used here because 30 cell types are the maximum a png figure of width 800*2.7 can take.
      # these values can be modified to tune to different figure sizes.
      png(paste0(dir_cellchat, "/cellchat/images/aggregate/", condition, "_outgoing_incoming_signal.png", sep=""), height = 600*1.8*(ceiling(pathway_num/50)), width = 800*2.7*ceiling(length(groupSize)/30), res = 300)
      ht1 <- netAnalysis_signalingRole_heatmap(X, pattern = "outgoing", height = 10*ceiling(pathway_num/50), width = 10*ceiling(length(groupSize)/30), font.size = 6)
      ht2 <- netAnalysis_signalingRole_heatmap(X, pattern = "incoming", height = 10*ceiling(pathway_num/50), width = 10*ceiling(length(groupSize)/30), font.size = 6)
      draw(ht1 + ht2)
      dev.off()
    },
    error=function(cond)
    {
      message(paste("\n**Check outgoing or incoming degree values in X@netP$centr for each pathway, at least one of the pathway need to have more than one value."))
      message(paste("**Outgoing/Incoming signaling role heatmap cannot be produced for this dataset."))
      message(paste("**Here's the original error message: ", cond))
      # Choose a return value in case of error
      return(NA)
    })
}

#' Run CellChat Visualization at the Pathway Level
#' 
#' Take as input a CellChat object and a pathway in which one
#' wants to focus on. Graph the communication network for that
#' specific pathway.
#' 
#' Could add one more parameter to choose specific figure layouts.
#' The CellChat object needs to have centrality scores calculated.
#' 
#' @param X a CellChat object.
#' @param Y a Seurat object which is corresponding to X.
#' @param pathway a signaling pathway in interest. 
#' @param condition the condition of the object used for naming files. 
#' 
#' @noRd

pathway_visu <- function(X, Y, pathway, condition, dir_cellchat, species){
  
  # interaction strength for the pathway
  png(paste0(dir_cellchat, "/cellchat/images/pathway/", pathway, "_", condition, "_signaling_strength_chord.png", sep=""), height = 600*2, width = 600*2, res = 300, pointsize = 8)
  netVisual_aggregate(X, signaling = pathway, title.space = 4, layout = "chord")
  dev.off()
  png(paste0(dir_cellchat, "/cellchat/images/pathway/",pathway,"_",condition,"_signaling_strength_circle.png", sep=""), height = 600*2.5, width = 600*2, res = 300)
  netVisual_aggregate(X, signaling = pathway, title.space=4, layout = "circle")
  dev.off()
  
  # need to load ComplexHeatmap
  png(paste0(dir_cellchat, "/cellchat/images/pathway/", pathway, "_", condition,"_signaling_strength_heatmap.png", sep=""),height = 600*3, width = 600*3, res = 300)
  ht3 <- netVisual_heatmap(X, signaling = pathway, color.heatmap = "Reds")
  draw(ht3)
  dev.off()
  
  # contribution of specific ligand/receptor pair to this pathway
  p1 <- netAnalysis_contribution(X, signaling = pathway)
  ggsave(file = paste0(dir_cellchat, "/cellchat/images/pathway/LR_gene/", pathway, "_", condition, "_LR_contribution.png", sep=""), plot = p1, height = 6, width = 8)
  
  # extract significant L-R pairs contributing to the pathway
  pairLR <- extractEnrichedLR(X, signaling = pathway, geneLR.return = FALSE)
  # cell-cell communication mediated by a single ligand-receptor pair
  for (eachLR in pairLR$interaction_name){
    png(paste0(dir_cellchat, "/cellchat/images/pathway/LR_gene/", pathway, "_", condition, "_", eachLR, ".png", sep=""), height = 600*2, width = 600*2, res = 300, pointsize = 8)
    netVisual_individual(X, signaling = pathway, pairLR.use = eachLR, layout = "chord")
    dev.off()
  }
  
  # plot signaling gene expression distribution related to the pathway
  pairLR <- extractEnrichedLR(X, signaling = pathway, geneLR.return = FALSE) # The extractEnrichedLR() function from CellChat returns ligand-receptor (LR) pairs in upper case by default, even if the CellChat object is based on mouse data.
  LRs_uni <- unique(unlist(strsplit(split = "_", x = pairLR$interaction_name)))
  if (species == "mouse") {
    genes <- rownames(X@data)
    indices <- match(LRs_uni, toupper(genes))
    LRs_uni <- genes[indices]
  }
  if (length(LRs_uni) == 1) {
    p2 <- VlnPlot(
      object = Y,
      features = LRs_uni, 
      pt.size = -1,
    )
    png(paste0(dir_cellchat, "/cellchat/images/pathway/LR_gene/", pathway, "_", condition, "_signaling_gene.png", sep=""), width = 300+150*length(levels(Y)), height = 1200, res = 300)
    print(p2)
    dev.off()
  } else {
    p2 <- VlnPlot(
      object = Y,
      features = LRs_uni, 
      pt.size = -1,
      stack = TRUE
    )
    png(paste0(dir_cellchat, "/cellchat/images/pathway/LR_gene/", pathway, "_", condition, "_signaling_gene.png", sep=""), width = 300+150*length(levels(Y)), height = 600+300*length(LRs_uni), res = 300)
    print(p2)
    dev.off()
  }

  # signaling role analysis on pathway of interest
  png(paste0(dir_cellchat, "/cellchat/images/pathway/", pathway, "_", condition, "_signaling_role_heatmap.png", sep=""),height = 600*1.2,width = 800*1.5, res=300)
  netAnalysis_signalingRole_network(X, signaling = pathway, font.size=6)
  dev.off()
  p3 <- netAnalysis_signalingRole_scatter(X, signaling = pathway)
  ggsave(file=paste0(dir_cellchat, "/cellchat/images/pathway/", pathway, "_", condition, "_signaling_role_scatter.png", sep=""), plot=p3, height = 6, width = 6)
  
  # Bubble plots for LR pairs
  p <- netVisual_bubble(X, signaling = pathway, remove.isolate = TRUE, font.size = 7)
  png(file=paste0(dir_cellchat, "/cellchat/images/pathway/LR_gene/", pathway, "_", condition, "_LR_bubble_plot.png"), res = 300, height = 600+120*length(unique(p$data$interaction_name)), width = 600+25*length(unique(p$data$source.target)))
  print(p)
  dev.off()
}

#' Takes as input a CellChat object with communication analysis results
#' and a vector of pathway names to show. Call aggregate_visu, pathway_visu,
#' to graph all visualizations.
#'
#' @param X a CellChat object.
#' @param Y a Seurat object which is corresponding to X.
#' @param pathways_to_show a vector of pathway names.
#' @param condition the condition of the object.
#' 
#' @noRd

doCellComVisu <- function(X, Y, pathways_to_show, condition, dir_cellchat, species){
  
  # communication at aggregated network level
  aggregate_visu(X, condition, dir_cellchat)
  
  # communication at signaling pathway level
  for (path in pathways_to_show) {
    pathway_visu(X, Y, path, condition, dir_cellchat, species)
  }
  
}

#' Take as input(s) one or two CellChat object(s) and output 
#' pathways with highest overeall communication probabilities.
#'
#' @param top_n the number of pathways.
#' @param X1 a CellChat object.
#' @param X2 a CellChat object.
#' @return the top "top_n" pathways in a vector.

top_pathways <- function(X1, X2=NULL, top_n=10){
  
  df.netP <- X1@netP$pathways[1:top_n]
  
  if (!(is.null(X2))){
    df.netP <- union(df.netP, X2@netP$pathways[1:top_n])
  }
  
  return(df.netP)
  
}


#' Take as inputs two CellChat objects with different cell type labels.
#' Prepare the two objects for merging by aligning cell type labels,
#' net values, and netP values.
#'
#' @param X1 a CellChat object.
#' @param X2 a CellChat object.
#' @return CellChat objects with aligned labels.
#' 
#' @noRd

align_cell_labels <- function(X1, X2){
  
  # get all unique cell types
  group.new <- unique(union(levels(X1@idents), levels(X2@idents)))
  
  # lift cell states for each object
  X1 <- liftCellChat(X1, group.new)
  X2 <- liftCellChat(X2, group.new)
  
  return(list(X1, X2))
  
}

#' Run CellChat visualization
#'
#' run CellChat analysis by comparing pair-wise conditions and generate 
#' figures, tables and rds files.
#'
#' @param dir_cellchat Path to the folder for storing CellChat results.
#' @param seurat_obj_cond1 A seurat object containing UMI counts and metadata for condition 1
#' @param cellchat_obj_cond1 A CellChat object corresponding to seurat_obj_cond1 for condition 1
#' @param seurat_obj_cond2 A seurat object containing UMI counts and metadata for condition 2
#' @param cellchat_obj_cond2 A CellChat object corresponding to seurat_obj_cond2 for condition 2
#' @param condition_col Name of the metadata column for conditions or groups in the Seurat object.
#' @param condition_1 Condition or group 1
#' @param condition_2 Condition or group 2
#' @param top_n the number of pathways.
#' 
#' @return NULL
#' 
#' @noRd

run_cellchatV2_cmp <- function(dir_cellchat, seurat_obj_cond1, cellchat_obj_cond1, seurat_obj_cond2, cellchat_obj_cond2, condition_col, condition_1, condition_2, top_n) {
  
  # if they do not have the same cell type labels
  if (!(identical(levels(cellchat_obj_cond1@idents), levels(cellchat_obj_cond2@idents)))){
    message("Aligning cell types between objects.")
    aligned <- align_cell_labels(cellchat_obj_cond1, cellchat_obj_cond2)
    cellchat_obj_cond1 <- unlist(aligned[1])
    cellchat_obj_cond2 <- unlist(aligned[2])
  }
  
  # merge cellchat objects from 2 biological conditions, the objects being merged need to have the same cell type annotations
  object_list <- list()
  object_list[condition_1] <- cellchat_obj_cond1
  object_list[condition_2] <- cellchat_obj_cond2
  cellchat <- mergeCellChat(object_list, add.names = names(object_list))
  
  # record conditions and pathways in comparison
  cond_in_compare <- levels(cellchat@meta$datasets)
  pathways_to_compare <- top_pathways(X1 = cellchat_obj_cond1, X2 = cellchat_obj_cond2, top_n = top_n)
  message("Top pathways to compare ", condition_1, " and ", condition_2,  " are calculated.")
  cat(pathways_to_compare, sep = ";\n")
  
  # Workflow and visualization for comparisons across conditions
  cellchat <- compareCellComVisu(dir_cellchat, cellchat, object_list, cond_in_compare, pathways_to_compare)
  
  # save merged cellchat object
  saveRDS(cellchat, file = paste0(dir_cellchat, "/cellchat/rds/", cond_in_compare[1], "_", cond_in_compare[2], "_CellChat.rds"))
  
  message("CellChat V2 Differential analysis completed.")
}

####### Functions below this line need to be organized

#' Take as input a merged CellChat object and a list of CellChat object
#' prior to the merge, and output the general comparison results such as
#' number of interactions and aggregated interaction strength.
#' 
#' @param X a merged CellChat object with two biological conditions.
#' @param object_list a list of CellChat objects before the merge.
#' @param cond_in_compare a vector of condition names in comparison.
#' 
#' @noRd

network_comparison <- function(dir_cellchat, X, object_list, cond_in_compare){

  # total number of interactions and strength between conditions
  gg1 <- compareInteractions(X, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(X, show.legend = F, group = c(1,2), measure = "weight")
  p1 <- gg1 + gg2
  ggsave(file=paste0(dir_cellchat, "/cellchat/images/comparison/Net/",cond_in_compare[1],"_",cond_in_compare[2],"_interactNum_histo", ".png", sep=""), plot=p1, height = 6, width = 8)

  # differential number of interactions and strength for each cell type in heatmap
  gg1 <- netVisual_heatmap(X)
  gg2 <- netVisual_heatmap(X, measure = "weight")
  png(paste0(dir_cellchat, "/cellchat/images/comparison/Net/",cond_in_compare[1],"_",cond_in_compare[2],"_diff_interaction", ".png", sep=""),height = 600*3,width = 800*4, res=300)
  draw(gg1 + gg2)
  dev.off()

  # compare the outgoing and incoming interaction strength in 2D space
  num.link <- sapply(object_list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(object_list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(object_list[[i]], title = names(object_list)[i], weight.MinMax = weight.MinMax)
  }
  p2 <- patchwork::wrap_plots(plots = gg)
  ggsave(file=paste0(dir_cellchat, "/cellchat/images/comparison/Net/",cond_in_compare[1],"_",cond_in_compare[2],"_diff_outgo_income_strength",".png", sep=""), plot=p2, height = 6, width = 10)

  # take cell type labels
  cell_groups <- levels(X@idents$joint)
  
  # identify specific signaling changes associated with each cell type
  for (i in unique(cell_groups)){
    # bug fixing: cellchat can have many dataset-specific errors due to the amount of analysis it supports. Here when drawing signaling changes
    # scatter plot, if a cell type has cell-cell communication that has been filtered out due to small number of cells in both conditions, i.e
    # cellchatObj@net$weight for that cell type is zero across both conditions, calling netAnalysis_signalingChanges_scatter() on it will cause
    # an error.
    tryCatch(
    {
      p <- netAnalysis_signalingChanges_scatter(X, idents.use = i)
      ggsave(file=paste0(dir_cellchat, "/cellchat/images/comparison/Net/",cond_in_compare[1],"_",cond_in_compare[2],"_signaling_change_",i,".png", sep=""), plot=p, height = 6, width = 6)
    },
    error=function(cond)
    {
      message(paste("\n**This cell type: ", i, " likely have zero weights in both conditions. Thus not revealing signaling changes"))
      message(paste("**To check: sum(object_list[[1]]@net$weight[, i]) + sum(object_list[[2]]@net$weight[, i]) = ", sum(object_list[[1]]@net$weight[, i])+sum(object_list[[2]]@net$weight[, i])))
      message(paste("**Here's the original error message: ", cond))
      # Choose a return value in case of error
      return(NA)
    })
  }

}


#' #' Take as input a CellChat object and quantify the similarity between
#' #' all significant signaling pathways and then group them based on their
#' #' cellular communication network similarity.
#' #' 
#' #' Need to use conda env to run python "umap-learn" package
#' #' 
#' #' @param X a CellChat object. 
#' #' @param cond_in_compare a vector of condition names in comparison.
#' #' @return X a CellChat object with manifold learning results. 
#' #' 
#' #' @noRd
#' 
#' manifold_learning <- function(X, cond_in_compare){
#' 
#'   # Manifold learning: functional
#'   X <- computeNetSimilarityPairwise(X, type = "functional")
#'   X <- netEmbedding(X, type = "functional")
#'   X <- netClustering(X, type = "functional", do.parallel=FALSE)
#' 
#'   # Visualization in 2D-space
#'   p1 <- netVisual_embeddingPairwise(X, type = "functional", label.size = 3.5)
#'   ggsave(file=paste0(dir_cellchat, "/cellchat/images/comparison/manifold/",cond_in_compare[1],"_",cond_in_compare[2],"_functional_similarity", ".png", sep=""), plot=p1,height = 6, width = 8)
#'   
#'   # Manifold learning: structural
#'   X <- computeNetSimilarityPairwise(X, type = "structural")
#'   X <- netEmbedding(X, type = "structural")
#'   X <- netClustering(X, type = "structural", do.parallel=FALSE)
#' 
#'   # Visualization in 2D-space
#'   p2 <- netVisual_embeddingPairwise(X, type = "structural", label.size = 3.5)
#'   ggsave(file=paste0(dir_cellchat, "/cellchat/images/comparison/manifold/",cond_in_compare[1],"_",cond_in_compare[2],"_structural_similarity", ".png", sep=""), plot=p2,height = 6, width = 8)
#' 
#'   # identify signaling networks with larger/less difference based on Euclidean distance
#'   p3 <- rankSimilarity(X, type = "functional")
#'   ggsave(file=paste0(dir_cellchat, "/cellchat/images/comparison/manifold/",cond_in_compare[1],"_",cond_in_compare[2],"_functional_rank", ".png", sep=""), plot=p3,height = 6, width = 6)
#'   p4 <- rankSimilarity(X, type = "structural")
#'   ggsave(file=paste0(dir_cellchat, "/cellchat/images/comparison/manifold/",cond_in_compare[1],"_",cond_in_compare[2],"_structural_rank", ".png", sep=""), plot=p4,height = 6, width = 6)
#' 
#'   return(X)
#' }


#' Take as input a CellChat object and output the information
#' flow comparison results.
#' 
#' @param X a merged CellChat object with two biological conditions.
#' @param object_list a list of CellChat objects before the merge.
#' @param cond_in_compare a vector of condition names in comparison.
#' 
#' @noRd

information_flow <- function(dir_cellchat, X, object_list,cond_in_compare){
 
  # significant signaling pathways based on differences in the overall information flow
  gg1 <- rankNet(X, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(X, mode = "comparison", stacked = F, do.stat = TRUE)
  p1 <- gg1 + gg2
  # use the number of pathways showing in plot to tune height
  pathway_in_plot <- length(levels(gg1$data$name))
  ggsave(file=paste0(dir_cellchat, "/cellchat/images/comparison/infoFlow/",cond_in_compare[1],"_",cond_in_compare[2],"_significant_pathway_rank", ".png", sep=""), plot=p1, height = 0.1*pathway_in_plot, width = 8)

  # compare outgoing signaling associated with each cell population 
  i = 1
  # use the number of union pathways to tune heatmap height, and number of cell types to tune heatmap width
  pathway_union <- union(object_list[[i]]@netP$pathways, object_list[[i+1]]@netP$pathways)
  pathway_union_length <- length(pathway_union)
  joint_cell_type <- length(levels(X@idents$joint))

  ht1 = netAnalysis_signalingRole_heatmap(object_list[[i]], pattern = "outgoing", signaling = pathway_union, title = names(object_list)[i], height=10*ceiling(pathway_union_length/50), width = 10*ceiling(joint_cell_type/30), font.size=6)
  ht2 = netAnalysis_signalingRole_heatmap(object_list[[i+1]], pattern = "outgoing", signaling = pathway_union, title = names(object_list)[i+1], height=10*ceiling(pathway_union_length/50), width = 10*ceiling(joint_cell_type/30), font.size=6)
  png(paste0(dir_cellchat, "/cellchat/images/comparison/infoFlow/",cond_in_compare[1],"_",cond_in_compare[2],"_diff_outgoing_interaction", ".png", sep=""),height = 600*1.8*(ceiling(pathway_union_length/50)), width = 800*2.7*ceiling(joint_cell_type/30), res=200)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()

  # compare incoming signaling associated with each cell population 
  ht1 = netAnalysis_signalingRole_heatmap(object_list[[i]], pattern = "incoming", signaling = pathway_union, title = names(object_list)[i], height=10*ceiling(pathway_union_length/50), width = 10*ceiling(joint_cell_type/30), font.size=6, color.heatmap = "GnBu")
  ht2 = netAnalysis_signalingRole_heatmap(object_list[[i+1]], pattern = "incoming", signaling = pathway_union, title = names(object_list)[i+1], height=10*ceiling(pathway_union_length/50), width = 10*ceiling(joint_cell_type/30), font.size=6, color.heatmap = "GnBu")
  png(paste0(dir_cellchat, "/cellchat/images/comparison/infoFlow/",cond_in_compare[1],"_",cond_in_compare[2],"_diff_incoming_interaction", ".png", sep=""),height = 600*1.8*(ceiling(pathway_union_length/50)), width = 800*2.7*ceiling(joint_cell_type/30), res=200)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()

}

#' Take as input a CellChat object and output the differential
#' ligand/receptor pair analysis result.
#' 
#' @param X a merged CellChat object with two biological conditions.
#' @param cond_in_compare a vector of condition names in comparison.
#' @return X a CellChat object with differential LR pair results.
#' 
#' @noRd

differential_ligand_receptor <- function(dir_cellchat, X, cond_in_compare){

  # DEG by communication probability: max.dataset = keep the communications with highest probability in max.dataset
  gg1 <- netVisual_bubble(X, comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in", cond_in_compare[2]), angle.x = 45, remove.isolate = T)
  gg2 <- netVisual_bubble(X, comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased signaling in", cond_in_compare[2]), angle.x = 45, remove.isolate = T)
  write.csv(gg1$data, file=paste0(dir_cellchat, "/cellchat/csv/",cond_in_compare[2],"_increased_signalingLR_commProb.csv", sep=""))
  write.csv(gg2$data, file=paste0(dir_cellchat, "/cellchat/csv/",cond_in_compare[2],"_decreased_signalingLR_commProb.csv", sep=""))

  # DEG by differential gene expression
  # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
  pos.dataset = cond_in_compare[2]
  features.name = pos.dataset
  # perform differential expression analysis
  X <- identifyOverExpressedGenes(X, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
  # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
  net <- netMappingDEG(X, features.name = features.name)
  # extract the ligand-receptor pairs with upregulated ligands in pos.dataset
  net.up <- subsetCommunication(X, net = net, datasets = cond_in_compare[2], ligand.logFC = 0.2, receptor.logFC = NULL)
  # extract the ligand-receptor pairs with upregulated ligands in the other dataset, i.e.,downregulated in pos.dataset
  net.down <- subsetCommunication(X, net = net, datasets = cond_in_compare[1], ligand.logFC = -0.1, receptor.logFC = -0.1)
  write.csv(net.up, file=paste0(dir_cellchat, "/cellchat/csv/",cond_in_compare[2],"_increased_signalingLR_diffExpession.csv", sep=""))
  write.csv(net.down, file=paste0(dir_cellchat, "/cellchat/csv/",cond_in_compare[2],"_decreased_signalingLR_diffExpession.csv", sep=""))

  return(X)

}

#' Take as input a CellChat object and output the side-by-side comparison
#' of a pathway's signaling strength in chord diagram.
#' 
#' @param X a merged CellChat object with two biological conditions.
#' @param object_list a list of CellChat object before the merge.
#' @param cond_in_compare a vector of condition names in comparison.
#' @param pathway a pathway of interest.
#' 
#' @noRd

side_by_side_path_compr <- function(dir_cellchat, X, object_list, cond_in_compare, pathway){
  
  png(paste0(dir_cellchat, "/cellchat/images/comparison/sidebyside/",cond_in_compare[1],"_",cond_in_compare[2],"_",pathway,"_sidebyside_strength", ".png", sep=""),height = 600*2,width = 800*3, res=200, pointsize = 10)
  par(mfrow = c(1,2), xpd=TRUE)
  par(mar = c(0.1, 1, 1, 1))
  for (i in 1:length(object_list)) {
    tryCatch(
    {
      netVisual_aggregate(object_list[[i]], signaling = pathway, layout = "chord", signaling.name = paste(pathway, names(object_list)[i]))
    },
    error=function(cond)
    {
      message(paste("Pathway does not exist", pathway))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    })
  }
  dev.off()

}

#' Do the workflow for cell-cell communication analysis on two
#' biological conditions and output relevant visualizations.
#'
#' @param X a merged CellChat object with two biological conditions.
#' @param object_list a list of CellChat object before the merge.
#' @param cond_in_compare a vector of condition names in comparison.
#' @param pathways_to_compare a vector of pathway names.
#' @return X a CellChat object after pairwise comparison workflow.
#' 
#' @noRd

compareCellComVisu <- function(dir_cellchat, X, object_list, cond_in_compare, pathways_to_compare){
  
  # general network inference and comparison
  network_comparison(dir_cellchat, X, object_list, cond_in_compare)
  # functional and structural similarity
  # X <- manifold_learning(X, cond_in_compare)
  # compare information flow
  information_flow(dir_cellchat, X, object_list,cond_in_compare)
  # find differential ligand-rceptor pairs
  X <- differential_ligand_receptor(dir_cellchat, X, cond_in_compare)
  # graph specific pathways of interests side by side for visual comparison
  for (pathway in pathways_to_compare) {
    side_by_side_path_compr(dir_cellchat, X, object_list, cond_in_compare, pathway)
  }

  return(X)

}

# modified original subsetCellChat() soure code so it will not raise dimension error
# See https://github.com/sqjin/CellChat/issues/210 for more information
# this issue only arises when using R 4.2.3, not R 4.1

subsetCellChatMod <- function(object, cells.use = NULL, idents.use = NULL, group.by = NULL, invert = FALSE, thresh = 0.05) {
    if (!is.null(idents.use)) {
      if (is.null(group.by)) {
        labels <- object@idents
        if (object@options$mode == "merged") {
          message("Use the joint cell labels from the merged CellChat object")
          labels <- object@idents$joint
        }
      } else {
        labels <- object@meta[[group.by]]
      }
      if (!is.factor(labels)) {
        labels <- factor(labels)
      }
      level.use0 <- levels(labels)
      level.use <- levels(labels)[levels(labels) %in% unique(labels)]
      
      if (invert) {
        level.use <- level.use[!(level.use %in% idents.use)]
      } else {
        level.use <- level.use[level.use %in% idents.use]
      }
      cells.use.index <- which(as.character(labels) %in% level.use)
      cells.use <- names(labels)[cells.use.index] # NULL
    } else if (!is.null(cells.use)) {
      labels <- object@idents
      if (object@options$mode == "merged") {
        message("Use the joint cell labels from the merged CellChat object")
        labels <- object@idents$joint
      }
      level.use0 <- levels(labels)
      level.use <- levels(labels)[levels(labels) %in% unique(as.character(labels[cells.use]))]
      cells.use.index <- which(names(labels) %in% cells.use)
    } else {
      stop("USER should define either `cells.use` or `idents.use`!")
    }
    cat("The subset of cell groups used for CellChat analysis are ", level.use, '\n')
    
    if (nrow(object@data) > 0) {
      data.subset <- object@data[, cells.use.index]
    } else {
      data.subset <- matrix(0, nrow = 0, ncol = 0)
    }
    if (nrow(object@data.project) > 0) {
      data.project.subset <- object@data.project[, cells.use.index]
    } else {
      data.project.subset <- matrix(0, nrow = 0, ncol = 0)
    }
    data.signaling.subset <- object@data.signaling[, cells.use.index]
    
    meta.subset <- object@meta[cells.use.index, , drop = FALSE]
    
    
    if (object@options$mode == "merged") {
      idents <- object@idents[1:(length(object@idents)-1)]
      group.existing <- level.use0[level.use0 %in% level.use]
      group.existing.index <- which(level.use0 %in% level.use)
      net.subset <- vector("list", length = length(object@net))
      netP.subset <- vector("list", length = length(object@netP))
      idents.subset <- vector("list", length = length(idents))
      names(net.subset) <- names(object@net)
      names(netP.subset) <- names(object@netP)
      names(idents.subset) <- names(object@idents[1:(length(object@idents)-1)])
      images.subset <- vector("list", length = length(idents))
      names(images.subset) <- names(object@idents[1:(length(object@idents)-1)])
      
      for (i in 1:length(idents)) {
        cat("Update slots object@images, object@net, object@netP, object@idents in dataset ", names(object@idents)[i],'\n')
        images <- object@images[[i]]
        for (images.j in names(images)) {
          values <- images[[images.j]]
          if (images.j %in% c("coordinates")) {
            values.new <- values[cells.use.index, ]
            images[[images.j]] <- values.new
          }
          if (images.j %in% c("distance")) {
            values.new <- values[group.existing.index, group.existing.index, drop = FALSE]
            images[[images.j]] <- values.new
          }
        }
        images.subset[[i]] <- images
        
        # cat("Update slot object@net...", '\n')
        net <- object@net[[i]]
        for (net.j in names(net)) {
          values <- net[[net.j]]
          if (net.j %in% c("prob","pval")) {
            values.new <- values[group.existing.index, group.existing.index, ]
            net[[net.j]] <- values.new
          }
          if (net.j %in% c("count","sum","weight")) {
            values.new <- values[group.existing.index, group.existing.index]
            net[[net.j]] <- values.new
          }
          # net[[net.j]] <- values.new
        }
        net.subset[[i]] <- net
        
        netP = computeCommunProbPathway(net = net.subset[[i]], pairLR.use = object@LR[[i]]$LRsig, thresh = thresh)
        netP$centr = netAnalysis_computeCentrality(net =  net.subset[[i]]$prob)
        netP.subset[[i]] <- netP
        idents.subset[[i]] <- idents[[i]][names(idents[[i]]) %in% cells.use]
        idents.subset[[i]] <- factor(idents.subset[[i]], levels = levels(idents[[i]])[levels(idents[[i]]) %in% level.use])
      }
      idents.subset$joint <- factor(object@idents$joint[cells.use.index], levels = level.use)
      
    } else {
      cat("Update slots object@images, object@net, object@netP in a single dataset...", '\n')
      
      group.existing <- level.use0[level.use0 %in% level.use]
      group.existing.index <- which(level.use0 %in% level.use)
      
      images <- object@images
      for (images.j in names(images)) {
        values <- images[[images.j]]
        if (images.j %in% c("coordinates")) {
          values.new <- values[cells.use.index, ]
          images[[images.j]] <- values.new
        }
        if (images.j %in% c("distance")) {
          values.new <- values[group.existing.index, group.existing.index, drop = FALSE]
          images[[images.j]] <- values.new
        }
      }
      images.subset <- images
      
      net <- object@net
      for (net.j in names(net)) {
        values <- net[[net.j]]
        if (net.j %in% c("prob","pval")) {
          ################## ISSUE FIXED ##################
          values.new <- values[group.existing.index, group.existing.index, ,drop = FALSE]
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("count","sum","weight")) {
          values.new <- values[group.existing.index, group.existing.index, drop = FALSE]
          net[[net.j]] <- values.new
        }
      }
      net.subset <- net
      
      netP = computeCommunProbPathway(net = net.subset, pairLR.use = object@LR$LRsig, thresh = thresh)
      netP$centr = netAnalysis_computeCentrality(net = net.subset$prob)
      netP.subset <- netP
      idents.subset <- object@idents[cells.use.index]
      idents.subset <- factor(idents.subset, levels = level.use)
    }
    
    object.subset <- methods::new(
      Class = "CellChat",
      data = data.subset,
      data.signaling = data.signaling.subset,
      data.project = data.project.subset,
      images = images.subset,
      net = net.subset,
      netP = netP.subset,
      meta = meta.subset,
      idents = idents.subset,
      var.features = object@var.features,
      LR = object@LR,
      DB = object@DB,
      options = object@options
    )
    return(object.subset)
}



# # ##################### NOT USED FOR NOW ######################
# #' Reorder the cell identities and net calculation results based on given input levels.
# #' https://github.com/sqjin/CellChat/issues/149#issuecomment-788102884
# #'
# #' @param object a CellChat object
# #' @param ident.use the name of the variable in object.meta
# #' @param levels set the levels of factor
# #' @return a CellChat object with updated cell type and net order.
# #' 
# #' @noRd

# reorder_ident <- function(object, ident.use, levels){

#   object@idents <- as.factor(object@meta[[ident.use]])
  
#   if (!is.null(levels)) {
#     object@idents <- factor(object@idents, levels = levels)
#   }
  
#   if (length(object@net) > 0) {
#     if (all(dimnames(object@net$prob)[[1]] %in% levels(object@idents) )) {
#       message("Reorder cell groups! ")
#       idx <- match(dimnames(object@net$prob)[[1]], levels(object@idents))
#       object@net$prob <- object@net$prob[idx, idx, ]
#       object@net$pval <- object@net$pval[idx, idx, ]
#       cat("The cell group order after reordering is ", dimnames(object@net$prob)[[1]],'\n')
#     } else {
#       message("Rename cell groups but do not change the order! ")
#       cat("The cell group order before renaming is ", dimnames(object@net$prob)[[1]],'\n')
#       dimnames(object@net$prob) <- list(levels(object@idents), levels(object@idents), dimnames(object@net$prob)[[3]])
#       dimnames(object@net$pval) <- dimnames(object@net$prob)
#       cat("The cell group order after renaming is ", dimnames(object@net$prob)[[1]],'\n')
#     }
#   }

#   return(object)

# }

