library(CellChat)
library(Seurat)
library(ComplexHeatmap)

#' Create Directories for Storing CellChat Results
#'
#' This function creates a set of subdirectories under the specified folder path 
#' to organize CellChat results, including directories for RDS files, figures, and tables.
#'
#' @param dir_cellchat The folder path where CellChat results will be stored, 
#' including subdirectories for RDS files, figures, and tables.
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
    dir.create(paste0(dir_cellchat, dir.i), showWarnings = FALSE, recursive = TRUE)
  }
}

#' Perform Cell-Cell Communication Analysis Using CellChat and Create a CellChat V2 Object
#'
#' This function performs cell-cell communication analysis using CellChat and generates a CellChat V2 object.
#' It takes a Seurat object as input, with cell annotation labels assigned as identities of the cells.
#' The Seurat object should have cell identities populated in `Idents(X)` and normalized count data in 
#' `X@assays$RNA@data`.
#'
#' @param X A Seurat object with cell-type identities assigned to `Idents(X)` and normalized counts 
#'          in `X@assays$RNA@data`.
#' @param species The species of the data, either 'human' or 'mouse'.
#' @return ccX, A CellChat object containing the results of the cell-cell communication analysis.
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

#' Run CellChat Visualization at the Aggregated Level
#' 
#' This function takes as input a CellChat object and generates a visualization of the aggregated 
#' cell-cell communication network.
#'
#'
#' @param X A CellChat object containing the results of the cell-cell communication analysis.
#' @param condition A character string representing the condition of the object, which is also used for naming output files.
#' @return NULL
#' 
#' @noRd

aggregate_visu <- function(X, condition, dir_cellchat){
  
  groupSize <- as.numeric(table(X@idents))
  numofcelltypes <- length(groupSize)

  # Circle plot: interaction strength and total interactions for all cell types
  # According to https://github.com/sqjin/CellChat/issues/499, position of vertex labels cannot be changed?
  png(paste0(dir_cellchat, "/cellchat/images/aggregate/", condition, "_net_interaction_and_weight.png", sep=""), height = 600*(numofcelltypes/4), width = 800*(numofcelltypes/4+1), res=300)
  par(mfrow = c(1, 2), xpd=TRUE)
  netVisual_circle(X@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(X@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
  # Circle plot: interaction strength for each individual cell type
  mat <- X@net$weight
  png(paste0(dir_cellchat, "/cellchat/images/aggregate/", condition, "_net_weight_per_celltype.png", sep=""), height = 600*3*ceiling(numofcelltypes/4), width = 600*4*3, res = 300)
  par(mfrow = c(ceiling(length(groupSize)/4),4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), vertex.label.cex = 1 + 4/numofcelltypes, title.name = rownames(mat)[i])
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
      png(paste0(dir_cellchat, "/cellchat/images/aggregate/", condition, "_outgoing_incoming_signal.png", sep=""), height = 600*3*(ceiling(pathway_num/50)), width = 800*3.5*ceiling(length(groupSize)/30), res = 300)
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

#' Run CellChat Visualization at the Aggregated Level
#' 
#' This function takes as input a CellChat object and generates a circle plot of the aggregated 
#' cell-cell communication network.
#'
#' @param X A CellChat object containing the results of the cell-cell communication analysis.
#' @param height User defined image height.
#' @param width User defined image width.
#' @param res User defined image resolution.
#' @return NULL
#' 
#' @noRd

aggregate_circleplot <- function(X, dir_cellchat, height, width, res) {
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  fname <- paste0(dir_cellchat, "/cellchat/images/aggregate/net_interaction_and_weight_", timestamp, ".png", sep="")
  cat("Image file saved as", fname, "\n")
  png(fname, height = height, width = width, res=res)
  par(mfrow = c(1, 2), xpd=TRUE)
  netVisual_circle(X@net$count, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(X@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
}

#' Run CellChat Visualization at the Aggregated Level
#' 
#' This function takes as input a CellChat object and generates a heatmap plot of the aggregated 
#' cell-cell communication network.
#'
#' @param X A CellChat object containing the results of the cell-cell communication analysis.
#' @param hp.height User defined heatmap height.
#' @param hp.width User defined heatmap width.
#' @param height User defined image height.
#' @param width User defined image width.
#' @param res User defined image resolution.
#' @return NULL
#' 
#' @noRd

aggregate_heatmap <- function(X, dir_cellchat, font.size = 6, hp.height, hp.width, height, width, res) {
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  fname <- paste0(dir_cellchat, "/cellchat/images/aggregate/outgoing_incoming_signal_", timestamp, ".png", sep="")
  cat("Image file saved as", fname, "\n")
  png(fname, height = height, width = width, res = res)
  ht1 <- netAnalysis_signalingRole_heatmap(X, pattern = "outgoing", height = hp.height, width = hp.width, font.size = font.size)
  ht2 <- netAnalysis_signalingRole_heatmap(X, pattern = "incoming", height = hp.height, width = hp.width, font.size = font.size)
  draw(ht1 + ht2)
  dev.off()
}

#' Run CellChat Visualization at the Aggregated Level
#' 
#' This function takes as input a CellChat object and generates a circle plot per cell type of the aggregated 
#' cell-cell communication network.
#'
#' @param X A CellChat object containing the results of the cell-cell communication analysis.
#' @param image.ncol User defined number of circle plots on each row.
#' @param vertex.label.cex User vertex label size.
#' @param height User defined image height.
#' @param width User defined image width.
#' @param res User defined image resolution.
#' @return NULL
#' 
#' @noRd

aggregate_circleplot_percelltype <- function(X, dir_cellchat, image.ncol, vertex.label.cex, height, width, res) {

  groupSize <- as.numeric(table(X@idents))
  numofcelltypes <- length(groupSize)
  image.nrow <- ceiling(numofcelltypes/image.ncol)

  # Circle plot: interaction strength for each individual cell type
  mat <- X@net$weight
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  fname = paste0(dir_cellchat, "/cellchat/images/aggregate/net_weight_per_celltype_", timestamp, ".png", sep="")
  cat("Image file saved as", fname, "\n")
  png(fname, height = height, width = width, res = 300)
  par(mfrow = c(image.nrow, image.ncol), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, weight.scale = T, edge.weight.max = max(mat), vertex.label.cex = vertex.label.cex, title.name = rownames(mat)[i])
  }
  dev.off()
}

#' Run CellChat Visualization at the Pathway Level
#' 
#' This function takes a CellChat object and a specified pathway to generate a visualization of the 
#' cell-cell communication network for that particular pathway. 
#' Additionally, a parameter for selecting figure layouts can be added. 
#' The CellChat object must have centrality scores calculated prior to running this function.
#' 
#' @param X A CellChat object containing the results of the cell-cell communication analysis.
#' @param Y A Seurat object corresponding to the CellChat object, providing additional context.
#' @param pathway A character string specifying the signaling pathway of interest.
#' @param condition A character string representing the condition of the object, used for naming output files.
#' @return NULL
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
  tryCatch({
    ht3 <- netVisual_heatmap(X, signaling = pathway, color.heatmap = "Reds")
    png(paste0(dir_cellchat, "/cellchat/images/pathway/", pathway, "_", condition,"_signaling_strength_heatmap.png", sep=""),height = 600*3, width = 600*3, res = 300)
    draw(ht3)
    dev.off()
  }, error = function(e) {
    # Print the error message (optional) and continue
    cat("Warning: ", e$message, "\n")
    cat("Creation of signaling_strength_heatmap.png file failed!\n")
  })
  
  # contribution of specific ligand/receptor pair to this pathway
  tryCatch({
    p1 <- netAnalysis_contribution(X, signaling = pathway)
    ggsave(file = paste0(dir_cellchat, "/cellchat/images/pathway/LR_gene/", pathway, "_", condition, "_LR_contribution.png", sep=""), plot = p1, height = 6, width = 8)
  }, error = function(e) {
    # Print the error message (optional) and continue
    cat("Warning: ", e$message, "\n")
    cat("Creation of LR_contribution.png file failed!\n")
  })

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
  # Pathway name correction: for some LR names it also contains pathway name which needs to be removed 
  # LRs_uni <- gsub("RetinoicAcid-RA-", "", LRs_uni)
  genes1 <- LRs_uni[LRs_uni %in% toupper(rownames(Y))]
  genes2 <- LRs_uni[!(LRs_uni %in% toupper(rownames(Y)))]
  genes22 <- sub(".?.?", "", genes2) 
  LRs_uni <- c(genes1, genes22[genes22 %in% toupper(rownames(Y))])
  
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
  p <- netVisual_bubble(X, signaling = pathway, remove.isolate = FALSE, font.size = 7)
  png(file=paste0(dir_cellchat, "/cellchat/images/pathway/LR_gene/", pathway, "_", condition, "_LR_bubble_plot.png"), res = 300, height = 600+120*length(unique(p$data$interaction_name)), width = 600+25*length(unique(p$data$source.target)))
  print(p)
  dev.off()
}

#' Generate Visualizations for Cell-Cell Communication Pathways
#' 
#' This function takes a CellChat object containing communication analysis results and a vector of 
#' pathway names to visualize. It calls `aggregate_visu` and `pathway_visu` to generate and display 
#' visualizations for each pathway in the provided list.
#'
#' @param X A CellChat object containing the results of the cell-cell communication analysis.
#' @param Y A Seurat object corresponding to the CellChat object, providing additional context.
#' @param pathways_to_show A character vector containing the names of the pathways to visualize.
#' @param condition A character string representing the condition of the object, used for naming output files.
#' @return NULL
#' 
#' @noRd

doCellComVisu <- function(X, Y, pathways_to_show, condition, dir_cellchat, species){
  
  # communication at signaling pathway level
  for (path in pathways_to_show) {
    pathway_visu(X, Y, path, condition, dir_cellchat, species)
  }
  
}


#' Calculate and Scale Information Flow for All Communication Pathways
#' 
#' This function calculates the information flow for communication pathways in a network object, 
#' scales the contribution of each pathway using a logarithmic transformation, and returns the 
#' results ordered by contribution.
#'
#' @param X A network object containing the communication probabilities, from which the `prob` matrix is extracted.
#' @param condition A character string representing the condition of the object, used for labeling in the output.
#' @return df_ordered, A data frame with the following columns:
#' 
#' @noRd

calc_infoflow <- function(X, condition) {
  prob <- methods::slot(X, "netP")$prob
  if (sum(prob) == 0) {
    stop("No inferred communications for the input!")
  }
  pSum <- apply(prob, 3, sum)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, length.out = length(idx1))
  position <- sort(pSum.original[idx1], index.return = TRUE)$ix
  pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  pair.name <- names(pSum)
  df <- data.frame(name = pair.name, contribution = pSum.original, contribution.scaled = pSum, group = condition)
  idx <- with(df, order(df$contribution))
  df <- df[idx, ]
  df$name <- factor(df$name, levels = as.character(df$name))
  df_ordered <- df[order(df$contribution, decreasing = TRUE), ]
  
  return(df_ordered)
}

#' Identify Top Pathways with the Highest Overall Communication Probabilities
#' 
#' This function takes one or two CellChat objects as input and returns the top pathways with the 
#' highest overall communication probabilities. The number of top pathways is specified by the 
#' `top_n` parameter.
#'
#' @param top_n An integer specifying the number of top pathways to return.
#' @param X1 A CellChat object containing the results of the first communication analysis.
#' @param X2 (Optional) A second CellChat object to compare communication probabilities. 
#'            If only one object is provided, the function will work with that single object.
#' @return df.netP, A character vector containing the names of the top "top_n" pathways.
#' 
#' @noRd

top_pathways <- function(X1, X2=NULL, top_n=10){
  
  df.netP <- X1@netP$pathways[1:top_n]
  
  if (!(is.null(X2))){
    df.netP <- union(df.netP, X2@netP$pathways[1:top_n])
  }
  
  return(df.netP)
  
}

#' Align Cell Type Labels Across Two CellChat Objects
#' 
#' This function takes two CellChat objects with different cell type labels and prepares them 
#' for merging by aligning their cell type labels, net values, and netP values.
#' The function ensures that both objects have consistent cell type labels.
#'
#' @param X1 A CellChat object with the first set of cell type labels.
#' @param X2 A CellChat object with the second set of cell type labels.
#' @return list(X1, X2), A list containing the two CellChat objects with aligned cell type labels.
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

#' Run CellChat Visualization with Pairwise Condition Comparison
#'
#' This function performs CellChat analysis by comparing cell-cell communication between two conditions,
#' generating figures, tables, and RDS files for differential analysis. The function compares pathways 
#' between two conditions and visualizes the results.
#'
#' @param dir_cellchat Path to the folder where CellChat results (figures, tables, and RDS files) will be stored.
#' @param seurat_obj_cond1 A Seurat object containing UMI counts and metadata for condition 1.
#' @param cellchat_obj_cond1 A CellChat object corresponding to `seurat_obj_cond1` for condition 1.
#' @param seurat_obj_cond2 A Seurat object containing UMI counts and metadata for condition 2.
#' @param cellchat_obj_cond2 A CellChat object corresponding to `seurat_obj_cond2` for condition 2.
#' @param condition_col Name of the metadata column used for conditions or groups in the Seurat object.
#' @param condition_1 The first condition or group for comparison.
#' @param condition_2 The second condition or group for comparison.
#' @param top_n The number of top pathways to compare between the two conditions.
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

#' Compare Cell-Cell Communication Networks Between Two Conditions
#' 
#' This function takes a merged CellChat object and a list of CellChat objects from 
#' two different biological conditions. It outputs general comparison results, 
#' such as the number of interactions, aggregated interaction strength, and 
#' signaling changes for each cell type.
#' 
#' @param dir_cellchat Path to the folder where CellChat results (figures, tables, and RDS files) will be stored.
#' @param X A merged CellChat object containing two biological conditions.
#' @param object_list A list of CellChat objects from each condition prior to merging.
#' @param cond_in_compare A vector of condition names being compared.
#' 
#' @return NULL
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

#' Compare Information Flow Between Two Conditions in CellChat
#' 
#' This function takes a merged CellChat object and a list of individual CellChat objects for each condition.
#' It outputs the comparison results of the information flow between the two biological conditions, 
#' including significant pathways and differential outgoing and incoming signaling associated with each cell population.
#' 
#' @param dir_cellchat Path to the folder where CellChat results (figures, tables, and RDS files) will be stored.
#' @param X A merged CellChat object containing data from two biological conditions.
#' @param object_list A list of CellChat objects from each condition before merging.
#' @param cond_in_compare A vector of condition names being compared.
#' 
#' @return NULL
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

#' Perform Differential Ligand-Receptor Pair Analysis
#' 
#' This function performs a differential analysis of ligand-receptor (LR) pairs using a merged 
#' CellChat object with two biological conditions. It outputs the differential signaling results 
#' based on communication probabilities and differential gene expression analysis. The function 
#' generates CSV files containing the increased and decreased signaling LR pairs.
#' 
#' @param dir_cellchat Path to the folder where CellChat results (figures, tables, and RDS files) will be stored.
#' @param X A merged CellChat object containing two biological conditions for comparison.
#' @param cond_in_compare A vector of condition names being compared.
#' 
#' @return A CellChat object with differential ligand-receptor pair analysis results.
#' 
#' @noRd

differential_ligand_receptor <- function(dir_cellchat, X, cond_in_compare){

  # DEG by communication probability: max.dataset = keep the communications with highest probability in max.dataset
  gg1 <- netVisual_bubble(X, comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in", cond_in_compare[2]), angle.x = 45, remove.isolate = T)
  gg2 <- netVisual_bubble(X, comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased signaling in", cond_in_compare[2]), angle.x = 45, remove.isolate = T)
  # write.csv(gg1$data, file=paste0(dir_cellchat, "/cellchat/csv/",cond_in_compare[2],"_increased_signalingLR_commProb.csv", sep=""))
  # write.csv(gg2$data, file=paste0(dir_cellchat, "/cellchat/csv/",cond_in_compare[2],"_decreased_signalingLR_commProb.csv", sep=""))

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
  # write.csv(net.up, file=paste0(dir_cellchat, "/cellchat/csv/",cond_in_compare[2],"_increased_signalingLR_diffExpession.csv", sep=""))
  # write.csv(net.down, file=paste0(dir_cellchat, "/cellchat/csv/",cond_in_compare[2],"_decreased_signalingLR_diffExpession.csv", sep=""))

  return(X)

}

#' Generate Side-by-Side Comparison of a Pathway's Signaling Strength in a Chord Diagram
#' 
#' This function takes a merged CellChat object and a list of individual CellChat objects to generate 
#' a side-by-side comparison of the signaling strength of a specified pathway. The comparison is visualized 
#' using chord diagrams for two biological conditions.
#' 
#' @param dir_cellchat Path to the folder where CellChat results (figures, tables, and RDS files) will be stored.
#' @param X A merged CellChat object containing two biological conditions.
#' @param object_list A list of CellChat objects from each condition before merging.
#' @param cond_in_compare A vector of condition names being compared.
#' @param pathway A character string representing the pathway of interest.
#' 
#' @return NULL
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

#' Perform Workflow for Cell-Cell Communication Analysis on Two Biological Conditions
#' 
#' This function performs the complete workflow for cell-cell communication analysis by comparing two biological conditions. 
#' It outputs relevant visualizations, including network comparisons, information flow analysis, differential ligand-receptor 
#' pair analysis, and side-by-side comparisons of selected pathways.
#' 
#' @param X A merged CellChat object containing data from two biological conditions.
#' @param object_list A list of CellChat objects from each condition prior to merging.
#' @param cond_in_compare A vector of condition names being compared.
#' @param pathways_to_compare A vector of pathway names to be compared.
#' 
#' @return X, A CellChat object after completing the pairwise comparison workflow, including updated communication results.
#' 
#' @noRd

compareCellComVisu <- function(dir_cellchat, X, object_list, cond_in_compare, pathways_to_compare){
  
  # general network inference and comparison
  network_comparison(dir_cellchat, X, object_list, cond_in_compare)
  # compare information flow
  information_flow(dir_cellchat, X, object_list, cond_in_compare)
  # find differential ligand-rceptor pairs
  X <- differential_ligand_receptor(dir_cellchat, X, cond_in_compare)
  # graph specific pathways of interests side by side for visual comparison
  for (pathway in pathways_to_compare) {
    side_by_side_path_compr(dir_cellchat, X, object_list, cond_in_compare, pathway)
  }

  return(X)

}

#' Subset a CellChat object based on user-defined cell identities or cells to analyze specific cell types or conditions.
#'
#' This function subsets a CellChat object based on specified cell identities (`idents.use`). 
#' It retains relevant data structures, such as cell-cell communication networks, images, and other metadata. 
#' The function is useful for focusing on a subset of cells or cell types for analysis in CellChat.
#'
#' @param object A CellChat object containing data and results of cell-cell communication analysis.
#' @param idents.use A vector of cell annotation identities (e.g., cell types) to subset, optional if `cells.use` is provided.
#' @param thresh A threshold for computing pathway probabilities. Default is 0.05.
#'
#' @return object.subset, A subsetted CellChat object with the same data structures but limited to the specified cells or identities.
#' 
#' @noRd

subsetCellChatMod <- function(object, idents.use, thresh = 0.05) {
  # Extract the labels as the idents which were defined as "annotation_column"
  labels <- object@idents
  if (object@options$mode == "merged") {
    message("Use the joint cell labels from the merged CellChat object")
    labels <- object@idents$joint
  }
  
  # Subsetting the cells based on provided idents.use
  if (!is.factor(labels)) {
    labels <- factor(labels)
  }
  level.use0 <- levels(labels)
  level.use <- levels(labels)[levels(labels) %in% unique(labels)]
  level.use <- level.use[level.use %in% idents.use]
  cells.use.index <- which(as.character(labels) %in% level.use)
  cells.use <- names(labels)[cells.use.index] # NULL
  cat("The subset of cell groups used for CellChat analysis are ", level.use, '\n')
  
  # Subsetting data for the selected cells
  data.subset <- object@data[, cells.use.index]
  data.signaling.subset <- object@data.signaling[, cells.use.index]
  meta.subset <- object@meta[cells.use.index, , drop = FALSE]
  
  # Handling for the merged CellChat object or single CellChat object
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
  
  # Return the subsetted CellChat object
  object.subset <- methods::new(
    Class = "CellChat",
    data = data.subset,
    data.signaling = data.signaling.subset,
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


