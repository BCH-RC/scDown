# load libraries that are installed in the docker container
# library(Seurat)
# library(monocle3)
# library(ggplot2)
# library(igraph)
# library(dplyr)
# library(stringr)
# library(tibble)
# library(pheatmap)
# library(SCENT)
# library(patchwork)

#' Learn trajectory graph for the inputted Seurat object.
#'
#' @param X a Seurat object
#' @param nDim an integer number of dimensions to calculate for PCA
#' @param batch a character string specifying the metadata used to correct for batch effect
#' @param transferUMAP boolean variable of whether to transfer Seurat umap coordinates
#' @param subset a character vector of cell types to subset
#' @param cond a character string specifying a condition, for output naming purpose
#' @return cds, a cell_data_set object with trajectory
#'
#' @importFrom grDevices dev.off jpeg pdf png
#' @importFrom graphics par
#' @importFrom checkmate expect_choice
#' @importFrom stats na.omit
#' @importFrom utils read.csv read.table write.csv
#' @import magrittr
#'
#'
#' @noRd

getTrajectory <- function(X, nDim=30, batch=NULL, transferUMAP=TRUE, subset=NULL, cond=NULL,outputDir="."){

  # when idents = NULL, WhichCells() will simply return all cells, thus no subsetting
  selected_cells <- Seurat::WhichCells(X, idents = subset)
  geneInfo <- data.frame(gene_short_name = rownames(X[,selected_cells]@assays$RNA))
  rownames(geneInfo) <- rownames(X[,selected_cells]@assays$RNA)
  cellInfo <- X[,selected_cells]@meta.data
  cds <- monocle3::new_cell_data_set(X[,selected_cells]@assays$RNA@counts,
                           cell_metadata = cellInfo,
                           gene_metadata = geneInfo)
  cds <- monocle3::preprocess_cds(cds, num_dim = nDim,method = "PCA")
  cds <- monocle3::reduce_dimension(cds,umap.fast_sgd=FALSE)

  if(!is.null(batch)){ # correct batch effect, if any
    p1 <- monocle3::plot_cells(cds, color_cells_by=batch, label_cell_groups=FALSE) +
      ggplot2::ggtitle(paste0("UMAP before batch correction", "\n","~",batch,sep=""))
    
    samples <- unique(SummarizedExperiment::colData(cds)[,batch])
    O <- lapply(samples, function(sample){
      d<-as.data.frame(SummarizedExperiment::colData(cds))
      cds_subset <- cds[,row.names(d[d[batch] == sample,])]
      p2 <- monocle3::plot_cells(cds_subset, show_trajectory_graph=FALSE,color_cells_by=batch, label_cell_groups=FALSE) +
        ggplot2::ggtitle(paste0(sample)) +
        ggplot2::theme(legend.position="none") + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      return(p2)
    })
    
    png(file.path(outputDir,"images","pseudotime",paste0('UMAP_before_batchCorrect',ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond),"",paste0("_",cond)),'.png',sep="")), width = 700*sqrt(length(samples))*2, height = 875*sqrt(length(samples)), res = 300)
    print(patchwork::wrap_plots(O,ncol = ceiling(sqrt(length(samples)))))
    dev.off()
    
    cds <- monocle3::align_cds(cds, alignment_group = batch)
    cds <- monocle3::reduce_dimension(cds,umap.fast_sgd=FALSE)
    
    samples <- unique(SummarizedExperiment::colData(cds)[,batch])
    O <- lapply(samples, function(sample){
      d<-as.data.frame(SummarizedExperiment::colData(cds))
      cds_subset <- cds[,row.names(d[d[batch] == sample,])]
      p2 <- monocle3::plot_cells(cds_subset, show_trajectory_graph=FALSE,color_cells_by=batch, label_cell_groups=FALSE) +
        ggplot2::ggtitle(paste0(sample)) +
        ggplot2::theme(legend.position="none") + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      return(p2)
    })
    
    png(file.path(outputDir,"images","pseudotime",paste0('UMAP_after_batchCorrect',ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond),"",paste0("_",cond)),'.png',sep="")), width = 700*sqrt(length(samples))*2, height = 875*sqrt(length(samples)), res = 300)
    print(patchwork::wrap_plots(O,ncol = ceiling(sqrt(length(samples)))))
    dev.off()
    
    p2 <- monocle3::plot_cells(cds, color_cells_by=batch, label_cell_groups=FALSE) +
      ggplot2::ggtitle(paste0("UMAP after batch correction", "\n","~",batch,sep="")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    combined_plot <- patchwork::wrap_plots(list(p1, p2), ncol = 2)
    ggplot2::ggsave(file=file.path(outputDir,"images","pseudotime",paste0('batchEffectCompare',ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond),"",paste0("_",cond)),'.png',sep="")), plot = combined_plot, width = 10, height = 5)
  }

  if(transferUMAP){     # transfer original Seurat umap coordinates
    if ("umap" %in% names(X@reductions)){
      SingleCellExperiment::reducedDims(cds)[['UMAP']] <- X[,selected_cells]@reductions$umap@cell.embeddings
    } else {
      stop("There is no UMAP coordinates to transfer! Check Seurat object for umap embedding.")
    }
  }
  ###Updated the method to "louvain" as the "leiden" method is not working
  #https://github.com/cole-trapnell-lab/monocle3/issues/666
  #https://github.com/cole-trapnell-lab/monocle3/issues/664
  ##Error in leidenbase::leiden_find_partition(graph_result[["g"]], partition_type = partition_type,  :
  ##REAL() can only be applied to a 'numeric', not a 'NULL'
  cds <- monocle3::cluster_cells(cds,cluster_method = "louvain")
  cds <- monocle3::learn_graph(cds, use_partition = FALSE) # learn a single graph for all partitions

  # plot umap by partitions
  png(filename = file.path(outputDir,"images","pseudotime",paste0("umap_partitions_","transferUMAP_",transferUMAP,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond),"",paste0("_",cond)),".png",sep="")), width = 2000*1.4, height = 1500*1.5, res = 400)
  p3 <- monocle3::plot_cells(cds, color_cells_by="partition", group_cells_by="partition", show_trajectory_graph=FALSE) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(paste0("UMAP by partitions", "\n", "transferUMAP=",transferUMAP, sep="")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p3)
  dev.off()

  # plot umap by cell types
  cds[["cell.type"]] <- unname(Seurat::Idents(X[,selected_cells])) # transfer labels
  png(filename = file.path(outputDir,"images","pseudotime",paste0("umap_celltypes_","transferUMAP_",transferUMAP,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond),"",paste0("_",cond)),".png",sep="")), width = 2000*1.4, height = 1500*1.5, res = 400)
  p4 <- monocle3::plot_cells(cds, color_cells_by="cell.type", show_trajectory_graph=FALSE,label_groups_by_cluster=FALSE) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(paste0("UMAP by cell types","\n", "transferUMAP=",transferUMAP, sep=""))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p4)
  dev.off()

  # plot umap by trajectory
  png(filename = file.path(outputDir,"images","pseudotime",paste0("umap_trajectory_","transferUMAP_",transferUMAP,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond),"",paste0("_",cond)),".png",sep="")), width = 2000*1.4, height = 1500*1.5, res = 400)
  p5 <- monocle3::plot_cells(cds, label_principal_points = TRUE,  color_cells_by = "cell.type", label_cell_groups=FALSE) +
    ggplot2::theme_void() +
    ggplot2::guides(color = ggplot2::guide_legend("Cell Type", override.aes = list(size=5))) +
    ggplot2::ggtitle(paste0("UMAP by trajectory","\n", "transferUMAP=",transferUMAP,sep=""))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p5)
  dev.off()
  return(cds)
}


#' Order the cells in pseudotime according to the chosen starting point.
#'
#' @param cds a cell_data_set object with trajectory
#' @param method a string specifying the method used to select a starting point for pseudotime scale
#' @param rootNodes a vector of strings of starting principal points
#' @param timePoint a character string specifying a category inside a metadata
#' @param timePointCol a character string specifying the name of metadata with the category @timePoint
#' @param species a character string specifying the species, for ppi calculation
#' @param subset a character vector of cell types to subset for, purely for output naming purpose
#' @param cond a character string specifying a condition, purely for output naming purpose
#' @return cds, a cell_data_set object with pseudotime values
#'
#' @noRd

orderCells <- function(cds, method, rootNodes=NULL, timePoint=NULL, timePointCol=NULL, species="mouse", subset=NULL, cond=NULL,outputDir="."){

  figure_title <- method

  if (method=="interactive"){

    cds <- monocle3::order_cells(cds)

  } else if (method == "rootNodes"){

    try(if(is.null(rootNodes) & is.null(timePoint) & is.null(timePointCol)) stop("No root node related information is provided. rootNodes cannot be a chosen method."))
    try(if((!is.null(timePoint) & is.null(timePointCol)) | (is.null(timePoint) & !is.null(timePointCol))) stop("Time point and time point column must be both specified."))

    if (is.null(rootNodes)){                    # based on time point information
      cds <- monocle3::order_cells(cds, root_pr_nodes=getRootPrincipalNodes(cds, timePoint, timePointCol))
      figure_title <- paste0(figure_title,"_",timePoint,sep="")
    } else {                                    # manual pick principal nodes
      cds <- monocle3::order_cells(cds, root_pr_nodes=rootNodes)
      figure_title <- paste0(figure_title,"_",rootNodes,sep="")
    }

  } else if (method == "potency"){

    cds <- getPotency(cds, species)

    # pick the cell with the highest potency score
    potency_df <- as.data.frame(cds$potency)
    rownames(potency_df) <- cds@colData@rownames
    colnames(potency_df)[1] <- "potency"
    max_potency <- which.max(potency_df$potency)
    root_cell <- rownames(potency_df)[max_potency]

    # define that cell as the root cell
    cds <- monocle3::order_cells(cds, root_cells = root_cell)

    # overlay potency scores
    png(filename = file.path(outputDir,"images","pseudotime",paste0('umap_potency_',figure_title,ifelse(!is.null(subset) ,paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond),"",paste0("_",cond)),'.png',sep="")), width = 2000*1.4, height = 1500*1.5, res = 400)
    P1 <- monocle3::plot_cells(cds, label_principal_points = TRUE, color_cells_by = "potency") +
      ggplot2::theme_void() +
      ggplot2::guides(colour = ggplot2::guide_colorbar("Potency"))
    print(P1)
    dev.off()

  }

  saveRDS(cds, file=file.path(outputDir,"rds",paste0("cds", ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond),"",paste0("_",cond)), ".rds",sep="")))

  png(filename = file.path(outputDir,"images","pseudotime",paste0("umap_pseudotime_",figure_title,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond),"",paste0("_",cond)),".png",sep="")), width = 2000*1.4, height = 1500*1.5, res = 400)
  p1 <- monocle3::plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(paste0("UMAP by pseudotime","\n", "method=",method,sep=""))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p1)
  dev.off()

  return(cds)
}


#' Computes the potency of cells based on their gene expression profiles
#' and a protein-protein interaction network.
#'
#' @param cds, a cell_data_set object containing gene expression profiles.
#' @param specie A character string specifying the species of the protein-protein interaction
#' network. Currently, only 'human' and 'mouse' are supported. Default is 'mouse'.
#' @param ppiThreshold An integer specifying the minimum confidence score for protein-protein
#' interactions to be included in the network. Default is 600.
#' @return cds, a cell_data_set object with a new slot 'potency', containing the potency scores of each cell.
#' @noRd

getPotency <- function(cds, specie = 'mouse', ppiThreshold = 600){
  # Getting species
  #specie <- c(human = 9606, mouse = 10090, ferret = 9669)[specie]
  specie <- c(human = 9606, mouse = 10090)[specie]

  # Downloading files
  if(!file.exists(paste0(specie,'PPI.txt.gz'))){
    utils::download.file(paste0('https://stringdb-static.org/download/protein.links.v11.5/',specie,'.protein.links.v11.5.txt.gz'), paste0(specie,'PPI.txt.gz'), method = 'wget',quiet = TRUE)
  }
  if(!file.exists(paste0(specie,'P2G.txt.gz'))){
    utils::download.file(paste0('https://stringdb-static.org/download/protein.info.v11.5/',specie,'.protein.info.v11.5.txt.gz'), paste0(specie,'P2G.txt.gz'), method = 'wget',quiet = TRUE)
  }

  # Reading files
  PPI <- read.table(paste0(specie,'PPI.txt.gz'), header = TRUE) # a data frame of PPI information
  P2G <- read.table(paste0(specie,'P2G.txt.gz'), sep = '\t')    # a data frame of each protein's metadata

  # Formating inputs: match Ensembl annotations to official gene symbols in PPI by using P2G
  geneDictionary <- P2G[,2]
  names(geneDictionary) <- P2G[,1]
  PPI <- PPI[PPI[,3] >= ppiThreshold,]
  PPI[,1] <- geneDictionary[PPI[,1]]
  PPI[,2] <- geneDictionary[PPI[,2]]
  PPI <- na.omit(PPI)
  PPI <- igraph::graph_from_data_frame(PPI, directed = FALSE)[,]  # creates a graph object from a data frame, now PPI is a sparse matrix

  # Computing potency: using the expression matrix and a protein-protein interaction (PPI) network
  # NOTE: this step requires **LARGE** memory because SCENT is trying to convert all count data from sparse to dense matrix
  # a count matrix with size ~0.4 GB would result in a dense matrix of size ~1.8 GB
  cds[["potency"]] <- SCENT::CompCCAT(as.matrix(cds@assays@data@listData$counts), PPI)
  return(cds)
}


#' Picks the node that is most heavily occupied by early timepoint cells and returns that as the root.
#'
#' @param cds a cell_data_set object with trajectory graph.
#' @param timePoint a character string specifying the earliest timepoint
#' @param timePointCol a character string specifying the metadata column with timepoint information
#' @return a vector of principle nodes to serve as roots
#' @noRd

getRootPrincipalNodes <- function(cds, timePoint, timePointCol){
  cell_ids <- which(SummarizedExperiment::colData(cds)[, timePointCol] == timePoint)

  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]

  return(root_pr_nodes)
}


#' Performs DEG analysis using regression method.
#'
#' This step takes some time to run if number of genes is large.
#'
#' @param cds a cell_data_set object.
#' @param model a character string of the name of a metadata column to conduct DEG analysis against.
#' @param batch a character string of the name of a metadata column with batch information
#' @param distribution a character string of the distribution to use in linear model
#' @param top_gene an integer number of top differentially expressed genes to plot
#' @param subset a character vector of cell types to subset for, for result naming purpose
#' @param cond1 a character string specifying a condition, for result naming purpose
#' @param cond2 a character string specifying a condition, for result naming purpose
#' @param outputDir outputDir
#' @param cores cores
#' @noRd

regressionAnalysis <- function(cds, model, batch, distribution, top_gene=10, subset=NULL, cond1=NULL, cond2=NULL,outputDir=".",cores=4){
  
  add_title <- ""
  if (!is.null(cond1) & !is.null(cond2)){
    add_title <- "+trajectory"
  }
  
  reduced_model <- monocle3::fit_models(cds, model_formula_str = paste0("~",model,sep=""), expression_family=distribution,cores=cores)
  fit_coefs <- monocle3::coefficient_table(reduced_model)
  
  fit_coefs <- fit_coefs[,!names(fit_coefs) %in% c("model", "model_summary")]
  deg_file_name <- file.path(outputDir,"csv",paste0("monocleDEG_by_",model,add_title,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond1),"",paste0("_",cond1)),ifelse(is.null(cond2),"",paste0("_",cond2)),".csv"))
  write.csv(fit_coefs, file=deg_file_name)
  fit_coefs_sig <- fit_coefs %>% dplyr::filter(term != "(Intercept)") %>% dplyr::filter (q_value < 0.05) %>% dplyr::select(gene_short_name, term, p_value, q_value, estimate)
  deg_sigfile_name <- file.path(outputDir,"csv",paste0("monocleDEG_significant_by_",model,add_title,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond1),"",paste0("_",cond1)),ifelse(is.null(cond2),"",paste0("_",cond2)),".csv"))
  write.csv(fit_coefs_sig, file=deg_sigfile_name)
  
  if (!is.null(batch)){
    full_model <- monocle3::fit_models(cds, model_formula_str = paste0("~",model,"+",batch, sep=""), expression_family=distribution)
    fit_coefs_full <- monocle3::coefficient_table(full_model)
    fit_coefs_full <- fit_coefs_full[,!names(fit_coefs_full) %in% c("model", "model_summary")]
    
    deg_file_name <- file.path(outputDir,"csv",paste0("monocleDEG_by_",model,"+",batch,add_title,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond1),"",paste0("_",cond1)),ifelse(is.null(cond2),"",paste0("_",cond2)),".csv"))
    write.csv(fit_coefs_full, file=deg_file_name)
    fit_coefs_sig_full <- fit_coefs_full %>% dplyr::filter(term != "(Intercept)") %>% dplyr::filter (q_value < 0.05) %>% dplyr::select(gene_short_name, term, q_value, estimate)
    deg_sigfile_name <- file.path(outputDir,"csv",paste0("monocleDEG_significant_by_",model,"+",batch,add_title,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond1),"",paste0("_",cond1)),ifelse(is.null(cond2),"",paste0("_",cond2)),".csv"))
    write.csv(fit_coefs_sig_full, file=deg_sigfile_name)
    
    model_compare <- monocle3::compare_models(full_model, reduced_model) %>% dplyr::select(gene_short_name, q_value)
    deg_modelfile_name <- file.path(outputDir,"csv",paste0("monocleDEG_modelCompare_",model,"_VS_",model,"+",batch,add_title,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond1),"",paste0("_",cond1)),ifelse(is.null(cond2),"",paste0("_",cond2)),".csv"))
    write.csv(model_compare, file=deg_modelfile_name)
    
    if (nrow(dplyr::filter(model_compare, q_value < 0.05)) >= (nrow(model_compare)/2)){
      warning("More than half of the genes' likelihood ratio tests are significant, indicating that there are substantial batch effects in the data. Consider using DEGs found by incorporating batch term!")
    }
    
    fit_coefs_sig_full <- fit_coefs_sig_full[order(fit_coefs_sig_full$q_value),]
    top_diff_genes <- head(fit_coefs_sig_full,n=top_gene)
    #top_diff_genes <- dplyr::slice_min(fit_coefs_sig_full, q_value, n=top_gene)
    
    if(nrow(top_diff_genes) > 0)
    {
      violinplot_filename<-file.path(outputDir,"images","DEG",paste0("significant_by_",model,"+",batch,add_title,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond1),"",paste0("_",cond1)),ifelse(is.null(cond2),"",paste0("_",cond2)),"_violinPlot",".png",sep=""))
      #png(violinplot_filename, width = 2000*5, height = 1500*7, res = 400)
      png(violinplot_filename, width = 1000*sqrt(top_gene)*2, height = 875*sqrt(top_gene), res = 300)
      p1 <- monocle3::plot_genes_violin(cds[SummarizedExperiment::rowData(cds)$gene_short_name %in% (top_diff_genes$gene_short_name), ], group_cells_by=model, ncol=ceiling(sqrt(top_gene))) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
      p2 <- monocle3::plot_genes_violin(cds[SummarizedExperiment::rowData(cds)$gene_short_name %in% (top_diff_genes$gene_short_name), ], group_cells_by=batch, ncol=ceiling(sqrt(top_gene))) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
      pfinal<-p1+p2
      print(pfinal)
      dev.off()
      
      featureplot_filename<-file.path(outputDir,"images","DEG",paste0("significant_by_",model,"+",batch,add_title,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond1),"",paste0("_",cond1)),ifelse(is.null(cond2),"",paste0("_",cond2)),"_featurePlot",".png",sep=""))
      png(featureplot_filename, width = 2000*1.4, height = 1500*1.5, res = 400)
      p2 <- monocle3::plot_cells(cds, genes=unique(top_diff_genes$gene_short_name),
                                 show_trajectory_graph=FALSE,
                                 label_cell_groups=FALSE,
                                 label_leaves=FALSE)
      print(p2)
      dev.off()
    }
  }
  
  fit_coefs_sig <- fit_coefs_sig[order(fit_coefs_sig$q_value),]
  top_diff_genes <- head(fit_coefs_sig,n=top_gene)
  #top_diff_genes <- dplyr::slice_min(fit_coefs_sig, q_value, n=top_gene)
  
  #print("Top genes DEGs number in regression function")
  #print(nrow(top_diff_genes))
  #print(top_gene)
  
  violinplot_filename<-file.path(outputDir,"images","DEG",paste0("significant_by_",model,add_title,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond1),"",paste0("_",cond1)),ifelse(is.null(cond2),"",paste0("_",cond2)),"_violinPlot",".png",sep=""))
  #png(violinplot_filename, width = 2000*1.5, height = 1500*3, res = 300)
  
  # png(violinplot_filename, width = 900*sqrt(top_gene), height = 875*sqrt(top_gene), res = 300)
  # p4 <- plot_genes_violin(cds[rowData(cds)$gene_short_name %in% (top_diff_genes$gene_short_name), ], group_cells_by=model, ncol=ceiling(sqrt(top_gene))) +
  #       theme(axis.text.x=element_text(angle=45, hjust=1))
  
  # print(p4)
  # dev.off()
  
  if(nrow(top_diff_genes) > 0)
  {
    featureplot_filename<-file.path(outputDir,"images","DEG",paste0("significant_by_",model,add_title,ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(is.null(cond1),"",paste0("_",cond1)),ifelse(is.null(cond2),"",paste0("_",cond2)),"_featurePlot",".png",sep=""))
    #png(featureplot_filename, width = 2000*1.4, height = 1500*1.5, res = 400)
    p5 <- monocle3::plot_cells(cds, genes=unique(top_diff_genes$gene_short_name),
                               show_trajectory_graph=FALSE,
                               label_cell_groups=FALSE,
                               label_leaves=FALSE)
    ggplot2::ggsave(file=featureplot_filename, plot = p5, width = 10, height = 8)
    
    O <- lapply(top_diff_genes$gene_short_name, function(gene){
      p4 <- monocle3::plot_genes_violin(cds[SummarizedExperiment::rowData(cds)$gene_short_name %in% (gene), ], group_cells_by=model, ncol=1) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)) + ggplot2::labs(x="")
    })
    print_violinPlot(O,violinplot_filename,top_gene)
  }
}

print_violinPlot<-function(O,filename,top_gene)
{
  #print("In Violin Plot function")
  png(filename, width = 700*sqrt(top_gene)*2, height = 875*sqrt(top_gene), res = 300)
  print(patchwork::wrap_plots(O,ncol = ceiling(sqrt(top_gene))))
  dev.off()
}

#' Finding genes that change as a function of pseudotime using graph auto-correlation method.
#'
#' @param cds a cell_data_set object with learned trajectory.
#' @param conditions_all conditions_all
#' @param colData_name a character string of the name of metadata that has the conditions.
#' @param top_gene an integer number of top differentially expressed genes to plot.
#' @param subset a character vector of cell types to subset for, for result naming purpose.
#' @param deg_method deg_method
#' @param batch batch
#' @param outputDir outputDir
#' @param cores cores
#'
#' @noRd

graphAutoCorrelation <- function(cds,conditions_all,colData_name, top_gene, subset=NULL,deg_method="quasipoisson",batch=NULL,outputDir=".",cores=6){

  # test whether cells at similar positions on the trajectory have correlated expression
  pr_graph_test_res <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = cores)

  print("Graph_test function finished successfully. Continue to csv saving...")
  #ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),"")
  write.csv(pr_graph_test_res, file=file.path(outputDir,"csv",paste0("monocleDEG_by_","trajectory",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),".csv", sep="")))

  # significant DEGs
  pr_graph_test_sig <- subset(pr_graph_test_res, q_value < 0.05 & morans_I > 0.1 & status == 'OK')
  write.csv(pr_graph_test_sig, file=file.path(outputDir,"csv",paste0("monocleDEG_significant_by_","trajectory",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),".csv", sep="")))

  # morans_I statistic: 0 indicates no effect, +1 indicates positive autocorrelation and nearby cells have similar values of a gene's expression
  pr_graph_test_sig_top <- pr_graph_test_sig[order(pr_graph_test_sig$morans_I,decreasing=TRUE),]
  DEG_ids <- row.names(pr_graph_test_sig_top[1:top_gene, ])

  # feature plot for top DEGs by trajectory
  png(filename = file.path(outputDir,"images","DEG",paste0("significant_by","_trajectory_","gene_featurePlot",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),".png",sep="")), width = 2000*1.4, height = 1500*1.5, res = 300)
  p1 <- monocle3::plot_cells(cds, genes=DEG_ids,
                   show_trajectory_graph=FALSE,
                   label_cell_groups=FALSE,
                   label_leaves=FALSE)
  print(p1)
  dev.off()

  print("CSV saving and feature plot plotting finished successfully. Continue to module finding...")

  # below few lines are repetitive and are just here to get around a bug caused by file too large (ex. object the size of p147)
  rm(pr_graph_test_sig)
  rm(pr_graph_test_sig_top)
  pr_graph_test_sig <- read.csv(file=file.path(outputDir,"csv",paste0("monocleDEG_significant_by_","trajectory",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),".csv", sep="")))
  rownames(pr_graph_test_sig) <- pr_graph_test_sig$X
  pr_graph_test_sig <- pr_graph_test_sig[,-1]

  print("Module finding start...")

  # organize significant trajectory-variable genes into modules
  # this function will give different results with every run if seed not set, see https://github.com/cole-trapnell-lab/monocle3/issues/494

  gene_module_df <- monocle3::find_gene_modules(cds[pr_graph_test_sig$gene_short_name, ], resolution=c(10^seq(-6,-1)), random_seed=1, cores=cores)
  #gene_module_df <- find_gene_modules_leiden(cds[pr_graph_test_sig$gene_short_name, ], resolution=c(10^seq(-6,-1)), random_seed=1, cores=cores)
  print("Module found successfully.")
  write.csv(gene_module_df, file=file.path(outputDir,"csv",paste0("monocleDEG_by_trajectory","_moduleInformation",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),".csv", sep="")))

  print("Module information saved to csv successfully.")

  cell_group_df <- tibble::tibble(cell=row.names(SummarizedExperiment::colData(cds)), cell_group=SummarizedExperiment::colData(cds)$cell.type)
  agg_mat <- monocle3::aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

  png(filename = file.path(outputDir,"images","DEG",paste0("significant_by_trajectory","_moduleHeatmap",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),".png",sep="")), width = 2000*1.4, height = 1800*1.5, res = 400)
  pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
  dev.off()

  png(filename = file.path(outputDir,"images","DEG",paste0("significant_by_trajectory","_moduleFeaturePlot",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),".png",sep="")), width = 2000*1.4, height = 1500*1.5, res = 400)
  p2 <- monocle3::plot_cells(cds,
                   genes=gene_module_df,
                   label_cell_groups=FALSE,
                   show_trajectory_graph=FALSE)
  print(p2)
  dev.off()

  # plot top genes' dynamics as a function of pseudotime
  png(filename = file.path(outputDir,"images","DEG",paste0("significant_by_trajectory_genesInPseudotime",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),".png",sep="")), width = 1000*sqrt(top_gene), height = 875*sqrt(top_gene), res = 300)
  #p3 <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name %in% DEG_ids , ], nrow = ceiling(top_gene/4), ncol = 4,color_cells_by="cell.type", min_expr=0.5)
  p3 <- monocle3::plot_genes_in_pseudotime(cds[SummarizedExperiment::rowData(cds)$gene_short_name %in% DEG_ids , ], ncol = ceiling(sqrt(top_gene)),color_cells_by="cell.type", min_expr=0.5)
  print(p3)
  dev.off()

  # Subset conditions pairwise and use identified genes that vary significantly across trajectory to do regression analysis
  seq <- length(conditions_all)
  if (seq > 1){
    for (i in 1:(seq-1)){
      for (j in (i+1):seq){
        # subset two conditions
        cds_for_compare <- cds[pr_graph_test_sig$gene_short_name, cds[[colData_name]] == c(conditions_all[i], conditions_all[j])]

        # run regression analysis
        regressionAnalysis(cds=cds_for_compare , model=colData_name, batch=batch, distribution=deg_method, top_gene=top_gene, subset=subset, cond1=conditions_all[i], cond2=conditions_all[j],outputDir=outputDir,cores=cores)
        # read in csv of significant differential gene along the trajectory AND between two conditions
        gene_csv <- read.csv(file=file.path(outputDir,"csv",paste0("monocleDEG_significant_by_",colData_name,ifelse(is.null(batch),"",paste0("+",batch)),"+trajectory",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""), "_",conditions_all[i],"_",conditions_all[j],".csv", sep="")))
        top_diff_genes <- dplyr::slice_min(gene_csv, q_value, n=top_gene)
        # plot expression dynamics as a function of pseudotime for top differential gene along the trajectory AND between two conditions
        png(filename = file.path(outputDir,"images","DEG",paste0("significant_by_",colData_name,ifelse(is.null(batch),"",paste0("+",batch)),"+trajectory",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""), "_",conditions_all[i],"_",conditions_all[j],"_genesInPseudotime",".png",sep="")), width = 1000*sqrt(top_gene), height = 875*sqrt(top_gene), res = 300)
        #p4 <- plot_genes_in_pseudotime(cds_for_compare[rowData(cds_for_compare)$gene_short_name %in% (top_diff_genes$gene_short_name) , ], nrow = ceiling(top_gene/4), ncol = 4,color_cells_by=colData_name, min_expr=0.5)
        p4 <-monocle3::plot_genes_in_pseudotime(cds_for_compare[SummarizedExperiment::rowData(cds_for_compare)$gene_short_name %in% (top_diff_genes$gene_short_name) , ], ncol = ceiling(sqrt(top_gene)),color_cells_by=colData_name, min_expr=0.5)
        print(p4)
        dev.off()
      }
    }
  }
}


#' Plot cell/cell type distribution density and histogram.
#'
#' @param cds a cell_data_set object.
#' @param colData_name a character string of the name of metadata that has the conditions.
#' @param subset a character vector of cell types to subset for, for result naming purpose.
#' @param cond a character string specifying a condition, for result naming purpose.
#' @noRd

cellTypeDistribution <- function(cds, colData_name, subset=NULL, cond=NULL,outputDir="."){

  # extract pseudotime, conditions, and cell type information into a single dataframe
  pseud_cells <- data.frame(
    pseudotime = cds@principal_graph_aux@listData$UMAP$pseudotime,
    cell_type = cds[["cell.type"]]
  )
  if (!is.null(colData_name)){
    pseud_cells$Condition = cds[[colData_name]]
  } else {
    pseud_cells$Condition = "ALL"
  }

  png(filename = file.path(outputDir,"images","cellDistribution",paste0("cell_distribution_density",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(!is.null(cond),paste0("_",cond),""),".png",sep="")), width = 2000*1.4, height = 1500*1.5, res = 400)
  p1 <- ggplot2::ggplot(pseud_cells, ggplot2::aes(x = pseudotime, color = Condition, fill = Condition)) +
    ggplot2::geom_density(alpha = .2) +
    ggplot2::facet_wrap(~Condition, ncol = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Pseudotime", y = "Density", title = "Distribution density of cells along the trajectory")
  print(p1)
  dev.off()

  # how many cell types, use this number for scaling figure sizes
  # "/20" is used here as an approximate for tuning figure widths, there is more flexibility in width since the column number is set to 2.
  # "/12" is used here as an approximate for tuning figure heights, this parameter might need to be adjusted with more tests.
  celltype_count <- length(unique(cds$cell.type))
  png(filename = file.path(outputDir,"images","cellDistribution",paste0("celltype_distribution_density",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(!is.null(cond),paste0("_",cond),""),".png",sep="")), width = 2000*1.4*ceiling(celltype_count/20), height = 1500*1.6*ceiling(celltype_count/12), res = 400)
  p2 <- ggplot2::ggplot(pseud_cells, ggplot2::aes(x = pseudotime, color = Condition, fill = Condition)) +
    ggplot2::geom_density(alpha = .2) +
    ggplot2::facet_wrap(~cell_type, ncol = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Pseudotime", y = "Density", title = "Distribution density of cell types along the trajectory")
  print(p2)
  dev.off()

  png(filename = file.path(outputDir,"images","cellDistribution",paste0("celltype_distribution_histogram",ifelse(!is.null(subset),paste0("_",paste(subset,collapse = '_')),""),ifelse(!is.null(cond),paste0("_",cond),""),".png",sep="")), width = 2000*1.4*ceiling(celltype_count/20), height = 1500*1.6*ceiling(celltype_count/12), res = 400)
  p3 <- ggplot2::ggplot(pseud_cells, ggplot2::aes(x = pseudotime, color = Condition, fill = Condition)) +
    ggplot2::geom_histogram(
      alpha = 0.2,
      position = "identity", bins = 50
    ) +
    ggplot2::facet_wrap(~cell_type, ncol = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Pseudotime", y = "Count", title = "Distribution of cell types along the trajectory")
  print(p3)
  dev.off()

}
