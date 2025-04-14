#' Function to run the monocle3 pipeline
#'
#' This function runs the moncole3 workflow to infer the development trajectory of the cells using the
#' single-cell RNASeq data. In this function, users can infer the trajectory of the whole data as well as
#' a subset of data and also identify the differential genes across the trajectory as well as between the
#' conditions.
#'
#' @param seurat_obj Seurat object containing the scRNASeq data (Required)
#' @param species A character variable specifying the species (Required, either human or mouse)
#' @param nDim A numeric variable specifying the PCA dimensions to be used for processing
#' @param conditions A character vector specifying the different conditions
#' @param annotation_column A character variable specifying the metdata column name of annotations. (Default: uses Idents)
#' @param group_column A character variable specifying the metdata column name containing conditions in the seurat object
#' @param transferUMAP A boolean variable specifying whether to transfer UMAP coordinates or not from Seurat object
#' @param rootNode_method A character variable specifying the method to identify root node (Either potency or rootNodes)
#' @param rootNode A character variable specifying the principal point on the trajectory, in the format of "Y_<node number>".
#' @param timePoint A character variable specifying the category inside a metadata to be used to automatically decide the root node.
#' @param timePoint_metadata A character variable specifying the metadata to be used to automatically decide the root node.
#' @param batch_metadata A character vector specifying the name of metadata column name with batch information to perform batch correction
#' @param celltype_groups A list of character vectors that represents a group of cell types to be subsetted and constructed trajectory by monocle
#' @param top_genes A numeric variable to set the number of top significant genes to visualize from DEG results
#' @param deg_method A character vector specifying the distribution used for regression method
#' @param metadata_deg_model A character vector specifying the model used for regression analysis, can be any term that exists in the metadata
#' @param graph_test A boolean variable to run graph auto-correlation to find differential genes along the trajectory
#' @param cores A numeric variable to set the number of cores to be used
#' @param output_dir A character vector specifying the output directory
#'
#' @return A list of monocle cell dataset objects
#'
#' @export
#'
#'

run_monocle3 <- function(seurat_obj,species,nDim=30,conditions=NULL,annotation_column=NULL,group_column=NULL,transferUMAP=TRUE,rootNode_method="potency",
                         rootNode = NULL, timePoint = NULL,timePoint_metadata=NULL,batch_metadata=NULL,celltype_groups=NULL,top_genes=10,
                         deg_method="quasipoisson",metadata_deg_model=NULL,graph_test=FALSE,cores=8,output_dir="."){

  `%dopar%` <- foreach::`%dopar%`
  check_required_variables(seurat_obj,species,output_dir,annotation_column,group_column)
  #checkmate::expect_class(seurat_obj,"Seurat",label="seurat_obj")
  checkmate::expect_choice(species,c("human","mouse"),label = "species")
  checkmate::expect_numeric(nDim, min.len = 1, max.len = 1, any.missing = FALSE,label="nDim")

  checkmate::expect_character(conditions, min.len = 1, any.missing = FALSE,label="conditions",null.ok = TRUE)
  if(checkmate::test_character(conditions, min.len = 1, any.missing = FALSE))
  {
    checkmate::expect_choice(group_column, colnames(seurat_obj@meta.data),label="group_column")
    #checkmate::expect_choice(conditions, names(table(seurat_obj@meta.data[,group_column])),label="group_column")
    if(!all(conditions %in% names(table(seurat_obj@meta.data[,group_column]))))
    {
      stop("The conditions should be the elements in the ",group_column," metadata. {",paste(names(table(seurat_obj@meta.data[,group_column])),collapse = ","),"}")
    }
  }

  library(dplyr)
  #checkmate::expect_choice(annotation_column, colnames(seurat_obj@meta.data),label="annotation_column",null.ok = TRUE)
  if(checkmate::test_character(annotation_column, min.len = 1, max.len = 1, any.missing = FALSE))
  {
    Seurat::Idents(seurat_obj) <- seurat_obj[[annotation_column]][,1]
  }
  checkmate::expect_flag(transferUMAP,label="transferUMAP")
  checkmate::expect_choice(rootNode_method,c("potency","rootNodes"),label = "rootNode_method")
  ##Use all conditions for potency method
  passed_conditions <- conditions
  
  if(rootNode_method == "rootNodes")
  {
    # if "rootNodes", need to supply either @rootNode, or a combination of @timePoint and @timePoint_metadata
    checkmate::expect_character(rootNode, min.len = 1, max.len = 1, any.missing = TRUE,label="rootNode",null.ok = TRUE)
    if(!checkmate::test_character(rootNode, min.len = 1, any.missing = FALSE))
    {
      # checkmate::expect_character(timePoint, min.len = 1, max.len = 1, any.missing = TRUE,label="timePoint",null.ok = TRUE)
      # if(checkmate::test_character(timePoint, min.len = 1, any.missing = FALSE))
      # {
      #   checkmate::expect_choice(timePoint_metadata,colnames(seurat_obj@meta.data),label="timePoint_metadata",null.ok = TRUE)
      # }
      checkmate::expect_choice(timePoint_metadata,colnames(seurat_obj@meta.data),label="timePoint_metadata",null.ok = TRUE)
      checkmate::expect_choice(timePoint,as.character(unique(seurat_obj@meta.data[,timePoint_metadata])),label="timePoint",null.ok = TRUE)
      if(checkmate::test_character(conditions, min.len = 1, any.missing = FALSE))
      {
        seurat_obj_subset <- seurat_obj[,row.names(seurat_obj@meta.data %>% dplyr::filter(seurat_obj@meta.data[,timePoint_metadata] == timePoint))]
        condList <- lapply(conditions, function(cond){  
          t<-table(seurat_obj_subset@meta.data[,group_column])
          if(t[cond] < 2)
          {
            msg<-paste0("There is less than 2 cells for the ",timePoint," timepoint for condition ",cond,". Trajectory analysis will not be carried out for this condition separately.")
            print(msg)
            return(NULL)
          }else{
            return(cond)
          }
        })
        passed_conditions<-condList[!sapply(condList, is.null)]
      }
    }
    if(!checkmate::test_character(rootNode, min.len = 1, any.missing = FALSE) && !checkmate::test_character(timePoint, min.len = 1, any.missing = FALSE))
    {
      stop("For rootNodes method, user need to supply either rootNode, or a combination of timePoint and timePoint_metadata")
    }
  }


  checkmate::expect_choice(batch_metadata, colnames(seurat_obj@meta.data),label="batch_metadata",null.ok = TRUE)
  checkmate::expect_class(celltype_groups,"list",label="celltype_groups",null.ok = TRUE)

  checkmate::expect_numeric(top_genes, min.len = 1, max.len = 1, any.missing = FALSE,label="top_genes")
  checkmate::expect_choice(deg_method,c("quasipoisson","negbinomial","poisson","binomial"),label="metadata_deg_model",null.ok = TRUE)

  checkmate::expect_choice(metadata_deg_model,colnames(seurat_obj@meta.data),label="metadata_deg_model",null.ok = TRUE)
  checkmate::expect_flag(graph_test,label="graph_test")
  checkmate::expect_numeric(cores, min.len = 1, max.len = 1, any.missing = FALSE,label="cores")
  #checkmate::expect_directory(output_dir,access="rw",label = "output_dir")

  ##########Create results subdirectories########
  subdirectories <- c(file.path("monocle","rds"),
                      file.path("monocle","csv"),
                      file.path("monocle","images","pseudotime"),
                      file.path("monocle","images","cellDistribution"),
                      file.path("monocle","images","DEG"))

  for(i in subdirectories){
    dir.create(file.path(output_dir,i), showWarnings = F, recursive = T)
  }
  output_dir <- file.path(output_dir,"monocle")

  ################################################
  subset <- ifelse(!is.null(celltype_groups) || length(celltype_groups) != 0, TRUE, FALSE)
  if(subset)
  {
    for (cell_group in celltype_groups){
      if(!all(cell_group %in% names(table(Seurat::Idents(seurat_obj)))))
      {
        stop("The cell types should be the elements in the data {",paste(names(table(Seurat::Idents(seurat_obj))),collapse = ","),"}")
      }
    }
  }

  regression <- ifelse(length(metadata_deg_model) != 0, TRUE, FALSE)
  conditions_to_compare <- ifelse(length(conditions) > 1, TRUE, FALSE)
  if (length(conditions) == 0){
    conditions <- c("ALL")
    group_column <- NULL
  }
  cds_by_group <- list()

  # Trajectory and pseudotime for the whole object
  cds <- getTrajectory(X=seurat_obj, nDim=nDim, batch=batch_metadata, transferUMAP=transferUMAP,outputDir=output_dir)
  cds <- orderCells(cds=cds, method=rootNode_method, rootNodes=rootNode, timePoint=timePoint, timePointCol=timePoint_metadata, species=species,outputDir=output_dir)
  cellTypeDistribution(cds=cds, colData_name=group_column,outputDir=output_dir)
  #cds_by_group <- c(cds_by_group, list(cds)) # add cds with completed trrajectory inference to list
  cds_by_group <- list("ALL"=cds) # add cds with completed trrajectory inference to list
  print("Trajectory and pseudotime for the entire object completed.")

  process_celltypegroup <- function(celltype_group)
  {
    cell_group <- celltype_groups[[celltype_group]]
    # trajectory and pseudotime for selected cell types
    cds.sub <- getTrajectory(X=seurat_obj, nDim=nDim, batch=batch_metadata, transferUMAP=transferUMAP, subset=unlist(cell_group),outputDir=output_dir)
    cds.sub <- orderCells(cds=cds.sub, method=rootNode_method, rootNodes=rootNode, timePoint=timePoint, timePointCol=timePoint_metadata, species=species, subset=unlist(cell_group),outputDir=output_dir)
    cellTypeDistribution(cds=cds.sub, colData_name=group_column, subset=unlist(cell_group),outputDir=output_dir)
    return(cds.sub)
    
  }
  # Individual trajectories and pseudotime by selected cell types
  if (subset){
    cl <- parallel::makeCluster(cores)
    cds_subset <- foreach::foreach(i=1:length(celltype_groups), .packages='monocle3') %dopar%
      process_celltypegroup(celltype_group=i)
    parallel::stopCluster(cl)
    names(cds_subset)<-unlist(lapply(celltype_groups,function(x){paste(unlist(x),collapse = '_')}))
    cds_by_group<-c(cds_by_group,cds_subset)
    print("Trajectory and pseudotime for all subsetted objects completed.")
  }

  # loop over each cell_data_set object (whole and subsetted, if any) and perform:
  #   - regression analysis on predefined variable (@DEG_models), if any
  #   - graph auto-correlation for finding genes that vary along pseudotime/trajectory
  #       - if there are multiple conditions, regression method is used to find genes that vary along pseudotime AND differentially expressed between any two conditions
  #   - split cds object by conditions and find trajectory and pseudotime for each condition
  print(length(cds_by_group))
  for (i in 1:length(cds_by_group))
  {

    cds.current <- cds_by_group[[i]]

    cell_type_use <- NULL
    if (i > 1){
      cell_type_use <- celltype_groups[[i-1]]
    }

    if (regression){
      for (model in metadata_deg_model){
        regressionAnalysis(cds=cds.current, model=model, batch=batch_metadata, distribution=deg_method, top_gene=top_genes, subset=cell_type_use,outputDir=output_dir,cores=cores)
      }
    }

    print(paste0("Regression analysis completed for cds object #",i," completed.", sep=""))

    if (graph_test){
      graphAutoCorrelation(cds=cds.current,conditions_all=conditions,colData_name=group_column, top_gene=top_genes, subset=cell_type_use,deg_method=deg_method,batch=batch_metadata,outputDir=output_dir,cores=cores)
      print(paste0("Graph autocorrelation analysis for cds object #",i," completed.", sep=""))
    }

    if (conditions_to_compare){
      runmonocle_percondition<-function(i){
        condition <- passed_conditions[i]
        cds.condition <- cds.current[ ,cds.current[[group_column]] == condition]
        cds.condition <- monocle3::cluster_cells(cds.condition,cluster_method = "louvain")
        cds.condition <- monocle3::learn_graph(cds.condition, use_partition = FALSE)

        # plot umap by cell types
        png(filename = file.path(output_dir,"images","pseudotime",paste0("umap_celltypes_","transferUMAP_",transferUMAP,transferUMAP,ifelse(!is.null(cell_type_use),paste0("_",paste(cell_type_use,collapse = '_')),""),"_",condition,".png",sep="")), width = 2000*1.4, height = 1500*1.5, res = 400)
        p4 <- monocle3::plot_cells(cds.condition, color_cells_by="cell.type", show_trajectory_graph=FALSE,label_groups_by_cluster=FALSE) +
          ggplot2::theme_void() +
          ggplot2::ggtitle(paste0("UMAP by cell types","\n", "transferUMAP=",transferUMAP, sep=""))+
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        print(p4)
        dev.off()

        # plot umap by trajectory
        png(filename = file.path(output_dir,"images","pseudotime",paste0("umap_trajectory_","transferUMAP_",transferUMAP,transferUMAP,ifelse(!is.null(cell_type_use),paste0("_",paste(cell_type_use,collapse = '_')),""),"_",condition,".png",sep="")), width = 2000*1.4, height = 1500*1.5, res = 400)
        p5 <- monocle3::plot_cells(cds.condition, label_principal_points = TRUE,  color_cells_by = "cell.type", label_cell_groups=FALSE) +
          ggplot2::theme_void() +
          ggplot2::guides(color= ggplot2::guide_legend("Cell Type", override.aes = list(size=5), ncol = 2)) +
          ggplot2::ggtitle(paste0("UMAP by trajectory","\n", "transferUMAP=",transferUMAP,sep=""))+
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        print(p5)
        dev.off()

        # order cells
        cds.condition <- orderCells(cds=cds.condition, method=rootNode_method, rootNodes=rootNode, timePoint=timePoint, timePointCol=timePoint_metadata, species=species, subset=cell_type_use, cond=condition,outputDir=output_dir)

        # cell/cell type distribution
        cellTypeDistribution(cds=cds.condition, colData_name=group_column, subset=cell_type_use, cond=condition,outputDir=output_dir)
      }
      
      cl <- parallel::makeCluster(cores)
      cds_subset <- foreach::foreach(i=1:length(passed_conditions), .packages=c('monocle3','magrittr')) %dopar%
        runmonocle_percondition(i=i)
      parallel::stopCluster(cl)
      
      print(paste0("Trajectory and pseudotime for all conditions subsetted from cds object #",i," completed.", sep=""))
    }
    #return(cds.current)
  }
  return(cds_by_group)
}

#remove.packages("scDown")
#devtools::check()
#devtools::build()
#install.packages("../scDown_0.1.0.tar.gz", repos = NULL, type="source")
#library("scDown")
#seurat<-readRDS("../cbl_sub_opc.rds")
#l<-run_monocle3(seurat,"human")
#vignette("scDown_monocle",package = "scDown")
#?run_monocle3

