# load needed library
library(symphony)

#' Read in a h5ad file and transfer the AnnData object into a Seurat object.
#'
#' @description This function uses the library(anndata) to transfer the normalized counts, metadata,
#' and lower dimensional embeddings (pca and umap) from a h5ad file into a Seurat object.
#' More data can be transferred between the two formats, and this function can be expanded
#' to accommodate that accordingly.
#'
#'
#' @param h5ad_file character string of full path to the h5ad file.
#' @return a Seurat object.
#'
#' @export
h5adToSeurat <- function(h5ad_file, annotation_column){
  miniconda_path <- "~/.local/share/r-miniconda"
  # Check if Miniconda is already installed at the specified path
  if (!file.exists(miniconda_path)) {
    reticulate::install_miniconda()
  }

  library(anndata)
  print("Convert from AnnData to Seurat...")
  adata <- read_h5ad(h5ad_file)
  # here normalized counts are used to create the Seurat object, raw counts are typically in adata$raw$X
  data <- CreateSeuratObject(counts = t(adata$X), meta.data = adata$obs)

  # transfer lower dimensional embeddings
  if ('X_pca' %in% names(adata$obsm)){
    umap_data <- adata$obsm$X_pca
    data <- AddDimReduc(data, reduction = umap_data, key = "pca")
  }
  if ('X_umap' %in% names(adata$obsm)){
    umap_data <- adata$obsm$X_umap
    data <- AddDimReduc(data, reduction = umap_data, key = "umap")
  }
  print("AnnData object successfully converted.")

  # save transferred Seurat object
  file_sub <- sub("\\.h5ad$", "", h5ad_file)
  saveRDS(data, file=paste0(file_sub,"_h5adTransferred.rds",sep=""))

  return(data)
}


#' Transfer cell labels from one Seurat object to another
#'
#' @description This function transfers cell labels from a reference Seurat object to a query Seurat object.
#' The labels are transferred based on k-NN classification of the query cells in the reference UMAP space.
#'
#' @param X Seurat object used as reference for the transfer of cell labels
#' @param Y Seurat object in which the cell labels will be transferred
#' @param varToHarmonize a character string specifying a variable to harmonize data over
#' @param transferCoordinates Logical. If TRUE, the UMAP coordinates of the query object will be
#' updated with the ones from the reference object.
#' @return A Seurat object with the transferred cell labels and UMAP coordinates (if transferCoordinates = TRUE)
#'
#' @export

doTransferLabel <- function(X, Y, varToHarmonize = NULL, transferCoordinates = FALSE){
  # Getting metadata from the reference
  MD <- data.frame(B = X@meta.data[,varToHarmonize], L = Idents(X))

  # Removing previous instances of the function
  if (file.exists('.UMAP')){
    file.remove('.UMAP')
  }

  # Adjust value so it will not generate error in either case when calling buildReferrence()
  varToHarmonize <- ifelse(is.null(varToHarmonize), NULL, 'B')

  # Generating a Symphony Reference
  S <- symphony::buildReference(X@assays$RNA@counts, MD,
                      vars = varToHarmonize,
                      do_umap = TRUE,
                      verbose = TRUE,
                      d = 30,
                      save_uwot_path = '.UMAP')

  # Replacing embeeding by the one in the Seurat object
  if(transferCoordinates){
    S$umap$embedding[,1] <- X@reductions$umap@cell.embeddings[,1]
    S$umap$embedding[,2] <- X@reductions$umap@cell.embeddings[,2]
    U <- uwot::load_uwot('.UMAP')
    U$embedding[,1] <- X@reductions$umap@cell.embeddings[,1]
    U$embedding[,2] <- X@reductions$umap@cell.embeddings[,2]
    file.remove('.UMAP')
    U <- uwot::save_uwot(U, file = '.UMAP', verbose = FALSE)
  }

  # Querying the new data into the generated reference
  set.seed(1)
  Q <- symphony::mapQuery(Y@assays$RNA@counts, Y@meta.data, ref_obj = S)

  # Transfering Labels
  Q <- symphony::knnPredict(Q, S, S$meta_data$L, k = 5)

  # Adding new metadata to the Y object
  Y$transferedCellType <- Q$meta_data$cell_type_pred_knn
  Y$labelProbability <- Q$meta_data$cell_type_pred_knn_prob
  Idents(Y) <- Y$transferedCellType

  if(transferCoordinates){
    Y@reductions$umap@cell.embeddings[,1] <- Q$umap[,1]
    Y@reductions$umap@cell.embeddings[,2] <- Q$umap[,2]
  }

  # Returning annotated object
  return(Y)
}
