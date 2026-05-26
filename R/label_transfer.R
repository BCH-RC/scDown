# load needed library
library(symphony)

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

  # Update object structure
  X <- SeuratObject::UpdateSeuratObject(X)
  Y <- SeuratObject::UpdateSeuratObject(X)

  # Adjust value so it will not generate error in either case when calling buildReferrence()
  varToHarmonize <- if (is.null(varToHarmonize)) NULL else 'B'

  # Generating a Symphony Reference
  cat("Building reference...\n")
  S <- symphony::buildReference(LayerData(X, layer="counts", assay="RNA"),
                                MD,
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
  cat("Mapping query...\n")
  set.seed(1)
  Q <- symphony::mapQuery(LayerData(Y, layer="counts", assay="RNA"),
                          Y@meta.data, 
                          ref_obj = S)

  # Transfering Labels
  cat("Predicting kNN...\n")
  Q <- symphony::knnPredict(Q, S, S$meta_data$L, k = 5)

  # Adding new metadata to the Y object
  Y$transferredCellType <- Q$meta_data$cell_type_pred_knn
  Y$labelProbability <- Q$meta_data$cell_type_pred_knn_prob
  Idents(Y) <- Y$transferredCellType

  if(transferCoordinates){
    Y@reductions$umap@cell.embeddings[,1] <- Q$umap[,1]
    Y@reductions$umap@cell.embeddings[,2] <- Q$umap[,2]
  }
  cat("Label transfer completed.\n")

  # Returning annotated object
  return(Y)
}
