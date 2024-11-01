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
