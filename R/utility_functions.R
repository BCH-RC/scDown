#' Read in a h5ad file and convert to a Seurat object.
#'
#' @description This function uses the library(zellkonverter) to convert the normalized counts, metadata,
#' lower dimensional embeddings (pca and umap), and spliced and unspliced counts if available from a h5ad 
#' file to a Seurat object.
#'
#' @param h5ad_file character string of full path to the h5ad file.
#' @param annotation_column character variable specifying the metdata column name of cell type annotations. 
#' @return a Seurat object.
#'
#' @export
#'
h5adToSeurat <- function(h5ad_file, annotation_column=NULL){
  # Check if python is available
  if (!reticulate::py_available(initialize = FALSE)) {
    reticulate::install_miniconda()
  }

  # Convert .h5ad with spliced and unspliced data to Seurat object
  library(zellkonverter)
  # Read the .h5ad file
  ad <- readH5AD(h5ad_file)
  ad
  file_sub <- sub("\\.h5ad$", "", h5ad_file)

  # Convert the main count matrix X to a Seurat object
  X <- as.Seurat(ad, counts = "X", data = NULL)
  X@assays[["RNA"]]<-X@assays$originalexp

  # Convert lower dimensional embeddings (pca and umap)
  if ('X_pca' %in% names(X@reductions)){
    X@reductions[["pca"]]<-X@reductions$X_pca
    colnames(X@reductions[["pca"]]@cell.embeddings)<-gsub("Xpca","PCA",colnames(X@reductions[["pca"]]@cell.embeddings))
  }
  if ('X_umap' %in% names(X@reductions)){
    X@reductions[["umap"]]<-X@reductions$X_umap
    colnames(X@reductions[["umap"]]@cell.embeddings)<-gsub("Xumap","UMAP",colnames(X@reductions[["umap"]]@cell.embeddings))
  }

  # Convert spliced and unspliced layers if available to Seurat objects
  if("spliced" %in% names(ad@assays)){
    s <- as.Seurat(ad, counts = "spliced", data = NULL)
    X@assays[["spliced"]]<-s@assays$originalexp
    file_sub<-paste0(file_sub,"_spliced")
  }
  if("unspliced" %in% names(ad@assays)){
    un <- as.Seurat(ad, counts = "unspliced", data = NULL)
    X@assays[["unspliced"]]<-un@assays$originalexp
    file_sub<-paste0(file_sub,"_unspliced")
  }

  # use cell type annotation column as identity
  if (!is.null(annotation_column) && (annotation_column %in% colnames(X@meta.data))) {
    Idents(X)<-X[[annotation_column]]
  } 

  # save converted Seurat object
  DefaultAssay(X) <- "RNA"
  saveRDS(X, file=paste0(file_sub,".rds"))

  return(X)
}


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


#' Check the required input objects for all the functions.
#'
#' @description This function checks the required input objects and variables for all the functions
#' in the scDown package using checkmate package
#'
#'
#' @param seurat_obj character string of full path to the h5ad file.
#' @param species species
#' @param output_dir output_dir
#' @param annotation_column annotation_column
#' @param group_column group_column
#'
#' @noRd

check_required_variables<-function(seurat_obj,species=NULL,output_dir,annotation_column,group_column)
{
  checkmate::expect_class(seurat_obj,"Seurat",label="seurat_obj")
  checkmate::expect_choice(species,c("human","mouse"),label = "species",null.ok = TRUE)
  checkmate::expect_choice(group_column, colnames(seurat_obj@meta.data),label="group_column",null.ok = TRUE)
  ###Be default we use seurat Idents
  checkmate::expect_choice(annotation_column, colnames(seurat_obj@meta.data),label="annotation_column",null.ok = TRUE)
  checkmate::expect_directory(output_dir,access="rw",label = "output_dir")
}

