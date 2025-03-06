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

  # Convert .h5ad to Seurat object
  library(zellkonverter)
  # Read the .h5ad file
  ad <- readH5AD(h5ad_file)
  ad
  file_sub <- sub("\\.h5ad$", "", h5ad_file)

  # Convert the main count matrix of h5ad to RNA assay of a Seurat object
  assays <- names(ad@assays)
  default_assay_name <- ifelse("X" %in% assays, "X", assays[1])
  X <- as.Seurat(ad, counts = default_assay_name, data = NULL)
  X@assays[["RNA"]]<-CreateAssayObject(counts = X@assays[["originalexp"]]@counts)
  X@assays[["RNA"]]@key<-"rna_"

  # Convert spliced and unspliced layers if available in original h5ad
  if("spliced" %in% names(ad@assays)){
    spliced_mat <- as.matrix(ad@assays@data$spliced)
    rownames(spliced_mat) <- rownames(ad)
    colnames(spliced_mat) <- colnames(ad)  
    X@assays[["spliced"]]<-CreateAssayObject(counts = spliced_mat)
    X@assays[["spliced"]]@key<-"spliced_"
    file_sub<-paste0(file_sub,"_spliced")
  }
  if("unspliced" %in% names(ad@assays)){
    unspliced_mat <- as.matrix(ad@assays@data$unspliced)
    rownames(unspliced_mat) <- rownames(ad)
    colnames(unspliced_mat) <- colnames(ad)  
    X@assays[["unspliced"]]<-CreateAssayObject(counts = unspliced_mat)
    X@assays[["unspliced"]]@key<-"unspliced_"
    file_sub<-paste0(file_sub,"_unspliced")
  }

  # Set Default Assay to RNA
  DefaultAssay(X) <- "RNA"

  # Convert lower dimensional embeddings (pca and umap) if available
  if ('X_pca' %in% names(X@reductions)){
    pca_coords<-X@reductions[["X_pca"]]@cell.embeddings
    X@reductions[["pca"]] <- CreateDimReducObject(embeddings = pca_coords, key = "PCA_", assay = "RNA")
  }
  if ('X_umap' %in% names(X@reductions)){
    umap_coords<-X@reductions[["X_umap"]]@cell.embeddings
    X@reductions[["umap"]] <- CreateDimReducObject(embeddings = umap_coords, key = "UMAP_", assay = "RNA")
  }

  # use cell type annotation column as identity if provided
  if (!is.null(annotation_column) && (annotation_column %in% colnames(X@meta.data))) {
    Idents(X)<-X[[annotation_column]]
  } 

  # save converted Seurat object
  saveRDS(X, file=paste0(file_sub,".rds"))

  return(X)
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

