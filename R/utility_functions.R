#' Read in a h5ad file and convert to a Seurat object.
#'
#' @description This function uses the library(anndataR) to convert the raw counts, metadata,
#' lower dimensional embeddings (pca and umap), and spliced and unspliced counts if available from a h5ad 
#' file to a Seurat object.
#'
#' @param h5ad_file character string of full path to the h5ad file.
#' @param annotation_column character variable specifying the metdata column name of cell type annotations. 
#' @param do_normalization boolean variable specifying whether to perform Seurat normalization after conversion
#' @return a Seurat object.
#'
#' @export
#'
h5adToSeurat <- function(h5ad_file, annotation_column=NULL, do_normalization=TRUE){

  # Convert .h5ad to Seurat object
  # library(zellkonverter)
  library(anndataR)
  library(Seurat)
  # Read the .h5ad file
  # ad <- readH5AD(h5ad_file)
  ad <- read_h5ad(h5ad_file)
  file_sub <- sub("\\.h5ad$", "", h5ad_file)

  # Convert the main count matrix of h5ad to RNA assay of a Seurat object
  # assays <- names(ad@assays)
  # default_assay_name <- ifelse("X" %in% assays, "X", assays[1])
  # X <- as.Seurat(ad, counts = default_assay_name, data = NULL)
  # X@assays[["RNA"]]<-CreateAssayObject(counts = LayerData(X, layer="counts", assay="originalexp"))
  # X@assays[["RNA"]]@key<-"rna_"
  X <- ad$as_Seurat(x_mapping="counts")

  # Save the raw counts as another assay for later use
  X[["originalexp"]] <- X[["RNA"]]

  # Convert spliced and unspliced layers if available in original h5ad
  # if("spliced" %in% names(ad@assays)){
  if("spliced" %in% Layers(X[["RNA"]])){
    # spliced_mat <- as.matrix(ad@assays@data$spliced)
    # rownames(spliced_mat) <- rownames(ad)
    # colnames(spliced_mat) <- colnames(ad)  
    # X@assays[["spliced"]]<-CreateAssayObject(counts = spliced_mat)
    # X@assays[["spliced"]]@key<-"spliced_"
    X[["spliced"]] <- CreateAssay5Object(counts = LayerData(X, layer="spliced", assay="RNA"))
    X[["RNA"]]$spliced <- NULL
    file_sub<-paste0(file_sub,"_spliced")
  }
  # if("unspliced" %in% names(ad@assays)){
  if("unspliced" %in% Layers(X[["RNA"]])){
    # unspliced_mat <- as.matrix(ad@assays@data$unspliced)
    # rownames(unspliced_mat) <- rownames(ad)
    # colnames(unspliced_mat) <- colnames(ad)  
    # X@assays[["unspliced"]]<-CreateAssayObject(counts = unspliced_mat)
    # X@assays[["unspliced"]]@key<-"unspliced_"
    X[["unspliced"]] <- CreateAssay5Object(counts = LayerData(X, layer="unspliced", assay="RNA"))
    X[["RNA"]]$unspliced <- NULL
    file_sub<-paste0(file_sub,"_unspliced")
  }
  if("ambiguous" %in% Layers(X[["RNA"]])){
    X[["RNA"]]$ambiguous <- NULL
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
    Idents(X)<-X@meta.data[[annotation_column]]
  }
  cat("Object converted.\n")

  if (do_normalization) {
    X = NormalizeData(X)
  }

  # Clean up the obs names in R
  colnames(X@meta.data) = make.names(colnames(X@meta.data), unique=TRUE)
  cat("Added meta data fields:\n")
  print(colnames(X@meta.data))

  # save converted Seurat object
  saveRDS(X, file=paste0(file_sub,".rds"))
  cat("Seurat object saved to ",file_sub,".rds\n", sep="")

  return(X)
}

#' Cleanup unsafe characters in the Seurat meta.data and Idents.
#'
#' @description This function converts unsafe characters in a file name to a underscore.
#'
#' @param seurat_obj Seurat object
#' @return a Seurat object.
#'
#' @export
#'
cleanSeuratMeta <- function(seurat_obj){
  library(Seurat)
  md <- seurat_obj@meta.data
  # Clean up the meta column names
  colnames(md) = make.names(colnames(md), unique=TRUE)

  # Loop through every column
  for (col in colnames(md)) {
    # If character, clean actual data
    if (is.character(md[[col]])) {
      md[[col]] <- gsub("[^[:alnum:]_()+-]", "_", md[[col]])
      # Clean multiple underscores
      md[[col]] <- gsub("_+", "_", md[[col]])

    } else if (is.factor(md[[col]])) {
      # If a factor, clean the levels
      clean_levels <- gsub("[^[:alnum:]_()+-]", "_", levels(md[[col]]))
      clean_levels <- gsub("_+", "_", clean_levels)
      levels(md[[col]]) <- clean_levels
    }
  }
  # Add the cleaned metadata back to object
  seurat_obj@meta.data <- md

  # Clean the Idents factor too
  active_idents <- Idents(seurat_obj)
  clean_levels <- gsub("[^[:alnum:]_()+-]", "_", levels(active_idents))
  clean_levels <- gsub("_+", "_", clean_levels)
  levels(active_idents) <- clean_levels
  Idents(seurat_obj) <- active_idents

  return(seurat_obj)
}

#' Check the required input objects for all the functions.
#'
#' @description This function checks the required input objects and variables for all the functions
#' in the scDown package using checkmate package
#'
#'
#' @param seurat_obj Seurat object
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

