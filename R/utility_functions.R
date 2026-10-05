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
  library(anndataR)
  library(Seurat)
  # Read the .h5ad file
  ad <- read_h5ad(h5ad_file)
  file_sub <- sub("\\.h5ad$", "", h5ad_file)

  # Convert the main count matrix of h5ad to RNA assay of a Seurat object
  X <- ad$as_Seurat(x_mapping="counts")

  # Save the raw counts as another assay for later use
  X[["originalexp"]] <- X[["RNA"]]

  # Convert spliced and unspliced layers if available in original h5ad
  if("spliced" %in% Layers(X[["RNA"]])){
    X[["spliced"]] <- CreateAssay5Object(counts = LayerData(X, layer="spliced", assay="RNA"))
    X[["RNA"]]$spliced <- NULL
    file_sub<-paste0(file_sub,"_spliced")
  }
  if("unspliced" %in% Layers(X[["RNA"]])){
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

#' Aligns genes (rows) of a matrix to a target gene list
#'
#' This function takes a matrix and a target gene list and consolidate the genes
#' at the same time padding missing genes with 0s.
#'
#' @param mat A input matrix.
#' @param target_genes A target gene list to conform to.
#' @return A sorted and padded matrix.
#'
#' @noRd
#'
alignGenes <- function(mat, target_genes) {
    # Subset to shared genes
    shared_genes <- intersect(rownames(mat), target_genes)
    mat <- mat[shared_genes, , drop = FALSE]

    # Pad missing genes with 0
    missing_genes <- setdiff(target_genes, shared_genes)
    if (length(missing_genes) > 0) {
        zero_mat <- Matrix::Matrix(0, nrow = length(missing_genes), ncol = ncol(mat), sparse = TRUE)
        rownames(zero_mat) <- missing_genes
        mat <- rbind(mat, zero_mat)
    }
    # Use all Seurat genes
    return(mat[target_genes, , drop = FALSE])
}

#' Extract a loom file, clean barcodes, and align it to the requested cells and genes
#'
#' This function takes a loom file and lists of cells and genes from Seurat object and
#' output spliced and unspliced matrices with shared barcodes and all target genes.
#'
#' @param loomFile A loom file path.
#' @param seurat_cells Seurat cell barcode list.
#' @param seurat_genes Seurat gene list.
#'
#' @return A named list of spliced and unspliced matrices extracted from loom.
#'
#' @noRd
#'
extractSUmatrices <- function(loomFile, seurat_cells, seurat_genes) {
    library(Matrix)
    library(velocyto.R)

    # Create a mapping between standardized cell barcodes and original Seurat barcodes
    clean_bcs <- gsub('^(?:.*?_)?([A-Z0-9]+)[_-].*$', '\\1', seurat_cells)
    names(clean_bcs) <- seurat_cells

    # Read spliced/unspliced matrices
    SU <- velocyto.R::read.loom.matrices(loomFile)
    spliced <- as(SU[[1]], "CsparseMatrix")
    unspliced <- as(SU[[2]], "CsparseMatrix")

    # Get standardized cell barcodes by removing prefix and suffix
    colnames(spliced) <- colnames(unspliced) <- gsub('^[[:print:]]+\\:|x$', '', colnames(spliced))

    # Identify shared barcodes and filter SU matrices
    shared_bcs <- intersect(colnames(spliced), clean_bcs)
    spliced <- spliced[, shared_bcs, drop = FALSE]
    unspliced <- unspliced[, shared_bcs, drop = FALSE]

    # Rename loom matrices to the original Seurat cell barcodes
    orig_bcs <- names(clean_bcs)[match(shared_bcs, clean_bcs)]
    colnames(spliced) <- colnames(unspliced) <- orig_bcs

    # Pad missing genes and align genes to Seurat RNA assay
    spliced <- alignGenes(spliced, seurat_genes)
    unspliced <- alignGenes(unspliced, seurat_genes)

    return(list(spliced = spliced, unspliced = unspliced))
}

#' Save a Seurat object to H5AD file.
#'
#' This function uses the library(anndataR) and saves a Seurat object to H5AD file. It extract and add
#' unspliced and spliced matrices from a list of loom files if provided.
#'
#' @param seurat_obj Seurat object
#' @param loom_files path and file names of Spliced and unspliced counts of the scRNA-seq data (Required if
#' not provided in `seurat_obj`). For multiple loom files, the file names in `loom_files` must correspond to
#' the values in the metadata column in the Seurat object specified by `loom_file_subset_column`
#' @param output_dir A character vector specifying the output directory
#' @param loom_file_subset_by A character variable specifying how the Seurat object should be subsetted in order
#' to match each indivdidual loom file name. This vector must correspond to a metadata column specified by
#' `loom_file_subset_column` in the Seurat object, and the order of `loom_file_subset_by` must match the order of
#' loom_files. If `loom_file_subset_by` is not provided (default `loom_file_subset_by <- NULL`), it will be
#' automatically extracted from file names in `loom_files` to ensure the correct order.
#' If there is only one loom file provided, this variable should be left blank: `loom_file_subset_by <- NULL`
#' @param loom_file_subset_column A string specifying the name of the metadata column in the Seurat object that
#' should be used for subsetting to match each of the loom files. This column must exist in the Seurat object
#' metadata. If there is only one loom file provided, this variable should be left blank:
#' `loom_file_subset_column <- NULL`; if there are multiple loom files provided, this variable needs to be provided.
#' @param annotation_column A character variable specifying which metadata column of the h5ad object contains
#' cell type annotations.
#' @param groups A list of character vectors representing groups of conditions or time points used to calculate
#' RNA velocity separately, default: NULL
#' @param group_column A string specifying the name of the metadata column in the h5ad object that contains groups.
#'
#' @export
#'
seuratToH5ad <- function(seurat_obj,
                         loom_files=NULL,
                         output_dir=".",
                         loom_file_subset_by=NULL,
                         loom_file_subset_column="orig.ident",
                         annotation_column='ID',
                         groups=NULL,
                         group_column=NULL){
  library(Matrix)
  library(Seurat)
  library(anndataR)

  # create subdirectories in the output directory
  subdirectories <- c("data")

  for(i in subdirectories){
    dir.create(file.path(output_dir,i), showWarnings = FALSE, recursive = TRUE)
  }
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)

  ### Input
  library(Seurat)
  library(SeuratObject)
  checkmate::test_class(seurat_obj, "Seurat")
  object_annotated <- seurat_obj

  # use cell type annotation column as identity
  if(checkmate::test_string(annotation_column, null.ok=FALSE)){
    checkmate::expect_choice(annotation_column, colnames(seurat_obj@meta.data), label = "annotation_column")
    Seurat::Idents(object_annotated) <- object_annotated[[annotation_column]][,1]
  }

  checkmate::expect_character(groups, min.len = 1, any.missing = FALSE, label="groups", null.ok = TRUE)
  if(!is.null(groups)){
    checkmate::assert_string(group_column, null.ok = FALSE)
    checkmate::expect_choice(group_column, colnames(seurat_obj@meta.data), label = "group_column")
  }

  # convert seurat v3/v4 to v5
  object_annotated <- UpdateSeuratObject(object_annotated)
  if (inherits(object_annotated[["RNA"]], "Assay")) {
    object_annotated[["RNA"]] <- as(object_annotated[["RNA"]], Class = "Assay5")
  }

  # check if spliced and unspliced data is already in seurat_obj
  if(!(("spliced" %in% names(object_annotated@assays) & ("unspliced" %in% names(object_annotated@assays))))){
    checkmate::assert_character(loom_files, min.len = 1, null.ok = FALSE, any.missing = FALSE)
    checkmate::assert_character(loom_file_subset_by, null.ok = TRUE, any.missing = FALSE)
    checkmate::assert_string(loom_file_subset_column, null.ok = FALSE)
    checkmate::expect_choice(loom_file_subset_column, colnames(seurat_obj@meta.data), label = "loom_file_subset_column")

    # add cell barcode as metadata
    object_annotated$orig.bc <- colnames(object_annotated)

    # add spliced and unspliced matrices as new assays
    if (length(loom_files) > 1){
      if(is.null(loom_file_subset_by)){
        loom_file_subset_by=gsub(".loom","",gsub(".*/","",loom_files))
      }

      expr <- FetchData(object = object_annotated, vars = loom_file_subset_column)

      matrices_list <- lapply(1:length(loom_files), function(i) {
        # Find the cells in this sample
        sample_cells <- colnames(object_annotated)[which(expr == loom_file_subset_by[i])]
        return(extractSUmatrices(loom_files[i], sample_cells, rownames(object_annotated)))
      })

      # Concatenate all split matrices
      cat("Merging split matrices into single matrices...\n")
      spliced_all <- do.call(cbind, lapply(matrices_list, `[[`, "spliced"))
      unspliced_all <- do.call(cbind, lapply(matrices_list, `[[`, "unspliced"))

    } else {
      # object_annotated <- addSUmatrices(object_annotated, loom_files[1])
      matrices_list <- extractSUmatrices(loom_files[1], colnames(object_annotated), rownames(object_annotated))
      spliced_all <- matrices_list$spliced
      unspliced_all <- matrices_list$unspliced
    }

    # Sync the barcodes in the Seurat object with looms
    shared_bcs <- intersect(colnames(object_annotated), colnames(spliced_all))
    if (length(shared_bcs) < ncol(object_annotated)) {
        cat("Filtering Seurat object to cells in loom files...\n")
        object_annotated <- object_annotated[, shared_bcs]
    }
    spliced_all <- spliced_all[, shared_bcs]
    unspliced_all <- unspliced_all[, shared_bcs]

    cat("Adding spliced and unspliced matrices to Seurat...\n")
    object_annotated[["spliced"]] <- CreateAssay5Object(counts = spliced_all)
    object_annotated[["unspliced"]] <- CreateAssay5Object(counts = unspliced_all)

    saveRDS(object_annotated, "data/obj_spliced_unspliced.rds")
  }

  # Add the originalexp assay if not there
  if (!"originalexp" %in% Assays(object_annotated)) {
    object_annotated[["originalexp"]] <- object_annotated[["RNA"]]
  }

  # Save to h5ad
  seurat_obj <- object_annotated
  print(seurat_obj)
  seurat_obj[[annotation_column]] <- as.character(seurat_obj[[annotation_column]][,1])
  Idents(seurat_obj)<-seurat_obj[[annotation_column]][,1]

  cat("Converting Seurat object to anndata...\n")
  adata <- as_AnnData(seurat_obj,
    assay_name = "originalexp",
    x_mapping = "counts",
    obs_mapping = TRUE,
    var_mapping = TRUE,
    obsm_mapping = list(X_umap = "umap"),
    output_class = "InMemory"
  )
  adata_unspliced <- as_AnnData(seurat_obj,
    assay_name = "unspliced",
    x_mapping = "counts",
    output_class = "InMemory"
  )
  adata_spliced <- as_AnnData(seurat_obj,
    assay_name = "spliced",
    x_mapping = "counts",
    output_class = "InMemory"
  )
  adata$layers['unspliced'] <- adata_unspliced$layers['counts']
  adata$layers['spliced'] <- adata_spliced$layers['counts']
  # Save the AnnData object
  h5ad_file <- file.path("data/obj_spliced_unspliced.h5ad")
  anndataR::write_h5ad(adata, path=h5ad_file, mode="w")

  # Split dataset by group and save to individual h5ad
  if (length(groups) != 0){
    file_base <- gsub(".h5ad", "", basename(h5ad_file))
    input_dir <- dirname(h5ad_file)
    for (group in groups){
      group_label=paste(group, collapse="_")
      file_name=paste0(file_base,"_",group_label,".h5ad")
      h5ad_group_file=file.path("data", file_name)
      # Create a subsetted h5ad file and save it to the output directory
      subset_adata <- adata[adata$obs[[group_column]] %in% group,]
      anndataR::write_h5ad(subset_adata, path=h5ad_group_file, mode="w")
    }
  }
}

#' Cleanup unsafe characters in the Seurat meta.data and Idents.
#'
#' @description This function converts unsafe characters in a file name to a underscore.
#'
#' @param seurat_obj Seurat object
#'
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
#'
check_required_variables<-function(seurat_obj,species=NULL,output_dir,annotation_column,group_column)
{
  checkmate::expect_class(seurat_obj,"Seurat",label="seurat_obj")
  checkmate::expect_choice(species,c("human","mouse"),label = "species",null.ok = TRUE)
  checkmate::expect_choice(group_column, colnames(seurat_obj@meta.data),label="group_column",null.ok = TRUE)
  ###Be default we use seurat Idents
  checkmate::expect_choice(annotation_column, colnames(seurat_obj@meta.data),label="annotation_column",null.ok = TRUE)
  checkmate::expect_directory(output_dir,access="rw",label = "output_dir")
}

