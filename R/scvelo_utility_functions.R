library(Seurat)
library(dplyr)
library(velocyto.R)


#' Add spliced and unspliced matrices from a Loom file to a Seurat object
#'
#' This function takes a Seurat object and a Loom file as inputs. It reads spliced and 
#' unspliced matrices from the Loom file and adds them as new assays to the Seurat object. 
#' The function first gets standardized barcodes for the Loom file and renames the cells in 
#' the Seurat object to match the barcodes. It then identifies the shared barcodes between 
#' the two datasets, filters the Seurat object to match the barcodes in the Loom file, and adds 
#' the spliced and unspliced matrices as new assays.
#'
#' @param X Seurat object to which spliced and unspliced matrices will be added
#' @param loomFile path to the Loom file containing spliced and unspliced matrices
#' @return Seurat object with spliced and unspliced matrices added as new assays
#'
#' @noRd

addSUmatrices <- function(X, loomFile){
    library(Matrix)

    # Reading spliced/unspliced matrices
    SU <- velocyto.R::read.loom.matrices(loomFile)

    # Getting standardized cell barcodes by removing prefix and suffix
    colnames(SU[[1]]) <- colnames(SU[[2]]) <- gsub('^[[:print:]]+\\:|x$', '', colnames(SU[[1]]))
    new.names <- gsub('^(?:.*?_)?([A-Z0-9]+)[_-].*$', '\\1', colnames(X))
    X <- Seurat::RenameCells(X, new.names= new.names)

    # Identifying shared barcodes
    sharedBarcodes <- intersect(colnames(SU[[1]]), colnames(X))
    sharedGenes <- intersect(rownames(SU[[1]]), rownames(X))

    # Filtering object to match with the barcodes available in the other matrices
    X <- X[sharedGenes,sharedBarcodes]

    # Adding spliced/unspliced matrices as assays
    X[['spliced']] <- SeuratObject::CreateAssay5Object(counts = as(SU[[1]][sharedGenes, sharedBarcodes], "CsparseMatrix"))
    X[['unspliced']] <- SeuratObject::CreateAssay5Object(counts = as(SU[[2]][sharedGenes, sharedBarcodes], "CsparseMatrix"))

    # Reverting to the original cell barcodes 
    X <- Seurat::RenameCells(X, new.names = X$orig.bc)

    # Returning object with new assays
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
#' @return A named list of spliced and unspliced matrices extracted from loom.
#'
#' @noRd

extractSUmatrices <- function(loomFile, seurat_cells, seurat_genes) {
    library(Matrix)

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