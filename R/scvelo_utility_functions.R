library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(scales)
library(velocyto.R)
library(ggrepel)


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
    X[['spliced']] <- Seurat::CreateAssayObject(SU[[1]][sharedGenes,sharedBarcodes])
    X[['unspliced']] <- Seurat::CreateAssayObject(SU[[2]][sharedGenes,sharedBarcodes])

    # Reverting to the original cell barcodes 
    X <- Seurat::RenameCells(X, new.names = X$orig.bc)
    X@reductions$umap <- as(X@reductions$umap, "DimReduc")

    # Returning object with new assays
    return(X)
}


