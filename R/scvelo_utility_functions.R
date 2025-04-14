library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(scales)
library(velocyto.R)
library(velociraptor)
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

    # Returning object with new assays
    return(X)
}


#' Compute single-cell velocity using scvelo package
#'
#' The function computes the single-cell velocity using the scvelo package.
#' The function takes in total, spliced and unspliced counts from a Seurat 
#' object and uses them to calculate single-cell velocity. The function 
#' returns a SingleCellExperiment object with the information required to 
#' generate the plots.
#'
#' @param X Seurat object containing the counts matrices.
#' @param mode Character mode specifying the type of velocity computation to use. Available modes are "steady_state" 
#'(original), "deterministic", "stochastic" (fastest:recommended), "dynamical".
#' @return A SingleCellExperiment object with the information required to generate the plots.
#'
#' @noRd

doVelocity <- function(X, mode = 'stochastic'){
    # Available Modes: "steady_state" (original), "deterministic", "stochastic" (fastest:recommended), "dynamical"
    # More info: https://scvelo.readthedocs.io/en/stable/scvelo.tl.velocity/
    # Takes the total/spliced/unspliced counts from the provided Seurat object
    countMatrices <- list(
        X = X@assays$RNA@counts, 
        spliced = X@assays$spliced@counts, 
        unspliced = X@assays$unspliced@counts
    )
    # Compute the single-cell velocity
    O <- velociraptor::scvelo(countMatrices, mode = mode)

    # Returns a SingleCellExperiment object with the information required to generate the plots
    return(O)
}


#' Compute vector field for plotting velocity arrows using scvelo package
#'
#' The function computes the vector field for plotting velocity arrows 
#' using the scvelo package. The function takes in a Seurat object 
#' containing the counts matrices, the scVeloOutput object returned by 
#' doVelocity function, the reduction method used for plotting (default: umap), 
#' dimensions to be plotted (default: 1:2) and the resolution of the grid (default: 50). 
#' The function returns a 4-column data.frame that can be used with 
#' geom_segment to plot velocity arrows.
#'
#' @param X Seurat object containing the counts matrices.
#' @param scVeloOutput scVeloOutput object returned by doVelocity function.
#' @param reduction Character specifying the dimensionality reduction method to use for plotting. Default is "umap".
#' @param dims Numeric vector specifying the dimensions to be plotted. Default is 1:2.
#' @param resolution Numeric value specifying the resolution of the grid to use. Default is 50.
#' @return A 4-column data.frame containing the x and y coordinates of the start and end points of the velocity arrows.
#' 
#'
#' @noRd

getVectorField <- function(X, scVeloOutput, reduction = 'umap', dims = 1:2, resolution = 50){
    # Getting the requested reduction
    E <- Seurat::Embeddings(X, reduction = reduction)

    # Making a compatible object
    E <- as.matrix(E[,dims])
    dimnames(E) <- NULL

    # Computing vectors
    O <- velociraptor::embedVelocity(E, scVeloOutput)
    O <- velociraptor::gridVectors(E,O,resolution = resolution)
    rownames(O) <- NULL

    # Returning a 4 columns data.frame to be used with geom_segment
    return(O)
}




#' Plot the vector field on lower-dimensional embeddings.
#'
#' This function takes in a Seurat object with dimensionality reduction and the four-column
#' data frame outputted by getVectorField() function, which contains velocity vector start
#' and end coordinates. Then the function uses geom_point to plot the lower-dimensional space
#' and geom_segment to overlay the velocity vectors.
#'
#' @param X a Seurat object with dimension reduction coordinates.
#' @param tpVF the 4-column dataframe outputted by getVectorField() with start and end point of each velocity vector.
#' @param group a character vector of conditions or timepoints in the datasets specifying cells from which time should 
#'be plotted with velocity arrows.
#' @param group_column a character string specifying name of the metadata that has timepoint information.
#' @param color_scale a character vector of colors to be used in plotting.
#' @param name_by a character string specifying which metadata in X is the color scale named by.
#' @param grid_res numeric value specifying the resolution of the grid, purely for figure naming purpose.
#' @param arrow_size a float or integer specifying the velocity vector arrow head size.
#' @param vector_width a float or integer specifying the velocity vector line width.
#'
#' @noRd

plotVectorField <- function(X, tpVF, group=NULL, group_column=NULL, color_scale=NULL, name_by=NULL, grid_res=NULL, arrow_size=0.5, vector_width=0.5){
    
    # get lower-dimensional coordinates
    D <- X@reductions$umap@cell.embeddings
    D <- as.data.frame.array(D)

    # add current cell labels as a metadata
    X[["cell.type"]] <- Idents(X)

    # add metadata from the Seurat object to the dataframe
    for (columns in colnames(X@meta.data)){
        D[[columns]] <- as.character(X@meta.data[[columns]])
    }

    # add color scale
    if (!is.null(color_scale) & !is.null(name_by)){
        labels <- unique(X[[name_by]][[name_by]])
        names(color_scale) <- labels
        D$Color <- color_scale[D[[name_by]]]
    } else {
        default_colors <- scales::hue_pal()(length(unique(Idents(X))))
        labels <- unique(Idents(X))
        names(default_colors) <- labels
        D$Color <- default_colors[D$cell.type]
    }
    
    # set uninterested cells to gray color
    if (!is.null(group)){
        D$Color[!X[[group_column]][[group_column]] %in% group] <- 'gray95'
    }

    png(file=paste0("scvelo/images/velocityField_",paste(group, collapse="_"),"_gridRes",grid_res,"_arrowSize",arrow_size,"_width",vector_width,".png",sep=""), width = 1800, height = 1800, res = 300)
    P <- ggplot2::ggplot(D, aes(UMAP_1, UMAP_2)) +
            geom_point(color = D$Color, size = 0.01) +
            theme_void() +
            theme(legend.position = 'None', 
                  panel.grid = element_blank(),
                  panel.border = element_blank()) +
            xlab('UMAP 1') +
            ylab('UMAP 2') +
            geom_segment(data=tpVF, 
            mapping=aes(
                x=start.1, 
                y=start.2, 
                xend=end.1, 
                yend=end.2), 
            linewidth = vector_width,
            lineend = 'round', 
            linejoin = 'round',
            arrow=arrow(length=unit(arrow_size, "mm")))
    print(P)
    dev.off()
}


#' Transfer computed pseudotime values using scvelo package to a Seurat object
#'
#' The function transfers the computed pseudotime values from the scvelo package
#' to the Seurat object. The function takes in a Seurat object containing the 
#' counts matrices and a scVeloOutput object returned by doVelocity function. 
#' The function returns the Seurat object with two additional columns in the metadata
#' containing the computed pseudotime values and the velocity confidence values.
#'
#' @param X Seurat object containing the counts matrices.
#' @param scVeloOutput scVeloOutput object returned by doVelocity function.
#' @return A Seurat object with two additional columns containing the computed pseudotime values and the velocity confidence values.
#' 
#'
#' @noRd

getVelocityPseudotime <- function(X, scVeloOutput){
    # Transfering computed pseudotime values to the Seurat object
    X$velocity_pseudotime <- scVeloOutput$velocity_pseudotime
    X$velocity_velocity_confidence <- scVeloOutput$velocity_confidence
    return(X)
}


# =============================== NOT USED FOR NOW ===============================
#'
#' Plot the pseudotime on lower-dimensional embeddings.
#'
#' This function takes in a Seurat object outputted by getVelocityPseudotime() function, which
#' contains pseudotime values in its metadata. Lower-dimensional space is plotted using geom_point
#' and pseudotime values are used as scales for coloring.
#'
#' @param X a Seurat object outputted by getVelocityPseudotime() function.
#' @param color_scale two colors to use to define color gradient when plotting pseudotime.
#'
#' @noRd

plotPseudotime <- function(X, color_scale=c("darkblue", "yellow")){
    
    # get lower-dimensional coordinates and pseudotime values
    D <- X@reductions$umap@cell.embeddings
    D <- as.data.frame.array(D)
    pseudo_df <- as.data.frame(X$velocity_pseudotime)
    D <- merge(D, pseudo_df, by='row.names', all=TRUE)
    colnames(D)[4] <- "Pseudotime"
    D <- D[order(D$Pseudotime, decreasing = TRUE), ]

    png('V-D0D3-pseudotime.png', width = 1800, height = 1800, res = 300)
    P <- ggplot2::ggplot(D, aes(UMAP_1, UMAP_2)) +
        geom_point(aes(color = Pseudotime), size = 0.01) +
        scale_color_gradient(low = color_scale[1], high = color_scale[2]) +
        scale_color_discrete(na.value="gray95") +
        theme_void() +
        theme(panel.grid = element_blank(),
              panel.border = element_blank()) +
        xlab('UMAP 1') +
        ylab('UMAP 2')
    print(P)
    dev.off()
}
