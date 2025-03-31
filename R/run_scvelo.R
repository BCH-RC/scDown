#' Function to run the scVelo pipeline using velociraptor 
#'
#' This function performs RNA velocity calculations from .loom files with the scVelo package.
#' In this function, users can calculate RNA velocity of the whole data as well as a subset of groups.
#'
#' @param seurat_obj Seurat object containing the scRNA-seq data (Required)
#' @param annotation_column A string specifying the metdata column name of cell type annotations. Default: 
#' `annotation_column <- NULL`, and `Idents(seurat_obj)` is used as the cell type annotation
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
#' @param mode Mode for scVelo velocity calculation, default: stochastic; can be one of four options: "steady_state" 
#' (original), "deterministic", "stochastic" (fastest: recommended, default), "dynamical"
#' @param grid_resolutions A vector of integers specifying the number of grids along each axis, essentially 
#' controlling the number of vectors on umap, default 50.
#' @param arrow_sizes A vector of integers or floats controlling velocity vector size (arrow head size), default 0.5
#' @param vector_widths A vector of integers or floats controlling velocity vector size (vector width), default 0.5
#' @param groups A list of character vectors representing groups of conditions or time points used to calculate RNA 
#' velocity separately, default: NULL
#' @param group_column A string specifying the name of the metadata column in the Seurat object that should be 
#' used for subsetting each group in `groups`
#' @param color_scale A character vector of colors to be used in plotting, must match number of unique values in 
#' the metadata column marked by @name_by
#' @param name_by A string specify the name of the metadata column in the Seurat object that should be used for colors
#' @param cores A numeric variable to set the number of cores to be used
#'
#' @return A list of scVelo data objects
#'
#' @export
#'
#' Estimate RNA velocity for spliced and unspliced counts of scRNA-seq data


run_scvelo <- function(seurat_obj,loom_files=NULL,output_dir=".",loom_file_subset_by=NULL,loom_file_subset_column="orig.ident",
                    annotation_column=NULL,mode='stochastic',grid_resolutions=c(50),arrow_sizes=c(0.5,1),vector_widths=c(0.25,0.5),
                    groups=NULL,group_column=NULL,color_scale=NULL,name_by=NULL,cores=8){

# create subdirectories in the output directory
setwd(output_dir)
subdirectories <- c("scvelo",
                    "scvelo/csv",
                    "scvelo/rds",
                    "scvelo/images")

for(i in subdirectories){
    dir.create(file.path(output_dir,i), showWarnings = F, recursive = T)
}

### Input
library(Seurat)
checkmate::test_class(seurat_obj, "Seurat")
object_annotated <- seurat_obj

# use cell type annotation column as identity
if(checkmate::test_string(annotation_column, null.ok=FALSE)){
  checkmate::expect_choice(annotation_column, colnames(seurat_obj@meta.data), label = "annotation_column")
  Seurat::Idents(object_annotated) <- object_annotated[[annotation_column]][,1]
}

checkmate::assert_list(groups, types = c("integer","numeric", "character"), null.ok = TRUE)
checkmate::expect_numeric(cores, min.len = 1, max.len = 1, any.missing = FALSE,label="cores")
if(!is.null(groups)){
  checkmate::assert_string(group_column, null.ok = FALSE)
  checkmate::expect_choice(group_column, colnames(seurat_obj@meta.data), label = "group_column")
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

    # empty list to store objects with spliced/unspliced matrices
    object_SU_list <- list()

    # subset to corresponding cells in loom file
    for (i in 1:length(loom_files)){
        expr <- FetchData(object = object_annotated, vars = loom_file_subset_column)
        object_subset <- object_annotated[, which(x = (expr == loom_file_subset_by[i]))]
    
        object_subset_SU <- addSUmatrices(object_subset, loom_files[i])
        object_SU_list <- append(object_SU_list, object_subset_SU)
    }

    # merge subsetted seurat objects
    object_annotated <- merge(object_SU_list[[1]], object_SU_list[-1], merge.dr = "umap")
    object_annotated <- RenameCells(object_annotated, new.names = object_annotated$orig.bc)

  } else {
    object_annotated <- addSUmatrices(object_annotated, loom_files[1])
  }

}

library(SeuratDisk)
# save to h5ad so if needed, can be used to conduct scvelo downstream analysis
#To prevent overwriting error, only saves if the file does not exist
seurat_obj<-object_annotated
seurat_obj[[annotation_column]] <- as.character(seurat_obj[[annotation_column]][,1])
Idents(seurat_obj)<-seurat_obj[[annotation_column]][,1]
if (!file.exists("scvelo/rds/obj_spliced_unspliced.h5Seurat")) {
    SeuratDisk::SaveH5Seurat(seurat_obj, filename = "scvelo/rds/obj_spliced_unspliced.h5Seurat")
}
if (!file.exists("scvelo/rds/obj_spliced_unspliced.h5ad")) {
    SeuratDisk::Convert("scvelo/rds/obj_spliced_unspliced.h5Seurat", dest = "h5ad")
}


### RNA velocity analysis
library(ggplot2)
# RNA velocity for the entire Seurat object
if(Sys.info()["machine"] != "aarch64") {

tpV <- doVelocity(object_annotated, mode=mode)
for (grid_resolution in grid_resolutions){
    tpVF <- getVectorField(object_annotated, tpV, reduction = 'umap', resolution = grid_resolution)
    save(tpVF, file = paste0('scvelo/rds/ALL_gridRes', grid_resolution,'.RData',sep=""))
    #load(paste0('scvelo/rds/ALL_gridRes', grid_resolution,'.RData',sep=""))
    for (arrow_size in arrow_sizes){
        for (vector_width in vector_widths){
            plotVectorField(object_annotated, tpVF, color_scale=color_scale, name_by=name_by, grid_res=grid_resolution, arrow_size=arrow_size, vector_width=vector_width)
        }
    }
}

print("RNA velocity done for the inputted, complete Seurat object.")

} else{
	print("Converting to h5ad... Use run_scvelo_full() function to complete RNA velocity")
}
  
# RNA velocity for specified conditions or time points, if any
# Function to process each condition
process_group <- function(group) {
  tpData <- object_annotated[ , object_annotated[[group_column]][[group_column]] %in% group]
  
  # Save as H5Seurat if it doesn't already exist
  seurat_obj<-tpData
  seurat_obj[[annotation_column]] <- as.character(seurat_obj[[annotation_column]][,1])
  Idents(seurat_obj)<-seurat_obj[[annotation_column]][,1]
  h5Seurat_file <- paste0("scvelo/rds/obj_spliced_unspliced_", paste(group, collapse = "_"), ".h5Seurat")
  if (!file.exists(h5Seurat_file)) {
    SeuratDisk::SaveH5Seurat(seurat_obj, filename = h5Seurat_file)
  }
  
  # Convert to h5ad if it doesn't already exist
  h5ad_file <- paste0("scvelo/rds/obj_spliced_unspliced_", paste(group, collapse = "_"), ".h5ad")
  if (!file.exists(h5ad_file)) {
    SeuratDisk::Convert(h5Seurat_file, dest = "h5ad")
  }

  if(Sys.info()["machine"] != "aarch64") {
  # Perform RNA velocity analysis
  tpV <- doVelocity(tpData, mode = mode)
  
  # Loop through grid resolutions
  for (grid_resolution in grid_resolutions) {
    tpVF <- getVectorField(tpData, tpV, reduction = "umap", resolution = grid_resolution)
    
    # Save the vector field object
    save(tpVF, file = paste0("scvelo/rds/", paste(group, collapse = "_"), "_gridRes", grid_resolution, ".RData"))
    
    # Generate plots for different arrow sizes and vector widths
    for (arrow_size in arrow_sizes) {
      for (vector_width in vector_widths) {
        plotVectorField(
          object_annotated, tpVF, group = group, group_column = group_column, 
          color_scale = color_scale, name_by = name_by, 
          grid_res = grid_resolution, arrow_size = arrow_size, vector_width = vector_width
        )
      }
    }
  }
  
  # Return status
  return(paste0("RNA velocity done for group ", paste(group, collapse = "_")))
  }
}

# Run in parallel
if (length(groups) != 0) {
  results <- parallel::mclapply(groups, process_group, mc.cores = cores)
  
  # Print results
  message("RNA velocity completed for all groups.")
  message(unlist(results))
}

sessionInfo()

}
