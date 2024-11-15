#' Function to run the scVelo pipeline using velociraptor 
#'
#' This function performs RNA velocity calculations from .loom files with the scVelo package.
#' In this function, users can calculate RNA velocity of the whole data as well as a subset of time points.
#'
#' @param seurat_object Seurat object containing the scRNA-seq data (Required)
#' @param loom_files Spliced and unspliced counts of the scRNA-seq data (Required)
#' @param output_dir A character vector specifying the output directory
#' @param loom_file_subset_by A character variable specifying how the Seurat object should be subsetted in order to match the loom files - the order of conditions must match the order of loom files for them to be matched
#' @param loom_file_subset_column A character variable specifying which metadata column of the Seurat object should be used for subsetting to match each of the loom files
#' @param mode Mode for scVelo velocity calculation, default stochastic - can be one of four options: "steady_state" (original), "deterministic", "stochastic" (fastest:recommended), "dynamical"
#' @param grid_resolutions A vector of integers specifying the number of grids along each axis, essentially controlling the number of vectors on umap, default 50.
#' @param arrow_sizes A vector of integer or float controlling velocity vector size (arrow head size), default 0.5
#' @param vector_widths A vector of integer or float controlling velocity vector size (vector width), default 0.5
#' @param time_point A list of character vectors representing a group of time points used to calculate RNA velocity together, can be left blank
#' @param time_point_column A character variable specify which metadata column of the Seurat object should be used for subsetting the group of time points
#' @param color_scale A character vector of colors to be used in plotting, must match number of unique values in the metadata column marked by @name_by
#' @param name_by A character variable specify which metadata column of the Seurat object should be used for colors
#'
#' @return A list of scVelo data objects
#'
#' @export
#'
#' Estimate RNA velocity for spliced and unspliced counts of scRNA-seq data


run_scvelo <- function(seurat_object,loom_files,output_dir=".",loom_file_subset_by=c(),loom_file_subset_column=NULL,
                    mode='stochastic',grid_resolutions=c(50),arrow_sizes=c(0.5),vector_widths=c(0.5),
                    time_point=list(),time_point_column=NULL,color_scale=NULL,name_by=NULL){

# create subdirectories in the output directory
setwd(output_dir)
subdirectories <- c("rds",
                    "csv",
                    "images",
                    "csv/RNA_velocity",
                    "rds/RNA_velocity",
                    "images/RNA_velocity")

for(i in subdirectories){
    dir.create(file.path(output_dir,i), showWarnings = F, recursive = T)
}

### Input
object_annotated <- readRDS(file = seurat_object)

# add cell barcode as metadata
object_annotated$orig.bc <- colnames(object_annotated)

# add spliced and unspliced matrices as new assays
if (length(loom_files) > 1){
    
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

# save to h5ad so if needed, can be used to conduct scvelo downstream analysis
#To prevent overwriting error, only saves if the file does not exist
if (!file.exists("rds/RNA_velocity/obj_spliced_unspliced.h5Seurat")) {
    SaveH5Seurat(object_annotated, filename = "rds/RNA_velocity/obj_spliced_unspliced.h5Seurat")
}
if (!file.exists("rds/RNA_velocity/obj_spliced_unspliced.h5ad")) {
    Convert("rds/RNA_velocity/obj_spliced_unspliced.h5Seurat", dest = "h5ad")
}


### RNA velocity analysis

# RNA velocity for the entire Seurat object
tpV <- doVelocity(object_annotated, mode=mode)
for (grid_resolution in grid_resolutions){
    tpVF <- getVectorField(object_annotated, tpV, reduction = 'umap', resolution = grid_resolution)
    save(tpVF, file = paste0('rds/RNA_velocity/ALL_gridRes', grid_resolution,'.RData',sep=""))

    for (arrow_size in arrow_sizes){
        for (vector_width in vector_widths){
            plotVectorField(object_annotated, tpVF, color_scale=color_scale, name_by=name_by, grid_res=grid_resolution, arrow_size=arrow_size, vector_width=vector_width)
        }
    }
}

print("RNA velocity done for the inputted, complete Seurat object.")

# RNA velocity for specified time points, if any
if (length(time_point) != 0){
    for (time in time_point){
        tpData <- object_annotated[ ,object_annotated[[time_point_column]][[time_point_column]] %in% time]
        tpV <- doVelocity(tpData, mode=mode)
        for (grid_resolution in grid_resolutions){
            tpVF <- getVectorField(tpData, tpV, reduction = 'umap', resolution = grid_resolution)
            save(tpVF, file = paste0('rds/RNA_velocity/',paste(time, collapse="_"),'_gridRes',grid_resolution,'.RData',sep=""))

            for (arrow_size in arrow_sizes){
                for (vector_width in vector_widths){
                    plotVectorField(object_annotated, tpVF, time_point=time, time_point_column=time_point_column, color_scale=color_scale, name_by=name_by, grid_res=grid_resolution, arrow_size=arrow_size, vector_width=vector_width)
                }
            }
        }

        print(paste0("RNA velocity done for timepoints ", time, sep=""))
    }
}

sessionInfo()

}
