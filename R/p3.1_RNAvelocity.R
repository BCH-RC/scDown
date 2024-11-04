## RNA Velocity Pipeline
# Estimate RNA velocity for inputted Seurat objects.

### read arguments from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=0) {
    data_dir_name <- args[1]
    code_dir_name <- args[2]
    source(paste(code_dir_name,"/universal_variables.R", sep=""))
    source(paste(code_dir_name,"/RNAvelocity_functions.R", sep=""))
    source(paste(code_dir_name,"/utility_functions.R", sep=""))


### Input
object_annotated <- readRDS(file = paste0(data_dir_name,"/",seurat_object, sep=""))

# if the current input object does not have cell type labels, run symphony to populate Idents(X) with cell type annotations
if (transfer_label) {
    # read in reference object
    object_reference <- readRDS(file = paste0(data_dir_name,"/",reference_object, sep=""))
    # umap before label transfer
    p1 <- DimPlot(object_annotated)
    ggsave(file="images/symphony/objectUMAP_before_label_transfer.png", plot=p1)
    # perform label transfer
    object_annotated <- doTransferLabel(object_reference, object_annotated, varToHarmonize = varToHarmonize, transferCoordinates = transferCoordinates)
    # umap after label transfer
    p2 <- DimPlot(object_annotated)
    ggsave(file="images/symphony/objectUMAP_after_label_transfer.png", plot=p2)
    # save annotated object
    saveRDS(object_annotated, file="rds/seurat_object_annotated_by_symphony.rds")
}

## RNA velocity analysis

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
