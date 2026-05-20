#' Function to run the scVelo pipeline using python script
#'
#' This function performs RNA velocity calculations from .h5ad file using scVelo python script.
#' Workflow of this function: 
#' 1. calculate RNA velocity using scVelo workflow
#' 2. cluster-specific differential velocity genes
#' 3. trajectory inference using PAGA
#' This function outputs basic figures such as spliced/unspliced count proportion, projected RNA velocity 
#' on umap, the phase portrait (ratio of spliced/unspliced RNA abundance) for top differential genes, and 
#' directed graphs of predicted lineages from PAGA trajectory inference 
#'
#' @param h5ad_file input h5ad file path and name, if running after `run_scvelo()`, this object has a fixed 
#' name and does not need to be changed
#' @param output_dir A character vector specifying the output directory
#' @param annotation_column A character variable specifying which metadata column of the h5ad object contains 
#' cell type annotations
#' @param mode Mode to conduct scvelo velocity calculation, either 'stochastic (default)', 'deterministic', 
#' or 'dynamical (slowest, recommended)'
#' @param top_gene The number of top differential velocity genes to plot phase portrait for
#' @param groups A list of character vectors representing groups of conditions or time points used to 
#' calculate RNA velocity separately, default: NULL
#' @param group_column A string specifying the name of the metadata column in the h5ad object that should be 
#' @param output_format Format of output figure: "png" or "pdf" (default: "png")
#' used for subsetting each group in `groups`
#'
#' @return A list of scVelo data objects
#'
#' @export
#'
#' Estimate RNA velocity for spliced and unspliced counts of scRNA-seq data


run_scvelo_full <- function(seurat_obj,
                            loom_files=NULL, 
                            output_dir=".", 
                            loom_file_subset_by=NULL,
                            loom_file_subset_column="orig.ident",
                            annotation_column='ID', 
                            mode='stochastic', 
                            top_gene=5,
                            groups=NULL, 
                            group_column=NULL,
                            output_format="png",
                            cores=4){

# Import scvelo python package now to avoid import error after write_h5ad
library(reticulate)
import('scvelo')

# create subdirectories in the output directory
subdirectories <- c("scvelo",
                    "scvelo/csv",
                    "scvelo/rds",
                    "scvelo/images")

for(i in subdirectories){
  dir.create(file.path(output_dir,i), showWarnings = FALSE, recursive = TRUE)
}
setwd(output_dir)

### Input
library(Seurat)
checkmate::test_class(seurat_obj, "Seurat")
object_annotated <- seurat_obj

# use cell type annotation column as identity
if(checkmate::test_string(annotation_column, null.ok=FALSE)){
  checkmate::expect_choice(annotation_column, colnames(seurat_obj@meta.data), label = "annotation_column")
  Seurat::Idents(object_annotated) <- object_annotated[[annotation_column]][,1]
}

checkmate::assert_list(groups, types = c("numeric", "integer","character"), null.ok = TRUE)
if(!is.null(groups)){
  groups <- lapply(groups, function(element) {
    if (is.integer(element) || is.numeric(element)) {
      return(as.character(element))
    }
    return(element)
  })
}
if(!is.null(groups)){
  checkmate::assert_string(group_column, null.ok = FALSE)
  checkmate::expect_choice(group_column, colnames(seurat_obj@meta.data), label = "group_column")
}

checkmate::assert_string(mode, null.ok = FALSE)
checkmate::assert_numeric(top_gene, null.ok = FALSE, any.missing=FALSE)
checkmate::expect_choice(output_format,c("png","pdf"),label = "output_format")

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

library(anndataR)
# save to h5ad so if needed, can be used to conduct scvelo downstream analysis
# To prevent overwriting error, only saves if the file does not exist
seurat_obj<-object_annotated
seurat_obj[[annotation_column]] <- as.character(seurat_obj[[annotation_column]][,1])
Idents(seurat_obj)<-seurat_obj[[annotation_column]][,1]

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
h5ad_file <- file.path(output_dir, "scvelo/rds/obj_spliced_unspliced.h5ad")
anndataR::write_h5ad(adata, path=h5ad_file, mode="w")


# Call the main python function from scvelo_workflow.py with parameters
py_script <- system.file("python", "scvelo_workflow.py", package = "scDown")
if (py_script == "") {
  stop("Python script not found in the installed scDown package")
}
source_python(py_script)

# RNA velocity for the entire object
run_scvelo_workflow(h5ad_file=h5ad_file, annotation_column=annotation_column, mode=mode, top_gene=top_gene, group_label="ALL", output_format=output_format, n_jobs=cores)
system("stty echo")

# RNA velocity for specified conditions or time points, if any
if (length(groups) != 0){
  library(anndata)
  file_base <- gsub(".h5ad", "", basename(h5ad_file))
  input_dir <- dirname(h5ad_file)
  for (group in groups){
      group_label=paste(group, collapse="_")
      file_name=paste0(file_base,"_",group_label,".h5ad")
      h5ad_group_file_input=file.path(input_dir, file_name)
      h5ad_group_file_output=file.path(paste0(output_dir,"/scvelo/rds"), file_name)
      # If the subsetted h5ad file already exists, use it. 
      # Otherwise, create a new subsetted h5ad file and save it to the output directory
      if(file.exists(h5ad_group_file_input)){
        h5ad_group_file <- h5ad_group_file_input
      } else {
        h5ad_group_file <- h5ad_group_file_output
        adata <- anndataR::read_h5ad(h5ad_file)
        subset_adata <- adata[adata$obs[[group_column]] %in% group,]
        anndataR::write_h5ad(subset_adata, path=h5ad_group_file, mode="w")
      }
      run_scvelo_workflow(h5ad_file=h5ad_group_file, annotation_column=annotation_column, mode=mode, top_gene=top_gene, group_label=group_label, output_format=output_format, n_jobs=cores)
  }
}
system("stty sane")
#system("stty echo")


sessionInfo()

}
