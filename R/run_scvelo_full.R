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
#' @param cores A numeric variable to set the number of cores to be used
#'
#' @return A list of scVelo data objects
#'
#' @export


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
  library(Matrix)
  import('scvelo')

  # create subdirectories in the output directory
  subdirectories <- c("scvelo",
                      "scvelo/csv",
                      "scvelo/rds",
                      "scvelo/images")

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

      # empty list to store objects with spliced/unspliced matrices
      # object_SU_list <- list()

      # # subset to corresponding cells in loom file
      # for (i in 1:length(loom_files)){
      #   object_subset <- object_annotated[, which(x = (expr == loom_file_subset_by[i]))]

      #   object_subset_SU <- addSUmatrices(object_subset, loom_files[i])
      #   DefaultAssay(object_subset_SU) <- "RNA"
      #   object_SU_list <- append(object_SU_list, object_subset_SU)
      # }

      # # merge subsetted seurat objects
      # cat("Merging Seurat objects with SU matrices...\n")
      # object_annotated <- suppressMessages(merge(object_SU_list[[1]], object_SU_list[-1], merge.dr=TRUE))
      # object_annotated <- RenameCells(object_annotated, new.names = object_annotated$orig.bc)
      # # join together individual layers
      # object_annotated <- JoinLayers(object_annotated)
      # object_annotated <- JoinLayers(object_annotated, assay = "spliced")
      # object_annotated <- JoinLayers(object_annotated, assay = "unspliced")

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

    saveRDS(object_annotated, "scvelo/rds/obj_spliced_unspliced.rds")

  }

  # # Match umap cell order to count matrix
  # cell_order <- colnames(object_annotated)
  # object_annotated@reductions$umap@cell.embeddings <- object_annotated@reductions$umap@cell.embeddings[cell_order, ]

  # # Align the gene orders in RNA, spliced, and unspliced assays
  # common_genes <- Reduce(intersect, list(
  #   rownames(object_annotated[["RNA"]]),
  #   rownames(object_annotated[["spliced"]]),
  #   rownames(object_annotated[["unspliced"]])
  # ))
  # for (assay_name in c("RNA", "spliced", "unspliced")) {
  #   object_annotated[[assay_name]] <- subset(object_annotated[[assay_name]], features = common_genes)
  # }
  # Add the originalexp assay if not there
  if (!"originalexp" %in% Assays(object_annotated)) {
    object_annotated[["originalexp"]] <- object_annotated[["RNA"]]
  }

  # Save to h5ad to conduct scvelo downstream analysis
  library(anndataR)
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
  h5ad_file <- file.path("scvelo/rds/obj_spliced_unspliced.h5ad")
  anndataR::write_h5ad(adata, path=h5ad_file, mode="w")


  # Call the main python function from scvelo_workflow.py with parameters
  cat("Starting velocity analysis...\n")
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
    file_base <- gsub(".h5ad", "", basename(h5ad_file))
    input_dir <- dirname(h5ad_file)
    for (group in groups){
      group_label=paste(group, collapse="_")
      file_name=paste0(file_base,"_",group_label,".h5ad")
      h5ad_group_file=file.path("scvelo/rds", file_name)
      # Create a subsetted h5ad file and save it to the output directory
      subset_adata <- adata[adata$obs[[group_column]] %in% group,]
      anndataR::write_h5ad(subset_adata, path=h5ad_group_file, mode="w")

      run_scvelo_workflow(h5ad_file=h5ad_group_file, annotation_column=annotation_column, mode=mode, top_gene=top_gene, group_label=group_label, output_format=output_format, n_jobs=cores)
    }
  }
  system("stty sane")
  #system("stty echo")

  sessionInfo()

}
