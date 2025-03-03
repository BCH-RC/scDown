## scDown: a pipeline for scRNASeq downstream analysis

### Installation
The **scDown** package can be installed using `remotes`: 
```r
install.packages("remotes")
remotes::install_github("BCH-RC/scDown")
```

The `velociraptor` package, along with several dependencies, needs to be installed via BiocManager before installing `scDown`, if not already available in your R library. 
```r
required_pkgs <- c("pcaMethods", "velociraptor", "Biobase", "BiocNeighbors", "BiocGenerics",
                   "DelayedArray", "DelayedMatrixStats", "limma", "lme4", "S4Vectors", 
                   "SingleCellExperiment", "SummarizedExperiment", "batchelor", 
                   "HDF5Array", "terra", "ggrastr")

missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]

if (length(missing_pkgs) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install(missing_pkgs, dependencies = TRUE)
}
```


### Tutorial 

The **scDown** package provides a single function for each purpose, integrating all necessary steps into one streamlined command, making the analysis more efficient and user-friendly. Below are the **key functions in scDown**, with links to their vignettes for detailed usage instructions and example outputs:
- [`run_scproportion`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/scProportionTest.html) - Implements scProportionTest to statistically assess the significance of differences in cell type proportions between all condition pairs. 
- [`run_cellchatV2`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/scDown_CellChatV2.html) - Utilizes CellChat V2 to perform comprehensive intercellular communications analysis based on ligand-recptor pair interactions across cell types. 
- [`run_monocle3`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/scDown_monocle.html) - Leverages Monocle3 to construct pseudotime trajectories to model the progression of cellular differentiation. 
- [`run_scvelo`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/run_scvelo.html) - Employs velocyto.R to incoporate spliced and unspliced counts to Seurat object and utilizes velociraptor to estimate RNA velocity by examining the ratio of unspliced and spliced mRNAs.
- [`run_scvelo_full`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/run_scvelo_full.html) - Calls the original scVelo for RNA velocity analysis from .h5ad files, providing enhanced visualizations and PAGA trajectory inference.
  
The latter 4 key functions in scDown can be applied to either entire data or selected groups of interest. 

To faciliate the scRNA-seq downstream analysis using scDown, scDown provides a function for cell type annotation when reference data with cell type annotation is available, which can be run beforehand if needed: 
- `doTransferLabel(ReferenceSeuratObject, SeuratObject)` - Transfers cell type labeling from a reference Seurat object to a query Seurat object, enabling automated annotation based on known cell identities prior to downstream analysis.

## Contact
- Please create an issue under our repository by clicking the issue tab, we will try to address it.

