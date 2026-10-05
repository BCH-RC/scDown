## scDown V2: a pipeline for scRNASeq downstream analysis

### 1 Run the Singularity image on E3 server
```r
# To run singularity image of scDown on HPC server
export TMPDIR="/path/to/tmp/directory/that/is/big/enough"
export SINGULARITY_CACHEDIR="/path/to/tmp/directory/that/is/big/enough"
export APPTAINER_CACHEDIR="/path/to/tmp/directory/that/is/big/enough"
cd /path/to/your/working/directory
singularity exec -B /path/to/your/input/data/directory:/input_dir -e /path/to/singularity/image/scdownv2.sif R --vanilla
library(scDownV2)
```
A local singularity image is at /lab-share/RC-DST-Bioinfo-e2/Public/DST_software/scdownv2/scdownv2.sif

### 2 Usage
Each key function in **scDown** is a wrap-up function of a workflow. Below are the main categories of key functions:
1. Preprocessing: doTransferLabel, cleanSeuratMeta, seuratToH5ad, h5adToSeurat
2. Compositional analysis: scProportionTest, scComp, scCODA
3. Cell-cell communications: CellChat V2, CellPhoneDB
4. Trajectory inference: Monocle3, scVelo and CellRank2
5. Gene set enrichment analysis: GSEApy

The test data used in the scDown vignettes is scRNA-seq data using 10X Genomics Chromium described in [Hochgerner et al. (2018)](https://www.nature.com/articles/s41593-017-0056-2). It is from dentate gyrus, a part of the hippocampus. The data consists of 25,919 genes across 2,930 cells with two time points. We converted the h5ad file of the dentate gyrus data (10X43_1.h5ad) to Seurat object (10X43_1_spliced_unspliced.rds) using our `h5adToSeurat` function in Preprocessing.  

#### 2.1 Preprocessing functions
- `h5adToSeurat` - Convert h5ad to Seurat rds as input for R functions.
- `seuratToH5ad` - Convert Seurat to anndata H5AD as input for python functions and optionally incorporate unspliced and spliced counts from loom files with velocyto.
- `cleanSeuratMeta` - Clean up unsafe characters in the Seurat meta.data and Idents.
- `doTransferLabel` - Transfers cell type annotation from a reference Seurat object to a query unannotated Seurat object, enabling automated annotation based on known cell types in reference scRNA-seq data.

#### 2.2 How to run each R function with the provided test data
- Load Seurat object and define variables
```r
library(scDownV2)
rds_path <- system.file("extdata", "10X43_1_spliced_unspliced.rds", package = "scDownV2")
seurat_obj <- readRDS(rds_path)
output_dir <- "."
annotation_column <- "clusters"
group_column <- "age.days."
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)  # Create result root folder
```

- `run_scproportion`- Uses scProportionTest to statistically assess the significance of differences in cell type proportions between all condition pairs. 
```r
run_scproportion(seurat_obj = seurat_obj,
                 annotation_column = annotation_column,
                 group_column = group_column, 
                 output_dir = output_dir,
                 n_jobs=2)
```
- `run_sccomp`- Uses scComp to perform robust cell-type compositional analysis based on outlier-aware statistical modeling across experimental conditions. (**For reproducibility, please use n_jobs=1.**)
```r
sample_column <- "age.days."

run_sccomp(seurat_obj = seurat_obj,
           annotation_column = annotation_column,
           sample_column = sample_column,
           group_column = group_column, 
           output_dir = output_dir,
           n_jobs=1)
```
- `run_cellchatV2` - Uses CellChat V2 to perform comprehensive intercellular communications analysis based on ligand-recptor pair interactions across cell types. 
```r
species <- "mouse"
celltypes_of_interest <- c("Granule immature", "Radial Glia-like", "Granule mature", "Neuroblast", "Microglia", "Cajal Retzius", "OPC", "Cck-Tox")
comparison_groups <- list(c("35", "12"))

run_cellchatV2(seurat_obj = seurat_obj,
               species = species,
               annotation_column = annotation_column,
               annotation_selected = celltypes_of_interest,
               group_column = group_column, 
               group_cmp = comparison_groups,
               output_dir = output_dir,
               n_jobs=8)
```
- `run_monocle3` - Uses Monocle3 to construct pseudotime trajectories to model the progression of cellular differentiation. 
```r
species <- "mouse"
nDim <- 30 
groups <- c("12","35") 

run_monocle3(seurat_obj = seurat_obj,
             species = species,
             nDim = nDim,
             annotation_column = annotation_column,
             graph_test = TRUE,
             groups = groups,
             group_column = group_column,
             deg_method = "quasipoisson",
             output_dir = output_dir,
             n_jobs = 8)
```
- `run_clusterProfiler` - Uses clusterProfiler to perform pathway enrichment analysis (KEGG, GO_BP, GO_CC, GO_MF). 
```r
species <- "mouse"

run_clusterProfiler(seurat_obj = seurat_obj,
                    species = species,
                    annotation_column = annotation_column,
                    group_column = group_column, 
                    output_dir = output_dir,
                    n_jobs = 2)
```

#### 2.3 How to run each python function with the provided test data
- `sccoda_workflow` - Calls scCODA to perform compositional data analysis based on Bayesian modeling to identify shifts in relative cell-type proportions across biological states.
```
singularity exec -e /pathto/scdownv2.sif python /app/sccoda_workflow.py --h5ad_file /opt/scdownv2/extdata/10X43_1.h5ad --out_path . --annotation_column clusters --sample_column "age(days)" --group_column "age(days)"
```
- `cellphonedb_workflow` - Calls CellPhoneDB to perform cell-cell communication analysis based on a curated repository of multi-subunit ligand-receptor complexes across interacting cell types. (**For reproducibility, please use n_jobs=1. To reproduce chord diagrams, use** `APPTAINERENV_PYTHONHASHSEED=42 singularity exec`)
```
singularity exec -e /pathto/scdownv2.sif python /app/cellphonedb_workflow.py --h5ad_file /opt/scdownv2/extdata/10X43_1.h5ad --cpdb_file_path /app/v5.0.0/cellphonedb.zip --out_path . --annotation_column clusters --species mouse --method 2 --n_jobs 8
```
- `scvelo_workflow` - Calls scVelo to perform RNA velocity analysis based on the dynamical modeling of spliced and unspliced transcripts to infer developmental trajectories. Provides PAGA trajectory inference.
```
singularity exec -e /pathto/scdownv2.sif python /app/scvelo_workflow.py --h5ad_file /opt/scdownv2/extdata/10X43_1.h5ad --out_path . --annotation_column clusters --mode dynamical --n_top_genes 10 --infer_per_group --group_column "age(days)" --n_jobs 8
```
 - `cellrank_workflow` - Calls CellRank2 to perform trajectory analysis based on cellular state transition probabilities to map fate decisions and identify driver genes.
 ```
singularity exec -e /pathto/scdownv2.sif python /app/cellrank_workflow.py --out_path . --annotation_column clusters --terminal_cell_types "Granule mature" "Astrocytes" "OL" "GABA" "Microglia" "Endothelial" "Cajal Retzius" --initial_cell_type "Radial Glia-like" --lineages_of_interest "Granule mature" "Astrocytes" "Microglia" --cluster_genes --n_jobs 8
 ```
 - `gseapy_workflow` - Calls GSEApy to perform Gene Set Enrichment Analysis (GSEA) and Over-Representation Analysis (ORA) in Python.
 ```
singularity exec -e /pathto/scdownv2.sif python /app/gseapy_workflow.py --h5ad_file /opt/scdownv2/extdata/10X43_1.h5ad --species mouse --out_path . --annotation_column clusters --group_column "age(days)" --control_cond 12 --output_format png --selected_cell_types "Astrocytes" "Granule immature" --cutoff_dotplot 0.05 --n_jobs 8
 ```