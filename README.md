## scDown: a pipeline for scRNASeq downstream analysis

### Table of Contents
- [1 Installation](#1-installation)
  - [1.1 Installation using Docker](#11-installation-using-docker)
  - [1.2 Installation using Singularity](#12-installation-using-singularity)
- [2 Tutorial](#2-tutorial)
  - [2.1 Preprocess Tutorial](#21-preprocess-tutorial)
  - [2.2 Functions Tutorial](#22-functions-tutorial)

### 1 Installation
#### 1.1 Installation using Docker
##### Docker images of scDown
We built the [docker image for scDown](https://hub.docker.com/repository/docker/rcbioinfo/scdown/general):
| Docker images | Platform | Supported Systems |
|----------|----------|----------|
| rcbioinfo/scdown::amd64 | linux/amd64 | Linux, Windows, Intel-based Mac (x86_64) |

##### Run the docker image of scDown on HPC server (amd64 platform)
```r
# To run docker image of scDown for amd64 (Linux)
cd /path/to/your/working/directory
docker run -it --platform linux/amd64 --rm -v /path/to/your/input/data/directory:/input_dir /path/to/your/working/directory:/workspace -w /workspace  --pid=rcbioinfo/scdown::amd64 R
library(scDown)
```

#### 1.2 Installation using Singularity 
##### Pull Docker image into Singularity
```r
# To pull docker image of scDown to be singularity image on HPC server
cd /path/to/singularity/image/
singularity pull docker://rcbioinfo/scdown::amd64
```

##### Run the Singularity image of scDown on HPC server
```r
# To run singularity image of scDown on HPC server
export TMPDIR="/path/to/tmp/directory/that/is/big/enough"
export SINGULARITY_CACHEDIR="/path/to/tmp/directory/that/is/big/enough"
export APPTAINER_CACHEDIR="/path/to/tmp/directory/that/is/big/enough"
cd /path/to/your/working/directory
singularity exec -B /path/to/your/input/data/directory:/input_dir /path/to/singularity/image/scdown_amd64.sif R
library(scDown)
```


### 2 Tutorial 
Each key function in **scDown** is a wrap-up function of a workflow. Below is a main flowchart of the key functions:
<img src="https://github.com/user-attachments/assets/51117ff1-dec9-47d0-9eca-794ca6e295ae" width="400" height="220">

The test data used in the scDown vignettes is from dentate gyrus, a part of the hippocampus. It is scRNA-seq data using 10X Genomics Chromium described in [Hochgerner et al. (2018)](https://www.nature.com/articles/s41593-017-0056-2). The data consists of 25,919 genes across 2,930 cells with two time points. We converted the original h5ad file of the dentate gyrus data (10X43_1.h5ad) to Seurat object (10X43_1_spliced_unspliced.rds) using our `h5adToSeurat` function in Preprocess tutorial below.  

#### 2.1 Preprocess Tutorial
First, we define universal variables that apply to all key functions in the following preprocess vignette:
- [Preprocess](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/refs/heads/main/vignettes/scDown_preProcess.html) - Set the universal variables for all key functions in scDown, and convert h5ad to Seurat rds or annotate cell type using reference scRNA-seq data if needed.
  - [`h5adToSeurat`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/refs/heads/main/vignettes/scDown_preProcess.html#annotated-seurat-object) - Convert h5ad to Seurat rds as input for key functions in **scDown**.
  - [`doTransferLabel`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/refs/heads/main/vignettes/scDown_preProcess.html#unannotated-seurat-object) - Transfers cell type annotation from a reference Seurat object to a query unannotated Seurat object, enabling automated annotation based on known cell types in reference scRNA-seq data.

#### 2.2 Functions Tutorial
The **scDown** package provides a single function for each purpose, integrating all necessary steps into one streamlined command, making the analysis more efficient and user-friendly. Below are the **key functions in scDown**, with links to their vignettes for detailed usage instructions and example outputs:
- [`run_scproportion`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/scProportionTest.html) - Implements scProportionTest to statistically assess the significance of differences in cell type proportions between all condition pairs. 
- [`run_cellchatV2`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/scDown_CellChatV2.html) - Utilizes CellChat V2 to perform comprehensive intercellular communications analysis based on ligand-recptor pair interactions across cell types. 
- [`run_monocle3`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/scDown_monocle.html) - Leverages Monocle3 to construct pseudotime trajectories to model the progression of cellular differentiation. 
- [`run_scvelo`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/run_scvelo.html) - Employs velocyto.R to incoporate spliced and unspliced counts to Seurat object and utilizes velociraptor to estimate RNA velocity by examining the ratio of unspliced and spliced mRNAs.
- [`run_scvelo_full`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/run_scvelo_full.html) - Calls the original scVelo for RNA velocity analysis from .h5ad files, providing enhanced visualizations and PAGA trajectory inference.
  
The latter 4 key functions in scDown can be applied to either entire data or selected conditions of interest. 

## Contact
- Please create an issue under our repository by clicking the issue tab, we will try to address it.

