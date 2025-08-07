## scDown: a pipeline for scRNASeq downstream analysis (version 1)

### Table of Contents
- [1 Installation](#1-installation)
  - [1.1 System requirement](#11-system-requirement)
  - [1.2 Installation using Docker](#12-installation-using-docker)
  - [1.3 Installation using Singularity](#13-installation-using-singularity)
- [2 Tutorial](#2-tutorial)
  - [2.1 Preprocess Tutorial](#21-preprocess-tutorial)
  - [2.2 Functions Tutorial](#22-functions-tutorial)
- [3 Citation](#3-citation)
- [4 Contact](#4-contact)  

### 1 Installation
#### 1.1 System requirement
We strongly recommend using an HPC (High-Performance Computing) Linux server for running most functions of scDown, especially for large datasets. HPC Linux server typically provides greater memory and more CPU cores, which significantly improve efficiency and performance. While small datasets (<10,000 cells) can be processed on a Mac, limited memory may result in slower performance. The only exception is the `run_scproportiontest` function, which can be run efficiently on a Mac, regardless of dataset size.

#### 1.2 Installation using Docker
##### Docker images of scDown
We built the [docker images for scDown](https://hub.docker.com/r/rcbioinfo/scdown/tags) supporting different system architectures:
| Docker images | Platform | Supported Systems |
|----------|----------|----------|
| rcbioinfo/scdown:amd64 | linux/amd64 | Linux, Windows, Intel-based Mac (x86_64) |
| rcbioinfo/scdown:arm64 | linux/arm64 | Apple Silicon Mac (M1/M2/M3) |

##### Run the docker image of scDown on HPC server (amd64 platform)
```r
# To run docker image of scDown for amd64 (Linux)
cd /path/to/your/working/directory
docker run -it --platform linux/amd64 --rm -v /path/to/your/input/data/directory:/input_dir -v /path/to/your/working/directory:/workspace -w /workspace  rcbioinfo/scdown:amd64 R
library(scDown)
```

#### 1.3 Installation using Singularity 
##### Pull Docker image into Singularity
```r
# To pull docker image of scDown to be singularity image on HPC server
cd /path/to/singularity/image/
singularity pull docker://rcbioinfo/scdown:amd64
```

##### Run the Singularity image of scDown on HPC server
```r
# To run singularity image of scDown on HPC server
export TMPDIR="/path/to/tmp/directory/that/is/big/enough"
export SINGULARITY_CACHEDIR="/path/to/tmp/directory/that/is/big/enough"
export APPTAINER_CACHEDIR="/path/to/tmp/directory/that/is/big/enough"
cd /path/to/your/working/directory
singularity exec -B /path/to/your/input/data/directory:/input_dir /path/to/singularity/image/scdown_amd64.sif R --vanilla
library(scDown)
```


### 2 Tutorial 
Each key function in **scDown** is a wrap-up function of a workflow. Below is a main flowchart of the key functions:

<img src="https://github.com/user-attachments/assets/79e2ac3c-694f-4da1-95e6-8e8001eadd52" width="600" height="330">

The test data used in the scDown vignettes is scRNA-seq data using 10X Genomics Chromium described in [Hochgerner et al. (2018)](https://www.nature.com/articles/s41593-017-0056-2). It is from dentate gyrus, a part of the hippocampus. The data consists of 25,919 genes across 2,930 cells with two time points. We converted the original h5ad file of the dentate gyrus data (10X43_1.h5ad) to Seurat object (10X43_1_spliced_unspliced.rds) using our `h5adToSeurat` function in Preprocess Tutorial below.  

#### 2.1 Preprocess Tutorial
First, we define universal variables that apply to all key functions in the following preprocess vignette:
- [Preprocess](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/scDown_preProcess.html) - Set the universal variables for all key functions in scDown, and convert h5ad to Seurat rds or annotate cell type using reference scRNA-seq data if needed.
  - [`h5adToSeurat`](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/scDown_preProcess.html#annotated-seurat-object) - Convert h5ad to Seurat rds as input for key functions in **scDown**.
  - [`doTransferLabel`](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/scDown_preProcess.html#unannotated-seurat-object) - Transfers cell type annotation from a reference Seurat object to a query unannotated Seurat object, enabling automated annotation based on known cell types in reference scRNA-seq data.

#### 2.2 Functions Tutorial
The **scDown** package provides a single function for each purpose, integrating all necessary steps into one streamlined command, making the analysis more efficient and user-friendly. Below are the **key functions in scDown**, with links to their vignettes for detailed usage instructions and example outputs:
- [`run_scproportion`](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/scProportionTest.html) - Implements scProportionTest to statistically assess the significance of differences in cell type proportions between all condition pairs. 
- [`run_cellchatV2`](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/scDown_CellChatV2.html) - Utilizes CellChat V2 to perform comprehensive intercellular communications analysis based on ligand-recptor pair interactions across cell types. 
- [`run_monocle3`](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/scDown_monocle.html) - Leverages Monocle3 to construct pseudotime trajectories to model the progression of cellular differentiation. 
- [`run_scvelo`](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/run_scvelo.html) - Employs velocyto.R to incoporate spliced and unspliced counts to Seurat object and utilizes velociraptor to estimate RNA velocity by examining the ratio of unspliced and spliced mRNAs.
- [`run_scvelo_full`](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/run_scvelo_full.html) - Calls the original scVelo for RNA velocity analysis from .h5ad files, providing enhanced visualizations and PAGA trajectory inference.
  
The latter 4 key functions in scDown can be applied to either entire data or selected conditions of interest. 

### 3 Citation
Sun, L.; Ma, Q.; Cai, C.; Labaf, M.; Jain, A.; Dias, C.; Rockowitz, S.; Sliz, P. scDown: A Pipeline for Single-Cell RNA-Seq Downstream Analysis. Int. J. Mol. Sci. 2025, 26, 5297. https://doi.org/10.3390/ijms26115297
### 4 Contact
Please submit an issue via our GitHub repository by clicking the 'Issues' tab. We'll do our best to address it promptly.
