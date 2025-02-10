## scDown: a pipeline for scRNASeq downstream analysis

### Installation
The `velociraptor` package, along with several dependencies, needs to be installed via BiocManager before installing scDown, if not already available in your R library. 
```r
if (!require("BiocManager", quietly = TRUE)){
   install.packages("BiocManager")
}

BiocManager::install(c("pcaMethods", "velociraptor"), dependencies = TRUE)
BiocManager::install(c("Biobase", "BiocNeighbors", "BiocGenerics"))
BiocManager::install(c('DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array','terra', 'ggrastr'))
```

Once all BiocManager dependencies are installed, scDown can be installed using remotes: 
```r
install.packages("remotes")
remotes::install_github("BCH-RC/scDown")
```


### Tutorial 

Below are the key functions in scDown, with links to their vignetts for detailed usage instructions and example outputs:
- [`run_scproportion`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/scProportionTest.html) - Implements scProportionTest to statistically assess the significance of differences in cell type proportions between all condition pairs. 
- [`run_cellchatV2`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/scDown_CellChatV2.html) - Utilizes CellChat V2 to perform comprehensive intercellular communications analysis based on ligand-recptor pair interactions across cell types. 
- [`run_monocle3`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/scDown_monocle.html) - Leverages Monocle3 to construct pseudotime trajectories to model the progression of cellular differentiation. 
- [`run_scvelo`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/run_scvelo.html) - Integrates velocyto.R to incoporate spliced and unspliced counts to Seurat object and utilizes velociraptor to estimate RNA velocity by examining the ratio of unspliced and spliced mRNAs.
- [`run_scvelo_full`](https://html-preview.github.io/?url=https://github.com/BCH-RC/scDown/blob/main/vignettes/run_scvelo_full.html) - Calls the original scVelo for RNA velocity analysis from .h5ad files, providing enhanced visualizations and PAGA trajectory inference.


## Example Commands

To run cellchat analysis in p1, modify the resources and email address in **control_p1_cellchat.sh** as needed. After defining all the required variables `cd` to the directory with the pipeline code and run:
```
# run cellchat analysis
sbatch control_p1_cellchat.sh
```

To run monocle analysis in p2, modify the resources, email address, **absolute path for the directory with input data**, and **absolute path for the directory with pipeline code** in **control_p2_monocle.sh** as needed. After defining all the required variables `cd` to the directory with the pipeline code and run:
```
# run monocle analysis
sbatch control_p2_monocle.sh
```

To run RNA velocity in p3.1, modify the resources and email address in **control_p3.1_RNAvelocity.sh** as needed. After defining all the required variables `cd` to the directory with the pipeline code and run:
```
# run RNA velocity
sbatch control_p3.1_RNAvelocity.sh
```

To run RNA velocity in p3.2, modify the resources and email address in **control_p3.2_RNAvelocity.sh** as needed. After defining all the required variables in **p3.2_RNAvelocity.py**, `cd` to the directory with the pipeline code and run:
```
# run complete RNA velocity in python
sbatch control_p3.2_RNAvelocity.sh
```

- As for now **p3.2_RNAvelocity.py** can only be run after running p3.1. The other scripts can be run in any order.
- An example input of all the variables using p147 data is in **universal_variables.R** for p1, p2, p3.1 and **p3.2_RNAvelocity.py** for p3.2.
- An example input of binding directories on E2 to the singularity container using p147 data is in **control_p2_monocle.sh**.
- In the folder **universal_variables_example**, there are additional examples of inputs using different project data.

## Outputs Explained

Google slides with all outputs from p147 test run: https://docs.google.com/presentation/d/15pUJRGZqxNDHssmXAaCgNS7FH72oK3qw/edit#slide=id.p1

All outputs of the pipeline will be in the **results** folder, which is divided into:
- csv
- images
- figures
- pipeline_output
- rds

The **csv** subfolder contains:
- **cellchat**: inferred cell-cell communications at the level of ligands/receptor, at the level of signaling pathways, and csv files of differentially expressed ligand/receptor pairs.
- **monocle**: full DEG analysis results and from them the significant DEGs per specified model and trajectory.
- **RNA_Velocity**: differential velocity genes and paga transition confidence matrix from scvelo python implementation.

The **pipeline_output** subfolder contains the standard output of the scripts.

The **rds** subfolder contains:
- **cellchat**: whole and subsetted CellChat objects. Merged CellChat object used for pairwise comparison between conditions.
- **monocle**: whole and subsetted cell_data_set objects, as well as whole and subsetted cell_data_set objects per condition.
- **RNA_Velocity**: The .RData files with velocity vector coordinates, h5ad file of seurat object with spliced/unspliced matrices, h5ad file of AnnData object (from python) with scvelo velocity calculation results.
- Seurat object with cell type labels transferred from the reference in case symphony is ran.

The **figures** subfolder contains all graphs outputted by scvelo in p3.2. This is named as 'figures' by itself instead of 'images' because scvelo package by default generates a 'figures' directory for storing all its outputs.

The **images** subfolder contains all figures produced by the pipline from p1, p2, and p3.1. It's divided into **images/cellchat**, **images/monocle**, **images/RNA_velocity**, and **images/symphony**:

In **images/cellchat** :
- For each single condition:
    - **cellchat/aggregate**: contains figures for signaling network by aggregating all L-R pairs:
        - number of interactions & the total interaction strength (weights) between any two cell groups.
        - signaling sent from each cell group to all other cell types.
        - predict if each cell type is the dominant sender or receiver of communication.
        - which signals contribute the most to the outgoing or incoming signaling of each cell group.
    - **cellchat/pathway**: contains figures for signaling network focused on specific signaling pathways:
        - signaling network focused on one pathway (circle plot, chord disgram, and heatmap).
        - predict if each cell type is a dominant sender or receiver of communication for this pathway.
    - **cellchat/pathway/LR_gene**: contains figures for specific ligand/receptor pairs of the pathway:
        - contribution of specific ligand/receptor pairs to this pathway.
        - signaling network for those ligand/receptor pairs.
        - gene expression of ligand/receptor pairs in all cell types.

- For comparisons between two conditions:
    - **cellchat/comparison/Net**:
        - total number of interactions and interaction strength of the communication networks from different biological conditions.
        - differential number of interactions and strength in heatmap.
        - compare the outgoing and incoming interaction strength in 2D space of two conditions.
        - specific signaling changes of each cell type between two conditions.
    - **cellchat/comparison/manifold**: 
        - functional and structural similarity clustering umaps.
        - pathway distance ranks in the joint manifold umaps.
    - **cellchat/comparison/infoFlow**: 
        - top signaling pathways enriched in each conditions.
        - compare the outgoing/incoming signaling patterns between two datasets to identify signaling pathways/ligand-receptors that exhibit different signaling patterns.
    - **cellchat/comparison/sidebyside**: 
        - visualize specific signaling pathways' communications side-by-side across conditions

In **images/monocle** :
- **monocle/cellDistribution**: cell and cell type distribution density plots plus histogram per (subsetted) object and condition.
- **monocle/DEG**: 
    - significant DEGs found by regression analysis on predefined models (defined in universal variable @DEG_models), visualized in feature plot and violin plot. _This is a general linear regression method, and might work better for models that have continuous values_.
    - significant DEGs found by graph auto-correlation analysis along the trajectory/pseudotime, visualized in feature plot, and plotted their expression level changes along pseudotime.
        - These trajectory-variable genes are subsetted to run DEG analysis with regression again, in an effort to identify genes that are differentially expressed along pseudotime AND between two conditions.
- **monocle/pseudotime**: umap by cell types, partitions, learned trajectory, and pseudotime.

In **images/RNA_velocity** :
- Umap overlayed with velocity vectors for the whole object, with specified grid resolution and arrow sizes.
- If any timepoints are specified, umap overlayed with velocity vectors for specific timepoints.

In **images/symphony** :
- UMAP of the Seurat object labeled by active idents before running symphony.
- UMAP of the Seurat object labeled by active idents after running symphony.

## Potential TODOS

- General:
    - Allow h5ad file input into the pipeline. Currently there is an incomplete h5adToSeurat() function inside **utility_functions.R**.
- CellChat: 
    - CellChat has multiple tutorials on Github and many visualization options. Only a secion is implemented here from "Full tutorial for CellChat analysis of a single dataset with detailed explanation of each function" and "Full tutorial for comparison analysis of multiple datasets."
- Monocle3: 
    - Currently the pipeline only generates one trajectory across all partitions. Trajectory per partition is not implemented.
    - TradeSeq library and its dependencies are also installed in the current monocle3 docker image. So can integrate TradeSeq into the pipeline if needed.
- RNA velocity: 
    - scVelo has other tutorials in python that are not implemented here. Specifically these are analysis in 'Dynamical Modeling' and 'Differential Kinetics' that are available if using the dynamical mode to calculate velocity.
