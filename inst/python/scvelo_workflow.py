import argparse
import warnings
from pathlib import Path
from typing import Literal

import pandas as pd
import scanpy as sc
import scvelo as scv

warnings.filterwarnings("ignore", category=DeprecationWarning)

def velocity_calculation(adata, annotation_column, out_path, mode='dynamical', group_label="ALL", n_jobs=2):
    """Velocity calculations by scVelo

    Args:
        adata (AnnData): anndata object with cell type annotation
        annotation_column (str): cell type column in adata.obs
        out_path (str): output path.
        mode (str, optional): velocity analysis mode ('deterministic', 'stochastic', or 'dynamical'). Defaults to 'dynamical'.
        group_label (str, optional): default "ALL" for entire data, or can be specific group. Defaults to "ALL".
        n_jobs (int, optional): number of parallel processes to use. Defaults to 2.
    """
    # change the column of metadata to categorical to comply with proportion plotting
    adata.obs[annotation_column] = adata.obs[annotation_column].astype('category')
    # visualize proportions of spliced/unspliced counts
    scv.pl.proportions(adata, groupby=annotation_column, fontsize=8, figsize=(12, 9), dpi=300, show=False, save=out_path+f'/plots/{group_label}_proportions')

    # preprocessing
    scv.pp.filter_genes(adata, min_shared_counts=20)
    adata.layers["counts"] = adata.X.copy()
    scv.pp.normalize_per_cell(adata, enforce=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, layer="counts", flavor='seurat_v3', n_top_genes=2000)

    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    if mode == 'dynamical':
        scv.tl.recover_dynamics(adata, n_jobs=n_jobs) # required if running dynamical model
    scv.tl.velocity(adata, mode=mode)
    scv.tl.velocity_graph(adata)
    if mode == 'dynamical':
        scv.tl.latent_time(adata)

    # save adata object after velocity calculation, since these results can take time to re-run.
    if adata._raw is not None:
        adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'}) # for getting around a bug
    adata.write_h5ad(out_path+f'/data/{group_label}_with_velocity_{mode}.h5ad')

    # set parameters for plotting
    kwargs = {"color": annotation_column, "dpi": 300, "figsize": (10, 10), "show": False}
    scv.pl.velocity_embedding_stream(adata, basis='umap', save=out_path+f'/plots/{group_label}_{mode}_embedding_stream', **kwargs)
    scv.pl.velocity_embedding_stream(adata, basis='umap', legend_loc='right margin', save=out_path+f'/plots/{group_label}_{mode}_embedding_stream_legend', **kwargs)
    scv.pl.velocity_embedding_grid(adata, basis='umap', save=out_path+f'/plots/{group_label}_{mode}_embedding_grid', **kwargs)
    scv.pl.velocity_embedding(adata, arrow_length=5, arrow_size=1, basis='umap', save=out_path+f'/plots/{group_label}_{mode}_embedding_arrow', **kwargs)
    if mode == 'dynamical':
        # also plot latent time for dynamical mode
        kwargs = {"figsize": (10, 10), "dpi": 300, "show": False}
        scv.pl.scatter(adata, basis='umap', color='latent_time', color_map='gnuplot', size=50, save=out_path+f'/plots/{group_label}_{mode}_latent_time', **kwargs)


def differential_velocity_genes(adata, annotation_column, out_path, top_gene=5, mode='dynamical', group_label="ALL"):
    """Performs differential velocity t-test to find genes that explain the directionality of velocity vectors

    Args:
        adata (AnnData): anndata object with cell type annotation
        annotation_column (str): cell type column in adata.obs
        out_path (str): output path.
        top_gene (int, optional): _description_. Defaults to 5.
        mode (str, optional): velocity analysis mode ('deterministic', 'stochastic', or 'dynamical'). Defaults to 'dynamical'.
        group_label (str, optional): default "ALL" for entire data, or can be specific group. Defaults to "ALL".
    """
    # perform differential velocity t-test
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scv.tl.rank_velocity_genes(adata, groupby=annotation_column, min_corr=.3)

    # extract top-ranking genes into pandas dataframe
    df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
    df.to_csv(out_path+f'/tables/{group_label}_{mode}_differential_velocity_genes_by_{annotation_column}.txt', sep="\t")

    # plot top 'top_gene' number of genes' phase portrait for each category in 'annotation_column'
    top_gene = int(top_gene)
    for item in adata.obs[annotation_column].unique():
        scv.pl.scatter(
            adata, 
            df[str(item)][:top_gene], 
            ylabel=item, 
            frameon=False, 
            linewidth=1.5, 
            fontsize=8, 
            save=out_path+f'/plots/{group_label}_{mode}_{item}_gene_phase', 
            color=annotation_column, 
            figsize=(2, 2), 
            dpi=500, 
            show=False
        )

        # convert from pandas series to list
        # Note: need to set colorbar=False below to bypass an error caused by matplotlib, in the generated figures darker colors indicate higher expression/velocity
        genes_to_plot = df[item][:top_gene].tolist()
        scv.pl.velocity(
            adata, 
            genes_to_plot, 
            colorbar=False, 
            ncols=2, 
            save=out_path+f'/plots/{group_label}_{mode}_{item}_gene_phase_complete_info', 
            color=annotation_column, 
            figsize=(6, 6), 
            dpi=500, 
            show=False
        )


def PAGA_trajectory_inference(adata, annotation_column, out_path, mode='dynamical', group_label="ALL"):
    """_summary_

    Args:
        adata (AnnData): anndata object with cell type annotation
        annotation_column (str): cell type column in adata.obs
        out_path (str): output path.
        mode (str, optional): velocity analysis mode ('deterministic', 'stochastic', or 'dynamical'). Defaults to 'dynamical'.
        group_label (str, optional): default "ALL" for entire data, or can be specific group. Defaults to "ALL".
    """
    # perform PAGA calculation
    scv.tl.paga(adata, groups=annotation_column)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
    df.to_csv(out_path+f'/tables/{group_label}_{mode}_paga_transition_confidence_matrix.txt', sep="\t")

    # generate a directed graph superimposed onto the UMAP embedding
    scv.pl.paga(
        adata, 
        basis='umap', 
        dashed_edges=None, 
        size=50, 
        alpha=0.05, 
        min_edge_width=2, 
        node_size_scale=1.5, 
        figsize=(10, 10), 
        dpi=300, 
        show=False, 
        save=out_path+f'/plots/{group_label}_{mode}_paga_graph'
    )


def run_scvelo_workflow(
        h5ad_file: str | None = None,
        out_path: str = ".",
        annotation_column: str = "clusters",
        mode: Literal["deterministic", "stochastic", "dynamical"] = "dynamical",
        n_top_genes: int = 5,
        infer_per_group: bool = False,
        group_column: str | None = None,
        output_format: str = "png",
        n_jobs: int = 2,
    ) -> None:
    """Main function of scVelo workflow

    Args:
        h5ad_file (str, optional): input h5ad file. Defaults to None.
        out_path (str, optional): output path. Defaults to '.'.
        annotation_column (str, optional): cell type column in adata.obs. Defaults to 'clusters'.
        mode (str, optional): velocity analysis mode ('deterministic', 'stochastic', or 'dynamical'). Defaults to 'dynamical'.
        n_top_genes (int, optional): number of top genes to show on plots. Defaults to 5.
        infer_per_group (bool, optional): whether to run analysis on each individual group. Defaults to False.
        group_column (str, optional): condition (group) column in adata.obs. Defaults to None.
        output_format (str, optional): output figure format. Defaults to "png".
        n_jobs (int, optional): number of parallel processes to use. Defaults to 2.
    """
    out_path = out_path+'/scvelo'
    for folder in [out_path, 
                   out_path+'/data',
                   out_path+'/tables',
                   out_path+'/plots']:
        Path(folder).mkdir(parents=True, exist_ok=True)

    scv.settings.verbosity = 3
    scv.settings.figdir = '.'
    scv.set_figure_params('scvelo', dpi=100, dpi_save=300, format=output_format, facecolor='white', transparent=False, frameon=False)

    # reading data
    if h5ad_file is None:
        h5ad_file = out_path+'/data/obj_spliced_unspliced.h5ad'
    adata = sc.read_h5ad(h5ad_file)
    print(adata)

    groups = ['ALL']
    if infer_per_group:
        if group_column is None:
            raise ValueError('Please specify "group_column" in adata.obs when infer_per_group=True.')
        if group_column not in adata.obs:
            raise ValueError('"group_column" cannot be found in adata.obs.')
        groups += list(adata.obs[group_column].unique())

    print(f'Performing analysis for {groups}')

    # group-wise velocity analysis
    for group_label in groups:
        if group_label == 'ALL':
            adata_sub = adata.copy()
        else:
            adata_sub = adata[adata.obs[group_column]==group_label].copy()

        # calculate RNA velocity using scVelo workflow
        velocity_calculation(adata_sub, annotation_column=annotation_column, out_path=out_path, mode=mode, group_label=group_label, n_jobs=n_jobs)

        # cluster-specific differential velocity genes
        differential_velocity_genes(adata_sub, annotation_column=annotation_column, out_path=out_path, top_gene=n_top_genes, mode=mode, group_label=group_label)

        # trajectory inference using PAGA
        PAGA_trajectory_inference(adata_sub, annotation_column=annotation_column, out_path=out_path, mode=mode, group_label=group_label)

    print('scVelo analysis completed.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scVelo workflow.")

    parser.add_argument(
        "--h5ad_file", 
        type=str, 
        default=None, 
        help="Input h5ad file path."
    )
    parser.add_argument(
        "--out_path", 
        type=str, 
        default=".", 
        help="Output directory path (default: '.')."
    )
    parser.add_argument(
        "--annotation_column", 
        type=str, 
        default="clusters", 
        help="Cell type column in adata.obs (default: 'clusters')."
    )
    parser.add_argument(
        "--mode", 
        type=str, 
        choices=["deterministic", "stochastic", "dynamical"], 
        default="dynamical", 
        help="Velocity analysis mode (default: 'dynamical')."
    )
    parser.add_argument(
        "--n_top_genes", 
        type=int, 
        default=5, 
        help="Number of top genes to show on plots (default: 5)."
    )
    parser.add_argument(
        "--infer_per_group", 
        action="store_true", 
        help="Flag to run analysis on each individual group. Include this flag to set to True."
    )
    parser.add_argument(
        "--group_column", 
        type=str, 
        default=None, 
        help="Condition (group) column in adata.obs (default: None)."
    )
    parser.add_argument(
        "--output_format", 
        type=str, 
        default="png", 
        help="Output figure format (default: 'png')."
    )
    parser.add_argument(
        "--n_jobs", 
        type=int, 
        default=2, 
        help="Number of parallel processes to use (default: 2)."
    )

    args = parser.parse_args()

    # Pass all arguments to the function as keyword arguments
    run_scvelo_workflow(**vars(args))
