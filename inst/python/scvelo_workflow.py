import os
import numpy as np
import pandas as pd
import scvelo as scv
import scanpy as sc
import warnings
scv.settings.figdir = '.'
warnings.filterwarnings("ignore", category=DeprecationWarning)

# ========================== SCVELO FUNCTIONS ==========================
# This function takes in an AnnData object and performs all basic velocity calculations
# enabled by scVelo. It also outputs basic figures such as spliced/unspliced count proportion
# and RNA velocity vectors on umap.
#
def velocity_calculation(adata, annotation_column, mode='stochastic', group_label="ALL", n_jobs=None):
    """
    :param adata: an AnnData object with cell type annotation
    :param annotation_column: a key inside adata.obs to use for calculating spliced/unspliced count proportion, normally cell types
    :param mode: can be 'stochastic (default)', 'deterministic', or 'dynamical (slowest)'
    :param group_label: default "ALL" for entire data, or can be specific group from groups in run_scvelo_full.R
    :return: an AnnData object with RNA velocity calculated
    """
    # Ensure raw is initialized
    if adata.raw is None:
        adata.raw = adata
    # set parameters for plotting
    kwargs = dict(color=annotation_column, figsize=(10, 10), dpi=500, show=False)
    # change the 'annotation_column' column of metadata from dtype: object into dtype: category to comply with proportion plotting
    adata.obs[annotation_column] = adata.obs[annotation_column].astype('category')
    # observe proportions of spliced/unspliced counts
    scv.pl.proportions(adata, groupby=annotation_column, fontsize=8, figsize=(14, 10), dpi=500, show=False, save=f'scvelo/images/{group_label}_proportions')
    # velocity calculation workflow
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
    adata.write(f'scvelo/rds/{group_label}_withVelocity_{mode}.h5ad', compression='gzip')
    # basic RNA velocity visualizations in various formats
    scv.pl.velocity_embedding_stream(adata, basis='umap', save=f'scvelo/images/{group_label}_{mode}_embedding_stream', **kwargs)
    scv.pl.velocity_embedding_stream(adata, basis='umap', legend_loc='right margin', save=f'scvelo/images/{group_label}_{mode}_embedding_stream_legend', **kwargs)
    scv.pl.velocity_embedding_grid(adata, basis='umap', save=f'scvelo/images/{group_label}_{mode}_embedding_grid', **kwargs)
    scv.pl.velocity_embedding(adata, arrow_length=5, arrow_size=1, basis='umap', save=f'scvelo/images/{group_label}_{mode}_embedding_arrow', **kwargs)
    if mode == 'dynamical':
        kwargs = dict(figsize=(10, 10), dpi=500, show=False)
        scv.pl.scatter(adata, basis='umap', color='latent_time', color_map='gnuplot', size=50, save=f'scvelo/images/{group_label}_{mode}_latent_time', **kwargs)
    return


# This function takes in an AnnData object and performs a differential velocity
# t-test to find genes that explain the directionality of calculated velocity vectors. 
# It tests which genes have group-specific differential velocity expression (definied in @annotation_column)
# i.e being siginificantly higher/lower compared to the remaining population, and visualizes
# the phase portrait (ratio of spliced/unspliced RNA abundance) for highly ranked genes.
#
def differential_velocity_genes(adata, annotation_column, top_gene=5, mode='stochastic', group_label="ALL"):
    """
    :param adata: an AnnData object with velocity calculated
    :param annotation_column: a key inside adata.obs that provides a grouping for cells, ex) cell type
    :param top_gene: an integer specifying number of top ranked genes to plot
    :param mode: can be 'stochastic (default)', 'deterministic', or 'dynamical (slowest)'
    :param group_label: default "ALL" for entire data, or can be specific group from groups in run_scvelo_full.R
    :return: an AnnData object with new data in adata.uns['rank_velocity_genes'] and adata.var['spearmans_score']
    """
    # perform differential velocity t-test
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scv.tl.rank_velocity_genes(adata, groupby=annotation_column, min_corr=.3)
    # extract top-ranking genes into pandas dataframe
    df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
    df.to_csv(f'scvelo/csv/{group_label}_{mode}_differential_velocity_genes_by_{annotation_column}.csv')
    # set parameters for plotting
    # plot top 'top_gene' number of genes' phase portrait for each category in 'annotation_column'
    top_gene = int(top_gene)
    for item in adata.obs[annotation_column].unique():
        kwargs = dict(color=annotation_column, figsize=(2, 2), dpi=500, show=False)
        scv.pl.scatter(adata, df[str(item)][:top_gene], ylabel=item, frameon=False, linewidth=1.5, save=f'scvelo/images/{group_label}_{mode}_{item}_genePhase', fontsize=8, **kwargs)
        # convert from pandas series to list
        # Note: need to set colorbar=False below to bypass an error caused by matplotlib, in the generated figures darker colors indicate higher expression/velocity
        genes_to_plot = df[item][:top_gene].tolist()
        kwargs = dict(color=annotation_column, figsize=(6, 6), dpi=500, show=False)
        scv.pl.velocity(adata, genes_to_plot, colorbar=False, ncols=2, save=f'scvelo/images/{group_label}_{mode}_{item}_genePhaseCompleteInfo', **kwargs)
    return


# This function takes an AnnData object and performs trajectory inference using the
# PAGA method (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1663-x).
# It provides a graph-like map of the data with solid edges corresponding to the transition confidence 
# between two groups (defined in @annotation_column). Here, PAGA is extended by velocity-inferred directionality
# and predicts transitions/lineages between groups.
#
# Note: there is a possible issue with PAGA with small transition probabilities (< 0.1):
# https://github.com/theislab/scvelo/issues/456
#
def PAGA_trajectory_inference(adata, annotation_column, mode='stochastic', group_label="ALL"):
    """
    :param adata: an AnnData object with velocity and velocity graph calculated
    :param annotation_column: a key inside adata.obs specifying how the cells should be grouped, normally cell types
    :param mode: can be 'stochastic (default)', 'deterministic', or 'dynamical (slowest)'
    :param group_label: default "ALL" for entire data, or can be specific group from groups in run_scvelo_full.R
    :return: an AnnData object with paga graph calculated and stored in adata.uns
    """
    # this is needed due to a current bug in scvelo that hasn't been fixed.
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
    # perform PAGA calculation
    scv.tl.paga(adata, groups=annotation_column)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
    df.to_csv(f'scvelo/csv/{group_label}_{mode}_paga_transition_confidence_matrix.csv')
    # generate a directed graph superimposed onto the UMAP embedding
    scv.pl.paga(adata, basis='umap', dashed_edges=None, size=50, alpha=.05, min_edge_width=2, node_size_scale=1.5, figsize=(10, 10), dpi=500, show=False, save=f'scvelo/images/{group_label}_{mode}_paga_graph')
    return


# Main function of scvelo workflow in python
# 
def run_scvelo_workflow(h5ad_file='scvelo/rds/obj_spliced_unspliced.h5ad', annotation_column='ID', mode='stochastic', top_gene=5, group_label="ALL", output_format = "png", n_jobs=None):
    """
    :param h5ad_file (str): Base h5ad file path.
    :param annotation_column (str): Annotation column name.
    :param mode (str): scvelo mode.
    :param top_gene (int): Number of top genes.
    :param groups (list): List of character vectors, where each vector defines a group.
    :param output_format (str): Format of output figure: "png" or "pdf" (default: "png")
    """
    # basic scvelo settings
    scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
    # set figure format for visualization
    if output_format == "png":
        scv.set_figure_params('scvelo', transparent=False, format='png')
    elif output_format == "pdf":
        scv.set_figure_params('scvelo', transparent=False, format='pdf')
    # reading data
    adata = sc.read(h5ad_file)
    # Workflow:
    # 1. calculate RNA velocity using scVelo workflow
    velocity_calculation(adata, annotation_column=annotation_column, mode=mode, group_label=group_label, n_jobs=n_jobs)
    # 2. cluster-specific differential velocity genes
    differential_velocity_genes(adata, annotation_column=annotation_column,top_gene=top_gene, mode=mode, group_label=group_label)
    # 3. trajectory inference using PAGA
    PAGA_trajectory_inference(adata, annotation_column=annotation_column, mode=mode, group_label=group_label)
    print('scVelo analysis completed.')

