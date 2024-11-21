import os
import numpy as np
import pandas as pd
import scvelo as scv

# ========================== SCVELO FUNCTIONS ==========================
# This function takes in an AnnData object and performs all basic velocity calculations
# enabled by scVelo. It also outputs basic figures such as spliced/unspliced count proportion
# and RNA velocity vectors on umap.
#
def velocity_calculation(adata, group_by, mode='stochastic'):
    """
    :param adata: an AnnData object with cell type annotation
    :param group_by: a key inside adata.obs to use for calculating spliced/unspliced count proportion, normally cell types
    :param mode: can be 'stochastic (default)', 'deterministic', or 'dynamical (slowest)'
    :return: an AnnData object with RNA velocity calculated
    """
    # set parameters for plotting
    kwargs = dict(color=group_by, figsize=(10, 10), dpi=500, show=False)
    # change the 'group_by' column of metadata from dtype: object into dtype: category to comply with proportion plotting
    adata.obs[group_by] = adata.obs[group_by].astype('category')
    # observe proportions of spliced/unspliced counts
    scv.pl.proportions(adata, groupby=group_by, fontsize=8, figsize=(10, 10), dpi=500, show=False, save=f'scvelo/images/scvelo_proportions')
    # velocity calculation workflow
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    if mode == 'dynamical':
        scv.tl.recover_dynamics(adata) # required if running dynamical model
    scv.tl.velocity(adata, mode=mode)
    scv.tl.velocity_graph(adata)
    # save adata object after velocity calculation, since these results can take time to re-run.
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'}) # for getting around a bug
    adata.write(f'scvelo/rds/scvelo_withVelocity_{mode}.h5ad', compression='gzip')
    # basic RNA velocity visualizations in various formats
    scv.pl.velocity_embedding_stream(adata, basis='umap', save=f'scvelo/images/scvelo_embedding_stream', **kwargs)
    scv.pl.velocity_embedding_grid(adata, basis='umap', save=f'scvelo/images/scvelo_embedding_grid', **kwargs)
    scv.pl.velocity_embedding(adata, arrow_length=5, arrow_size=1, basis='umap', save=f'scvelo/images/scvelo_embedding_arrow', **kwargs)
    return


# This function takes in an AnnData object and performs a differential velocity
# t-test to find genes that explain the directionality of calculated velocity vectors. 
# It tests which genes have group-specific differential velocity expression (definied in @group_by)
# i.e being siginificantly higher/lower compared to the remaining population, and visualizes
# the phase portrait (ratio of spliced/unspliced RNA abundance) for highly ranked genes.
#
def differential_velocity_genes(adata, group_by, top_gene=5):
    """
    :param adata: an AnnData object with velocity calculated
    :param group_by: a key inside adata.obs that provides a grouping for cells, ex) cell type
    :param top_gene: an integer specifying number of top ranked genes to plot
    :return: an AnnData object with new data in adata.uns['rank_velocity_genes'] and adata.var['spearmans_score']
    """
    # perform differential velocity t-test
    scv.tl.rank_velocity_genes(adata, groupby=group_by, min_corr=.3)
    # extract top-ranking genes into pandas dataframe
    df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
    df.to_csv(f'scvelo/csv/scvelo_differential_velocity_genes_by_{group_by}.csv')
    # set parameters for plotting
    # plot top 'top_gene' number of genes' phase portrait for each category in 'group_by'
    for item in adata.obs[group_by].unique():
        kwargs = dict(color=group_by, figsize=(2, 2), dpi=500, show=False)
        scv.pl.scatter(adata, df[item][:top_gene], ylabel=item, frameon=False, linewidth=1.5, save=f'scvelo/images/{item}_genePhase', fontsize=8, **kwargs)
        # convert from pandas series to list
        # Note: need to set colorbar=False below to bypass an error caused by matplotlib, in the generated figures darker colors indicate higher expression/velocity
        genes_to_plot = df[item][:top_gene].tolist()
        kwargs = dict(color=group_by, figsize=(10, 10), dpi=500, show=False)
        scv.pl.velocity(adata, genes_to_plot, colorbar=False, ncols=2, save=f'scvelo/images/{item}_genePhaseCompleteInfo', **kwargs)
    return


# This function takes an AnnData object and performs trajectory inference using the
# PAGA method (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1663-x).
# It provides a graph-like map of the data with solid edges corresponding to the transition confidence 
# between two groups (defined in @group_by). Here, PAGA is extended by velocity-inferred directionality
# and predicts transitions/lineages between groups.
#
# Note: there is a possible issue with PAGA with small transition probabilities (< 0.1):
# https://github.com/theislab/scvelo/issues/456
#
def PAGA_trajectory_inference(adata, group_by):
    """
    :param adata: an AnnData object with velocity and velocity graph calculated
    :param group_by: a key inside adata.obs specifying how the cells should be grouped, normally cell types
    :return: an AnnData object with paga graph calculated and stored in adata.uns
    """
    # this is needed due to a current bug in scvelo that hasn't been fixed.
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
    # perform PAGA calculation
    scv.tl.paga(adata, groups=group_by)
    df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
    df.to_csv('scvelo/csv/scvelo_paga_transition_confidence_matrix.csv')
    # generate a directed graph superimposed onto the UMAP embedding
    scv.pl.paga(adata, basis='umap', dashed_edges=None, size=50, alpha=.05, min_edge_width=2, node_size_scale=1.5, figsize=(10, 10), dpi=500, show=False, save=f'scvelo/images/paga_graph')
    return


# Main function of scvelo workflow in python
# 
def run_scvelo_workflow(h5ad_file='scvelo/rds/obj_spliced_unspliced.h5ad', annotation_column='ID', mode='stochastic', top_gene=5):
    # basic scvelo settings
    scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
    scv.set_figure_params('scvelo', transparent=False, format='png')  # set figure format for visualization
    # reading data
    adata = scv.read(h5ad_file)
    # Workflow:
    # 1. calculate RNA velocity using scVelo workflow
    velocity_calculation(adata, group_by=annotation_column, mode=mode)
    # 2. cluster-specific differential velocity genes
    differential_velocity_genes(adata, group_by=annotation_column,top_gene=top_gene)
    # 3. trajectory inference using PAGA
    PAGA_trajectory_inference(adata, group_by=annotation_column)
    print('scVelo analysis completed.')



