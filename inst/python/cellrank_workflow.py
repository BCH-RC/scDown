import os

os.environ["PETSC_OPTIONS"] = "-no_signal_handler"

import argparse
import re
import warnings
from pathlib import Path

import cellrank as cr
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scvelo as scv

warnings.simplefilter("ignore", category=DeprecationWarning)

def save_fig(name, output_format='png'):
    """Save a figure

    Args:
        name (str): prefix of the figure to save
        output_format (str, optional): Output figure format. Defaults to 'png'.
    """
    # saves the currently active Matplotlib plot to disk and clears memory.
    plt.savefig(
        name,
        format=output_format,
        dpi=300,
        bbox_inches='tight',
        facecolor='white'
    )
    plt.close('all')


def set_robust_states(estimator, expected_states, annotation_column, state_type='terminal'):
    """Set initial or terminal macrostates using cell type labels

    Args:
        estimator (GPCCA): estimator used by CellRank
        expected_states (list[str]): expected cell types to be either terminal or initial states.
        annotation_column (str): cell type column in adata.obs
        state_type (str, optional): type of state (either initial or terminal). Defaults to 'terminal'.
    """
    adata = estimator.adata

    # extract the names of the automatically discovered macrostates
    computed_macrostates = estimator.macrostates.cat.categories.tolist()

    state_dict = {}
    for expected in expected_states:
        # regex to match the exact name or the name followed by _1, _2, etc.
        pattern = re.compile(rf"^{re.escape(expected)}(_\d+)?$")
        matches = [m for m in computed_macrostates if pattern.match(m)]

        if matches:
            # if matched, grab the cells from the computed macrostates
            for match in matches:
                # Filter for cells assigned to this specific macrostate
                cells = estimator.macrostates[estimator.macrostates == match].index.tolist()
                state_dict[match] = cells
                print(f"Matched '{expected}' to macrostate: {match}")
        else:
            # state missing from macrostates
            cells = adata.obs.index[adata.obs[annotation_column] == expected].tolist()
            if not cells:
                print(f"Error: '{expected}' not found in macrostates or {annotation_column}.")
                continue
            print(f"Using cells from '{expected}': Missing from macrostates.")
            state_dict[expected] = cells

    # use the custom dictionary to set the states
    if state_type == 'terminal':
        estimator.set_terminal_states(state_dict, allow_overlap=True)
    elif state_type == 'initial':
        estimator.set_initial_states(state_dict, allow_overlap=True)


def get_lineage_clusters(estimator, fates_df, lineage, annotation_column):
    """Find relevant cell type clusters associated with a lineage

    Args:
        estimator (GPCCA): Estimator used by CellRank
        fates_df (DataFrame): precomputed dataframe with cell-wise fate probabilities
        lineage (str): lineage (terminal state) to focus on
        annotation_column (str): cell type column in adata.obs

    Returns:
        list[str]: list of clusters along the specified lineage
    """
    # aggregate split clusters (e.g., Astrocyte_1 + Astrocyte_2 -> Astrocyte)
    rename_map = {}
    for col in fates_df.columns:
        base_name = re.sub(r'_\d+$', '', col)
        rename_map[col] = base_name

    aggregated_fates_df = fates_df.groupby(rename_map, axis=1).sum()

    # check if lineage is in the fate dataframe
    lineage = re.sub(r'_\d+$', '', lineage)
    if lineage not in aggregated_fates_df.columns:
        raise ValueError(f"Error: '{lineage}' not found.")

    # add the cluster annotations
    aggregated_fates_df['cluster'] = estimator.adata.obs[annotation_column].values

    # calculate the mean fate probability for each cluster towards every lineage
    mean_fates_per_cluster = aggregated_fates_df.groupby('cluster').mean()

    # find the lineage with the highest mean probability for each cluster
    vote_lineage = mean_fates_per_cluster.idxmax(axis=1)

    # filter for target lineage
    relevant_clusters = vote_lineage[vote_lineage == lineage].index.tolist()

    print(f"Assigned these clusters to '{lineage}': {relevant_clusters}")
    return relevant_clusters


def lineage_analysis(
        adata,
        estimator,
        fates_df,
        lineage,
        model,
        model_small,
        annotation_column,
        out_path,
        cluster_genes,
        clustering_kwargs,
        neighbors_kwargs,
        output_format='png',
        n_jobs=1
    ):
    """Perform lineage-specific analysis and plotting

    Args:
        adata (AnnData): input adata object
        estimator (GPCCA): estimator used by CellRank
        fates_df (DataFrame): precomputed dataframe with cell-wise fate probabilities
        lineage (str): lineage (terminal state) to focus on
        model (GAM): large GAM model
        model_small (GAM): large GAM model
        annotation_column (str): cell type column in adata.obs
        out_path (str): output_path
        cluster_genes (bool): whether to cluster genes by gene trends
        clustering_kwargs (dict): keyword arguments for clustering genes
        neighbors_kwargs (dict): keyword arguments for computing gene neighbors for UMAP
        output_format (str, optional): output figure format. Defaults to 'png'.
        n_jobs (int, optional): number of parallel processes to use. Defaults to 1.
    """
    fig_path = out_path+'/plots'

    # find clusters involved in a lineage
    clusters_relevant = get_lineage_clusters(estimator, fates_df, lineage, annotation_column)

    # uncover driver genes (highly correlated with expression) for a lineage
    driver_df = estimator.compute_lineage_drivers(
        lineages=[lineage],
        cluster_key=annotation_column,
        clusters=clusters_relevant
    )
    driver_df.to_csv(out_path + f'/tables/driver_genes_{lineage}.txt', sep='\t', index_label="Gene")

    adata.obs[f"fate_probabilities_{lineage}"] = estimator.fate_probabilities[lineage].X.flatten()

    sc.pl.embedding(
        adata,
        basis="umap",
        color=[f"fate_probabilities_{lineage}"] + list(driver_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        show=False
    )
    save_fig(f"{fig_path}/driver_genes_fate_umap_{lineage}.{output_format}", output_format)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # visualize gene trends
        cr.pl.gene_trends(
            adata,
            model=model,
            genes=driver_df.index[:8],
            same_plot=True,
            ncols=2,
            time_key="latent_time",
            hide_cells=True,
            data_key="X"
        )
        save_fig(f"{fig_path}/driver_genes_gene_trend_{lineage}.{output_format}", output_format)

        cr.pl.heatmap(
            adata,
            model=model,
            lineages=lineage,
            cluster_key=annotation_column,
            show_fate_probabilities=True,
            data_key="X",
            genes=driver_df.index[:40],
            time_key="latent_time",
            figsize=(12, 10),
            show_all_genes=True
        )
        save_fig(f"{fig_path}/driver_genes_heatmap_{lineage}.{output_format}", output_format)

        # cluster gene expression trends. Can be slow.
        if cluster_genes:
            cr.pl.cluster_trends(
                adata,
                model=model_small,
                lineage=lineage,
                data_key="X",
                genes=adata[:, adata.var["highly_variable"]].var_names,
                time_key="latent_time",
                n_jobs=n_jobs,
                backend="threading",
                random_state=0,
                recompute=True,
                clustering_kwargs=clustering_kwargs,
                neighbors_kwargs=neighbors_kwargs
            )
            save_fig(f"{fig_path}/gene_cluster_trends_{lineage}.{output_format}", output_format)

            gdata = adata.uns[f"lineage_{lineage}_trend"].copy()
            sc.tl.umap(gdata, random_state=0, init_pos='random')
            sc.pl.embedding(gdata, basis="umap", color="clusters", show=False)
            save_fig(f"{fig_path}/gene_clusters_{lineage}.{output_format}", output_format)

            gdata.write_h5ad(f'{out_path}/data/gdata_{lineage}.h5ad')


def run_cellrank_workflow(
        h5ad_file: str | None = None,
        out_path: str = ".",
        annotation_column: str = "clusters",
        n_states: int | None = None,
        terminal_cell_types: str | list[str] | None = None,
        initial_cell_types: str | list[str] | None = None,
        lineages_of_interest: str | list[str] | None = None,
        velocity_kernel_weight: float = 0.9,
        cluster_genes: bool = False,
        cluster_resolution: float | None = None,
        n_neighbors: int | None = None,
        output_format: str = "png",
        n_jobs: int = 2,
    ) -> None:
    """Main function of CellRank workflow

    Args:
        h5ad_file (str, optional): input h5ad file. Defaults to None.
        out_path (str, optional): output path. Defaults to '.'.
        annotation_column (str, optional): cell type column in adata.obs. Defaults to 'clusters'.
        n_states (int, optional): number of macrostates to compute. Defaults to None.
        terminal_cell_types (str | list[str], optional): expected list of terminal cell types. Defaults to None.
        initial_cell_types (str | list[str], optional): expected list of initial cell types. Defaults to None.
        lineages_of_interest (str | list[str], optional): list of lineages to perform lineage-specific analysis. Defaults to None.
        velocity_kernel_weight (float, optional): weight of velocity_kernel. connectivity kernel weight is 1 - this value. Defaults to 0.9.
        cluster_genes (bool, optional): whether to cluster genes by gene trends. Defaults to False.
        cluster_resolution (float, optional): resolution for clustering genes with leiden. Defaults to None.
        n_neighbors (int, optional): number of neighbors for computing gene neighbors. Defaults to None.
        output_format (str, optional): output figure format. Defaults to "png".
        n_jobs (int, optional): number of parallel processes to use. Defaults to 2.
    """
    out_path = out_path+'/cellrank'
    for folder in [out_path,
                   out_path+'/data',
                   out_path+'/tables',
                   out_path+'/plots']:
        Path(folder).mkdir(parents=True, exist_ok=True)
    fig_path = out_path+'/plots'

    scv.settings.verbosity = 3
    cr.settings.logging_level = "INFO"
    scv.set_figure_params('scvelo', dpi=100, dpi_save=300, format=output_format, facecolor='white', transparent=False, frameon=False)

    if h5ad_file is None:
        h5ad_file = out_path+'/../scvelo/data/ALL_with_velocity_dynamical.h5ad'
    adata = sc.read_h5ad(h5ad_file)
    if 'latent_time' not in adata.obs.columns:
        scv.tl.latent_time(adata)
    print(adata)

    # set up kernel
    vk = cr.kernels.VelocityKernel(adata)
    vk.compute_transition_matrix()
    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()
    combined_kernel = velocity_kernel_weight * vk + (1-velocity_kernel_weight) * ck
    print(combined_kernel)

    # visualize cell-type transitions
    combined_kernel.plot_projection(show=False)
    save_fig(f"{fig_path}/projection.{output_format}", output_format)

    # initialize an estimator
    g = cr.estimators.GPCCA(combined_kernel)
    print(g)
    g.compute_schur()
    _, _ = plt.subplots(figsize=(6, 4))
    g.plot_spectrum(real_only=True)
    save_fig(f"{fig_path}/eigenvalues.{output_format}", output_format)

    # compute macrostates depending on n_states
    if n_states is None:
        n_states = adata.obs[annotation_column].nunique()
    g.compute_macrostates(n_states=n_states, cluster_key="clusters")
    g.plot_macrostates(which="all", legend_loc="right", s=100, show=False)
    save_fig(f"{fig_path}/macrostates_umap.{output_format}", output_format)

    g.plot_macrostate_composition(key=annotation_column, figsize=(7, 4), show=False)
    save_fig(f"{fig_path}/state_compositions.{output_format}", output_format)

    plt.figure(figsize=(6, 6))
    g.plot_coarse_T(annotate=False)
    save_fig(f"{fig_path}/transition_matrix.{output_format}", output_format)

    # classify macrostates
    if terminal_cell_types is None:
        g.predict_terminal_states()
        print(f"Found terminal states: {g.terminal_states.cat.categories.tolist()}")
    else:
        if isinstance(terminal_cell_types, str):
            terminal_cell_types = [terminal_cell_types]
        set_robust_states(g, terminal_cell_types, annotation_column, 'terminal')
    g.plot_macrostates(which="terminal", legend_loc="right", s=100, show=False)
    save_fig(f"{fig_path}/terminal_states.{output_format}", output_format)

    if initial_cell_types is None:
        g.predict_initial_states()
        print(f"Found initial states: {g.initial_states.cat.categories.tolist()}")
    else:
        if isinstance(initial_cell_types, str):
            initial_cell_types = [initial_cell_types]
        set_robust_states(g, initial_cell_types, annotation_column, 'initial')
    g.plot_macrostates(which="initial", s=100, show=False)
    save_fig(f"{fig_path}/initial_states.{output_format}", output_format)

    # estimate fate probabilities
    n_cells = adata.n_obs
    print(f"Dataset contains {n_cells} cells.")
    if n_cells <= 20000:
        print("Using 'direct' solver for fate probabilities...")
        g.compute_fate_probabilities(solver="direct", use_petsc=False)
    else:
        print("Using 'iterative' solver with ILU preconditioner...")
        g.compute_fate_probabilities(solver="gmres", n_jobs=n_jobs, use_petsc=False)

    g.plot_fate_probabilities(same_plot=False, show=False)
    save_fig(f"{fig_path}/fate_probabilities_split.{output_format}", output_format)
    g.plot_fate_probabilities(same_plot=True, legend_loc='right margin', show=False)
    save_fig(f"{fig_path}/fate_probabilities_shared.{output_format}", output_format)

    cr.pl.circular_projection(
        adata,
        keys=[annotation_column],
        legend_loc="right"
    )
    save_fig(f"{fig_path}/fate_probabilities_circular.{output_format}", output_format)

    sc.tl.paga(adata, groups=annotation_column)
    cr.pl.aggregate_fate_probabilities(
        adata,
        mode="paga_pie",
        cluster_key=annotation_column,
        legend_kwargs={"title_fontsize": 10, "fontsize": 8, "loc": "upper right out"},
        fontsize=6,
        edge_width_scale=0.3,
        show=False
    )
    save_fig(f"{fig_path}/fate_probabilities_paga_pie.{output_format}", output_format)

    # save fate probabilities to table
    fates_matrix = g.fate_probabilities.X
    if hasattr(fates_matrix, "todense"):
        fates_matrix = fates_matrix.todense()
    fates_df = pd.DataFrame(
        fates_matrix,
        index=g.adata.obs_names,
        columns=g.fate_probabilities.names
    )
    fates_df.to_csv(out_path + '/tables/fate_probabilities.txt', sep="\t", index_label="Cell_ID")

    # precompute gene trends
    model = cr.models.GAM(adata)
    model_small = None
    if cluster_genes:
        model_small = cr.models.GAM(adata, max_iter=100)  # instead of default 2000
    computed_macrostates = g.macrostates.cat.categories.tolist()
    if lineages_of_interest is not None:
        if isinstance(lineages_of_interest, str):
            lineages_of_interest = [lineages_of_interest]

        lineages_matched = []
        for lineage in lineages_of_interest:
            # regex to match the exact name or the name followed by _1, _2, etc.
            pattern = re.compile(rf"^{re.escape(lineage)}(_\d+)?$")
            matches = [m for m in computed_macrostates if pattern.match(m)]
            if not matches:
                print(f"Warning: No macrostates found matching '{lineage}'.")
            else:
                lineages_matched.extend(matches)

        # default clustering and neighborhood computation arguments
        clustering_kwargs = {"resolution": 0.2 if cluster_resolution is None else cluster_resolution, "random_state": 0}
        neighbors_kwargs = {"n_neighbors": 15 if n_neighbors is None else n_neighbors, "random_state": 0}

        for lineage in lineages_matched:
            lineage_analysis(adata,
                             g,
                             fates_df,
                             lineage,
                             model,
                             model_small,
                             annotation_column,
                             out_path,
                             cluster_genes,
                             clustering_kwargs,
                             neighbors_kwargs,
                             output_format,
                             1)

    print('CellRank2 analysis completed.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CellRank workflow.")

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
        "--n_states",
        type=int,
        default=None,
        help="Number of macrostates to compute (default: None)."
    )
    parser.add_argument(
        "--terminal_cell_types",
        type=str,
        nargs='+',
        default=None,
        help="Expected terminal cell types. Can pass multiple space-separated values."
    )
    parser.add_argument(
        "--initial_cell_types",
        type=str,
        nargs='+',
        default=None,
        help="Expected initial cell types. Can pass multiple space-separated values."
    )
    parser.add_argument(
        "--lineages_of_interest",
        type=str,
        nargs='+',
        default=None,
        help="List of lineages for lineage-specific analysis. Can pass multiple space-separated values."
    )
    parser.add_argument(
        "--velocity_kernel_weight",
        type=float,
        default=0.9,
        help="Weight of velocity_kernel. Connectivity kernel weight is 1 - this value (default: 0.9)."
    )
    parser.add_argument(
        "--cluster_genes",
        action="store_true",
        help="Include this flag to cluster genes by gene trends (sets to True)."
    )
    parser.add_argument(
        "--cluster_resolution",
        type=float,
        default=None,
        help="Resolution for clustering genes with leiden (default: None)."
    )
    parser.add_argument(
        "--n_neighbors",
        type=int,
        default=None,
        help="Number of neighbors for computing gene neighbors (default: None)."
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

    # Pass all parsed arguments to the function as keyword arguments
    run_cellrank_workflow(**vars(args))
