import argparse
from pathlib import Path
from typing import Literal

import ktplotspy as kpy
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from plotnine import facet_wrap
from plotnine.options import set_option

set_option('limitsize', False)

def process_h5ad(h5ad_file, out_path):
    """Preprocess h5ad data

    Args:
        h5ad_file (str): input h5ad path
        out_path (str): processed h5ad path

    Returns:
        AnnData: processed adata object
    """
    adata = sc.read_h5ad(h5ad_file)
    print(adata)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.write_h5ad(out_path + '/data/lognorm.h5ad')
    return adata


def convert_mouse_genes(adata, out_path):
    """Convert mouse genes to human orthologs

    Args:
        adata (AnnData): input adata object
        out_path (str): h5ad to save the converted adata

    Returns:
        AnnData: converted adata object
    """
    from mousipy import translate
    humanized_adata = translate(adata)
    print(humanized_adata)
    humanized_adata.write_h5ad(out_path + '/data/lognorm_humanized.h5ad')
    return humanized_adata


def save_fig(p, name, output_format='png', plot_type='dotplot'):
    """Save a figure

    Args:
        p (Figure): figure object from plotnine, matplotlib, or seaborn
        name (str): prefix of the figure to save
        output_format (str, optional): output figure format. Defaults to 'png'.
        plot_type (str, optional): type of the plot. Defaults to 'dotplot'.
    """
    if plot_type == 'dotplot':
        if output_format == 'png':
            p.save(name + '.png', dpi=300)
        elif output_format == 'pdf':
            p.save(name + '.pdf', )
        elif output_format == 'jpg' or output_format == 'jpeg':
            p.save(name + '.jpg', dpi=300)
    elif plot_type == 'heatmap':
        if output_format == 'png':
            p.savefig(name + '.png', bbox_inches="tight", dpi=300)
        elif output_format == 'pdf':
            p.savefig(name + '.pdf', bbox_inches="tight")
        elif output_format == 'jpg' or output_format == 'jpeg':
            p.savefig(name + '.jpg', bbox_inches="tight", dpi=300)
    elif plot_type == 'chord':
        if output_format == 'png':
            p.savefig(name + '.png', bbox_inches="tight", pad_inches=2.5, dpi=300)
        elif output_format == 'pdf':
            p.savefig(name + '.pdf', bbox_inches="tight", pad_inches=2.5)
        elif output_format == 'jpg' or output_format == 'jpeg':
            p.savefig(name + '.jpg', bbox_inches="tight", pad_inches=2.5, dpi=300)


def cellphonedb_method1(cpdb_file_path, meta_file_path, counts_file_path, out_path, threads=2):
    """CellPhoneDB Method 1: simple analysis

    Args:
        adata (AnnData): input adata object
        cpdb_file_path (str): cpdb database file location
        meta_file_path (str): metadata file location that contains cell to cell type mapping
        counts_file_path (str): in memory adata object
        out_path (str): output path
        threads (int, optional): number of parallel processes to use. Defaults to 5.

    Returns:
        dict: cpdb result tables
    """
    from cellphonedb.src.core.methods import cpdb_analysis_method

    cpdb_results = cpdb_analysis_method.call(
        cpdb_file_path = cpdb_file_path,           # mandatory: CellphoneDB database zip file.
        meta_file_path = meta_file_path,           # mandatory: tsv file defining barcodes to cell label.
        counts_file_path = counts_file_path,       # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
        counts_data = 'hgnc_symbol',               # defines the gene annotation in counts matrix.
        microenvs_file_path = None,                # optional (default: None): defines cells per microenvironment.
        score_interactions = True,                 # optional: whether to score interactions or not.
        threshold = 0.1,                           # defines the min % of cells expressing a gene for this to be employed in the analysis.
        threads = threads,                         # number of threads to use in the analysis.
        result_precision = 3,                      # Sets the rounding for the mean values in significan_means.
        separator = '|',                           # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
        debug = False,                             # Saves all intermediate tables emplyed during the analysis in pkl format.
        output_path = out_path+'/tables',          # Path to save results.
        output_suffix = ''                         # Replaces the timestamp in the output files by a user defined string (default: None).
    )

    print(cpdb_results.keys())
    return cpdb_results


def cellphonedb_method2(
        adata,
        cpdb_file_path,
        meta_file_path,
        counts_file_path,
        out_path,
        cell_type_col,
        n_top_interactions=10,
        n_top_genes=3,
        threads=2,
        significant_threshold=0.05,
        output_format='png'
    ):
    """CellPhoneDB Method 2: statistical analysis

    Args:
        adata (AnnData): input adata object
        cpdb_file_path (str): cpdb database file location
        meta_file_path (str): metadata file location that contains cell to cell type mapping
        counts_file_path (str): in memory adata object
        out_path (str): output path
        cell_type_col (str): cell type column in adata.obs
        n_top_interactions (int, optional): number of top interactions to show on plots. Defaults to 10.
        n_top_genes (int, optional): number of top genes to show on plots. Defaults to 3.
        threads (int, optional): number of parallel processes to use. Defaults to 5.
        significant_threshold (float, optional): threshold for significance. Defaults to 0.05.
        output_format (str, optional): output figure format. Defaults to 'png'.

    Returns:
        dict: cpdb result tables
    """
    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

    cpdb_results = cpdb_statistical_analysis_method.call(
        cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
        meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
        counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
        counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
        active_tfs_file_path = None,                     # optional: defines cell types and their active TFs.
        microenvs_file_path = None,                      # optional (default: None): defines cells per microenvironment.
        score_interactions = True,                       # optional: whether to score interactions or not.
        iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
        threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
        threads = threads,                               # number of threads to use in the analysis.
        debug_seed = 42,                                 # debug randome seed. To disable >=0.
        result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
        pvalue = significant_threshold,                  # P-value threshold to employ for significance.
        subsampling = False,                             # To enable subsampling the data (geometri sketching).
        subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
        subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
        subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
        separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
        debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
        output_path = out_path+'/tables',                # Path to save results.
        output_suffix = ''                               # Replaces the timestamp in the output files by a user defined string (default: None).
        )

    # Extract the top pairs of interactions based on counts of significance
    sig_means = cpdb_results['significant_means']
    cell_pair_cols = sig_means.columns[14:]
    sig_means['significant_count'] = sig_means[cell_pair_cols].notna().sum(axis=1)

    unique_interactions = sig_means[['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b', 'gene_a', 'gene_b', 'significant_count']].drop_duplicates()
    top_unique_interactions = unique_interactions.sort_values(by='significant_count', ascending=False).head(n_top_interactions)
    top_pairs = [str(x) for x in top_unique_interactions['interacting_pair'].values]

    # Filter unique gene list
    top_interaction_genes = pd.concat([top_unique_interactions['gene_a'], top_unique_interactions['gene_b']])
    unique_genes = top_interaction_genes.value_counts().index.to_numpy()

    p = kpy.plot_cpdb_heatmap(
        pvals = cpdb_results['pvalues'],
        degs_analysis = False,
        alpha = significant_threshold,
        figsize = (6, 6),
        title = "Sum of significant interactions (sym)"
    )
    save_fig(p, out_path+'/plots/interaction_heatmap_symmetrical', output_format=output_format, plot_type='heatmap')
    plt.close(p.fig)

    table = kpy.plot_cpdb_heatmap(
        pvals = cpdb_results['pvalues'],
        degs_analysis = False,
        alpha = significant_threshold,
        return_tables = True,
        symmetrical = False # X: target, Y: source
    )
    table["count_network"].to_csv(out_path + '/tables/count_network_directional.txt', sep='\t')
    table["interaction_edges"].to_csv(out_path + '/tables/interaction_edges_directional.txt', sep='\t', index=False)
    p = kpy.plot_cpdb_heatmap(
        pvals = cpdb_results['pvalues'],
        degs_analysis = False,
        alpha = significant_threshold,
        figsize = (6, 6),
        title = "Sum of significant interactions (dir)",
        symmetrical = False # X: target, Y: source
    )
    p.ax_heatmap.set_ylabel("Senders", fontsize=13, labelpad=12)
    p.ax_heatmap.set_xlabel("Receivers", fontsize=13, labelpad=12)
    save_fig(p, out_path+'/plots/interaction_heatmap_directional', output_format=output_format, plot_type='heatmap')
    plt.close(p.fig)

    p = kpy.plot_cpdb(
        adata = adata,
        cell_type1 = ".",
        cell_type2 = ".",
        means = cpdb_results['means'],
        pvals = cpdb_results['pvalues'],
        celltype_key = cell_type_col,
        interacting_pairs = top_pairs,
        figsize = (9+len(cell_pair_cols)/10, 5+len(top_pairs)/2),
        title = "Interactions between all cell types",
        max_size = 5,
        highlight_size = 0.75,
        degs_analysis = False,
        standard_scale = True,
        interaction_scores = cpdb_results['interaction_scores'],
        scale_alpha_by_interaction_scores = True
    )
    save_fig(p, out_path+'/plots/interaction_dotplot_top_pairs', output_format=output_format, plot_type='dotplot')

    p = kpy.plot_cpdb(
        adata = adata,
        cell_type1 = ".",
        cell_type2 = ".",
        means = cpdb_results['means'],
        pvals = cpdb_results['pvalues'],
        celltype_key = cell_type_col,
        interacting_pairs = top_pairs,
        figsize = (9+len(cell_pair_cols)/10, 10+len(top_pairs)),
        title = "Interactions between all cell types\n grouped by pathway",
        max_size = 5,
        highlight_size = 0.75,
        degs_analysis = False,
        standard_scale = True,
        interaction_scores = cpdb_results['interaction_scores'],
        scale_alpha_by_interaction_scores = True
    )
    p = p + facet_wrap("~ classification", ncol = 1)
    save_fig(p, out_path+'/plots/interaction_dotplot_top_pairs_by_pathway', output_format=output_format, plot_type='dotplot')

    for i in range(n_top_genes):
        gi = unique_genes[i]
        interaction_count = cpdb_results['means']['interacting_pair'].str.contains(gi).sum()
        p = kpy.plot_cpdb(
            adata = adata,
            cell_type1 = ".",
            cell_type2 = ".",
            means = cpdb_results['means'],
            pvals = cpdb_results['pvalues'],
            celltype_key = cell_type_col,
            genes = [gi],
            figsize = (9+len(cell_pair_cols)/10, 3+interaction_count/2),
            title = "Interactions between all cell types",
            max_size = 5,
            highlight_size = 0.75,
            degs_analysis = False,
            standard_scale = True,
            interaction_scores = cpdb_results['interaction_scores'],
            scale_alpha_by_interaction_scores = True
        )
        save_fig(p, out_path+f'/plots/interaction_dotplot_{gi}', output_format=output_format, plot_type='dotplot')

    # Filter the DataFrames to contain top pairs
    means_subset = cpdb_results['means'][cpdb_results['means']['interacting_pair'].isin(top_pairs)].copy()
    pvals_subset = cpdb_results['pvalues'][cpdb_results['pvalues']['interacting_pair'].isin(top_pairs)].copy()
    means_subset = means_subset.sort_values('interacting_pair')
    pvals_subset = pvals_subset.sort_values('interacting_pair')
    p = kpy.plot_cpdb_chord(
        adata=adata,
        cell_type1=".",
        cell_type2=".",
        means=means_subset,
        pvals=pvals_subset,
        deconvoluted=cpdb_results['deconvoluted'],
        celltype_key=cell_type_col,
        interaction=None,
        link_kwargs={"direction": 1, "allow_twist": True, "r1": 90, "r2": 90},
        sector_text_kwargs={"color": "black", "size": 12, "r": 105, "adjust_rotation": True},
        legend_kwargs={"loc": "upper left", "bbox_to_anchor": (1.05, 1), "fontsize": 8},
        link_offset=1
    )
    p = p.ax.get_figure()
    save_fig(p, out_path+'/plots/interaction_chord_diagram', output_format=output_format, plot_type='chord')
    plt.close(p)

    print(cpdb_results.keys())
    return cpdb_results


def calculate_deg(
        adata,
        cell_type_col,
        condition_col,
        disease_label,
        control_label,
        out_path,
        pval_cutoff=0.05,
        logfc_cutoff=0.5
    ):
    """Calcualte differentially expressed genes (DEGs) in each cell type comparing conditions

    Args:
        adata (AnnData): input adata object
        cell_type_col (str): cell type column in adata.obs
        condition_col (str): condition (group) column in adata.obs
        disease_label (str or list[str]): the treatment or disease group labels in conditions
        control_label (str): the control or background group label in conditions
        out_path (str): output path
        pval_cutoff (float, optional): P adjusted value cutoff. Defaults to 0.05.
        logfc_cutoff (float, optional): log2 fold change cutoff. Defaults to 0.5.

    Returns:
        str: path to the output DEG table
    """
    deg_list = []
    unique_cell_types = adata.obs[cell_type_col].unique()

    # Make disease label a list if it's not
    if isinstance(disease_label, str):
        disease_label = [disease_label]

    # Loop through each cell type
    for cell_type in unique_cell_types:
        # Subset the adata to this cell type
        adata_sub = adata[adata.obs[cell_type_col] == cell_type].copy()

        # Ensure this cell type exists in both conditions
        if (adata_sub.obs[condition_col].isin(disease_label).sum()<2) or ((adata_sub.obs[condition_col].values==control_label).sum()<2):
            print(f"Skipping {cell_type}: missing or having too few cells in one of the conditions.")
            continue

        print(f"Calculating DEGs for {cell_type}")

        # Compare disease vs control for this cell type
        sc.tl.rank_genes_groups(
            adata_sub,
            groupby=condition_col,
            groups=disease_label,
            reference=control_label,
            method='wilcoxon'
        )
        result = adata_sub.uns['rank_genes_groups']

        # Build a df for this cell type from all disease labels
        df_list = []
        for label in disease_label:
            df_dis = pd.DataFrame({
                'gene': result['names'][label],
                'pvals_adj': result['pvals_adj'][label],
                'logfoldchanges': result['logfoldchanges'][label]
            })
            df_list.append(df_dis)

        df = pd.concat(df_list, ignore_index=True)

        # Filter for significance and up-regulated in any disease state
        filtered_df = df[(df['pvals_adj'] < pval_cutoff) & (df['logfoldchanges'] > logfc_cutoff)].copy()
        filtered_df['cluster'] = cell_type
        filtered_df = filtered_df.drop_duplicates()
        deg_list.append(filtered_df[['cluster', 'gene']])

    deg_df = pd.concat(deg_list, ignore_index=True)
    deg_df.to_csv(out_path+'/tables/DEG_significant.txt', sep='\t', index=False)
    return out_path+'/tables/DEG_significant.txt'


def cellphonedb_method3(
        adata,
        cpdb_file_path,
        meta_file_path,
        counts_file_path,
        degs_file_path,
        out_path,
        cell_type_col,
        n_top_interactions=10,
        n_top_genes=3,
        threads=2,
        output_format='png'
    ):
    """CellPhoneDB Method 3: DEGs analysis

    Args:
        adata (AnnData): input adata object
        cpdb_file_path (str): cpdb database file location
        meta_file_path (str): metadata file location that contains cell to cell type mapping
        counts_file_path (str): in memory adata object
        degs_file_path (str): file that contains the significant DEG results
        out_path (str): ouptut path
        cell_type_col (str): cell type column in adata.obs
        n_top_interactions (int, optional): number of top interactions to show on plots. Defaults to 10.
        threads (int, optional): number of parallel processes to use. Defaults to 5.
        output_format (str, optional): output figure format. Defaults to 'png'.

    Returns:
        dict: cpdb result tables
    """
    import seaborn as sns
    from cellphonedb.src.core.methods import cpdb_degs_analysis_method
    from scipy.stats import zscore

    cpdb_results = cpdb_degs_analysis_method.call(
        cpdb_file_path = cpdb_file_path,                            # mandatory: CellphoneDB database zip file.
        meta_file_path = meta_file_path,                            # mandatory: tsv file defining barcodes to cell label.
        counts_file_path = counts_file_path,                        # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
        degs_file_path = degs_file_path,                            # mandatory: tsv file with DEG to account.
        counts_data = 'hgnc_symbol',                                # defines the gene annotation in counts matrix.
        active_tfs_file_path = None,                                # optional: defines cell types and their active TFs.
        microenvs_file_path = None,                                 # optional (default: None): defines cells per microenvironment.
        score_interactions = True,                                  # optional: whether to score interactions or not.
        threshold = 0.1,                                            # defines the min % of cells expressing a gene for this to be employed in the analysis.
        threads = threads,                                          # number of threads to use in the analysis.
        result_precision = 3,                                       # Sets the rounding for the mean values in significan_means.
        separator = '|',                                            # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
        debug = False,                                              # Saves all intermediate tables emplyed during the analysis in pkl format.
        output_path = out_path+'/tables',                           # Path to save results.
        output_suffix = ''                                          # Replaces the timestamp in the output files by a user defined string (default: None).
        )

    sig_means = cpdb_results['significant_means']
    cell_pair_cols = sig_means.columns[14:]

    # Get annotation columns and interactions
    annotation = list(cpdb_results['relevant_interactions'].columns[:13])
    interaction = list(cpdb_results['relevant_interactions'].columns[13:])

    # Convert relevant_interactions from wide to long
    relevant_interactions_long = pd.melt(
        cpdb_results['relevant_interactions'],
        id_vars = annotation,
        var_name = 'Interacting_cell',
        value_vars = interaction,
        value_name = 'Relevance'
    )

    relevant_interactions_long[['Cell_a', 'Cell_b']] = relevant_interactions_long['Interacting_cell'].str.split('|', expand = True).rename(columns={0: 'Cell_a', 1: 'Cell_b'})

    # Convert means file from wide to long
    means_long = pd.melt(
        cpdb_results['means'],
        id_vars = annotation,
        var_name = 'Interacting_cell',
        value_vars = interaction,
        value_name = 'Mean'
    )

    # Create a dictionary with recurrence of relevance of any given interaction
    id_cp_dict = relevant_interactions_long.groupby('id_cp_interaction')['Relevance'].sum().to_dict()

    # Add new column to indicate the recurrence of an interaction
    relevant_interactions_long['Recurrence'] = relevant_interactions_long['id_cp_interaction'].map(id_cp_dict)

    # Sort according to the recurrence
    relevant_interactions_long = relevant_interactions_long.sort_values(['Recurrence'], ascending = True)

    # Find the top interactions
    unique_interactions = relevant_interactions_long[['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b', 'gene_a', 'gene_b', 'Recurrence']].drop_duplicates()
    top_unique_interactions = unique_interactions.sort_values(by='Recurrence', ascending=False).head(n_top_interactions)
    top_interactions_id = top_unique_interactions["id_cp_interaction"].values
    top_pairs = [str(x) for x in top_unique_interactions["interacting_pair"].values]

    # Filter unique gene list
    top_interaction_genes = pd.concat([top_unique_interactions['gene_a'], top_unique_interactions['gene_b']])
    unique_genes = top_interaction_genes.value_counts().index.to_numpy()

    # Filter to top interactions
    idx = [id in top_interactions_id for id in relevant_interactions_long['id_cp_interaction']]
    relevant_interactions_plot = relevant_interactions_long[idx].copy()

    # Add mean value of the interacting partners
    relevant_interactions_plot = relevant_interactions_plot.merge(
        means_long[['id_cp_interaction', 'Interacting_cell', 'Mean']],
        on = ['id_cp_interaction', 'Interacting_cell'],
        how = 'inner'
    )

    # Scale within interaction
    relevant_interactions_plot['Mean_scaled'] = relevant_interactions_plot.groupby('id_cp_interaction', group_keys = False)['Mean'].transform(lambda x : zscore(x, ddof = 1))
    relevant_interactions_plot = relevant_interactions_plot.sort_values('Interacting_cell')

    p = sns.relplot(
        data = relevant_interactions_plot,
        x = "Interacting_cell",
        y = "interacting_pair",
        hue = "Relevance",
        size = "Mean_scaled",
        palette = "vlag",
        hue_norm = (-1, 1),
        height = 4,
        aspect = 3,
        sizes = (0, 200)
    )
    p.set_xticklabels(rotation = 90)
    save_fig(p, out_path+'/plots/interaction_bubble', output_format=output_format, plot_type='heatmap')
    plt.close(p.figure)

    p = kpy.plot_cpdb_heatmap(
        pvals = cpdb_results['relevant_interactions'],
        degs_analysis = True,
        figsize = (6, 6),
        title = "Sum of significant interactions"
    )
    save_fig(p, out_path+'/plots/interaction_heatmap_symmetrical', output_format=output_format, plot_type='heatmap')
    plt.close(p.fig)

    table = kpy.plot_cpdb_heatmap(
        pvals = cpdb_results['relevant_interactions'],
        degs_analysis = True,
        return_tables = True,
        symmetrical = False # X: target, Y: source
    )
    table["count_network"].to_csv(out_path + '/tables/count_network_directional.txt', sep='\t')
    table["interaction_edges"].to_csv(out_path + '/tables/interaction_edges_directional.txt', sep='\t', index=False)
    p = kpy.plot_cpdb_heatmap(
        pvals = cpdb_results['relevant_interactions'],
        degs_analysis = True,
        figsize = (6, 6),
        title = "Sum of significant interactions (dir)",
        symmetrical = False # X: target, Y: source
    )
    p.ax_heatmap.set_ylabel("Senders", fontsize=13, labelpad=12)
    p.ax_heatmap.set_xlabel("Receivers", fontsize=13, labelpad=12)
    save_fig(p, out_path+'/plots/interaction_heatmap_directional', output_format=output_format, plot_type='heatmap')
    plt.close(p.fig)

    p = kpy.plot_cpdb(
        adata = adata,
        cell_type1 = ".",
        cell_type2 = ".",
        means = cpdb_results['means'],
        pvals = cpdb_results['relevant_interactions'],
        celltype_key = cell_type_col,
        interacting_pairs = top_pairs,
        figsize = (9+len(cell_pair_cols)/10, 5+len(top_pairs)/2),
        title = "Interactions between all cell types",
        max_size = 5,
        highlight_size = 0.75,
        degs_analysis = True,
        standard_scale = True,
        interaction_scores = cpdb_results['interaction_scores'],
        scale_alpha_by_interaction_scores = True
    )
    save_fig(p, out_path+'/plots/interaction_dotplot_top_pairs', output_format=output_format, plot_type='dotplot')

    p = kpy.plot_cpdb(
        adata = adata,
        cell_type1 = ".",
        cell_type2 = ".",
        means = cpdb_results['means'],
        pvals = cpdb_results['relevant_interactions'],
        celltype_key = cell_type_col,
        interacting_pairs = top_pairs,
        figsize = (9+len(cell_pair_cols)/10, 10+len(top_pairs)),
        title = "Interactions between all cell types\n grouped by pathway",
        max_size = 5,
        highlight_size = 0.75,
        degs_analysis = True,
        standard_scale = True,
        interaction_scores = cpdb_results['interaction_scores'],
        scale_alpha_by_interaction_scores = True
    )
    p = p + facet_wrap("~ classification", ncol = 1)
    save_fig(p, out_path+'/plots/interaction_dotplot_top_pairs_by_pathway', output_format=output_format, plot_type='dotplot')

    for i in range(n_top_genes):
        gi = unique_genes[i]
        interaction_count = cpdb_results['means']['interacting_pair'].str.contains(gi).sum()
        p = kpy.plot_cpdb(
            adata = adata,
            cell_type1 = ".",
            cell_type2 = ".",
            means = cpdb_results['means'],
            pvals = cpdb_results['relevant_interactions'],
            celltype_key = cell_type_col,
            genes = [gi],
            figsize = (9+len(cell_pair_cols)/10, 3+interaction_count/2),
            title = "Interactions between all cell types",
            max_size = 5,
            highlight_size = 0.75,
            degs_analysis = True,
            standard_scale = True,
            interaction_scores = cpdb_results['interaction_scores'],
            scale_alpha_by_interaction_scores = True
        )
        save_fig(p, out_path+f'/plots/interaction_dotplot_{gi}', output_format=output_format, plot_type='dotplot')

    # Filter the DataFrames to contain top pairs
    means_subset = cpdb_results['means'][cpdb_results['means']['interacting_pair'].isin(top_pairs)].copy()
    pvals_subset = cpdb_results['relevant_interactions'][cpdb_results['relevant_interactions']['interacting_pair'].isin(top_pairs)].copy()
    means_subset = means_subset.sort_values('interacting_pair')
    pvals_subset = pvals_subset.sort_values('interacting_pair')
    p = kpy.plot_cpdb_chord(
        adata = adata,
        cell_type1 = ".",
        cell_type2 = ".",
        means = means_subset,
        pvals = pvals_subset,
        deconvoluted = cpdb_results['deconvoluted'],
        celltype_key = cell_type_col,
        interaction = None,
        link_kwargs = {"direction": 1, "allow_twist": True, "r1": 90, "r2": 90},
        sector_text_kwargs = {"color": "black", "size": 12, "r": 105, "adjust_rotation": True},
        legend_kwargs = {"loc": "upper left", "bbox_to_anchor": (0.95, 1), "fontsize": 8},
        link_offset = 1
    )
    p = p.ax.get_figure()
    save_fig(p, out_path+'/plots/interaction_chord_diagram', output_format=output_format, plot_type='chord')
    plt.close(p)
    return cpdb_results


def run_cellphonedb_workflow(
        h5ad_file: str | None = None,
        cpdb_file_path: str = "/app/v5.0.0/cellphonedb.zip",
        out_path: str = ".",
        annotation_column: str = "clusters",
        species: Literal["human", "mouse"] = "human",
        method: Literal[1, 2, 3] = 2,
        n_top_interactions: int = 10,
        n_top_genes: int = 3,
        significant_threshold: float = 0.05,
        degs_file_path: str | None = None,
        group_column: str | None = None,
        disease_label: str | list[str] | None = None,
        control_label: str | None = None,
        deg_pval_cutoff: float = 0.05,
        deg_logfc_cutoff: float = 0.5,
        output_format: str = "png",
        n_jobs: int = 2,
    ) -> None:
    """Main function of CellPhoneDB workflow

    Args:
        h5ad_file (str, optional): input h5ad file. Defaults to None.
        cpdb_file_path (str, optional): cpdb database file location. Defaults to '/app/v5.0.0/cellphonedb.zip'.
        out_path (str, optional): output path. Defaults to '.'.
        annotation_column (str, optional): cell type column in adata.obs. Defaults to 'clusters'.
        species (str, optional): species, supporting human or mouse. Defaults to 'human'.
        method (int, optional): cpdb method to run (1,2,3). Defaults to 2.
        n_top_interactions (int, optional): number of top interactions to show on plots. Defaults to 10.
        n_top_genes (int, optional): number of top genes to show on plots. Defaults to 3.
        significant_threshold (float, optional): threshold for significance (method 2). Defaults to 0.05.
        degs_file_path (str, optional): user provided DEG file path. Defaults to None.
        group_column (str, optional): condition (group) column in adata.obs. Defaults to None.
        disease_label (str | list[str], optional): the treatment or disease group labels in conditions. Defaults to None.
        control_label (str, optional): the control or background group label in conditions. Defaults to None.
        deg_pval_cutoff (float, optional): P adjusted value cutoff. Defaults to 0.05.
        deg_logfc_cutoff (float, optional): log2 fold change cutoff. Defaults to 0.5.
        output_format (str, optional): output figure format. Defaults to "png".
        n_jobs (int, optional): number of parallel processes to use. Defaults to 2.
    """
    # Check if folders exist, if not create them
    out_path = out_path+'/cellphonedb'
    for folder in [out_path,
                   out_path+'/data',
                   out_path+'/method'+str(method)+'/tables',
                   out_path+'/method'+str(method)+'/plots']:
        Path(folder).mkdir(parents=True, exist_ok=True)

    # Process anndata
    if h5ad_file is None:
        h5ad_file = out_path+'/data/obj_spliced_unspliced.h5ad'
    adata = process_h5ad(h5ad_file, out_path)
    if species == 'mouse':
        adata = convert_mouse_genes(adata, out_path)
        counts_file_path = out_path+'/data/lognorm_humanized.h5ad'
    elif species == 'human':
        adata = adata.copy()
        counts_file_path = out_path+'/data/lognorm.h5ad'
    else:
        raise ValueError('This species is not supported.')

    # Make a cell type map
    metadata = pd.DataFrame({"Cell": adata.obs.index, "cell_type": adata.obs[annotation_column]})
    meta_file_path = out_path + '/celltype_metadata.txt'
    metadata.to_csv(meta_file_path, sep='\t', index=False)

    # Run the specified method
    if method == 1:
        _ = cellphonedb_method1(
            cpdb_file_path,
            meta_file_path,
            counts_file_path,
            out_path+'/method'+str(method)
        )
    elif method == 2:
        _ = cellphonedb_method2(
            adata,
            cpdb_file_path,
            meta_file_path,
            counts_file_path,
            out_path+'/method'+str(method),
            cell_type_col=annotation_column,
            n_top_interactions=n_top_interactions,
            n_top_genes=n_top_genes,
            threads=n_jobs,
            significant_threshold=significant_threshold,
            output_format=output_format
        )
    elif method == 3:
        if degs_file_path is None:
            if disease_label is None or control_label is None:
                raise ValueError('For DEG analysis, please specify disease and control groups.')
            # Compute DEGs if not supplied
            degs_file_path = calculate_deg(
                adata,
                cell_type_col=annotation_column,
                condition_col=group_column,
                disease_label=disease_label,
                control_label=control_label,
                pval_cutoff=deg_pval_cutoff,
                logfc_cutoff=deg_logfc_cutoff,
                out_path=out_path+'/method'+str(method)
            )
        _ = cellphonedb_method3(
            adata,
            cpdb_file_path,
            meta_file_path,
            counts_file_path,
            degs_file_path,
            out_path+'/method'+str(method),
            cell_type_col=annotation_column,
            n_top_interactions=n_top_interactions,
            n_top_genes=n_top_genes,
            threads=n_jobs,
            output_format=output_format
        )
    print(f'CellPhoneDB method {method} completed.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CellPhoneDB workflow.")

    parser.add_argument(
        "--h5ad_file",
        type=str,
        default=None,
        help="Input h5ad file path."
    )
    parser.add_argument(
        "--cpdb_file_path",
        type=str,
        default="/app/v5.0.0/cellphonedb.zip",
        help="CellPhoneDB database file location (default: '/app/v5.0.0/cellphonedb.zip')."
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
        "--species",
        type=str,
        choices=["human", "mouse"],
        default="human",
        help="Species, supporting 'human' or 'mouse' (default: 'human')."
    )
    parser.add_argument(
        "--method",
        type=int,
        choices=[1, 2, 3],
        default=2,
        help="CellPhoneDB method to run (1, 2, or 3) (default: 2)."
    )
    parser.add_argument(
        "--n_top_interactions",
        type=int,
        default=10,
        help="Number of top interactions to show on plots (default: 10)."
    )
    parser.add_argument(
        "--n_top_genes",
        type=int,
        default=3,
        help="Number of top genes to show on plots (default: 3)."
    )
    parser.add_argument(
        "--significant_threshold",
        type=float,
        default=0.05,
        help="Threshold for significance (method 2) (default: 0.05)."
    )
    parser.add_argument(
        "--degs_file_path",
        type=str,
        default=None,
        help="User provided DEG file path (default: None)."
    )
    parser.add_argument(
        "--group_column",
        type=str,
        default=None,
        help="Condition (group) column in adata.obs (default: None)."
    )
    parser.add_argument(
        "--disease_label",
        type=str,
        nargs='+',
        default=None,
        help="The treatment or disease group labels in conditions. Can pass multiple space-separated values."
    )
    parser.add_argument(
        "--control_label",
        type=str,
        default=None,
        help="The control or background group label in conditions (default: None)."
    )
    parser.add_argument(
        "--deg_pval_cutoff",
        type=float,
        default=0.05,
        help="P adjusted value cutoff (default: 0.05)."
    )
    parser.add_argument(
        "--deg_logfc_cutoff",
        type=float,
        default=0.5,
        help="log2 fold change cutoff (default: 0.5)."
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
    run_cellphonedb_workflow(**vars(args))
