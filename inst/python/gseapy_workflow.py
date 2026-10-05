import argparse
import shutil
import tempfile
import warnings
from pathlib import Path
from typing import Literal

import gseapy as gp
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from gseapy.scipalette import SciPalette

warnings.simplefilter("ignore", category=DeprecationWarning)

def draw_network(res, name, column, cutoff):
    # Network Visualization
    nodes, edges = gp.enrichment_map(res, column=column, cutoff=cutoff)

    G = nx.from_pandas_edgelist(edges,
            source='src_idx',
            target='targ_idx',
            edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])
    # Add missing node if there is any
    for node in nodes.index:
        if node not in G.nodes():
            G.add_node(node)

    # Draw network
    _, _ = plt.subplots(figsize=(10, 10))
    pos=nx.layout.spiral_layout(G)
    nx.draw_networkx_nodes(G,
        pos=pos,
        cmap=plt.cm.viridis,
        node_color=list(nodes.NES),
        node_size=list(nodes.Hits_ratio *1000))
    nx.draw_networkx_labels(G,
        pos=pos,
        labels=nodes.Term.to_dict(),
        font_size=10)
    edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
    nx.draw_networkx_edges(G,
        pos=pos,
        width=[x*10 for x in edge_weight],
        edge_color='#CDDBD4')
    plt.margins(0.3)
    plt.savefig(name)
    plt.close('all')


def analyze_result(
        res,
        out_path,
        gene_set_label,
        method,
        ct,
        exp_group,
        control_cond,
        cutoff_dotplot=0.05,
        n_top_terms=10,
        selected_terms=None,
        output_format="png"
    ):
    res_sorted = res.res2d.copy()
    res_sorted = res_sorted.sort_values(["FDR q-val", "NES"], ascending=[True, False])
    term = res_sorted.Term
    _ = res.plot(terms=term[:5], figsize=(8, 6), ofname=out_path+f'/plots/{gene_set_label}/{method}/{ct}_{exp_group}_vs_{control_cond}_gseaESn5.{output_format}')
    _ = res.plot(terms=term[:10], figsize=(8, 6), ofname=out_path+f'/plots/{gene_set_label}/{method}/{ct}_{exp_group}_vs_{control_cond}_gseaESn10.{output_format}')
    _ = res.plot(terms=term[:20], figsize=(8, 6), ofname=out_path+f'/plots/{gene_set_label}/{method}/{ct}_{exp_group}_vs_{control_cond}_gseaESn20.{output_format}')

    sig = res_sorted[res_sorted['FDR q-val'] < 0.25]
    top_terms = sig.Term.tolist()  # convert to list
    if len(top_terms) == 0:
        print("No pathways to plot with FDR < 0.25.")
    else:
        # Plot top n terms
        _ = res.plot(terms=top_terms[:n_top_terms], figsize=(8, 6), ofname=out_path+f'/plots/{gene_set_label}/{method}/{ct}_{exp_group}_vs_{control_cond}_gseaES0.25_top{n_top_terms}.{output_format}')

    if selected_terms is not None:
        if isinstance(selected_terms, str):
            search_strings = [selected_terms.lower()]
        else:
            search_strings = [x.lower() for x in selected_terms]
        terms_to_plot = [full_term for full_term in term if any(sub in full_term.lower() for sub in search_strings)]
        if terms_to_plot:
            _ = res.plot(terms=terms_to_plot, figsize=(8, 6), ofname=out_path+f'/plots/{gene_set_label}/{method}/{ct}_{exp_group}_vs_{control_cond}_gseaES_selected_terms.{output_format}')

    if (res_sorted['FDR q-val'] < cutoff_dotplot).sum() == 0:
        print('GSEA dotplot cutoff set to 1.')
        cutoff_dotplot_tmp = 1
        title_suffix = ' (no significance)'
    else:
        cutoff_dotplot_tmp = cutoff_dotplot
        title_suffix = ''
    _ = gp.dotplot(res_sorted[:n_top_terms],
        column="FDR q-val",
        title=gene_set_label + title_suffix,
        cmap=plt.cm.viridis,
        size=5,
        figsize=(5,5),
        cutoff=cutoff_dotplot_tmp,
        top_term=n_top_terms,
        ofname=out_path+f'/plots/{gene_set_label}/{method}/{ct}_{exp_group}_vs_{control_cond}_dotplotFDR_top{n_top_terms}.{output_format}')

    _ = gp.dotplot(res.res2d,
        column="NES",
        title=gene_set_label + title_suffix,
        cmap=plt.cm.viridis,
        size=5,
        figsize=(5,5),
        top_term=n_top_terms,
        ofname=out_path+f'/plots/{gene_set_label}/{method}/{ct}_{exp_group}_vs_{control_cond}_dotplotNES_top{n_top_terms}.{output_format}')

    draw_network(res_sorted[:n_top_terms],
        out_path+f'/plots/{gene_set_label}/{method}/{ct}_{exp_group}_vs_{control_cond}_network.{output_format}',
        column='FDR q-val',
        cutoff=cutoff_dotplot_tmp)


def run_gseapy_workflow(
        h5ad_file: str | None = None,
        species: Literal["human", "mouse"] = "human",
        out_path: str = ".",
        annotation_column: str = "clusters",
        group_column: str | None = None,
        control_cond: str | None = None,
        selected_cell_types: str | list[str] | None = None,
        n_top_terms: int = 10,
        gsea_graph_num: int = 20,
        significance_threshold: float = 0.05,
        cutoff_dotplot: float = 0.05,
        selected_terms: str | list[str] | None = None,
        output_format: str = "png",
        n_jobs: int = 2,
    ) -> None:
    """Main function of GSEApy workflow

    Args:
        h5ad_file (str, optional): input h5ad file. Defaults to None.
        species (str, optional): species to use for enrichment analysis ('human' or 'mouse'). Defaults to 'human'.
        out_path (str, optional): output path. Defaults to '.'.
        annotation_column (str, optional): cell type column in adata.obs. Defaults to 'clusters'.
        group_column (str, optional): condition (group) column in adata.obs. Defaults to None.
        control_cond (str, optional): the control or background group label in conditions. Defaults to None.
        selected_cell_types (str | list[str], optional): specific cell types to run the analysis on. Defaults to None.
        n_top_terms (int, optional): number of top enriched terms to show on plots. Defaults to 10.
        gsea_graph_num (int, optional): number of GSEA enrichment plots for individual pathway to generate. Defaults to 20.
        significance_threshold (float, optional): threshold for significance (pvals_adj of scanpy.rank_genes_groups). Defaults to 0.05.
        cutoff_dotplot (float, optional): cutoff value for dotplots. Defaults to 0.05.
        selected_terms (str | list[str], optional): specific pathway terms of interest to plot. Defaults to None.
        output_format (str, optional): output figure format. Defaults to "png".
        n_jobs (int, optional): number of parallel processes to use. Defaults to 2.
    """
    out_path = out_path+'/gseapy'
    for folder in [out_path,
                   out_path+'/tables/DEGs',
                   out_path+'/tables/GO_BP/GSEA',
                   out_path+'/tables/GO_BP/prerank',
                   out_path+'/tables/GO_BP/enrichr',
                   out_path+'/tables/Hallmark/GSEA',
                   out_path+'/tables/Hallmark/prerank',
                   out_path+'/tables/Hallmark/enrichr',
                   out_path+'/plots/GO_BP/GSEA',
                   out_path+'/plots/GO_BP/prerank',
                   out_path+'/plots/GO_BP/enrichr',
                   out_path+'/plots/Hallmark/GSEA',
                   out_path+'/plots/Hallmark/prerank',
                   out_path+'/plots/Hallmark/enrichr']:
        Path(folder).mkdir(parents=True, exist_ok=True)

    # Reading data
    if h5ad_file is None:
        h5ad_file = out_path+'/data/obj_spliced_unspliced.h5ad'
    adata = sc.read_h5ad(h5ad_file)
    print(adata)

    # Preprocess data
    adata.layers['counts'] = adata.X
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X
    background_genes = adata.var_names.tolist()
    if len(background_genes) < 10000:
        print("WARNING: Does the input h5ad include all genes?")

    if group_column is None:
        raise ValueError('Please specify the group column name in adata.obs.')
    if species == 'mouse':
        gene_sets = ['/app/gmt_files/m5.go.bp.v2026.1.Mm.symbols.gmt','/app/gmt_files/mh.all.v2026.1.Mm.symbols.gmt']
        # gene_sets = ['GO_Biological_Process_2026','MSigDB_Hallmark_2020']
        gene_set_labels = ["GO_BP", "Hallmark"]
    elif species == 'human':
        gene_sets=['/app/gmt_files/c5.go.bp.v2026.1.Hs.symbols.gmt','/app/gmt_files/h.all.v2026.1.Hs.symbols.gmt']
        # gene_sets = ['GO_Biological_Process_2026', 'MSigDB_Hallmark_2020']
        gene_set_labels = ["GO_BP", "Hallmark"]
    else:
        raise ValueError(f'This species {species} is not supported.')

    if selected_cell_types is not None:
        if isinstance(selected_cell_types, str):
            cell_types = [selected_cell_types]
        else:
            cell_types = selected_cell_types
    else:
        cell_types = adata.obs[annotation_column].unique()
    all_groups = adata.obs[group_column].unique()

    # Remove the control group from the list and only loop through contrasting groups
    if control_cond is None:
        raise ValueError("Please specify the control group with 'control_cond'.")
    experimental_groups = [g for g in all_groups if g != control_cond]

    # Configure matplotlib save settings
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.facecolor'] = 'white'
    plt.rcParams['savefig.bbox'] = 'tight'

    for ct in cell_types:
        print(f'======Analyzing {ct}...')
        for exp_group in experimental_groups:
            print(f'------Comparing {exp_group} and {control_cond}...')
            # Create boolean masks for the subset
            mask_celltype = adata.obs[annotation_column] == ct
            mask_groups = adata.obs[group_column].isin([control_cond, exp_group])

            bdata = adata[mask_celltype & mask_groups].copy()
            # Set categorical orders
            bdata.obs[group_column] = pd.Categorical(bdata.obs[group_column], categories=[exp_group, control_cond], ordered=True)
            indices = bdata.obs.sort_values([annotation_column, group_column]).index
            bdata = bdata[indices,:].copy()
            print(f'Selected {bdata.n_obs} cells.')

            # Check number of cells in each condition
            group_counts = bdata.obs[group_column].value_counts()
            num_control = group_counts.get(control_cond, 0)
            num_exp = group_counts.get(exp_group, 0)

            print(f"{control_cond} has {num_control} cells.")
            print(f"{exp_group} has {num_exp} cells.")
            if num_control < 3 or num_exp < 3:
                print(f"Skipping {ct} - {exp_group} vs {control_cond}. Not enough cells.")
                continue

            # DEG Analysis
            sc.tl.rank_genes_groups(bdata,
                groupby=group_column,
                use_raw=False,
                method='wilcoxon',
                groups=[exp_group],
                reference=control_cond)
            result = bdata.uns['rank_genes_groups']  # ranked by Wilcoxon scores (z-score)
            degs = pd.DataFrame(
                {group + '_' + key: result[key][group]
                for group in result['names'].dtype.names for key in ['names','scores', 'pvals','pvals_adj','logfoldchanges']})
            degs.to_csv(out_path+f'/tables/DEGs/{ct}_{exp_group}_vs_{control_cond}_degs.csv', index=False)

            for gene_set, gene_set_label in zip(gene_sets, gene_set_labels):
                print(f'Using {gene_set_label}: {gene_set}')
                with tempfile.TemporaryDirectory() as tmpdirname:
                    tmp_path = Path(tmpdirname)

                    # Standard GSEA
                    res = gp.gsea(data=bdata.to_df().T, # row -> genes, column-> samples
                        gene_sets=gene_set,
                        cls=bdata.obs[group_column],
                        organism=species,
                        permutation_num=1000,
                        permutation_type='phenotype',
                        method='s2n', # signal_to_noise
                        threads=n_jobs,
                        outdir=tmpdirname,
                        format=output_format,
                        graph_num=gsea_graph_num,
                        seed=7)

                    # Copy specific items to the final directory
                    source_csv = tmp_path / 'gseapy.phenotype.gsea.report.csv'
                    source_plots = tmp_path / 'gsea'
                    if source_csv.exists():
                        shutil.copy(source_csv, out_path+f'/tables/{gene_set_label}/GSEA/{ct}_{exp_group}_vs_{control_cond}.csv')
                    if source_plots.exists() and source_plots.is_dir():
                        dest_plots = Path(out_path+f'/plots/{gene_set_label}/GSEA/{ct}_{exp_group}_vs_{control_cond}')
                        if dest_plots.exists():
                            shutil.rmtree(dest_plots)
                        shutil.copytree(source_plots, dest_plots)

                    analyze_result(
                        res,
                        out_path,
                        gene_set_label,
                        "GSEA",
                        ct,
                        exp_group,
                        control_cond,
                        cutoff_dotplot,
                        n_top_terms,
                        selected_terms,
                        output_format)

                    # Prerank
                    rnk_df = degs.loc[:, [f'{exp_group}_names', f'{exp_group}_logfoldchanges']].dropna()
                    pre_res = gp.prerank(rnk=rnk_df,
                        gene_sets=gene_set,
                        organism=species,
                        method="multilevel",            # use fgsea multilevel p-values
                        sample_size=101,                # multilevel sampling depth (default 101)
                        eps=1e-50,                      # smallest p-value resolved; 0 => machine precision
                        min_size=5,
                        max_size=1000,
                        threads=n_jobs,
                        outdir=tmpdirname,
                        format=output_format,
                        graph_num=gsea_graph_num,
                        seed=6)

                    # Copy specific items to the final directory
                    source_csv = tmp_path / 'gseapy.gene_set.prerank.report.csv'
                    source_plots = tmp_path / 'prerank'
                    if source_csv.exists():
                        shutil.copy(source_csv, out_path+f'/tables/{gene_set_label}/prerank/{ct}_{exp_group}_vs_{control_cond}.csv')
                    if source_plots.exists() and source_plots.is_dir():
                        dest_plots = Path(out_path+f'/plots/{gene_set_label}/prerank/{ct}_{exp_group}_vs_{control_cond}')
                        if dest_plots.exists():
                            shutil.rmtree(dest_plots)
                        shutil.copytree(source_plots, dest_plots)

                    analyze_result(
                        pre_res,
                        out_path,
                        gene_set_label,
                        "prerank",
                        ct,
                        exp_group,
                        control_cond,
                        cutoff_dotplot,
                        n_top_terms,
                        selected_terms,
                        output_format)

                    # Over-representation analysis using Enrichr
                    # subset up or down regulated genes
                    degs_sig = degs[degs[degs.columns[3]] < significance_threshold]
                    degs_up = degs_sig[degs_sig[degs_sig.columns[4]] > 0]
                    degs_dw = degs_sig[degs_sig[degs_sig.columns[4]] < 0]
                    print(f"There are {degs_up.shape[0]} significantly upregulated, and {degs_dw.shape[0]} significantly downregulated genes.")

                    # Up result
                    genes = degs_up[degs_up.columns[0]].tolist()
                    if len(genes) < 5:
                        print("Skipping Enrichr for upregulated genes. Need at least 5.")
                        enr_up = None
                    else:
                        enr_up = gp.enrichr(genes,
                            gene_sets=gene_set,
                            organism=species,
                            background=background_genes,
                            outdir=None)

                        if (enr_up.res2d['Adjusted P-value'] < cutoff_dotplot).sum() == 0:
                            print('Enrichr UP dotplot cutoff set to 1.')
                            cutoff_dotplot_tmp = 1
                        else:
                            cutoff_dotplot_tmp = cutoff_dotplot
                        _ = gp.dotplot(enr_up.res2d,
                            title=f'{gene_set_label} UP',
                            cmap = plt.cm.autumn_r,
                            size=5,
                            figsize=(5,5),
                            cutoff=cutoff_dotplot_tmp,
                            top_term=n_top_terms,
                            ofname=out_path+f'/plots/{gene_set_label}/enrichr/{ct}_{exp_group}_vs_{control_cond}_dotplotUP_top{n_top_terms}.{output_format}')

                        enr_up.res2d['-log10(FDR)']=-np.log10(enr_up.res2d['Adjusted P-value'])
                        enr_up.res2d = enr_up.res2d.sort_values(by=['Adjusted P-value'])
                        _ = gp.dotplot(enr_up.res2d,
                            column='Combined Score',
                            x='-log10(FDR)',
                            title=f'{gene_set_label} Up',
                            cmap = plt.cm.autumn_r,
                            size=5,
                            figsize=(5,5),
                            top_term=n_top_terms,
                            ofname=out_path+f'/plots/{gene_set_label}/enrichr/{ct}_{exp_group}_vs_{control_cond}_dotplotUP2_top{n_top_terms}.{output_format}')

                    # Down result
                    genes = degs_dw[degs_dw.columns[0]].tolist()
                    if len(genes) < 5:
                        print("Skipping Enrichr for downregulated genes. Need at least 5.")
                        enr_dw = None
                    else:
                        enr_dw = gp.enrichr(degs_dw[degs_dw.columns[0]],
                            gene_sets=gene_set,
                            organism=species,
                            background=background_genes,
                            outdir=None)

                        if (enr_dw.res2d['Adjusted P-value'] < cutoff_dotplot).sum() == 0:
                            print('Enrichr Down dotplot cutoff set to 1.')
                            cutoff_dotplot_tmp = 1
                        else:
                            cutoff_dotplot_tmp = cutoff_dotplot
                        _ = gp.dotplot(enr_dw.res2d,
                            title=f'{gene_set_label} DOWN',
                            cmap = plt.cm.winter_r,
                            size=5,
                            figsize=(5,5),
                            cutoff=cutoff_dotplot_tmp,
                            top_term=n_top_terms,
                            ofname=out_path+f'/plots/{gene_set_label}/enrichr/{ct}_{exp_group}_vs_{control_cond}_dotplotDW_top{n_top_terms}.{output_format}')

                        enr_dw.res2d['-log10(FDR)']=-np.log10(enr_dw.res2d['Adjusted P-value'])
                        enr_dw.res2d = enr_dw.res2d.sort_values(by=['Adjusted P-value'])
                        _ = gp.dotplot(enr_dw.res2d,
                            column='Combined Score',
                            x='-log10(FDR)',
                            title=f'{gene_set_label} Up',
                            cmap = plt.cm.winter_r,
                            size=5,
                            figsize=(5,5),
                            top_term=n_top_terms,
                            ofname=out_path+f'/plots/{gene_set_label}/enrichr/{ct}_{exp_group}_vs_{control_cond}_dotplotDW2_top{n_top_terms}.{output_format}')

                    # Combine result
                    valid_res2d = []
                    color_dict = {'UP': 'r', 'DOWN': 'b'}

                    if enr_up is not None and hasattr(enr_up, 'res2d') and not enr_up.res2d.empty:
                        up_df = enr_up.res2d.copy()
                        up_df['Direction'] = "UP"
                        up_df = up_df.sort_values(by=['Adjusted P-value'])
                        valid_res2d.append(up_df)

                    if enr_dw is not None and hasattr(enr_dw, 'res2d') and not enr_dw.res2d.empty:
                        dw_df = enr_dw.res2d.copy()
                        dw_df['Direction'] = "DOWN"
                        dw_df = dw_df.sort_values(by=['Adjusted P-value'])
                        valid_res2d.append(dw_df)

                    # Proceed if at least one direction had valid enrichment
                    if valid_res2d:
                        enr_res_full = pd.concat(valid_res2d)
                        enr_res_full.to_csv(out_path+f'/tables/{gene_set_label}/enrichr/{ct}_{exp_group}_vs_{control_cond}.csv', index=False)

                        valid_top_terms = [df.head(n_top_terms) for df in valid_res2d]
                        enr_res = pd.concat(valid_top_terms)

                        plot_directions = sorted(enr_res['Direction'].unique())
                        plot_colors = [color_dict[d] for d in plot_directions]

                        sci = SciPalette()
                        NbDr = sci.create_colormap()
                        # display multi-datasets
                        _ = gp.dotplot(enr_res,
                            x='Direction',
                            x_order=plot_directions,
                            title=gene_set_label,
                            cmap = NbDr.reversed(),
                            size=3,
                            figsize=(5,5),
                            cutoff=cutoff_dotplot,
                            show_ring=True,
                            top_term=n_top_terms,
                            ofname=out_path+f'/plots/{gene_set_label}/enrichr/{ct}_{exp_group}_vs_{control_cond}_dotplot_top{n_top_terms}.{output_format}')

                        _ = gp.barplot(enr_res,
                            group='Direction',
                            title =gene_set_label,
                            figsize=(5,5),
                            cutoff=cutoff_dotplot,
                            color=plot_colors,
                            top_term=n_top_terms,
                            ofname=out_path+f'/plots/{gene_set_label}/enrichr/{ct}_{exp_group}_vs_{control_cond}_barplot_top{n_top_terms}.{output_format}')
                    else:
                        print("Skipping plotting: No valid enrichr results.")

    print('GSEApy analysis completed.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run GSEApy workflow.")

    parser.add_argument(
        "--h5ad_file",
        type=str,
        default=None,
        help="Input h5ad file path."
    )
    parser.add_argument(
        "--species",
        type=str,
        choices=["human", "mouse"],
        default="human",
        help="Species to use for enrichment analysis ('human' or 'mouse') (default: 'human')."
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
        "--group_column",
        type=str,
        default=None,
        help="Condition (group) column in adata.obs (default: None)."
    )
    parser.add_argument(
        "--control_cond",
        type=str,
        default=None,
        help="The control or background group label in conditions (default: None)."
    )
    parser.add_argument(
        "--selected_cell_types",
        type=str,
        nargs='+',
        default=None,
        help="Specific cell types to run the analysis on. Can pass multiple space-separated values (default: None)."
    )
    parser.add_argument(
        "--n_top_terms",
        type=int,
        default=10,
        help="Number of top enriched terms to show on plots (default: 10)."
    )
    parser.add_argument(
        "--gsea_graph_num",
        type=int,
        default=20,
        help="Number of GSEA enrichment plots for individual pathway to generate (default: 20)."
    )
    parser.add_argument(
        "--significance_threshold",
        type=float,
        default=0.05,
        help="Threshold for significance (pvals_adj of scanpy.rank_genes_groups) (default: 0.05)."
    )
    parser.add_argument(
        "--cutoff_dotplot",
        type=float,
        default=0.05,
        help="Cutoff value for dotplots (default: 0.05)."
    )
    parser.add_argument(
        "--selected_terms",
        type=str,
        nargs='+',
        default=None,
        help="Specific pathway terms of interest to plot (default: None)."
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
    run_gseapy_workflow(**vars(args))
