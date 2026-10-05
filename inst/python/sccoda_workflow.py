import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pertpy as pt
import scanpy as sc
import seaborn as sns
from adjustText import adjust_text

warnings.simplefilter("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")

def save_fig(p, name, output_format='png'):
    """Save a figure

    Args:
        p (Figure): figure object from matplotlib or seaborn
        name (str): prefix of the figure to save
        output_format (str, optional): Output figure format. Defaults to 'png'.
    """
    if output_format == 'png':
        p.savefig(name + '.png', bbox_inches="tight", dpi=300)
    elif output_format == 'pdf':
        p.savefig(name + '.pdf', bbox_inches="tight")
    elif output_format == 'jpg' or output_format == 'jpeg':
        p.savefig(name + '.jpg', bbox_inches="tight", dpi=300)
    plt.close('all')


def run_sccoda_workflow(
        h5ad_file: str | None = None,
        out_path: str = ".",
        annotation_column: str = "clusters",
        sample_column: str | None = None,
        group_column: str | None = None,
        control_cond: str | None = None,
        reference_cell_type: str | None = None,
        additional_covariates: str | list[str] | None = None,
        significant_threshold: float = 0.05,
        output_format: str = "png",
    ) -> None:
    """Main function of scCODA workflow

    Args:
        h5ad_file (str, optional): input h5ad file. Defaults to None.
        out_path (str, optional): output path. Defaults to '.'.
        annotation_column (str, optional): cell type column in adata.obs. Defaults to 'clusters'.
        sample_column (str, optional): sample column in adata.obs. Defaults to None.
        group_column (str, optional): condition (group) column in adata.obs. Defaults to None.
        control_cond (str, optional): the control or background group label in conditions. Defaults to None.
        reference_cell_type (str, optional): the low-variance cell-type label in annotations. Defaults to None.
        additional_covariates (str | list[str], optional): additional covariates to include in statistical model. Defaults to None.
        significant_threshold (float): threshold for expected FDR to calculate credible effects. Defaults to 0.05.
        output_format (str, optional): output figure format. Defaults to "png".
    """
    out_path = out_path+'/sccoda'
    for folder in [out_path, 
                   out_path+'/tables',
                   out_path+'/plots']:
        Path(folder).mkdir(parents=True, exist_ok=True)

    # Reading data
    if h5ad_file is None:
        h5ad_file = out_path+'/data/obj_spliced_unspliced.h5ad'
    adata = sc.read_h5ad(h5ad_file)
    print(adata)

    if sample_column is None:
        raise ValueError('Please specify the sample column name in adata.obs.')
    if group_column is None or group_column == sample_column:
        print('group and sample columns are the same -- running with one replicate per condition.')
        adata.obs['sccoda'] = adata.obs[sample_column]
        group_column = 'sccoda'

    # Get unique conditions
    conditions = adata.obs[group_column].dropna().unique().tolist()
    if len(conditions) < 2:
        raise ValueError('Need at least two conditions to run comparisons.')

    # Determine which controls
    if control_cond is None:
        controls_to_run = conditions
    else:
        if control_cond not in conditions:
            raise ValueError(f"control_cond '{control_cond}' not found in group_column '{group_column}'.")
        controls_to_run = [control_cond]

    # Covariate setup
    if additional_covariates is None:
        additional_covariates = []
    elif isinstance(additional_covariates, str):
        additional_covariates = [additional_covariates]
    covariates = [group_column] + additional_covariates

    # Initialize and load the base model data once
    sccoda_model = pt.tl.Sccoda()

    general_plots = False

    sccoda_data_base = sccoda_model.load(
        adata,
        type="cell_level",
        generate_sample_level=True,
        cell_type_identifier=annotation_column,
        sample_identifier=sample_column,
        covariate_obs=covariates
    )

    if reference_cell_type is None:
        print("Calculating optimal global reference cell type...")
        sccoda_data = sccoda_model.prepare(
            sccoda_data_base.copy(),
            modality_key="coda",
            formula=group_column,
            reference_cell_type="automatic",
            automatic_reference_absence_threshold=0.1
        )
        reference_cell_type = sccoda_data["coda"].uns["scCODA_params"]["reference_cell_type"]

    for current_control in controls_to_run:
        print(f"\nRunning scCODA with baseline control: {current_control}.")

        # Use current_control as the Intercept
        formula = f"C({group_column}, Treatment('{current_control}'))"
        if additional_covariates:
            formula += " + " + " + ".join(additional_covariates)
        print(f"Formula: {formula}")

        # Prepare data
        sccoda_data = sccoda_model.prepare(
            sccoda_data_base.copy(),
            modality_key="coda",
            formula=formula,
            reference_cell_type="automatic" if reference_cell_type is None else reference_cell_type,
            automatic_reference_absence_threshold=0.1,
        )

        # Generate basic plots without statistical analysis
        if not general_plots:
            num_cell_types = sccoda_data["coda"].var.shape[0]
            dynamic_width_cells = max(20, num_cell_types * 0.4)

            fig = sccoda_model.plot_boxplots(sccoda_data, modality_key="coda", feature_name=group_column, plot_facets=False, add_dots=False, return_fig=True)
            fig.set_size_inches(dynamic_width_cells, 6)
            save_fig(fig, out_path+'/plots/boxplot', output_format)

            fig = sccoda_model.plot_boxplots(sccoda_data, modality_key="coda", feature_name=group_column, plot_facets=True, add_dots=True, return_fig=True)
            save_fig(fig, out_path+'/plots/boxplot_facets', output_format)

            num_samples = sccoda_data["coda"].obs.shape[0]
            dynamic_width_samples = max(10, num_samples * 0.5)

            fig = sccoda_model.plot_stacked_barplot(sccoda_data, feature_name="samples", return_fig=True)
            fig.set_size_inches(dynamic_width_samples, 6)
            save_fig(fig, out_path+'/plots/stacked_barplot', output_format)

            fig = sccoda_model.plot_rel_abundance_dispersion_plot(sccoda_data, return_fig=True)
            ax = fig.axes[0]
            ax.set_yscale("symlog", linthresh=0.01)
            ax.set_ylim(bottom=-0.001)
            if ax.collections:
                # Shrink the scatter dots
                scatter_plot = ax.collections[0]
                scatter_plot.set_sizes([15])
                # Extract the un-jittered locations of the red dots
                true_dot_coords = scatter_plot.get_offsets()
            texts = []
            for i, txt in enumerate(ax.texts):
                # Make the font significantly smaller
                txt.set_fontsize(6)
                # Force the text anchor to match the real red dot
                true_x, true_y = true_dot_coords[i]
                txt.set_position((true_x, true_y))
                texts.append(txt)
            # Automatically repel text and draw connecting lines
            adjust_text(texts, ax=ax, arrowprops={"arrowstyle": "-", "color": "grey", "lw": 0.5})

            save_fig(fig, out_path+'/plots/abundance_dispersion', output_format)

            general_plots = True

        # Run MCMC
        sccoda_model.run_nuts(sccoda_data, modality_key="coda", rng_key=0)
        sccoda_model.set_fdr(sccoda_data, modality_key="coda", est_fdr=significant_threshold)

        # Collect results for each comparison group against the current control
        comparison_groups = [g for g in conditions if g != current_control]
        effect_list = []
        plot_list = []
        safe_ctrl_name = current_control.replace(" ", "_").replace(",", "")

        for comp_group in comparison_groups:
            varm_keys = list(sccoda_data["coda"].varm.keys())
            target_substr = f"[T.{comp_group}]"
            possible_keys = [k for k in varm_keys if k.startswith("effect_df_") and target_substr in k]
            if not possible_keys:
                print(f"Warning: could not find effects for {comp_group} vs {current_control} in model output.")
                continue

            group_string = possible_keys[0]
            group_effects = sccoda_data["coda"].varm[group_string].copy()

            # Save raw table
            safe_comp_name = comp_group.replace(" ", "_").replace(",", "")
            group_effects.to_csv(f"{out_path}/tables/effects_{safe_comp_name}_vs_{safe_ctrl_name}.txt", sep="\t")

            # Prepare unfiltered data for plotting
            plot_group = group_effects[["log2-fold change", "Final Parameter"]].copy()
            plot_group["Cell Type"] = plot_group.index
            plot_group["Reference"] = current_control
            plot_group["Comp. Group"] = comp_group
            plot_list.append(plot_group)

            # Prepare filtered Data for summary table
            sig_group = plot_group[plot_group["Final Parameter"] != 0].copy()
            if not sig_group.empty:
                effect_list.append(sig_group)

                max_abs_fc = plot_group["log2-fold change"].abs().max()
                if max_abs_fc == 0:
                    max_abs_fc = 0.1

                sig_cell_types = sig_group["Cell Type"].tolist()
                # Create a temporary boolean sorting  column: True if significant, False if not
                rna_adata = sccoda_data.mod["rna"]
                rna_adata.obs["_temp_sig_sort"] = rna_adata.obs[annotation_column].isin(sig_cell_types)
                # Sort the object so False (grey dots) go first, True (colored dots) go last
                # Reassign them back into the MuData container
                sccoda_data.mod["rna"] = rna_adata[rna_adata.obs["_temp_sig_sort"].argsort()].copy()
                sccoda_data.update()
                if "X_umap" in sccoda_data["rna"].obsm:
                    fig_umap = sccoda_model.plot_effects_umap(
                        sccoda_data,
                        effect_name=group_string,
                        cluster_key=annotation_column,
                        return_fig=True,
                        title=f"Significant log2FC: {comp_group} vs {current_control}",
                        color_map="vlag",   # vlag is Blue(negative) -> Grey(zero) -> Red(positive)
                        vcenter=0,          # Forces neutral grey to 0.0
                        vmin=-max_abs_fc,
                        vmax=max_abs_fc,
                        size=3.5,
                        sort_order=False
                    )
                    if len(fig_umap.axes) > 1:
                        cbar_ax = fig_umap.axes[-1]
                        cbar_ax.set_ylabel("log2-fold change", rotation=270, labelpad=15, fontsize=12)
                    save_fig(fig_umap, f"{out_path}/plots/umap_effects.{safe_comp_name}_vs_{safe_ctrl_name}", output_format)

        plot_df = pd.concat(plot_list, ignore_index=True) if plot_list else pd.DataFrame()
        effect_df = pd.concat(effect_list, ignore_index=True) if effect_list else pd.DataFrame()

        if not effect_df.empty:
            # Save plots
            num_comp = len(comparison_groups)
            dynamic_width = max(8, 8 * num_comp)
            fig = sccoda_model.plot_effects_barplot(sccoda_data, return_fig=True)
            fig.fig.set_size_inches(dynamic_width, 4)
            save_fig(fig, f"{out_path}/plots/effects_barplot_all_comp.vs_{safe_ctrl_name}", output_format)

            # Identify which groups actually have significant effects
            valid_groups = effect_df["Comp. Group"].unique()

            # Filter plot_df to drop completely empty panels
            plot_df_filtered = plot_df[plot_df["Comp. Group"].isin(valid_groups)].copy()
            # If significant, keep its name. If not, label as 'Non-significant'
            plot_df_filtered["Color Group"] = plot_df_filtered.apply(
                lambda row: row["Cell Type"] if row["Final Parameter"] != 0 else "Non-significant", 
                axis=1
            )
            unique_cells = plot_df_filtered["Cell Type"].unique()
            base_colors = sns.color_palette("tab20", len(unique_cells)) # Standard categorical palette
            custom_pal = dict(zip(unique_cells, base_colors))
            custom_pal["Non-significant"] = "lightgrey"

            g = sns.catplot(
                data=plot_df_filtered,
                x="Cell Type",
                y="log2-fold change",
                row="Comp. Group",   # This forces vertical stacking
                hue="Color Group",
                palette=custom_pal,  # Apply the custom color
                kind="bar",
                height=3,            # Height of each panel (inches)
                aspect=3,            # Width multiplier (3 * 3 = 9 inches wide)
                sharex=True,         # Keep all cell types
                dodge=False,
                legend=False
            )

            # Clean up titles and labels
            g.set_titles("{row_name}", size=11)
            g.set_axis_labels("Cell Type", "log2-fold change")
            g.set_xticklabels(rotation=90, ha="center")
            g.fig.suptitle(f"Significant Effects (Reference: {current_control})", y=1.05, fontweight="bold", fontsize=12)
            save_fig(g, f"{out_path}/plots/effects_barplot.vs_{safe_ctrl_name}", output_format)

            # Save summary table
            effect_df.to_csv(f"{out_path}/tables/summary_credible_effects.vs_{safe_ctrl_name}.txt", sep="\t", index=False)
        else:
            print(f"No significant effects for reference: {current_control}")

    print('scCODA analysis completed.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scCODA workflow.")

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
        "--sample_column", 
        type=str, 
        default=None, 
        help="Sample column in adata.obs (default: None)."
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
        "--reference_cell_type", 
        type=str, 
        default=None, 
        help="The low-variance cell-type label in annotations (default: None)."
    )
    parser.add_argument(
        "--additional_covariates", 
        type=str, 
        nargs='+', 
        default=None, 
        help="Additional covariates to include in the statistical model. Can pass multiple space-separated values."
    )
    parser.add_argument(
        "--significant_threshold", 
        type=float, 
        default=0.05, 
        help="Threshold for expected FDR to calculate credible effects (default: 0.05)."
    )
    parser.add_argument(
        "--output_format", 
        type=str, 
        default="png", 
        help="Output figure format (default: 'png')."
    )

    args = parser.parse_args()

    # Pass all parsed arguments to the function as keyword arguments
    run_sccoda_workflow(**vars(args))
