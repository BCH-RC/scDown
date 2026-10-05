#!/usr/bin/env python3

import os
import re
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread


def clean_barcode(barcode: str) -> str:
    """Remove trailing '-1' if present."""
    return barcode[:-2] if barcode.endswith("-1") else barcode


def discover_sample_dirs(starsolo_root: str) -> dict:
    """
    Discover STARsolo velocyto raw directories.

    Expected structure somewhere under starsolo_root:
      <sample>_Solo.out/Velocyto/raw

    Returns:
      dict: {sample_name: raw_dir}
    """
    starsolo_root = Path(starsolo_root)
    sample_dirs = {}

    for raw_dir in starsolo_root.rglob("Velocyto/raw"):
        raw_dir = raw_dir.resolve()

        required = [
            raw_dir / "spliced.mtx",
            raw_dir / "unspliced.mtx",
            raw_dir / "barcodes.tsv",
            raw_dir / "features.tsv",
        ]
        if not all(f.exists() for f in required):
            continue

        parent_name = raw_dir.parents[1].name  # e.g. Dupco1_Solo.out
        sample_name = re.sub(r"_Solo\.out$", "", parent_name)

        sample_dirs[sample_name] = str(raw_dir)

    if not sample_dirs:
        raise FileNotFoundError(
            f"No STARsolo Velocyto/raw directories found under: {starsolo_root}"
        )

    return dict(sorted(sample_dirs.items()))


def insert_sparse(full_matrix, small_matrix, rows, cols):
    """
    Insert a smaller sparse matrix into a full-sized sparse matrix using row/col mapping.
    """
    coo = small_matrix.tocoo()
    global_rows = rows[coo.row]
    global_cols = cols[coo.col]

    update = sparse.coo_matrix(
        (coo.data, (global_rows, global_cols)),
        shape=full_matrix.shape,
        dtype=full_matrix.dtype,
    ).tocsr()

    return full_matrix + update


def add_velocity_layers(
    input_h5ad: str,
    starsolo_root: str,
    output_dir: str,
    obs_separator: str = "_",
    output_name: str | None = None,
):
    """
    Add spliced/unspliced layers from STARsolo Velocyto outputs to an h5ad file.

    Parameters
    ----------
    input_h5ad : str
        Input h5ad file.
    starsolo_root : str
        Root directory containing STARsolo outputs.
    output_dir : str
        Folder where output h5ad will be written.
    obs_separator : str
        Separator between sample name and barcode in adata.obs_names.
        Example: Dupco1_AAAC... -> separator is "_"
    output_name : str | None
        Optional output filename. If None, derived from input name.
    """
    os.makedirs(output_dir, exist_ok=True)

    if output_name is None:
        base = Path(input_h5ad).stem
        output_name = f"{base}_spliced_unspliced.h5ad"

    output_file = str(Path(output_dir) / output_name)

    print("\n========================================")
    print("Reading h5ad")
    print("========================================")
    adata = sc.read_h5ad(input_h5ad)

    obs_names = np.asarray(adata.obs_names.astype(str))
    adata_genes = np.asarray(adata.var_names.astype(str))

    n_cells = len(obs_names)
    n_genes = len(adata_genes)

    print(f"H5AD dimensions: {n_cells} x {n_genes}")

    print("\n========================================")
    print("Discovering STARsolo sample directories")
    print("========================================")
    sample_dirs = discover_sample_dirs(starsolo_root)

    print(f"Found {len(sample_dirs)} samples:")
    for sample_name, sample_dir in sample_dirs.items():
        print(f"  {sample_name}: {sample_dir}")

    print("\n========================================")
    print("Extracting barcodes from obs_names")
    print("========================================")

    adata_barcodes = np.array([
        x.split(obs_separator, 1)[1] if obs_separator in x else x
        for x in obs_names
    ])
    adata_barcodes = np.array([clean_barcode(x) for x in adata_barcodes])

    print("\nCreating full-size matrices...")
    spliced_all = sparse.csr_matrix((n_cells, n_genes), dtype=np.float32)
    unspliced_all = sparse.csr_matrix((n_cells, n_genes), dtype=np.float32)

    for sample_name, sample_dir in sample_dirs.items():
        print("\n\n########################################")
        print(f"# SAMPLE: {sample_name}")
        print("########################################")

        splice_file = os.path.join(sample_dir, "spliced.mtx")
        unsplice_file = os.path.join(sample_dir, "unspliced.mtx")
        barcode_file = os.path.join(sample_dir, "barcodes.tsv")
        features_file = os.path.join(sample_dir, "features.tsv")

        for f in [splice_file, unsplice_file, barcode_file, features_file]:
            if not os.path.exists(f):
                raise FileNotFoundError(f"Missing file for {sample_name}: {f}")

        sample_prefix = sample_name + obs_separator
        sample_idx = np.array([x.startswith(sample_prefix) for x in obs_names])
        sample_cell_idx = np.where(sample_idx)[0]
        sample_barcodes = adata_barcodes[sample_cell_idx]

        print(f"H5AD cells: {len(sample_cell_idx)}")

        if len(sample_cell_idx) == 0:
            print(f"WARNING: No cells found for {sample_name}")
            continue

        print("Reading spliced.mtx...")
        splice = mmread(splice_file).tocsr()

        print("Reading unspliced.mtx...")
        unsplice = mmread(unsplice_file).tocsr()

        print(f"Spliced: {splice.shape[0]} x {splice.shape[1]}")
        print(f"Unspliced: {unsplice.shape[0]} x {unsplice.shape[1]}")

        mtx_barcodes = pd.read_csv(
            barcode_file,
            header=None,
            sep="\t",
            dtype=str
        ).iloc[:, 0].to_numpy()

        print(f"MTX barcodes: {len(mtx_barcodes)}")

        mtx_barcodes_clean = np.array([clean_barcode(x) for x in mtx_barcodes])

        if splice.shape[1] != len(mtx_barcodes):
            raise ValueError(
                f"{sample_name}: spliced.mtx columns ({splice.shape[1]}) != barcode count ({len(mtx_barcodes)})"
            )
        if unsplice.shape[1] != len(mtx_barcodes):
            raise ValueError(
                f"{sample_name}: unspliced.mtx columns ({unsplice.shape[1]}) != barcode count ({len(mtx_barcodes)})"
            )

        barcode_to_idx = {barcode: i for i, barcode in enumerate(mtx_barcodes_clean)}

        cell_idx = np.array([barcode_to_idx.get(barcode, -1) for barcode in sample_barcodes])
        valid_cells = cell_idx >= 0

        print(f"Matched cells: {valid_cells.sum()} / {len(sample_barcodes)}")

        if valid_cells.sum() == 0:
            raise ValueError(f"No matching barcodes found for {sample_name}")

        if (~valid_cells).any():
            print(f"Missing cells: {(~valid_cells).sum()}")
            print("First missing barcodes:")
            for barcode in sample_barcodes[~valid_cells][:10]:
                print(f"  {barcode}")

        features = pd.read_csv(
            features_file,
            header=None,
            sep="\t",
            dtype=str
        )

        mtx_gene_ids = features.iloc[:, 0].astype(str).to_numpy()
        mtx_gene_ids = np.array([x.split(".", 1)[0] for x in mtx_gene_ids])

        mtx_gene_names = features.iloc[:, 1].astype(str).to_numpy()

        print(f"MTX genes: {len(mtx_gene_ids)}")

        gene_id_to_idx = {}
        for i, gene_id in enumerate(mtx_gene_ids):
            if gene_id not in gene_id_to_idx:
                gene_id_to_idx[gene_id] = i

        gene_name_to_idx = {}
        for i, gene_name in enumerate(mtx_gene_names):
            if gene_name not in gene_name_to_idx:
                gene_name_to_idx[gene_name] = i

        gene_idx_id = np.array([gene_id_to_idx.get(gene, -1) for gene in adata_genes])
        gene_idx_name = np.array([gene_name_to_idx.get(gene, -1) for gene in adata_genes])

        n_id = np.sum(gene_idx_id >= 0)
        n_name = np.sum(gene_idx_name >= 0)

        print(f"Genes matched by ID: {n_id}")
        print(f"Genes matched by name: {n_name}")

        gene_idx = gene_idx_name.copy()
        use_id = (gene_idx_name < 0) & (gene_idx_id >= 0)
        gene_idx[use_id] = gene_idx_id[use_id]
        valid_genes = gene_idx >= 0

        print(f"Genes matched using name: {np.sum(gene_idx_name >= 0)}")
        print(f"Genes matched using ID fallback: {use_id.sum()}")
        print(f"Total matched genes: {valid_genes.sum()} / {n_genes}")

        valid_cell_indices = cell_idx[valid_cells]

        splice_cells = splice[:, valid_cell_indices]
        unsplice_cells = unsplice[:, valid_cell_indices]

        valid_gene_indices = gene_idx[valid_genes]

        splice_matched = splice_cells[valid_gene_indices, :]
        unsplice_matched = unsplice_cells[valid_gene_indices, :]

        h5ad_rows = sample_cell_idx[valid_cells]
        h5ad_cols = np.where(valid_genes)[0]

        splice_insert = splice_matched.T.tocsr()
        unsplice_insert = unsplice_matched.T.tocsr()

        spliced_all = insert_sparse(spliced_all, splice_insert, h5ad_rows, h5ad_cols)
        unspliced_all = insert_sparse(unspliced_all, unsplice_insert, h5ad_rows, h5ad_cols)

        print("\n----------------------------------------")
        print(f"Sample: {sample_name}")
        print(f"Cells matched: {valid_cells.sum()}")
        print(f"Genes matched: {valid_genes.sum()}")
        print("----------------------------------------")

        print("\n----------------------------------------")
        print(f"COUNT CHECK: {sample_name}")
        print("----------------------------------------")
        print(f"Spliced nonzero: {splice_matched.nnz}")
        print(f"Unspliced nonzero: {unsplice_matched.nnz}")
        print(f"Spliced total: {splice_matched.sum()}")
        print(f"Unspliced total: {unsplice_matched.sum()}")

    print("\n\n========================================")
    print("FINAL CHECK")
    print("========================================")
    print(f"Expected: {n_cells} x {n_genes}")
    print(f"Spliced: {spliced_all.shape}")
    print(f"Unspliced: {unspliced_all.shape}")

    expected_shape = (n_cells, n_genes)

    if spliced_all.shape != expected_shape:
        raise ValueError("Spliced dimensions do not match h5ad.")
    if unspliced_all.shape != expected_shape:
        raise ValueError("Unspliced dimensions do not match h5ad.")

    print("\n========================================")
    print("Adding layers")
    print("========================================")

    adata.layers["spliced"] = spliced_all
    adata.layers["unspliced"] = unspliced_all

    print("\nAvailable layers:")
    print(list(adata.layers.keys()))

    print("\nLayer shapes:")
    for layer in ["spliced", "unspliced"]:
        x = adata.layers[layer]
        print(
            f"{layer}: shape={x.shape}, "
            f"nnz={x.nnz if sparse.issparse(x) else np.count_nonzero(x)}, "
            f"total={x.sum()}"
        )

    print("\n========================================")
    print("Writing h5ad")
    print("========================================")

    adata.write_h5ad(output_file, compression="gzip")

    print("\nDONE!")
    print("Output:")
    print(output_file)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add STARsolo Velocyto spliced/unspliced layers to an h5ad file."
    )
    parser.add_argument(
        "-i", "--input-h5ad",
        required=True,
        help="Input h5ad file"
    )
    parser.add_argument(
        "-s", "--starsolo-root",
        required=True,
        help="Root folder containing STARsolo outputs"
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--obs-separator",
        default="_",
        help="Separator between sample and barcode in adata.obs_names (default: _)"
    )
    parser.add_argument(
        "--output-name",
        default=None,
        help="Optional output h5ad filename"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    add_velocity_layers(
        input_h5ad=args.input_h5ad,
        starsolo_root=args.starsolo_root,
        output_dir=args.output_dir,
        obs_separator=args.obs_separator,
        output_name=args.output_name,
    )


if __name__ == "__main__":
    main()