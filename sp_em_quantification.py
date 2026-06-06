"""
sp_em_quantification.py
-----------------------
Single-cell / spatial transcriptomics isoform quantification built on
miniQuant's community-level EM framework.

Entry point: run_sc_em_quantification()
  - SC mode  (output_dir_name='SC_output', barcode_label='cell_barcode')
  - ST mode  (output_dir_name='ST_output', barcode_label='spot_barcode')
    Additional ST parameters: barcode_whitelist_path, tissue_positions_path

EM mode: per-unit EM
  For each gene community:
    1. Build molecule cond_prob via UMI dedup (element-wise max across reads
       sharing the same barcode+UMI).
    2. For each cell/spot: run unit-level EM independently.
    3. Hard-assign each molecule to its best isoform (argmax posterior).
"""

import glob
import os
import pickle
import re
import time
from collections import defaultdict
from pathlib import Path

import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.sparse


# ---------------------------------------------------------------------------
# Spatial helpers
# ---------------------------------------------------------------------------

def _load_barcode_whitelist(path):
    """Load spatial barcode whitelist; return frozenset or None."""
    if path is None:
        return None
    valid = set()
    with open(path) as fh:
        for line in fh:
            bc = line.strip()
            if bc:
                valid.add(bc)
    print(f'[ST-EM] Barcode whitelist: {len(valid)} valid barcodes loaded', flush=True)
    return frozenset(valid)


def _load_tissue_positions(path):
    """
    Load Visium tissue_positions.csv.

    Supports both classic Visium (no header, 6 columns) and Visium HD / Space Ranger
    v2+ (with header row). Returns a DataFrame indexed by barcode with columns:
    in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres.
    Returns None if path is None.
    """
    if path is None:
        return None
    # Peek at the first line to decide whether a header is present
    with open(path) as fh:
        first = fh.readline().strip()
    if first.startswith('barcode'):
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(path, header=None,
                         names=['barcode', 'in_tissue', 'array_row', 'array_col',
                                'pxl_row_in_fullres', 'pxl_col_in_fullres'])
    df = df.set_index('barcode')
    n_in = int((df['in_tissue'] == 1).sum()) if 'in_tissue' in df.columns else len(df)
    print(f'[ST-EM] Tissue positions: {len(df)} spots total, {n_in} in-tissue', flush=True)
    return df


# ---------------------------------------------------------------------------
# EM primitives (same E/M logic as bulk miniQuant)
# ---------------------------------------------------------------------------

def _e_step(cond_prob, theta):
    """
    E-step: posterior assignment probability for each molecule.

    Parameters
    ----------
    cond_prob : ndarray (n_mol × n_iso)
    theta     : ndarray (n_iso,)

    Returns
    -------
    q : ndarray (n_mol × n_iso), rows sum to 1
    """
    weighted = cond_prob * theta
    row_sum = weighted.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0] = 1.0
    return weighted / row_sum


def _m_step(q):
    """M-step: new theta = normalised column sums of q."""
    counts = q.sum(axis=0)
    s = counts.sum()
    return counts / s if s > 0 else np.ones(len(counts)) / len(counts)


def _run_em(cond_prob, theta_init=None, max_iter=200, min_diff=1e-6):
    """
    Run EM on a molecule × isoform cond_prob array.

    Parameters
    ----------
    cond_prob  : ndarray (n_mol × n_iso)
    theta_init : ndarray (n_iso,) or None — None → uniform initialisation
    max_iter   : int
    min_diff   : float — convergence threshold

    Returns
    -------
    theta : ndarray (n_iso,)
    """
    n_mol, n_iso = cond_prob.shape
    if n_mol == 0 or n_iso == 0:
        return np.ones(n_iso) / n_iso if n_iso > 0 else np.array([])

    if theta_init is not None:
        theta = np.array(theta_init, dtype=np.float64)
        s = theta.sum()
        theta = theta / s if s > 0 else np.ones(n_iso) / n_iso
    else:
        theta = np.ones(n_iso, dtype=np.float64) / n_iso

    for _ in range(max_iter):
        q = _e_step(cond_prob, theta)
        new_theta = _m_step(q)
        if np.max(np.abs(new_theta - theta)) < min_diff:
            theta = new_theta
            break
        theta = new_theta

    return theta


# ---------------------------------------------------------------------------
# UMI aggregation
# ---------------------------------------------------------------------------

def _build_mol_matrix(umi_to_rows, cond_prob_dense):
    """
    Aggregate reads sharing the same (barcode, UMI) into one molecule row.

    Strategy: element-wise maximum cond_prob across reads in the group.

    Parameters
    ----------
    umi_to_rows    : dict {umi_key: [local_row_indices]}
    cond_prob_dense: ndarray (n_community_reads × n_community_iso)

    Returns
    -------
    mol_matrix : ndarray (n_molecules × n_community_iso)
    """
    rows = []
    for local_rows in umi_to_rows.values():
        if len(local_rows) == 1:
            rows.append(cond_prob_dense[local_rows[0]])
        else:
            rows.append(cond_prob_dense[np.array(local_rows)].max(axis=0))
    return (np.array(rows, dtype=np.float64) if rows
            else np.zeros((0, cond_prob_dense.shape[1]), dtype=np.float64))


# ---------------------------------------------------------------------------
# CSR memmap helpers — allow workers to share cond_prob without copying
# ---------------------------------------------------------------------------

def _save_csr_memmap(csr, path_prefix):
    """Save CSR matrix components as .npy files for mmap-based sharing."""
    np.save(f'{path_prefix}_data.npy',    csr.data)
    np.save(f'{path_prefix}_indices.npy', csr.indices)
    np.save(f'{path_prefix}_indptr.npy',  csr.indptr)
    np.save(f'{path_prefix}_shape.npy',   np.array(csr.shape))


def _load_csr_memmap(path_prefix):
    """Load CSR matrix from mmap .npy files (read-only, no extra copy)."""
    data    = np.load(f'{path_prefix}_data.npy',    mmap_mode='r')
    indices = np.load(f'{path_prefix}_indices.npy', mmap_mode='r')
    indptr  = np.load(f'{path_prefix}_indptr.npy',  mmap_mode='r')
    shape   = tuple(np.load(f'{path_prefix}_shape.npy'))
    return scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)


# ---------------------------------------------------------------------------
# Initialisation helpers
# ---------------------------------------------------------------------------

def _unique_mapping_init(cond_prob):
    """
    Initialise theta from uniquely-mapping molecules.

    A molecule is unique-mapping if it has a non-zero probability for
    exactly one isoform. Returns None if no unique-mapping molecules exist
    (caller falls back to uniform initialisation).
    """
    unique_mask = (cond_prob > 0).sum(axis=1) == 1
    if unique_mask.sum() == 0:
        return None
    theta = cond_prob[unique_mask].sum(axis=0)
    s = theta.sum()
    return theta / s if s > 0 else None


# ---------------------------------------------------------------------------
# Per-community processing (runs inside each worker)
# ---------------------------------------------------------------------------

def _process_community(community_id, lr_rows, iso_cols, sub_cp,
                        row_barcode, row_umi, col_to_isoform,
                        max_iter):
    """
    Process one community (gene): run per-unit EM, hard-assign, return UMI counts.

    Returns
    -------
    dict  {barcode: {isoform_id: count}}
    """
    # Group reads by (barcode, umi)
    barcode_umi = defaultdict(lambda: defaultdict(list))

    for local_row, g_row in enumerate(lr_rows):
        barcode = row_barcode.get(g_row)
        umi_key = row_umi.get(g_row, g_row)
        if barcode is not None:
            barcode_umi[barcode][umi_key].append(local_row)

    if not barcode_umi:
        return {}

    counts = defaultdict(lambda: defaultdict(int))

    for barcode, umi_dict in barcode_umi.items():
        cell_mol_cp = _build_mol_matrix(umi_dict, sub_cp)
        n_mols = len(cell_mol_cp)
        if n_mols == 0:
            continue
        theta_init = _unique_mapping_init(cell_mol_cp)
        theta_cell = _run_em(cell_mol_cp,
                             theta_init=theta_init,
                             max_iter=max_iter)

        # Largest remainder method: distribute n_mols counts proportionally
        # to theta, rounding to integers while preserving the total.
        fractional   = theta_cell * n_mols
        floor_counts = np.floor(fractional).astype(int)
        deficit      = n_mols - floor_counts.sum()
        if deficit > 0:
            remainders = fractional - floor_counts
            top_indices = np.argsort(remainders, kind='stable')[::-1][:deficit]
            floor_counts[top_indices] += 1

        for local_idx, cnt in enumerate(floor_counts):
            if cnt > 0:
                isoform_id = col_to_isoform[iso_cols[local_idx]]
                counts[barcode][isoform_id] += cnt

    return {bc: dict(iso_cnt) for bc, iso_cnt in counts.items()}


# ---------------------------------------------------------------------------
# Worker function (one process per worker)
# ---------------------------------------------------------------------------

def _sc_worker(args):
    """Process a subset of communities in a worker process."""
    os.nice(10)
    (worker_id, community_ids, cond_prob_prefix, temp_dir,
     LR_labels, isoform_labels, col_to_isoform, max_iter) = args

    cond_prob_global = _load_csr_memmap(cond_prob_prefix)

    with open(os.path.join(temp_dir, 'sc_row_barcode.pkl'), 'rb') as f:
        row_barcode = pickle.load(f)
    with open(os.path.join(temp_dir, 'sc_row_umi.pkl'), 'rb') as f:
        row_umi = pickle.load(f)

    cell_counts = defaultdict(lambda: defaultdict(int))

    for community_id in community_ids:
        lr_rows  = np.where(LR_labels      == community_id)[0]
        iso_cols = np.where(isoform_labels == community_id)[0]

        if len(lr_rows) == 0 or len(iso_cols) == 0:
            continue

        sub_cp = np.asarray(
            cond_prob_global[lr_rows, :][:, iso_cols].todense(),
            dtype=np.float64
        )

        community_counts = _process_community(
            community_id, lr_rows, iso_cols, sub_cp,
            row_barcode, row_umi, col_to_isoform,
            max_iter
        )

        for barcode, iso_cnt in community_counts.items():
            for iso, cnt in iso_cnt.items():
                cell_counts[barcode][iso] += cnt

    return {bc: dict(iso_cnt) for bc, iso_cnt in cell_counts.items()}


def _callback_error(e):
    print(f'[SC-EM] Worker error: {e}', flush=True)


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _write_outputs(output_path, cell_isoform_counts, isoform_gene_dict,
                   isoform_len_dict, quant_mode,
                   output_dir_name='SC_output', barcode_label='cell_barcode',
                   tissue_positions=None):
    """Build sparse count matrix and write all output files (SC or ST mode)."""
    all_barcodes = sorted(cell_isoform_counts.keys())
    all_isoforms = sorted(isoform_gene_dict.keys())
    barcode_idx  = {bc: i for i, bc in enumerate(all_barcodes)}
    isoform_idx  = {iso: i for i, iso in enumerate(all_isoforms)}
    n_cells    = len(all_barcodes)
    n_isoforms = len(all_isoforms)

    rows_list, cols_list, data_list = [], [], []
    for bc, iso_counts in cell_isoform_counts.items():
        bc_i = barcode_idx[bc]
        for iso, cnt in iso_counts.items():
            if iso in isoform_idx and cnt > 0:
                rows_list.append(bc_i)
                cols_list.append(isoform_idx[iso])
                data_list.append(int(cnt))

    count_matrix = scipy.sparse.coo_matrix(
        (data_list, (rows_list, cols_list)), shape=(n_cells, n_isoforms)
    )

    # Output directory structure:  <output_dir_name>/isoform/  and  <output_dir_name>/gene/
    sc_out_dir   = os.path.join(output_path, output_dir_name)
    iso_dir      = os.path.join(sc_out_dir, 'isoform')
    gene_dir     = os.path.join(sc_out_dir, 'gene')
    Path(iso_dir).mkdir(parents=True, exist_ok=True)
    Path(gene_dir).mkdir(parents=True, exist_ok=True)

    barcodes_content = '\n'.join(all_barcodes) + '\n'

    # ------------------------------------------------------------------
    # isoform/ : barcodes.tsv, features.tsv, matrix.mtx, cpm_matrix.mtx
    # ------------------------------------------------------------------
    with open(os.path.join(iso_dir, 'barcodes.tsv'), 'w') as f:
        f.write(barcodes_content)

    with open(os.path.join(iso_dir, 'features.tsv'), 'w') as f:
        for iso in all_isoforms:
            gene = isoform_gene_dict.get(iso, 'NA')
            f.write(f'{iso}\t{gene}\tGene Expression\n')

    mtx_csc = count_matrix.T.tocsc()
    with open(os.path.join(iso_dir, 'matrix.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate integer general\n')
        f.write('%metadata_json: {"software_version": "miniQuant-SC"}\n')
        f.write(f'{n_isoforms} {n_cells} {mtx_csc.nnz}\n')
        cx = mtx_csc.tocoo()
        for r, c, v in zip(cx.row, cx.col, cx.data):
            f.write(f'{r + 1} {c + 1} {int(v)}\n')

    # Isoform CPM matrix
    count_csr   = count_matrix.tocsr()
    cell_totals = np.asarray(count_csr.sum(axis=1)).flatten()
    cell_totals[cell_totals == 0] = 1.0
    cpm_rows, cpm_cols, cpm_vals = [], [], []
    for bc_i in range(n_cells):
        row = count_csr[bc_i]
        if row.nnz == 0:
            continue
        iso_indices = row.indices
        iso_counts  = np.array(row.data, dtype=np.float64)
        cpm_values  = iso_counts / cell_totals[bc_i] * 1e6
        for iso_i, cpm in zip(iso_indices, cpm_values):
            if cpm > 0:
                cpm_rows.append(int(iso_i))
                cpm_cols.append(bc_i)
                cpm_vals.append(float(cpm))

    cpm_matrix_coo = scipy.sparse.coo_matrix(
        (cpm_vals, (cpm_rows, cpm_cols)), shape=(n_isoforms, n_cells)
    )
    with open(os.path.join(iso_dir, 'cpm_matrix.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write('%metadata_json: {"software_version": "miniQuant-SC"}\n')
        f.write(f'{n_isoforms} {n_cells} {cpm_matrix_coo.nnz}\n')
        for r, c, v in zip(cpm_matrix_coo.row,
                           cpm_matrix_coo.col,
                           cpm_matrix_coo.data):
            f.write(f'{r + 1} {c + 1} {v:.4f}\n')

    # ------------------------------------------------------------------
    # gene/ : barcodes.tsv, features.tsv, matrix.mtx, cpm_matrix.mtx
    # ------------------------------------------------------------------
    all_genes    = sorted(set(isoform_gene_dict.values()))
    gene_idx_map = {gene: i for i, gene in enumerate(all_genes)}
    n_genes      = len(all_genes)

    agg_rows = [i for i, iso in enumerate(all_isoforms)
                if isoform_gene_dict.get(iso, 'NA') in gene_idx_map]
    agg_cols = [gene_idx_map[isoform_gene_dict[all_isoforms[i]]] for i in agg_rows]
    iso_to_gene = scipy.sparse.coo_matrix(
        (np.ones(len(agg_rows)), (agg_rows, agg_cols)),
        shape=(n_isoforms, n_genes)
    ).tocsr()

    gene_count_matrix_T = (count_matrix @ iso_to_gene).T.tocoo()

    with open(os.path.join(gene_dir, 'barcodes.tsv'), 'w') as f:
        f.write(barcodes_content)

    with open(os.path.join(gene_dir, 'features.tsv'), 'w') as f:
        for gene in all_genes:
            f.write(f'{gene}\t{gene}\tGene Expression\n')

    with open(os.path.join(gene_dir, 'matrix.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate integer general\n')
        f.write('%metadata_json: {"software_version": "miniQuant-SC"}\n')
        f.write(f'{n_genes} {n_cells} {gene_count_matrix_T.nnz}\n')
        for r, c, v in zip(gene_count_matrix_T.row,
                           gene_count_matrix_T.col,
                           gene_count_matrix_T.data):
            f.write(f'{r + 1} {c + 1} {int(v)}\n')

    gene_count_csr = (count_matrix @ iso_to_gene).tocsr()
    gene_cpm_coo = (gene_count_csr
                    .multiply(1.0 / cell_totals[:, np.newaxis])
                    .multiply(1e6)).T.tocoo()

    with open(os.path.join(gene_dir, 'cpm_matrix.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write('%metadata_json: {"software_version": "miniQuant-SC"}\n')
        f.write(f'{n_genes} {n_cells} {gene_cpm_coo.nnz}\n')
        for r, c, v in zip(gene_cpm_coo.row,
                           gene_cpm_coo.col,
                           gene_cpm_coo.data):
            f.write(f'{r + 1} {c + 1} {v:.4f}\n')

    # ------------------------------------------------------------------
    # README.txt (one per subdirectory)
    # ------------------------------------------------------------------
    n_nonzero_cells = int((count_matrix.sum(axis=1) > 0).sum())
    n_nonzero_iso   = int((count_matrix.sum(axis=0) > 0).sum())
    total_counts    = int(count_matrix.sum())

    is_spatial   = (barcode_label == 'spot_barcode')
    tag          = 'ST-EM' if is_spatial else 'SC-EM'
    unit         = 'spot'  if is_spatial else 'cell'
    unit_cap     = unit.capitalize()
    software_ver = 'miniQuant-ST' if is_spatial else 'miniQuant-SC'

    common_header = [
        f'Mode                : {"spatial transcriptomics" if is_spatial else "single-cell"}',
        f'Quantification mode : {quant_mode}',
        f'{unit_cap}s{" " * max(0, 16 - len(unit))}: {n_cells}',
        f'Total UMI counts    : {total_counts}',
    ]
    mtx_format_lines = [
        'File format (MEX — Market Exchange Format)',
        '------------------------------------------',
        'barcodes.tsv',
        f'  One {unit} barcode per line, no header.',
        f'  Line number (1-based) = {unit} index = column index in matrix.mtx.',
        '',
        'features.tsv',
        '  Tab-separated, three columns, no header: feature_id  feature_name  feature_type',
        '  Line number (1-based) = feature index = row index in matrix.mtx.',
        '',
        'matrix.mtx',
        '  Sparse coordinate format (integer UMI counts).',
        '  Line 1 : %%MatrixMarket matrix coordinate integer general',
        '  Line 2 : comment (%...)',
        '  Line 3 : n_features  n_barcodes  n_nonzero',
        '  Line 4+: row_index  col_index  value  (1-based, non-zero entries only)',
        '',
        'cpm_matrix.mtx',
        '  Sparse coordinate format (real CPM values).',
        f'  CPM normalised within each {unit} independently (UMI counts / total {unit} UMIs × 10^6).',
        '  Same structure as matrix.mtx but with "real" type and 4 decimal places.',
        '',
        f'  Compatible with Seurat::ReadMtx() and scanpy.read_mtx().',
    ]

    iso_readme = [
        f'miniQuant {output_dir_name}/isoform — isoform-level quantification',
        '=' * 55,
        *common_header,
        f'Isoforms            : {n_isoforms}',
        f'Non-zero {unit}s     : {n_nonzero_cells}',
        f'Non-zero isoforms   : {n_nonzero_iso}',
        '',
        'features.tsv columns: isoform_id  gene_id  Gene Expression',
        '',
        *mtx_format_lines,
        '',
        'README.txt',
        '  This file.',
    ]
    with open(os.path.join(iso_dir, 'README.txt'), 'w') as f:
        f.write('\n'.join(iso_readme) + '\n')

    gene_readme = [
        f'miniQuant {output_dir_name}/gene — gene-level quantification',
        '=' * 55,
        *common_header,
        f'Genes               : {n_genes}',
        '',
        'features.tsv columns: gene_id  gene_id  Gene Expression',
        '  (gene_id appears twice to match the standard MEX features.tsv convention)',
        '  matrix.mtx value = sum of UMI counts across all isoforms of that gene.',
        '',
        *mtx_format_lines,
        '',
        'README.txt',
        '  This file.',
    ]
    with open(os.path.join(gene_dir, 'README.txt'), 'w') as f:
        f.write('\n'.join(gene_readme) + '\n')

    print(f'[{tag}] Results written to: {sc_out_dir}', flush=True)
    print(f'[{tag}]   isoform/matrix.mtx     — isoform UMI counts  ({n_isoforms} × {n_cells})',
          flush=True)
    print(f'[{tag}]   isoform/cpm_matrix.mtx — isoform CPM         ({n_isoforms} × {n_cells})',
          flush=True)
    print(f'[{tag}]   gene/matrix.mtx        — gene UMI counts     ({n_genes} × {n_cells})',
          flush=True)
    print(f'[{tag}]   gene/cpm_matrix.mtx    — gene CPM            ({n_genes} × {n_cells})',
          flush=True)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_sc_em_quantification(output_path, isoform_gene_dict, isoform_len_dict,
                              alpha=10.0, max_iter=200, umi_dedup_hamming=0,
                              threads=1,
                              barcode_whitelist_path=None, tissue_positions_path=None,
                              output_dir_name='SC_output', barcode_label='cell_barcode'):
    """
    Single-cell / spatial-transcriptomics isoform quantification.

    Parameters
    ----------
    output_path           : str   — run output directory
    isoform_gene_dict     : dict  {isoform_id: gene_id}
    isoform_len_dict      : dict  {isoform_id: length_in_bp}
    alpha                 : float — reserved
    max_iter              : int   — max EM iterations per community
    umi_dedup_hamming     : int   — Hamming distance for UMI dedup (0 = exact match)
    threads               : int   — parallel worker processes
    barcode_whitelist_path: str   — path to whitelist file (one barcode per line); ST mode
    tissue_positions_path : str   — path to tissue_positions.csv (Visium format); ST mode
    output_dir_name       : str   — output subdirectory name ('SC_output' or 'ST_output')
    barcode_label         : str   — index column label ('cell_barcode' or 'spot_barcode')
    """
    is_spatial = (barcode_label == 'spot_barcode')
    tag  = 'ST-EM' if is_spatial else 'SC-EM'
    unit = 'spot'  if is_spatial else 'cell'

    mode_label = f'per-{unit} EM'
    data_type  = 'spatial transcriptomics' if is_spatial else 'single-cell'
    print(f'[{tag}] Starting {data_type} quantification '
          f'({mode_label}, threads={threads})...', flush=True)
    _start_time = time.time()

    barcode_whitelist = _load_barcode_whitelist(barcode_whitelist_path)
    tissue_positions  = _load_tissue_positions(tissue_positions_path)

    lr_align_dir  = os.path.join(output_path, 'temp', 'LR_alignments')
    cond_prob_dir = os.path.join(output_path, 'temp', 'cond_prob')
    community_dir = os.path.join(output_path, 'temp', 'community')
    temp_dir      = os.path.join(output_path, 'temp')

    # ------------------------------------------------------------------
    # 1. Load global cond_prob and save as memmap for worker sharing
    # ------------------------------------------------------------------
    cond_prob_files = sorted(glob.glob(os.path.join(cond_prob_dir, '*_cond_prob.npz')))
    if not cond_prob_files:
        print(f'[{tag}] No cond_prob files found. Aborting.', flush=True)
        return

    worker_ids = sorted(
        int(re.match(r'(\d+)_cond_prob', os.path.basename(f)).group(1))
        for f in cond_prob_files
    )
    print(f'[{tag}] Loading cond_prob from {len(worker_ids)} worker(s)...', flush=True)
    cond_prob_global = scipy.sparse.vstack(
        [scipy.sparse.load_npz(os.path.join(cond_prob_dir, f'{wid}_cond_prob.npz'))
         for wid in worker_ids]
    ).tocsr()
    n_reads_global, n_isoforms_global = cond_prob_global.shape
    print(f'[{tag}] cond_prob: {n_reads_global} reads × {n_isoforms_global} isoforms',
          flush=True)

    # Save cond_prob as memmap so workers can load without copying
    cond_prob_prefix = os.path.join(temp_dir, 'sc_cond_prob')
    _save_csr_memmap(cond_prob_global, cond_prob_prefix)

    # ------------------------------------------------------------------
    # 2. Isoform list — must match column order used in bulk EM
    # ------------------------------------------------------------------
    isoform_list   = sorted(isoform_len_dict.keys())
    col_to_isoform = {i: iso for i, iso in enumerate(isoform_list)}

    # ------------------------------------------------------------------
    # 3. Build global_row → (barcode, umi) mapping and save to disk
    # ------------------------------------------------------------------
    cb_umi_map = {}
    for fpath in sorted(glob.glob(os.path.join(lr_align_dir, 'cb_umi_*'))):
        with open(fpath, 'rb') as f:
            cb_umi_map.update(pickle.load(f))

    if not cb_umi_map:
        flag = '--st_mode' if is_spatial else '--sc_mode'
        print(f'[{tag}] No CB/UMI data found. Did you run with {flag}? Aborting.',
              flush=True)
        return

    row_barcode = {}
    row_umi     = {}
    global_row  = 0
    n_whitelist_filtered = 0
    for wid in worker_ids:
        for fpath in sorted(glob.glob(os.path.join(lr_align_dir, f'reads_{wid}_*'))):
            with open(fpath, 'rb') as f:
                reads_isoform_info, _ = pickle.load(f)
            for read_name, isoform_info in reads_isoform_info.items():
                if len(isoform_info) > 0:
                    cb_umi = cb_umi_map.get(read_name)
                    if cb_umi is not None:
                        barcode, umi = cb_umi
                        if barcode is not None:
                            if barcode_whitelist is None or barcode in barcode_whitelist:
                                row_barcode[global_row] = barcode
                                row_umi[global_row] = (umi if umi is not None
                                                       else read_name)
                            else:
                                n_whitelist_filtered += 1
                    global_row += 1

    print(f'[{tag}] Reads with valid barcode: {len(row_barcode)} / {global_row}', flush=True)
    if barcode_whitelist is not None:
        print(f'[{tag}] Reads filtered (barcode not in whitelist): {n_whitelist_filtered}',
              flush=True)

    # Save to disk for workers to load
    with open(os.path.join(temp_dir, 'sc_row_barcode.pkl'), 'wb') as f:
        pickle.dump(row_barcode, f)
    with open(os.path.join(temp_dir, 'sc_row_umi.pkl'), 'wb') as f:
        pickle.dump(row_umi, f)

    # ------------------------------------------------------------------
    # 4. Load community labels
    # ------------------------------------------------------------------
    lr_labels_path  = os.path.join(community_dir, 'LR_labels.npy')
    iso_labels_path = os.path.join(community_dir, 'isoform_labels.npy')
    if not os.path.exists(lr_labels_path):
        print(f'[{tag}] LR_labels.npy not found. Aborting.', flush=True)
        return

    LR_labels      = np.load(lr_labels_path)
    isoform_labels = np.load(iso_labels_path)
    unique_communities = np.unique(LR_labels)
    print(f'[{tag}] Communities (genes): {len(unique_communities)}', flush=True)

    # ------------------------------------------------------------------
    # 5. Distribute communities across workers and run EM
    # ------------------------------------------------------------------
    community_chunks = np.array_split(unique_communities, threads)

    worker_args = [
        (wid, chunk.tolist(), cond_prob_prefix, temp_dir,
         LR_labels, isoform_labels, col_to_isoform, max_iter)
        for wid, chunk in enumerate(community_chunks)
        if len(chunk) > 0
    ]

    if threads == 1:
        all_results = [_sc_worker(worker_args[0])]
    else:
        pool = mp.Pool(threads)
        futures = [
            pool.apply_async(_sc_worker, (args,), error_callback=_callback_error)
            for args in worker_args
        ]
        all_results = [f.get() for f in futures]
        pool.close()
        pool.join()

    # ------------------------------------------------------------------
    # 6. Merge results from all workers
    # ------------------------------------------------------------------
    cell_isoform_counts = defaultdict(lambda: defaultdict(int))
    for partial in all_results:
        for barcode, iso_cnt in partial.items():
            for iso, cnt in iso_cnt.items():
                cell_isoform_counts[barcode][iso] += cnt

    print(f'[{tag}] EM complete. {unit.capitalize()}s with counts: {len(cell_isoform_counts)}',
          flush=True)

    # ------------------------------------------------------------------
    # 7. Write output files
    # ------------------------------------------------------------------
    _write_outputs(output_path, cell_isoform_counts, isoform_gene_dict,
                   isoform_len_dict, mode_label,
                   output_dir_name=output_dir_name,
                   barcode_label=barcode_label,
                   tissue_positions=tissue_positions)

    elapsed = time.time() - _start_time
    print(f'[{tag}] Done in {elapsed:.1f} s.', flush=True)
