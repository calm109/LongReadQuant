"""
sp_output.py
------------
Single-cell post-processing step for miniQuant-SC (legacy global-theta path).

After bulk EM gives final isoform proportions (theta), this module:
1. Loads reads_isoform_info and CB/UMI maps saved by parse_alignment_EM
2. For each read: assigns it to isoforms probabilistically using theta
3. Per cell: deduplicates UMIs per isoform (exact match or hamming)
4. Outputs — mirroring the sp_em_quantification structure:

   SC_output/
   ├── isoform/
   │   ├── barcodes.tsv        one cell barcode per line
   │   ├── features.tsv        isoform_id  gene_id  Gene Expression
   │   ├── matrix.mtx          sparse UMI count matrix  (isoforms x cells)
   │   ├── cpm_matrix.mtx      sparse CPM matrix        (isoforms x cells)
   │   └── README.txt
   └── gene/
       ├── barcodes.tsv
       ├── features.tsv        gene_id  gene_id  Gene Expression
       ├── matrix.mtx          gene-level UMI count matrix (genes x cells)
       ├── cpm_matrix.mtx      gene-level CPM matrix
       └── README.txt

Note: this function is the *legacy* path used when the bulk EM result
(Isoform_abundance.out) is the source of theta.  The newer path —
run_sc_em_quantification() in sp_em_quantification.py — runs its own
per-cell / per-spot EM and writes the same directory layout directly.
"""

import glob
import pickle
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse


# ---------------------------------------------------------------------------
# UMI deduplication
# ---------------------------------------------------------------------------

def _hamming(a, b):
    """Return hamming distance between two equal-length strings."""
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def _deduplicate_umis(umi_list, max_hamming=0):
    """
    Cluster UMIs by hamming distance and return the number of unique clusters.
    max_hamming=0  → exact match (fastest)
    max_hamming=1  → standard UMI dedup (1-mismatch allowed)
    """
    if max_hamming == 0:
        return len(set(umi_list))

    clusters = []
    for umi in umi_list:
        merged = False
        for cluster in clusters:
            if any(_hamming(umi, c) <= max_hamming for c in cluster):
                cluster.append(umi)
                merged = True
                break
        if not merged:
            clusters.append([umi])
    return len(clusters)


# ---------------------------------------------------------------------------
# Output writer  (isoform/ and gene/ subdirectories)
# ---------------------------------------------------------------------------

def _write_sc_outputs(output_path, cell_isoform_counts, isoform_gene_dict,
                      isoform_len_dict):
    """
    Build sparse count matrices and write all output files in the
    gene/ + isoform/ subdirectory layout.

    Parameters
    ----------
    output_path        : str   — base output directory (e.g. <run_dir>/SC_output)
    cell_isoform_counts: dict  {cell_barcode: {isoform_id: umi_count}}
    isoform_gene_dict  : dict  {isoform_id: gene_id}
    isoform_len_dict   : dict  {isoform_id: length_in_bp}
    """
    all_barcodes = sorted(cell_isoform_counts.keys())
    all_isoforms = sorted(isoform_gene_dict.keys())
    barcode_idx  = {bc: i for i, bc in enumerate(all_barcodes)}
    isoform_idx  = {iso: i for i, iso in enumerate(all_isoforms)}
    n_cells    = len(all_barcodes)
    n_isoforms = len(all_isoforms)

    # Build isoform-level count matrix  (cells x isoforms)
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

    # Output subdirectory paths
    iso_dir  = os.path.join(output_path, 'isoform')
    gene_dir = os.path.join(output_path, 'gene')
    Path(iso_dir).mkdir(parents=True, exist_ok=True)
    Path(gene_dir).mkdir(parents=True, exist_ok=True)

    barcodes_content = '\n'.join(all_barcodes) + '\n'

    # ------------------------------------------------------------------
    # isoform/
    # ------------------------------------------------------------------
    with open(os.path.join(iso_dir, 'barcodes.tsv'), 'w') as f:
        f.write(barcodes_content)

    with open(os.path.join(iso_dir, 'features.tsv'), 'w') as f:
        for iso in all_isoforms:
            gene = isoform_gene_dict.get(iso, 'NA')
            f.write(f'{iso}\t{gene}\tGene Expression\n')

    # UMI count matrix  (isoforms x cells, MEX convention)
    mtx_csc = count_matrix.T.tocsc()
    with open(os.path.join(iso_dir, 'matrix.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate integer general\n')
        f.write('%metadata_json: {"software_version": "LongReadQuant-SC"}\n')
        f.write(f'{n_isoforms} {n_cells} {mtx_csc.nnz}\n')
        cx = mtx_csc.tocoo()
        for r, c, v in zip(cx.row, cx.col, cx.data):
            f.write(f'{r + 1} {c + 1} {int(v)}\n')

    # CPM matrix  (isoforms x cells)
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
        f.write('%metadata_json: {"software_version": "LongReadQuant-SC"}\n')
        f.write(f'{n_isoforms} {n_cells} {cpm_matrix_coo.nnz}\n')
        for r, c, v in zip(cpm_matrix_coo.row,
                           cpm_matrix_coo.col,
                           cpm_matrix_coo.data):
            f.write(f'{r + 1} {c + 1} {v:.4f}\n')

    # ------------------------------------------------------------------
    # gene/  — aggregate isoform counts to gene level
    # ------------------------------------------------------------------
    all_genes    = sorted(set(isoform_gene_dict.values()))
    gene_idx_map = {gene: i for i, gene in enumerate(all_genes)}
    n_genes      = len(all_genes)

    # Isoform → gene aggregation matrix  (n_isoforms x n_genes)
    agg_rows = [i for i, iso in enumerate(all_isoforms)
                if isoform_gene_dict.get(iso, 'NA') in gene_idx_map]
    agg_cols = [gene_idx_map[isoform_gene_dict[all_isoforms[i]]] for i in agg_rows]
    iso_to_gene = scipy.sparse.coo_matrix(
        (np.ones(len(agg_rows)), (agg_rows, agg_cols)),
        shape=(n_isoforms, n_genes)
    ).tocsr()

    # Gene count matrix  (genes x cells, MEX)
    gene_count_matrix_T = (count_matrix @ iso_to_gene).T.tocoo()

    with open(os.path.join(gene_dir, 'barcodes.tsv'), 'w') as f:
        f.write(barcodes_content)

    with open(os.path.join(gene_dir, 'features.tsv'), 'w') as f:
        for gene in all_genes:
            f.write(f'{gene}\t{gene}\tGene Expression\n')

    with open(os.path.join(gene_dir, 'matrix.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate integer general\n')
        f.write('%metadata_json: {"software_version": "LongReadQuant-SC"}\n')
        f.write(f'{n_genes} {n_cells} {gene_count_matrix_T.nnz}\n')
        for r, c, v in zip(gene_count_matrix_T.row,
                           gene_count_matrix_T.col,
                           gene_count_matrix_T.data):
            f.write(f'{r + 1} {c + 1} {int(v)}\n')

    # Gene CPM matrix  (genes x cells)
    gene_count_csr = (count_matrix @ iso_to_gene).tocsr()
    gene_cpm_coo = (gene_count_csr
                    .multiply(1.0 / cell_totals[:, np.newaxis])
                    .multiply(1e6)).T.tocoo()

    with open(os.path.join(gene_dir, 'cpm_matrix.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write('%metadata_json: {"software_version": "LongReadQuant-SC"}\n')
        f.write(f'{n_genes} {n_cells} {gene_cpm_coo.nnz}\n')
        for r, c, v in zip(gene_cpm_coo.row,
                           gene_cpm_coo.col,
                           gene_cpm_coo.data):
            f.write(f'{r + 1} {c + 1} {v:.4f}\n')

    # ------------------------------------------------------------------
    # README files
    # ------------------------------------------------------------------
    n_nonzero_cells = int((count_matrix.sum(axis=1) > 0).sum())
    n_nonzero_iso   = int((count_matrix.sum(axis=0) > 0).sum())
    total_counts    = int(count_matrix.sum())

    common_header = [
        'Mode                : single-cell (legacy global-theta path)',
        f'Cells               : {n_cells}',
        f'Total UMI counts    : {total_counts}',
    ]
    mtx_format_lines = [
        'File format (MEX — Market Exchange Format)',
        '------------------------------------------',
        'barcodes.tsv',
        '  One cell barcode per line, no header.',
        '  Line number (1-based) = cell index = column index in matrix.mtx.',
        '',
        'features.tsv',
        '  Tab-separated, three columns, no header: feature_id  feature_name  feature_type',
        '  Line number (1-based) = feature index = row index in matrix.mtx.',
        '',
        'matrix.mtx',
        '  Sparse coordinate format (integer UMI counts).',
        '  Line 1 : %%MatrixMarket matrix coordinate integer general',
        '  Line 2 : comment (%...)',
        '  Line 3 : n_features  n_cells  n_nonzero',
        '  Line 4+: row_index  col_index  value  (1-based, non-zero entries only)',
        '',
        'cpm_matrix.mtx',
        '  Sparse coordinate format (real CPM values).',
        '  CPM normalised within each cell independently (UMI counts / total cell UMIs x 10^6).',
        '  Same structure as matrix.mtx but with "real" type and 4 decimal places.',
        '',
        '  Compatible with Seurat::ReadMtx() and scanpy.read_mtx().',
    ]

    iso_readme = [
        'miniQuant SC_output/isoform — isoform-level quantification',
        '=' * 56,
        *common_header,
        f'Isoforms            : {n_isoforms}',
        f'Non-zero cells      : {n_nonzero_cells}',
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
        'miniQuant SC_output/gene — gene-level quantification',
        '=' * 52,
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

    print(f'[SC] Output written to: {output_path}', flush=True)
    print(f'[SC]   isoform/matrix.mtx     — isoform UMI counts  ({n_isoforms} x {n_cells})',
          flush=True)
    print(f'[SC]   isoform/cpm_matrix.mtx — isoform CPM         ({n_isoforms} x {n_cells})',
          flush=True)
    print(f'[SC]   gene/matrix.mtx        — gene UMI counts     ({n_genes} x {n_cells})',
          flush=True)
    print(f'[SC]   gene/cpm_matrix.mtx    — gene CPM            ({n_genes} x {n_cells})',
          flush=True)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def generate_sc_matrix(output_path, isoform_gene_dict, isoform_len_dict,
                       umi_dedup_hamming=0):
    """
    Main entry point called from EM.py after bulk EM completes.

    Parameters
    ----------
    output_path      : str
        The run output directory (same as passed to EM_hybrid).
        Output files are written to <output_path>/SC_output/isoform/
        and <output_path>/SC_output/gene/.
    isoform_gene_dict : dict
        {isoform_id: gene_id}
    isoform_len_dict  : dict
        {isoform_id: length_in_bp}  — used for TPM normalisation.
    umi_dedup_hamming : int
        Hamming distance threshold for UMI deduplication (0 = exact match).
    """
    print('[SC] Starting single-cell isoform matrix generation...', flush=True)

    # ------------------------------------------------------------------
    # 1. Load final theta from Isoform_abundance.out
    # ------------------------------------------------------------------
    abundance_path = os.path.join(output_path, 'Isoform_abundance.out')
    if not os.path.exists(abundance_path):
        print('[SC][ERROR] Isoform_abundance.out not found. Skipping SC output.')
        return

    abundance_df = pd.read_csv(abundance_path, sep='\t')
    theta_dict = dict(zip(abundance_df['Isoform'], abundance_df['num_expected_LRs']))
    total_lr = abundance_df['num_expected_LRs'].sum()
    if total_lr > 0:
        theta_dict = {iso: cnt / total_lr for iso, cnt in theta_dict.items()}
    else:
        print('[SC][WARNING] No expected LRs found in Isoform_abundance.out.')
        return

    all_isoforms = sorted(theta_dict.keys())
    isoform_idx  = {iso: i for i, iso in enumerate(all_isoforms)}

    # ------------------------------------------------------------------
    # 2. Load CB/UMI maps
    # ------------------------------------------------------------------
    lr_align_dir = os.path.join(output_path, 'temp', 'LR_alignments')

    cb_umi_map = {}
    for fpath in glob.glob(os.path.join(lr_align_dir, 'cb_umi_*')):
        with open(fpath, 'rb') as f:
            batch_map = pickle.load(f)
        cb_umi_map.update(batch_map)

    if not cb_umi_map:
        print('[SC][WARNING] No CB/UMI data found. '
              'Did you run with --sc_mode? Skipping SC output.')
        return

    print(f'[SC] Loaded CB/UMI info for {len(cb_umi_map)} reads.', flush=True)

    # ------------------------------------------------------------------
    # 3. Assign reads to isoforms using theta, record (cell, isoform, UMI)
    # ------------------------------------------------------------------
    cell_isoform_umis = defaultdict(lambda: defaultdict(set))

    num_reads_processed = 0
    num_reads_no_cb     = 0
    num_reads_no_iso    = 0

    for fpath in sorted(glob.glob(os.path.join(lr_align_dir, 'reads_*'))):
        with open(fpath, 'rb') as f:
            reads_isoform_info, _ = pickle.load(f)

        for read_name, isoform_info in reads_isoform_info.items():
            if not isoform_info:
                continue

            cb_umi = cb_umi_map.get(read_name)
            if cb_umi is None:
                num_reads_no_cb += 1
                continue
            cell_barcode, umi = cb_umi
            if cell_barcode is None:
                num_reads_no_cb += 1
                continue

            compatible = [iso for iso in isoform_info if iso in theta_dict]
            if not compatible:
                num_reads_no_iso += 1
                continue

            weights    = np.array([theta_dict[iso] for iso in compatible])
            weight_sum = weights.sum()
            if weight_sum == 0:
                weights    = np.ones(len(compatible))
                weight_sum = weights.sum()

            best_isoform = compatible[int(np.argmax(weights))]
            iso_i        = isoform_idx[best_isoform]

            umi_key = umi if umi is not None else read_name
            cell_isoform_umis[cell_barcode][iso_i].add(umi_key)
            num_reads_processed += 1

    print(f'[SC] Reads processed              : {num_reads_processed}', flush=True)
    print(f'[SC] Reads skipped (no CB)        : {num_reads_no_cb}',     flush=True)
    print(f'[SC] Reads skipped (no iso match) : {num_reads_no_iso}',    flush=True)

    # ------------------------------------------------------------------
    # 4. UMI deduplication → cell x isoform count dict
    # ------------------------------------------------------------------
    cell_isoform_counts = defaultdict(dict)
    for barcode, iso_umis in cell_isoform_umis.items():
        for iso_i, umi_set in iso_umis.items():
            count = _deduplicate_umis(list(umi_set), umi_dedup_hamming)
            if count > 0:
                iso_id = all_isoforms[iso_i]
                cell_isoform_counts[barcode][iso_id] = count

    n_cells = len(cell_isoform_counts)
    print(f'[SC] Cells with counts: {n_cells}', flush=True)

    # ------------------------------------------------------------------
    # 5. Write outputs (isoform/ and gene/ subdirectories)
    # ------------------------------------------------------------------
    sc_out_dir = os.path.join(output_path, 'SC_output')
    _write_sc_outputs(
        output_path=sc_out_dir,
        cell_isoform_counts=cell_isoform_counts,
        isoform_gene_dict=isoform_gene_dict,
        isoform_len_dict=isoform_len_dict,
    )

    n_nonzero_cells = sum(1 for v in cell_isoform_counts.values() if v)
    print(f'[SC] Matrix generation complete. '
          f'Cells with ≥1 count: {n_nonzero_cells}', flush=True)
