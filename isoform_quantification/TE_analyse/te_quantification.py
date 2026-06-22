"""
te_quantification.py
--------------------
Per-spot / per-cell TE quantification for miniQuant spatial / single-cell data.

After running  cal_TE --skip-quantification  (which produces
transcript_quantification_with_TE.tsv for every transcript in the GTF),
this module merges the TE annotation with the spot x isoform UMI count
matrix from a previous  quantify --st_mode  (or  --sc_mode) run to
compute per-spot TE metrics.

Entry point
-----------
    run_sc_te_analysis(te_table_path, quant_dir, output_dir,
                       percent_threshold=0.5, output_loci=False)

Inputs
------
  te_table_path : str
      Path to  transcript_quantification_with_TE.tsv  produced by cal_TE.
  quant_dir     : str
      Directory that contains the isoform-level MEX output from a previous
      miniQuant spatial / SC run.  The function searches for barcodes.tsv,
      features.tsv, and matrix.mtx in the following sub-paths (in order):
        <quant_dir>/ST_output/isoform/
        <quant_dir>/SC_output/isoform/
        <quant_dir>/isoform/
        <quant_dir>/            (direct path)
  output_dir    : str
      Directory where per-spot TE output files are written.
  percent_threshold : float
      Minimum TE-overlap proportion (overlap_length / TE_length) required
      to count a TE as "associated" with an isoform.  Default 0.5 (50 %).
  output_loci   : bool
      If True, also write spot_te_loci_counts.tsv (spot x TE-locus matrix).
      Can be very large; disabled by default.

Outputs (all written to <output_dir>/)
---------------------------------------
  spot_transcript_classification_counts.tsv
      Per-spot UMI count breakdown by transcript classification
      (Gene-alone, TE-Gene, TE-alone and sub-types).

  spot_te_class_counts.tsv
      Per-spot UMI counts aggregated by TE class (LINE, SINE, LTR, …).

  spot_te_family_counts.tsv
      Per-spot UMI counts aggregated by TE family.

  spot_te_subfamily_counts.tsv
      Per-spot UMI counts aggregated by TE subfamily.

  spot_te_loci_counts.tsv   [only if output_loci=True]
      Per-spot UMI counts aggregated by individual TE locus.

  README.txt
      Description of output files.

Aggregation logic
-----------------
  For each isoform in the count matrix:
    - Parse the semicolon-separated (locus / subfamily / family / class /
      length / overlap_length) columns from the TE annotation table.
    - For each TE entry:  prop = overlap_length / TE_length
    - If prop >= percent_threshold, the isoform is associated with that TE.
  Per spot:
    - spot_te_class[spot, class] += UMI_count for isoforms associated with class.
    (Each isoform may contribute its UMI count to multiple TE-class columns if
     it overlaps TEs from different classes — "attribute-to-all" policy.)
"""

import gzip
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse
import scipy.io


# ---------------------------------------------------------------------------
# MTX / MEX loading helpers
# ---------------------------------------------------------------------------

def _resolve_mex_file(directory: str, base_name: str):
    """
    Resolve a MEX component file to its actual on-disk path, accepting either
    the plain filename or a gzip-compressed '<base_name>.gz' variant (the
    latter is the default CellRanger/10x output layout).

    Returns the full path if either form exists, else None.
    """
    plain_path = os.path.join(directory, base_name)
    gz_path = plain_path + ".gz"
    if os.path.exists(plain_path):
        return plain_path
    if os.path.exists(gz_path):
        return gz_path
    return None


def _open_text(path: str):
    """Open a text file transparently, decompressing it if it ends in .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _find_isoform_mex_dir(quant_dir: str) -> str:
    """
    Locate the isoform MEX directory inside quant_dir.
    Tries several common sub-paths produced by miniQuant spatial/SC mode.
    Each of barcodes.tsv / features.tsv / matrix.mtx may be plain text or
    gzip-compressed (e.g. barcodes.tsv.gz).
    """
    candidates = [
        os.path.join(quant_dir, "ST_output", "isoform"),
        os.path.join(quant_dir, "SC_output", "isoform"),
        os.path.join(quant_dir, "isoform"),
        quant_dir,
    ]
    for d in candidates:
        if (_resolve_mex_file(d, "barcodes.tsv") is not None
                and _resolve_mex_file(d, "features.tsv") is not None
                and _resolve_mex_file(d, "matrix.mtx") is not None):
            print(f"[TE-SC] Found MEX isoform dir: {d}", flush=True)
            return d
    raise FileNotFoundError(
        f"Cannot find MEX isoform directory under '{quant_dir}'.\n"
        "Expected barcodes.tsv(.gz), features.tsv(.gz), matrix.mtx(.gz) in one of:\n"
        + "\n".join(f"  {d}" for d in candidates)
    )


def load_count_matrix(quant_dir: str):
    """
    Load isoform UMI count matrix from MEX format.
    barcodes.tsv, features.tsv and matrix.mtx may each be plain text or
    gzip-compressed (e.g. barcodes.tsv.gz — the default CellRanger layout).

    Returns
    -------
    count_csr : scipy.sparse.csr_matrix  (n_spots x n_isoforms)
    barcodes  : list[str]
    isoforms  : list[str]
    """
    mex_dir = _find_isoform_mex_dir(quant_dir)
    barcodes_path = _resolve_mex_file(mex_dir, "barcodes.tsv")
    features_path = _resolve_mex_file(mex_dir, "features.tsv")
    matrix_path = _resolve_mex_file(mex_dir, "matrix.mtx")

    # --- barcodes ---
    barcodes = []
    with _open_text(barcodes_path) as f:
        for line in f:
            bc = line.strip()
            if bc:
                barcodes.append(bc)

    # --- features (isoform IDs are the first column) ---
    isoforms = []
    with _open_text(features_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                isoforms.append(parts[0])

    # --- matrix.mtx ---
    # MEX convention: rows = features (isoforms), cols = barcodes (spots), 1-indexed
    # scipy.io.mmread honours the %%MatrixMarket header and returns a coo_matrix
    # in the same orientation (features x barcodes). It also auto-detects gzip
    # compression from the '.gz' suffix, so matrix_path is passed through as-is.
    # We transpose to spots x isoforms.
    n_iso = len(isoforms)
    n_spots = len(barcodes)

    try:
        mat_raw = scipy.io.mmread(matrix_path)
        # mat_raw shape: (n_iso x n_spots) — transpose to (n_spots x n_iso)
        count_csr = mat_raw.T.tocsr()
        # Sanity check: sizes must match
        if count_csr.shape != (n_spots, n_iso):
            raise ValueError(
                f"Matrix shape {count_csr.shape} does not match "
                f"({n_spots} spots x {n_iso} isoforms)"
            )
    except Exception as e:
        raise RuntimeError(f"Failed to load matrix.mtx: {e}") from e

    print(f"[TE-SC] Count matrix: {n_spots} spots x {n_iso} isoforms, "
          f"{count_csr.nnz} non-zero entries, total counts = {int(count_csr.sum())}",
          flush=True)
    return count_csr, barcodes, isoforms


# ---------------------------------------------------------------------------
# TE annotation helpers
# ---------------------------------------------------------------------------

def _parse_te_entries(row, percent_threshold: float):
    """
    Parse semicolon-separated TE fields from a transcript row and return
    a list of entries whose TE-overlap proportion >= percent_threshold.

    Each returned entry is a dict:
      {locus, subfamily, family, class_name, prop}
    """
    entries = []

    def _split(val):
        s = str(val) if not pd.isna(val) else ""
        return s.split(";") if s else []

    loci       = _split(row.get("locus", ""))
    lengths    = _split(row.get("length", ""))
    overlaps   = _split(row.get("overlap_length", ""))
    subfams    = _split(row.get("subfamily", ""))
    families   = _split(row.get("family", ""))
    classes    = _split(row.get("class", ""))

    for i, locus in enumerate(loci):
        if not locus or locus == "nan":
            continue
        try:
            length  = float(lengths[i])  if i < len(lengths)  else 0.0
            overlap = float(overlaps[i]) if i < len(overlaps) else 0.0
            prop    = overlap / length if length > 0 else 0.0
        except (ValueError, ZeroDivisionError):
            prop = 0.0

        if prop < percent_threshold:
            continue

        entries.append({
            "locus":      locus,
            "subfamily":  subfams[i]  if i < len(subfams)  else "Unknown",
            "family":     families[i] if i < len(families) else "Unknown",
            "class_name": classes[i]  if i < len(classes)  else "Unknown",
            "prop":       prop,
        })
    return entries


def _strip_version(transcript_id) -> str:
    """Strip a trailing Ensembl-style version suffix (e.g. '.2') from an ID."""
    return re.sub(r"\.\d+$", "", str(transcript_id))


def build_te_id_resolver(isoforms, te_df):
    """
    Resolve each count-matrix isoform ID to the matching te_df index label.

    Tries an exact match first; if that fails, falls back to comparing
    version-stripped IDs (e.g. 'ENST00000448958' matches 'ENST00000448958.2').
    This handles the common case where the GTF used by cal_TE carries
    versioned transcript IDs but the quantification result (from
    LongReadQuant or an external tool) does not, or vice versa.

    Parameters
    ----------
    isoforms : list[str]   isoform IDs from the count matrix
    te_df    : DataFrame   indexed by transcript_id (TE annotation table)

    Returns
    -------
    dict  {isoform_id: te_df_index_label_or_None}
    """
    exact_ids = set(te_df.index)
    stripped_to_exact = {}
    for tid in te_df.index:
        stripped_to_exact.setdefault(_strip_version(tid), tid)

    resolved = {}
    for iso in isoforms:
        if iso in exact_ids:
            resolved[iso] = iso
        else:
            resolved[iso] = stripped_to_exact.get(_strip_version(iso))
    return resolved


def build_aggregation_matrices(isoforms, te_df, percent_threshold: float):
    """
    Build isoform → TE-level aggregation matrices for efficient sparse
    matrix multiplication.

    For each level (class / family / subfamily / locus), creates a binary
    isoform x TE sparse matrix:  A[i, j] = 1  if isoform i is associated
    with TE-level j (proportion >= threshold).

    Also builds isoform → transcript_classification mapping.

    Isoform IDs are matched against te_df's transcript_id exactly first,
    falling back to a version-stripped comparison (see _strip_version) so
    that a version-number mismatch alone does not cause every isoform to be
    treated as unmatched.

    Returns
    -------
    agg : dict
        Keys: 'classification', 'class', 'family', 'subfamily', 'locus'
        Each value: (scipy.sparse.csr_matrix, list_of_column_labels)
          shape: (n_isoforms x n_levels)
    """
    n_iso = len(isoforms)
    isoform_to_idx = {iso: i for i, iso in enumerate(isoforms)}

    # Index the TE table by transcript_id for fast lookup
    te_df = te_df.copy()
    # drop duplicates, keep first
    te_df = te_df.drop_duplicates(subset="transcript_id", keep="first")
    te_df = te_df.set_index("transcript_id")

    resolved_ids = build_te_id_resolver(isoforms, te_df)

    print(f"[TE-SC] Building aggregation matrices (threshold={percent_threshold:.2f})...",
          flush=True)

    # Collect unique labels per level
    loci_set    = set()
    subfam_set  = set()
    family_set  = set()
    class_set   = set()

    iso_entries = {}   # isoform_idx -> list of TE entry dicts

    for iso, iso_idx in isoform_to_idx.items():
        te_id = resolved_ids[iso]
        if te_id is None:
            continue
        row = te_df.loc[te_id]
        entries = _parse_te_entries(row, percent_threshold)
        if entries:
            iso_entries[iso_idx] = entries
            for e in entries:
                if e["locus"]:      loci_set.add(e["locus"])
                if e["subfamily"] not in ("", "Unknown", "nan"):
                    subfam_set.add(e["subfamily"])
                if e["family"] not in ("", "Unknown", "nan"):
                    family_set.add(e["family"])
                if e["class_name"] not in ("", "Unknown", "nan"):
                    class_set.add(e["class_name"])

    # --- classification mapping ---
    all_classifications = sorted(set(
        str(te_df.loc[resolved_ids[iso], "transcript_classification"])
        if resolved_ids[iso] is not None else "Gene-alone"
        for iso in isoforms
    ) | {"Gene-alone"})
    class_label_to_idx = {c: i for i, c in enumerate(all_classifications)}

    class_rows, class_cols = [], []
    for iso_idx, iso in enumerate(isoforms):
        te_id = resolved_ids[iso]
        if te_id is not None:
            cls = str(te_df.loc[te_id, "transcript_classification"])
        else:
            cls = "Gene-alone"
        c_idx = class_label_to_idx.get(cls, class_label_to_idx.get("Gene-alone", 0))
        class_rows.append(iso_idx)
        class_cols.append(c_idx)

    clf_mat = scipy.sparse.csr_matrix(
        (np.ones(n_iso), (class_rows, class_cols)),
        shape=(n_iso, len(all_classifications)),
        dtype=np.float32,
    )

    def _build_mat(id_set, key):
        id_list = sorted(id_set)
        id_to_idx = {t: i for i, t in enumerate(id_list)}
        rows, cols = [], []
        for iso_idx, entries in iso_entries.items():
            seen = set()
            for e in entries:
                te_id = e[key]
                if te_id in id_to_idx and te_id not in seen:
                    rows.append(iso_idx)
                    cols.append(id_to_idx[te_id])
                    seen.add(te_id)
        if not rows:
            return scipy.sparse.csr_matrix((n_iso, len(id_list)), dtype=np.float32), id_list
        mat = scipy.sparse.csr_matrix(
            (np.ones(len(rows), dtype=np.float32), (rows, cols)),
            shape=(n_iso, len(id_list)),
        )
        return mat, id_list

    loci_mat,   loci_ids   = _build_mat(loci_set,   "locus")
    subfam_mat, subfam_ids = _build_mat(subfam_set,  "subfamily")
    family_mat, family_ids = _build_mat(family_set,  "family")
    class_mat,  class_ids  = _build_mat(class_set,   "class_name")

    print(f"[TE-SC] TE classification labels : {len(all_classifications)}", flush=True)
    print(f"[TE-SC] Unique TE classes        : {len(class_ids)}",          flush=True)
    print(f"[TE-SC] Unique TE families       : {len(family_ids)}",         flush=True)
    print(f"[TE-SC] Unique TE subfamilies    : {len(subfam_ids)}",         flush=True)
    print(f"[TE-SC] Unique TE loci           : {len(loci_ids)}",           flush=True)
    print(f"[TE-SC] Isoforms with TE hits    : {len(iso_entries)}",        flush=True)

    return {
        "classification": (clf_mat,   all_classifications),
        "class":          (class_mat, class_ids),
        "family":         (family_mat, family_ids),
        "subfamily":      (subfam_mat, subfam_ids),
        "locus":          (loci_mat,   loci_ids),
    }


# ---------------------------------------------------------------------------
# Spot x TE matrix computation
# ---------------------------------------------------------------------------

def compute_spot_te_matrices(count_csr, agg, output_loci: bool = False):
    """
    Compute per-spot TE count matrices via sparse matrix multiplication.

    spot_level_matrix = count_csr  @  iso_level_mat
      (n_spots x n_iso)            (n_iso x n_levels)
    = (n_spots x n_levels)

    Returns dict:  level_name → (scipy.sparse.csr_matrix, column_labels)
    """
    results = {}
    for level, (iso_mat, labels) in agg.items():
        if level == "locus" and not output_loci:
            continue
        if iso_mat.shape[1] == 0:
            continue
        spot_mat = count_csr.dot(iso_mat).tocsr()
        results[level] = (spot_mat, labels)
        print(f"[TE-SC] Computed spot x {level}: "
              f"{spot_mat.shape[0]} x {spot_mat.shape[1]}", flush=True)
    return results


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

# Matrices with more columns than this threshold are written in MEX / MTX
# format instead of TSV.  TSV is fine for small wide tables (class=18,
# family=59, subfamily=~4K) but becomes prohibitively slow and large at
# locus scale (79K+ columns x 77K+ cells = billions of cells).
_MTX_COL_THRESHOLD = 500


def _write_mtx_mex(output_dir: str, base: str, mat, barcodes: list,
                   features: list, barcode_label: str) -> list:
    """Write a sparse matrix in MEX (Market Exchange) format.

    Creates three files:
      <base>_barcodes.tsv   — one barcode per line
      <base>_features.tsv   — one feature ID per line
      <base>_matrix.mtx     — sparse coordinate format (feature x barcode)

    Returns list of (filename, description) tuples for the README.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    bc_path  = os.path.join(output_dir, f"{base}_barcodes.tsv")
    feat_path = os.path.join(output_dir, f"{base}_features.tsv")
    mtx_path  = os.path.join(output_dir, f"{base}_matrix.mtx")

    with open(bc_path, "w") as f:
        f.write("\n".join(barcodes) + "\n")
    with open(feat_path, "w") as f:
        f.write("\n".join(str(x) for x in features) + "\n")

    # Standard MEX convention: rows = features, cols = barcodes
    coo = mat.T.tocoo().astype(int)
    n_feat = len(features)
    n_bc   = len(barcodes)
    with open(mtx_path, "w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"%metadata_json: {{\"software_version\": \"LongReadQuant-TE\"}}\n")
        f.write(f"{n_feat} {n_bc} {coo.nnz}\n")
        for r, c, v in zip(coo.row, coo.col, coo.data):
            f.write(f"{r + 1} {c + 1} {int(v)}\n")

    return [
        (f"{base}_barcodes.tsv", f"{n_bc} {barcode_label}s, one per line."),
        (f"{base}_features.tsv", f"{n_feat} feature IDs, one per line."),
        (f"{base}_matrix.mtx",
         f"Sparse MTX: {n_feat} features x {n_bc} {barcode_label}s, "
         f"{coo.nnz} non-zero entries.  "
         "Load with: ReadMtx() in Seurat or scanpy.read_mtx()."),
    ]


def _level_filenames(unit: str) -> dict:
    return {
        "classification": f"{unit}_transcript_classification_counts.tsv",
        "class":          f"{unit}_te_class_counts.tsv",
        "family":         f"{unit}_te_family_counts.tsv",
        "subfamily":      f"{unit}_te_subfamily_counts.tsv",
        "locus":          f"{unit}_te_loci_counts.tsv",
    }

def _level_descriptions(unit: str) -> dict:
    return {
        "classification": f"Per-{unit} UMI counts by transcript classification (Gene-alone / TE-Gene / TE-alone).",
        "class":          f"Per-{unit} UMI counts aggregated to TE class level (LINE, SINE, LTR, DNA, ...).",
        "family":         f"Per-{unit} UMI counts aggregated to TE family level.",
        "subfamily":      f"Per-{unit} UMI counts aggregated to TE subfamily level.",
        "locus":          f"Per-{unit} UMI counts aggregated to individual TE locus (may be large).",
    }


def write_spot_te_outputs(output_dir: str, spot_matrices: dict, barcodes: list,
                          barcode_label: str = "barcode"):
    """Write per-spot/per-cell TE matrices to TSV files and a README."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    unit = "spot" if barcode_label == "spot_barcode" else "cell"
    filenames    = _level_filenames(unit)
    descriptions = _level_descriptions(unit)
    written = []

    for level, (mat, labels) in spot_matrices.items():
        if not labels:
            continue

        n_units, n_cols = mat.shape

        if n_cols > _MTX_COL_THRESHOLD:
            # Large matrix: write as MEX/MTX sparse format to avoid writing
            # billions of cells as text (e.g. 77K cells x 79K loci = 6 B cells).
            base = os.path.splitext(
                filenames.get(level, f"{unit}_{level}_counts.tsv")
            )[0]
            mtx_entries = _write_mtx_mex(
                output_dir, base, mat, barcodes, labels, barcode_label
            )
            for mfname, mdesc in mtx_entries:
                written.append((mfname, n_units, n_cols, mdesc))
            print(f"[TE-SC] Written MTX: {output_dir}/{base}_*  "
                  f"({n_units} {unit}s x {n_cols} {level})", flush=True)
        else:
            fname = filenames.get(level, f"{unit}_{level}_counts.tsv")
            fpath = os.path.join(output_dir, fname)
            if n_cols > 500:
                df = pd.DataFrame.sparse.from_spmatrix(
                    mat.tocsr(), index=barcodes, columns=labels
                )
            else:
                df = pd.DataFrame(
                    mat.toarray().astype(int), index=barcodes, columns=labels
                )
            df.index.name = barcode_label
            try:
                df = df.astype(pd.SparseDtype(int, fill_value=0))
            except Exception:
                pass
            df.to_csv(fpath, sep="\t")
            written.append((fname, n_units, n_cols, descriptions.get(level, "")))
            print(f"[TE-SC] Written: {fpath}  ({n_units} {unit}s x {n_cols} {level})",
                  flush=True)

    # README
    readme_lines = [
        f"miniQuant TE {unit} output — per-{unit} TE quantification",
        "=" * 54,
        "",
        "Generated by: miniQuant cal_TE  (te_quantification.py)",
        "",
        "Files",
        "-----",
    ]
    for fname, n_units, n_cols, desc in written:
        readme_lines += [
            f"  {fname}",
            f"    {n_units} {unit}s x {n_cols} columns.",
            f"    {desc}",
            f"    Row index: {barcode_label}.  Values: integer UMI counts.",
            "",
        ]
    readme_lines += [
        "Notes",
        "-----",
        "  - UMI counts attributed to all overlapping TE entries of each isoform",
        f"    (attribute-to-all policy; counts may therefore exceed total per-{unit}",
        "    isoform counts when isoforms overlap multiple TE classes/families).",
        "  - Only TEs with overlap_proportion >= percent_threshold are considered.",
        "  - Isoforms present in the GTF but absent from the count matrix are ignored.",
        "",
        "Downstream analysis",
        "-------------------",
        "  Compatible with Seurat, Scanpy, Squidpy, etc.",
        "  Example (Seurat):",
        f"    te_mat <- read.table('{unit}_te_class_counts.tsv', header=TRUE, row.names=1, sep='\\t')",
        "    seu[['TE_class']] <- CreateAssayObject(counts = t(te_mat))",
    ]
    with open(os.path.join(output_dir, "README.txt"), "w") as f:
        f.write("\n".join(readme_lines) + "\n")

    print(f"[TE-SC] README.txt written to {output_dir}", flush=True)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_sc_te_analysis(
    te_table_path: str,
    quant_dir: str,
    output_dir: str,
    percent_threshold: float = 0.5,
    output_loci: bool = True,
    barcode_label: str = "barcode",
):
    """
    Main entry point for per-spot TE analysis.

    Parameters
    ----------
    te_table_path : str
        Path to  transcript_quantification_with_TE.tsv  from a cal_TE run.
    quant_dir : str
        Directory containing the isoform-level MEX output from a prior
        LongReadQuant quantify --st_mode (or --sc_mode) run.
    output_dir : str
        Output directory for per-spot TE metric files.
    percent_threshold : float
        Minimum TE-overlap proportion for associating an isoform with a TE.
    output_loci : bool
        If True, also output spot_te_loci_counts.tsv (may be very large).
    barcode_label : str
        Row index label in output TSVs ('barcode', 'spot_barcode', etc.).
    """
    print("\n" + "=" * 60, flush=True)
    print("[TE-SC] Per-spot TE quantification", flush=True)
    print("=" * 60, flush=True)

    # 1. Load TE annotation table
    print(f"[TE-SC] Loading TE annotation: {te_table_path}", flush=True)
    te_df = pd.read_csv(te_table_path, sep="\t", low_memory=False)
    print(f"[TE-SC] {len(te_df)} transcripts in TE table", flush=True)

    # 2. Load count matrix
    count_csr, barcodes, isoforms = load_count_matrix(quant_dir)

    # 2b. Sanity check: isoform ID overlap between the count matrix and the TE
    # annotation table. A low match rate usually means -gtf does not match the
    # annotation that produced this quantification result (e.g. an external
    # tool's novel/discovered transcripts are absent from -gtf). Matching uses
    # the same exact-or-version-stripped resolver as build_aggregation_matrices
    # (see _strip_version), so a version-number-only mismatch (e.g.
    # 'ENST00000448958' vs 'ENST00000448958.2') no longer reads as 0% here.
    resolved_ids = build_te_id_resolver(isoforms, te_df.set_index("transcript_id"))
    matched = sum(1 for te_id in resolved_ids.values() if te_id is not None)
    match_rate = matched / len(isoforms) if isoforms else 0.0
    print(f"[TE-SC] Isoform ID match rate vs TE table: {match_rate:.1%} "
          f"({matched}/{len(isoforms)})", flush=True)
    if match_rate < 0.5:
        print("[TE-SC][WARNING] Less than half of the isoform IDs in the count "
              "matrix were found in the TE annotation table. TE counts will be "
              "incomplete for the unmatched isoforms. Check that -gtf matches "
              "the annotation used to produce this quantification result "
              "(e.g. novel transcripts from an external tool's discovery mode "
              "will not be found unless -gtf includes them).", flush=True)

    # 3. Build isoform → TE aggregation matrices
    agg = build_aggregation_matrices(isoforms, te_df, percent_threshold)

    # 4. Compute spot x TE matrices
    spot_matrices = compute_spot_te_matrices(count_csr, agg, output_loci=output_loci)

    # 5. Write outputs
    write_spot_te_outputs(output_dir, spot_matrices, barcodes, barcode_label=barcode_label)

    print(f"\n[TE-SC] Done. Results saved in: {output_dir}", flush=True)
