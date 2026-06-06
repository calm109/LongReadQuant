#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Replace TPM in "skip quantification" result files with real quantification results (e.g., isoQuant)

Example usage:
    python replace_tpm_from_quant.py \
        -i TE_initiated_transcripts.tsv \
        -q isoquant_isoform_expression.tsv \
        -o TE_initiated_transcripts.isoquant_TPM.tsv

If column names cannot be automatically detected, they can be specified manually:
    python replace_tpm_from_quant.py \
        -i TE_initiated_transcripts.tsv \
        -q isoquant_isoform_expression.tsv \
        -o TE_initiated_transcripts.isoquant_TPM.tsv \
        --te-id-col transcript_id \
        --te-tpm-col transcript_TPM \
        --quant-id-col transcript_id \
        --quant-tpm-col TPM
"""

import argparse
from pathlib import Path
import sys

import pandas as pd


def detect_col(df, candidates, kind):
    """
    Automatically identify column name in DataFrame:
    - First try exact match with candidate list
    - Then case-insensitive match
    - For TPM column, if still not found, look for any column containing 'tpm'
    """
    # Exact match
    for c in candidates:
        if c in df.columns:
            return c

    # Case-insensitive match
    lower2col = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower2col:
            return lower2col[c.lower()]

    # Extra fallback for TPM: contains 'tpm'
    if kind.lower() == "tpm":
        tpm_like = [c for c in df.columns if "tpm" in c.lower()]
        if len(tpm_like) == 1:
            return tpm_like[0]

    raise SystemExit(
        f"Unable to automatically identify {kind} column in table."
        f"\n  Current columns: {list(df.columns)}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Replace TPM column in 'skip quantification' result files with real quantification results (e.g., isoQuant)"
    )
    parser.add_argument(
        "-i", "--input-te",
        required=True,
        help="Path to the result file produced by skip quantification (e.g., TE_initiated_transcripts.tsv or transcript_quantification_with_TE.tsv)"
    )
    parser.add_argument(
        "-q", "--quant",
        required=True,
        help="Path to real quantification result file (e.g., isoQuant isoform quantification result)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to output file (TPM replaced)"
    )
    parser.add_argument(
        "--te-id-col",
        help="Transcript ID column name in TE result file (default auto-detected, e.g., transcript_id / Isoform / isoform_id)"
    )
    parser.add_argument(
        "--te-tpm-col",
        help="TPM column name in TE result file to be replaced (default auto-detected, e.g., transcript_TPM / TPM / tpm)"
    )
    parser.add_argument(
        "--quant-id-col",
        help="Transcript ID column name in quantification result file (default auto-detected, e.g., transcript_id / Isoform / isoform_id / target_id etc.)"
    )
    parser.add_argument(
        "--quant-tpm-col",
        help="TPM column name in quantification result file (default auto-detected, column containing TPM/tpm)"
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help="Input/output file separator, default tab ('\\t')"
    )

    args = parser.parse_args()

    te_path = Path(args.input_te)
    quant_path = Path(args.quant)
    out_path = Path(args.output)

    if not te_path.exists():
        sys.exit(f"Error: Skip quantification result file not found: {te_path}")
    if not quant_path.exists():
        sys.exit(f"Error: Quantification result file not found: {quant_path}")

    print(f"Reading skip quantification result file: {te_path}")
    te_df = pd.read_csv(te_path, sep=args.sep)

    print(f"Reading quantification result file: {quant_path}")
    quant_df = pd.read_csv(quant_path, sep=args.sep)

    # Auto/manual column identification
    te_id_col = args.te_id_col or detect_col(
        te_df,
        ["transcript_id", "Isoform", "isoform_id", "transcript", "target_id", "Name"],
        kind="transcript ID (TE file)"
    )
    te_tpm_col = args.te_tpm_col or detect_col(
        te_df,
        ["transcript_TPM", "TPM", "tpm", "Tpm"],
        kind="TPM (TE file)"
    )
    quant_id_col = args.quant_id_col or detect_col(
        quant_df,
        ["transcript_id", "Isoform", "isoform_id", "transcript", "target_id", "Name"],
        kind="transcript ID (quant file)"
    )
    quant_tpm_col = args.quant_tpm_col or detect_col(
        quant_df,
        ["TPM", "tpm", "Tpm", "IsoTPM"],
        kind="TPM (quant file)"
    )

    print(f"TE file transcript ID column: {te_id_col}")
    print(f"TE file TPM column: {te_tpm_col}")
    print(f"Quant file transcript ID column: {quant_id_col}")
    print(f"Quant file TPM column: {quant_tpm_col}")

    # Build ID -> TPM mapping
    # If a transcript appears multiple times in quant file, the last TPM is used by default
    quant_map = quant_df.set_index(quant_id_col)[quant_tpm_col]

    # Replace using map
    new_tpm_series = te_df[te_id_col].map(quant_map)

    updated_mask = new_tpm_series.notna()
    replaced = int(updated_mask.sum())
    total = len(te_df)
    missing = total - replaced

    # Only replace TPM for IDs found in quantification results, keep original for others
    te_df.loc[updated_mask, te_tpm_col] = new_tpm_series[updated_mask].astype(float)

    print(f"Successfully updated TPM for transcripts found in quantification: {replaced}/{total}")
    print(f"Transcripts not found in quantification file: {missing} (TPM kept unchanged)")

    # Output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    te_df.to_csv(out_path, sep=args.sep, index=False)
    print(f"Saved new file (only TPM replaced, all other columns preserved): {out_path}")


if __name__ == "__main__":
    main()