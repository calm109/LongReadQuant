#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create TE GTF annotation from RepeatMasker .out (revised version)
Field conventions:
  - gene_id       = subfamily
  - transcript_id = locus_name (unique name for this TE copy; subfamily.continuous number)
  - family_id     = family
  - class_id      = class
  - locus         = "1:268690-268968" (without chr prefix)
Also retain consensus_length / truncated_proportion as optional attributes.
"""

import argparse
import os
from collections import defaultdict

def _nochr(chrom: str) -> str:
    """Remove chr prefix and normalize: chrM→MT."""
    c = chrom
    if c.lower().startswith("chr"):
        c = c[3:]
    c = c.upper()
    if c == "M":
        c = "MT"
    return c

def convert_repeatmasker_to_gtf(repeatmasker_file: str, output_gtf: str):
    """
    Convert RepeatMasker .out to GTF (gene_id=subfamily, transcript_id=locus_name, locus=1:beg-end)
    """
    print(f"Reading RepeatMasker file: {repeatmasker_file}")
    with open(repeatmasker_file, "r") as f:
        lines = f.readlines()

    # Skip first three header lines
    data_lines = lines[3:]

    gtf_lines = []
    # Count occurrences per subfamily to generate locus_name (subfamily.sequence number)
    subfamily_counter = defaultdict(int)
    skipped = 0

    print("Starting conversion …")
    for ln, line in enumerate(data_lines, 4):  # line number from actual file line count
        if not line.strip():
            continue

        # RepeatMasker .out uses column-width format, split() works safely
        fields = [x for x in line.split() if x]
        # Classic .out has at least 15 columns: see UCSC/RepeatMasker specification
        if len(fields) < 15:
            skipped += 1
            continue

        try:
            # Index alignment (according to UCSC/RepeatMasker standard .out):
            # 0=score 1=div 2=del 3=ins 4=chrom 5=start 6=end 7=left
            # 8=strand 9=repeat_name 10=class/family 11=rep_start 12=rep_end 13=rep_left 14=ID
            chrom_raw = fields[4]
            beg = int(fields[5])
            end = int(fields[6])
            strand = fields[8]  # '+' | 'C' | '?'
            subfamily = fields[9]  # matching repeat = subfamily
            class_family = fields[10]

            if not chrom_raw.startswith("chr"):
                chrom = f"chr{chrom_raw}"
            else:
                chrom = chrom_raw

            # Normalize strand
            if strand == "C":
                strand = "-"
            elif strand not in ["+", "-"]:
                strand = "+"

            # Split class/family
            if "/" in class_family:
                cls, fam = class_family.split("/", 1)
            else:
                cls, fam = class_family, class_family

            # Generate locus_name: incrementing number within same subfamily
            subfamily_counter[subfamily] += 1
            locus_name = f"{subfamily}.{subfamily_counter[subfamily]}"

            # locus (without chr prefix)
            locus = f"{_nochr(chrom)}:{beg}-{end}"

            # Estimate consensus length and truncation proportion (more accurate if rep_start/rep_end/rep_left present)
            consensus_length = end - beg + 1
            truncated_proportion = 0.0
            try:
                rep_start = int(fields[11].strip("()"))
                rep_end = int(fields[12])
                rep_left = int(fields[13].strip("()"))
                if strand == "+":
                    total_consensus = rep_end + rep_left
                    actual = rep_end - rep_start
                else:
                    total_consensus = rep_start + rep_end
                    actual = rep_end - rep_left
                if total_consensus > 0:
                    truncated_proportion = round(actual / total_consensus, 3)
                    consensus_length = total_consensus
            except Exception:
                pass  # if fields missing, fall back to genomic interval length

            # Assemble GTF attributes — core five-level information
            attributes = (
                f'gene_id "{subfamily}"; '               # subfamily
                f'transcript_id "{locus_name}"; '        # locus_name
                f'family_id "{fam}"; '                   # family
                f'class_id "{cls}"; '                    # class
                f'locus "{locus}"; '                     # specific location (without chr)
                f'consensus_length "{consensus_length}"; '
                f'truncated_proportion "{truncated_proportion}";'
            )

            gtf_line = f"{chrom}\tRepeatMasker\texon\t{beg}\t{end}\t.\t{strand}\t.\t{attributes}\n"
            gtf_lines.append(gtf_line)

            if len(gtf_lines) % 10000 == 0:
                print(f"  Processed {len(gtf_lines)} entries")

        except Exception as e:
            skipped += 1
            if skipped <= 10:
                print(f"  Skipping line {ln} (format issue): {e}")

    print(f"\nWriting GTF: {output_gtf}")
    with open(output_gtf, "w") as fout:
        fout.writelines(gtf_lines)

    print("Conversion completed!")
    print(f"  Success: {len(gtf_lines)} entries")
    print(f"  Skipped: {skipped} entries")
    return gtf_lines

def _cli():
    p = argparse.ArgumentParser(description="RepeatMasker .out → GTF (gene_id=subfamily, transcript_id=locus_name, locus=1:beg-end)")
    p.add_argument("-i", "--input", required=True, help="Path to .out file")
    p.add_argument("-o", "--output", required=True, help="Path to output GTF")
    args = p.parse_args()
    convert_repeatmasker_to_gtf(args.input, args.output)

if __name__ == "__main__":
    _cli()