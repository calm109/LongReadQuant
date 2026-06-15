"""
cal_te_main.py
--------------
Entry point for the cal_TE subcommand of miniQuant.

Wraps the transcript-TE overlap analysis from integrated_quantification_sam_only.py
and provides a clean interface for main.py, with **optional multi-process parallelism**.

Parallelism strategy
--------------------
The two computationally heavy steps are:
  • find_first_exon_TE_overlaps()          — O(n_transcripts x log n_TE)
  • calculate_enhanced_all_exon_TE_overlaps() — same order

Both iterate over transcripts independently, so they can be parallelized by
grouping transcripts by chromosome and distributing chromosome groups across
worker processes.  Each worker only receives the TE intervals for its own
chromosomes, keeping pickling overhead small.

When threads=1 (default), everything runs sequentially in the main process —
no subprocess overhead.

Workflow
--------
1. Parse transcript GTF → per-transcript exon structures
2. Parse TE GTF        → TE interval trees per chromosome/strand
3. [parallel] Find first-exon TE overlaps (TSS/TES window analysis)
4. [parallel] Find all-exon TE overlaps
5. Classify each transcript (Gene-alone / TE-Gene / TE-alone)
6. Write  transcript_quantification_with_TE.tsv  and  TE_derived_transcripts.tsv
"""

import os
import sys
import math
import multiprocessing as mp
from collections import defaultdict
from pathlib import Path


# ===========================================================================
# Multiprocessing worker  (must be at module level for spawn-safe pickling)
# ===========================================================================

def _chrom_overlap_worker(packed):
    """
    Worker function for parallel TE overlap computation.

    Receives a subset of transcripts (belonging to specific chromosomes)
    and the TE interval data for those same chromosomes.
    Returns (first_exon_overlaps, all_exon_overlaps) dicts for the subset.

    Parameters
    ----------
    packed : tuple
        (worker_id, transcript_exons_chunk, te_intervals_plain_chunk)
        • transcript_exons_chunk : dict {tid -> exon_list}
        • te_intervals_plain_chunk : dict {chrom -> {strand -> IntervalTree}}
          (plain dict, NOT defaultdict — lambdas are not picklable)
    """
    worker_id, transcript_exons_chunk, te_plain_chunk = packed

    # Ensure TE_analyse directory is on sys.path in the worker process
    _te_dir = os.path.dirname(os.path.abspath(__file__))
    if _te_dir not in sys.path:
        sys.path.insert(0, _te_dir)

    from integrated_quantification_sam_only import (
        find_first_exon_TE_overlaps,
        calculate_enhanced_all_exon_TE_overlaps,
    )

    fe = find_first_exon_TE_overlaps(transcript_exons_chunk, te_plain_chunk)
    ae = calculate_enhanced_all_exon_TE_overlaps(transcript_exons_chunk, te_plain_chunk)
    return fe, ae


# ===========================================================================
# Helpers
# ===========================================================================

def _te_intervals_to_plain(te_intervals):
    """
    Convert  defaultdict(lambda: defaultdict(IntervalTree))  to a regular
    nested dict so that it can be pickled by multiprocessing workers.

    The IntervalTree objects themselves are left as-is (they are picklable).
    Only the two outer defaultdict wrappers — which use un-picklable lambdas
    as their default_factory — are converted to plain dicts.
    """
    return {
        chrom: dict(strand_dict)
        for chrom, strand_dict in te_intervals.items()
    }


def _split_transcripts_by_chrom(transcript_exons):
    """
    Group transcripts by the chromosome of their first exon.

    Returns
    -------
    dict : {chrom -> {tid -> exon_list}}
    """
    by_chrom = defaultdict(dict)
    for tid, exons in transcript_exons.items():
        if exons:
            chrom = exons[0]['chrom']
            by_chrom[chrom][tid] = exons
    return dict(by_chrom)


def _distribute_chroms(chroms, n_workers):
    """
    Assign chromosomes to workers using round-robin so load is roughly
    balanced (human genome chromosomes vary considerably in size, so
    strict round-robin won't be perfectly balanced — but it's fast to
    compute and avoids a dependency on chrom size data).

    Returns
    -------
    list of lists : [[chrom, ...], ...]  — one inner list per worker
    """
    buckets = [[] for _ in range(n_workers)]
    for i, chrom in enumerate(sorted(chroms)):
        buckets[i % n_workers].append(chrom)
    return buckets


def _run_parallel_overlaps(transcript_exons, te_intervals, threads):
    """
    Run find_first_exon_TE_overlaps and calculate_enhanced_all_exon_TE_overlaps
    in parallel across  `threads`  worker processes.

    Groups transcripts by chromosome; each worker handles a subset of
    chromosomes and only receives the TE data for those chromosomes.

    Returns
    -------
    first_exon_overlaps, all_exon_overlaps : merged dicts
    """
    te_plain = _te_intervals_to_plain(te_intervals)
    by_chrom = _split_transcripts_by_chrom(transcript_exons)
    chroms   = sorted(by_chrom.keys())

    effective_workers = min(threads, len(chroms))
    print(f"[cal_TE] Parallel mode: {effective_workers} worker(s) "
          f"for {len(chroms)} chromosomes", flush=True)

    buckets = _distribute_chroms(chroms, effective_workers)

    worker_args = []
    for wid, chrom_list in enumerate(buckets):
        if not chrom_list:
            continue
        # Gather transcripts and TE data for this worker's chromosomes
        chunk_exons = {}
        for chrom in chrom_list:
            chunk_exons.update(by_chrom.get(chrom, {}))
        chunk_te = {chrom: te_plain[chrom]
                    for chrom in chrom_list if chrom in te_plain}
        worker_args.append((wid, chunk_exons, chunk_te))

    if effective_workers == 1:
        # Avoid subprocess overhead when only one worker is needed
        results = [_chrom_overlap_worker(worker_args[0])]
    else:
        # Use 'spawn' context (Windows-safe; also safe on Linux/Mac)
        ctx = mp.get_context('spawn')
        with ctx.Pool(effective_workers) as pool:
            results = pool.map(_chrom_overlap_worker, worker_args)

    # Merge results from all workers
    first_exon_overlaps = {}
    all_exon_overlaps   = {}
    for fe, ae in results:
        first_exon_overlaps.update(fe)
        all_exon_overlaps.update(ae)

    return first_exon_overlaps, all_exon_overlaps


# ===========================================================================
# Public API
# ===========================================================================

def run_te_annotation(
    gtf: str,
    te_gtf: str,
    output_dir: str,
    skip_quantification: bool = True,
    threads: int = 1,
    first_exon_threshold: float = 50.0,
    total_threshold: float = 50.0,
    te_overlap_threshold: int = 10,
    te_ratio_threshold: float = 80.0,
    te_feature_threshold: int = 50,
):
    """
    Run transcript-TE overlap analysis.

    Parameters
    ----------
    gtf : str
        Path to transcript GTF annotation file.
    te_gtf : str
        Path to TE GTF annotation file (created by create_simple_TE_gtf.py).
    output_dir : str
        Output directory.
    skip_quantification : bool
        If True (default for cal_TE), use dummy TPM=0 for all transcripts.
        Real TPM values can be merged later with replace_tpm_from_quant.py.
    threads : int
        Number of parallel worker processes for the TE overlap computation
        (steps 3 & 4).  Transcripts are grouped by chromosome and distributed
        across workers.  threads=1 runs sequentially (no subprocess overhead).
    first_exon_threshold : float
        First-exon TE proportion threshold (%) for TE-derived transcript detection.
    total_threshold : float
        Full-transcript TE proportion threshold (%) for TE-derived transcript detection.
    te_overlap_threshold : int
        Minimum TE overlap length (bp) to distinguish Gene-alone from others.
    te_ratio_threshold : float
        TE proportion in full transcript (%) to identify TE-alone transcripts.
    te_feature_threshold : int
        TE overlap length (bp) within TSS/TES 200 bp window for sub-classification.

    Outputs
    -------
    <output_dir>/transcript_quantification_with_TE.tsv
        Per-transcript TE annotation table (all transcripts in GTF).
    <output_dir>/TE_derived_transcripts.tsv
        Subset of transcripts meeting TE-derived thresholds, with dominant TE info.

    Returns
    -------
    str — path to transcript_quantification_with_TE.tsv
    """
    # Make TE_analyse functions importable
    _te_dir = os.path.dirname(os.path.abspath(__file__))
    if _te_dir not in sys.path:
        sys.path.insert(0, _te_dir)

    from integrated_quantification_sam_only import (
        parse_transcript_annotation,
        parse_te_annotation,
        find_first_exon_TE_overlaps,
        calculate_enhanced_all_exon_TE_overlaps,
        generate_enhanced_integrated_output,
    )

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    print("=" * 60, flush=True)
    print("[cal_TE] Transcript TE annotation", flush=True)
    mode_str = "--skip-quantification (dummy TPM=0)" if skip_quantification else "full"
    print(f"[cal_TE] Mode    : {mode_str}", flush=True)
    print(f"[cal_TE] Threads : {threads}", flush=True)
    print("=" * 60, flush=True)

    # ------------------------------------------------------------------
    # Step 1: Parse transcript annotation
    # ------------------------------------------------------------------
    print("\n[cal_TE] Step 1: Parsing transcript annotation...", flush=True)
    transcript_exons, transcript_info = parse_transcript_annotation(gtf)
    print(f"[cal_TE] Transcripts parsed: {len(transcript_exons)}", flush=True)

    # ------------------------------------------------------------------
    # Step 2: Parse TE annotation
    # ------------------------------------------------------------------
    print("\n[cal_TE] Step 2: Parsing TE annotation...", flush=True)
    te_intervals = parse_te_annotation(te_gtf)

    # ------------------------------------------------------------------
    # Steps 3 & 4: TE overlap analysis  (sequential or parallel)
    # ------------------------------------------------------------------
    if threads > 1:
        print(f"\n[cal_TE] Step 3 & 4: TE overlap analysis "
              f"(parallel, {threads} threads)...", flush=True)
        first_exon_overlaps, all_exon_overlaps = _run_parallel_overlaps(
            transcript_exons, te_intervals, threads
        )
    else:
        print("\n[cal_TE] Step 3: Finding first-exon TE overlaps...", flush=True)
        first_exon_overlaps = find_first_exon_TE_overlaps(
            transcript_exons, te_intervals
        )
        print("\n[cal_TE] Step 4: Calculating all-exon TE overlaps...", flush=True)
        all_exon_overlaps = calculate_enhanced_all_exon_TE_overlaps(
            transcript_exons, te_intervals
        )

    print(f"[cal_TE] Transcripts with TE overlaps: "
          f"{len(all_exon_overlaps)}", flush=True)

    # ------------------------------------------------------------------
    # Step 5: Generate integrated output
    # ------------------------------------------------------------------
    print("\n[cal_TE] Step 5: Generating integrated output...", flush=True)
    generate_enhanced_integrated_output(
        transcript_output=output_dir,      # unused in skip mode
        first_exon_overlaps=first_exon_overlaps,
        all_exon_overlaps=all_exon_overlaps,
        output_path=output_dir,
        transcript_info=transcript_info,
        first_exon_threshold=first_exon_threshold,
        total_threshold=total_threshold,
        skip_quantification=skip_quantification,
        te_overlap_threshold=te_overlap_threshold,
        te_ratio_threshold=te_ratio_threshold,
        te_feature_threshold=te_feature_threshold,
    )

    te_table = os.path.join(output_dir, "transcript_quantification_with_TE.tsv")
    print("\n[cal_TE] Annotation complete.", flush=True)
    print(f"[cal_TE] Main output : {te_table}", flush=True)
    return te_table
