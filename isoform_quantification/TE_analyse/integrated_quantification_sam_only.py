#!/usr/bin/env python3
"""
Integrated transcript quantification and TE detection - fixed version
Fix miniQuant command format issue
"""

import argparse
import os
import sys
import pandas as pd
from pathlib import Path
from collections import defaultdict
from intervaltree import IntervalTree
import re
import subprocess

# Path configuration
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MINIQUANT_PATH = os.path.join(SCRIPT_DIR, 'miniQuant_old', 'isoform_quantification')
MINIQUANT_MAIN = os.path.abspath(os.path.join(MINIQUANT_PATH, 'main.py'))

def parse_GTF_feature(feature_string):
    """Parse GTF feature string"""
    feature_dict = {}
    key_val_pairs = [s for s in feature_string.split(";") if (s!= "" and not s.isspace())]
    for key_val_pair_string in key_val_pairs:
        key_val_pair = [s for s in key_val_pair_string.split(" ") if s != ""]
        if len(key_val_pair) >= 2:
            key, val = key_val_pair[0:2]
            feature_dict[key] = val.replace('"','')
    return feature_dict

def sync_reference_name(ref_name):
    """Synchronize reference sequence name"""
    ref_name = ref_name.upper()
    match = re.search("(?<=CHR).*", ref_name)
    if match:
        ref_name = match.group(0)
    if ref_name == "M":
        ref_name = "MT"
    return ref_name

def parse_transcript_annotation(gtf_path):
    """Parse transcript annotation"""
    transcript_exons = defaultdict(list)
    transcript_info = {}
    
    print(f"Parsing transcript annotation: {gtf_path}")
    
    exon_count = 0
    with open(gtf_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue
                
            chrom, source, feature, start, end, score, strand, frame, attributes = parts
            
            if feature == 'exon':
                feature_dict = parse_GTF_feature(attributes)
                transcript_id = feature_dict.get('transcript_id')
                gene_id = feature_dict.get('gene_id', '')
                
                if transcript_id:
                    chrom = sync_reference_name(chrom)
                    start, end = int(start), int(end)
                    
                    transcript_exons[transcript_id].append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'length': end - start +1
                    })
                    
                    if transcript_id not in transcript_info:
                        transcript_info[transcript_id] = {
                            'gene_id': gene_id,
                            'length': 0
                        }
                    
                    transcript_info[transcript_id]['length'] += (end - start +1)
                    exon_count += 1
    
    # Sort exons for each transcript
    # For positive strand (+) sort by start ascending, for negative strand (-) sort by start descending
    for transcript_id, exons in transcript_exons.items():
        if exons:
            strand = exons[0]['strand']  # get strand information from first exon
            if strand == '+':
                transcript_exons[transcript_id].sort(key=lambda x: x['start'])
            else:  # strand == '-'
                transcript_exons[transcript_id].sort(key=lambda x: x['start'], reverse=True)
    
    print(f"Parsing completed: {len(transcript_exons)} transcripts, {exon_count} exons")
    return transcript_exons, transcript_info

def parse_te_annotation(te_gtf_path):
    """Parse new format TE annotation - adapted for updated GTF format
    New format field mapping:
    - gene_id = subfamily
    - transcript_id = locus_name (e.g., L1HS.1)
    - family_id = family
    - class_id = class
    - locus = specific location (e.g., 1:268690-268968)
    """
    te_intervals = defaultdict(lambda: defaultdict(IntervalTree))
    
    print(f"Parsing new format TE annotation: {te_gtf_path}")
    
    # Statistics
    parsed_lines = 0
    class_stats = defaultdict(int)
    family_stats = defaultdict(int)
    subfamily_stats = defaultdict(int)
    
    with open(te_gtf_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue
                
            chrom, source, feature, start, end, score, strand, frame, attributes = parts
            
            # Parse attributes
            feature_dict = parse_GTF_feature(attributes)
            
            # New format field mapping
            subfamily = feature_dict.get('gene_id', 'Unknown')  # gene_id = subfamily
            locus_name = feature_dict.get('transcript_id', 'Unknown')  # transcript_id = locus_name
            family = feature_dict.get('family_id', 'Unknown')  # family_id = family
            te_class = feature_dict.get('class_id', 'Unknown')  # class_id = class
            locus_position = feature_dict.get('locus', 'Unknown')  # locus = position
            
            # Clean field values
            if subfamily and subfamily != 'Unknown':
                subfamily = subfamily.strip('"').strip("'")
            if locus_name and locus_name != 'Unknown':
                locus_name = locus_name.strip('"').strip("'")
            if family and family != 'Unknown':
                family = family.strip('"').strip("'")
            if te_class and te_class != 'Unknown':
                te_class = te_class.strip('"').strip("'")
            if locus_position and locus_position != 'Unknown':
                locus_position = locus_position.strip('"').strip("'")
            
            # repeat_name uses locus_name, if not available then subfamily
            repeat_name = locus_name if locus_name != 'Unknown' else subfamily
            
            chrom = sync_reference_name(chrom)
            start, end = int(start), int(end)
            
            te_info = {
                'repeat_name': repeat_name,
                'subfamily': subfamily,
                'family': family,
                'class': te_class,
                'locus_name': locus_name,
                'locus_position': locus_position,
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'length': end - start +1
            }
            
            te_intervals[chrom][strand][start:end] = te_info
            parsed_lines += 1
            
            # Statistics
            class_stats[te_class] += 1
            family_stats[family] += 1
            subfamily_stats[subfamily] += 1
            
            # Debug: show parsing results for first 5 lines
            if line_num <= 5:
                print(f"  Line{line_num}: {repeat_name} | subfamily:{subfamily} | family:{family} | class:{te_class}")
    
    print(f"New format TE parsing completed: {parsed_lines} elements")
    print(f"Main TE classes (top5): {dict(list(sorted(class_stats.items(), key=lambda x: x[1], reverse=True))[:5])}")
    print(f"Main TE families (top5): {dict(list(sorted(family_stats.items(), key=lambda x: x[1], reverse=True))[:5])}")
    print(f"Main TE subfamilies (top5): {dict(list(sorted(subfamily_stats.items(), key=lambda x: x[1], reverse=True))[:5])}")
    
    return te_intervals

def find_first_exon_TE_overlaps(transcript_exons, te_intervals):
    """Find overlaps between first exon and TE (fixed overlap calculation)"""
    print("Finding first exon TE overlaps (fixed)...")
    
    first_exon_overlaps = {}
    
    for i, (transcript_id, exons) in enumerate(transcript_exons.items(), 1):
        if i % 5000 == 0:
            print(f"Processed {i}/{len(transcript_exons)}")
        
        if not exons:
            continue
        
        first_exon = exons[0]
        chrom = first_exon['chrom']
        strand = first_exon['strand']
        
        if chrom not in te_intervals or strand not in te_intervals[chrom]:
            continue
        
        overlapping_tes_iterator = te_intervals[chrom][strand].overlap(first_exon['start'], first_exon['end'])
        
        # Convert iterator to list for multiple uses
        overlapping_tes = list(overlapping_tes_iterator)
        
        if not overlapping_tes:
            continue

        # 1. Calculate total unique overlap length after deduplication (for exact proportion calculation)
        # Create a temporary IntervalTree and merge overlapping TEs to avoid double counting
        merged_tes = IntervalTree(overlapping_tes)
        merged_tes.merge_overlaps()
        
        total_unique_overlap_length = 0
        for merged_interval in merged_tes:
            # Calculate intersection between merged interval and exon
            intersection_start = max(first_exon['start'], merged_interval.begin)
            intersection_end = min(first_exon['end'], merged_interval.end)
            intersection_length = intersection_end - intersection_start +1
            
            if intersection_length > 0:
                total_unique_overlap_length += intersection_length

        # 2. Collect detailed information for each TE (for finding dominant TE)
        te_overlaps_details = []
        for te_interval in overlapping_tes:
            overlap_start = max(first_exon['start'], te_interval.begin)
            overlap_end = min(first_exon['end'], te_interval.end)
            overlap_length = overlap_end - overlap_start +1
            
            if overlap_length > 0:
                te_info = te_interval.data
                te_overlaps_details.append({
                    'repeat_name': te_info['repeat_name'],
                    'subfamily': te_info['subfamily'],
                    'family': te_info['family'],
                    'class': te_info['class'],
                    'overlap_length': overlap_length,
                    'te_start': te_info['start'],
                    'te_end': te_info['end'],
                    'te_length': te_info['length'],
                    'locus': f"{chrom}:{te_info['start']}-{te_info['end']}"
                })

        # 3. Calculate TE overlap length within 200bp downstream of TSS and 200bp upstream of TES (considering strand)
        # Determine transcript start and end coordinates (all exons)
        last_exon = exons[-1]
        transcript_start = min(e['start'] for e in exons)
        transcript_end = max(e['end'] for e in exons)

        # Determine TSS and TES windows based on strand
        if strand == '+':
            tss = transcript_start
            tss_win_start = tss
            tss_win_end = min(tss + 200 - 1, first_exon['end'])

            tes = transcript_end
            tes_win_start = max(last_exon['start'], tes - 200 + 1)
            tes_win_end = tes
        else:
            # '-' strand: start site is the highest coordinate, "downstream" means decreasing coordinate
            tss = transcript_end
            tss_win_start = max(first_exon['start'], tss - 200 + 1)
            tss_win_end = tss

            tes = transcript_start
            tes_win_start = tes
            tes_win_end = min(tes + 200 - 1, last_exon['end'])

        # Calculate deduplicated overlap length within a window (query only the same chromosome and strand)
        def window_overlap_length(win_start, win_end):
            if chrom not in te_intervals or strand not in te_intervals[chrom]:
                return 0
            overlapping = list(te_intervals[chrom][strand].overlap(win_start, win_end))
            if not overlapping:
                return 0
            merged = IntervalTree(overlapping)
            merged.merge_overlaps()
            total = 0
            for mi in merged:
                inter_start = max(win_start, mi.begin)
                inter_end = min(win_end, mi.end)
                if inter_end >= inter_start:
                    total += (inter_end - inter_start + 1)
            return total

        tss_200_overlap = window_overlap_length(tss_win_start, tss_win_end)
        tes_200_overlap = window_overlap_length(tes_win_start, tes_win_end)

        # Store results
        if te_overlaps_details:
            first_exon_overlaps[transcript_id] = {
                'first_exon_length': first_exon['length'],
                'te_overlaps': te_overlaps_details, # original list, for finding dominant TE
                'total_unique_overlap': total_unique_overlap_length, # corrected total length
                'tss_200_overlap': tss_200_overlap,  # overlap length between TSS 200bp downstream window and TE
                'tes_200_overlap': tes_200_overlap   # overlap length between TES 200bp upstream window and TE
            }
    
    print(f"Completed: {len(first_exon_overlaps)} transcripts contain TE")
    return first_exon_overlaps

def calculate_enhanced_all_exon_TE_overlaps(transcript_exons, te_intervals):
    """Calculate enhanced all-exon TE overlaps, including detailed information (fixed overlap calculation)"""
    print("Calculating enhanced all-exon TE overlaps (fixed)...")
    
    all_exon_overlaps = {}
    
    for i, (transcript_id, exons) in enumerate(transcript_exons.items(), 1):
        if i % 5000 == 0:
            print(f"Processed {i}/{len(transcript_exons)}")
        
        if not exons:
            continue
        
        chrom = exons[0]['chrom']
        strand = exons[0]['strand']
        
        if chrom not in te_intervals or strand not in te_intervals[chrom]:
            continue
        
        total_exon_length = sum(e['length'] for e in exons)
        
        # --- Collect detailed information of all TEs overlapping exons and raw Interval objects ---
        raw_overlapping_tes = []
        te_overlap_dict = defaultdict(int)
        te_detail_dict = defaultdict(lambda: {
            'overlap_length': 0, 'family': 'Unknown', 'subfamily': 'Unknown', 
            'locus': 'Unknown', 'te_length': 0, 'te_start': 0, 'te_end': 0
        })

        for exon in exons:
            overlapping_tes_per_exon = te_intervals[chrom][strand].overlap(exon['start'], exon['end'])
            
            for te_interval in overlapping_tes_per_exon:
                overlap_start = max(exon['start'], te_interval.begin)
                overlap_end = min(exon['end'], te_interval.end)
                overlap_length = overlap_end - overlap_start +1
                
                if overlap_length > 0:
                    raw_overlapping_tes.append(te_interval)
                    
                    # Accumulate information for finding dominant TE etc. (this sum is not used for final proportion)
                    te_info = te_interval.data
                    te_key = (te_info['repeat_name'], te_info['class'])
                    te_overlap_dict[te_key] += overlap_length
                    
                    te_detail_dict[te_key]['overlap_length'] += overlap_length
                    te_detail_dict[te_key]['family'] = te_info.get('family', 'Unknown')
                    te_detail_dict[te_key]['subfamily'] = te_info.get('subfamily', 'Unknown')
                    te_detail_dict[te_key]['locus'] = f"{chrom}:{te_info['start']}-{te_info['end']}"
                    te_detail_dict[te_key]['te_length'] = te_info['length']
                    te_detail_dict[te_key]['te_start'] = te_info['start']
                    te_detail_dict[te_key]['te_end'] = te_info['end']
        
        if not raw_overlapping_tes:
            continue

        # --- Calculate deduplicated total overlap length ---
        # 1. Merge all overlapping TE intervals
        merged_tes = IntervalTree(raw_overlapping_tes)
        merged_tes.merge_overlaps()
        
        # 2. Calculate total length of intersections between merged TE intervals and all exons
        total_unique_overlap_length = 0
        # Create an IntervalTree for exons for efficient query
        exon_tree = IntervalTree.from_tuples([(e['start'], e['end']) for e in exons if e['start'] < e['end']])

        # Iterate over merged, non-overlapping TE intervals
        for te_interval in merged_tes:
            # Find all exons overlapping this TE interval
            overlapping_exons = exon_tree.overlap(te_interval.begin, te_interval.end)
            for exon_interval in overlapping_exons:
                # Calculate intersection and accumulate its length
                intersection_start = max(te_interval.begin, exon_interval.begin)
                intersection_end = min(te_interval.end, exon_interval.end)
                total_unique_overlap_length += (intersection_end - intersection_start +1)
        
        # --- Store results ---
        te_classes = list(set(te_key[1] for te_key in te_overlap_dict.keys()))
        
        all_exon_overlaps[transcript_id] = {
            'total_exon_length': total_exon_length,
            'total_overlap_length': total_unique_overlap_length, # use corrected value
            'num_tes': len(te_overlap_dict),
            'te_classes': te_classes,
            'te_overlap_details': dict(te_detail_dict) # keep original accumulated information
        }
            
    print("Enhanced all-exon analysis completed")
    return all_exon_overlaps

def find_dominant_te_in_transcript(transcript_id, all_exon_overlaps):
    """Find the TE with the largest overlap proportion in the transcript"""
    if transcript_id not in all_exon_overlaps:
        return None
    
    te_overlap_details = all_exon_overlaps[transcript_id].get('te_overlap_details', {})
    if not te_overlap_details:
        return None
    
    # Find TE with maximum overlap length
    max_overlap_te = max(te_overlap_details.items(), key=lambda x: x[1]['overlap_length'])
    te_key, te_data = max_overlap_te
    
    repeat_name, te_class = te_key
    
    return {
        'repeat_name': repeat_name,
        'class': te_class,
        'family': te_data.get('family', 'Unknown'),
        'subfamily': te_data.get('subfamily', 'Unknown'),
        'locus': te_data.get('locus', 'Unknown'),
        'overlap_length': te_data['overlap_length'],
        'te_length': te_data.get('te_length', 0)
    }

def find_dominant_te_in_first_exon(transcript_id, first_exon_overlaps):
    """Find the TE with the largest overlap proportion in the first exon, including proportion relative to TE total length"""
    if transcript_id not in first_exon_overlaps:
        return None
    
    te_overlaps = first_exon_overlaps[transcript_id].get('te_overlaps', [])
    if not te_overlaps:
        return None
    
    # Find TE with maximum overlap length
    dominant_te = max(te_overlaps, key=lambda x: x['overlap_length'])
    
    # Calculate proportion of overlap length relative to TE total length
    te_total_length = dominant_te.get('te_length', 0)
    overlap_length = dominant_te['overlap_length']
    
    # Calculate overlap proportion relative to TE total length (percentage)
    te_overlap_proportion = (overlap_length / te_total_length * 100) if te_total_length > 0 else 0.0
    
    # Create enhanced result
    result = dict(dominant_te)  # copy original data
    result['te_overlap_proportion'] = te_overlap_proportion  # new: proportion of TE length overlapped
    
    return result

def run_miniquant(gtf_path, long_sam, short_sam, output_dir, **kwargs):
    """Run miniQuant - fixed version"""
    print("\nRunning miniQuant...")
    
    # Convert to absolute paths
    gtf_path = os.path.abspath(gtf_path)
    long_sam = os.path.abspath(long_sam)
    output_dir = os.path.abspath(output_dir)
    if short_sam:
        short_sam = os.path.abspath(short_sam)
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Validate input files
    print("\nValidating input files...")
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")
    if not os.path.exists(long_sam):
        raise FileNotFoundError(f"Long-read SAM file not found: {long_sam}")
    if short_sam and not os.path.exists(short_sam):
        raise FileNotFoundError(f"Short-read SAM file not found: {short_sam}")
    
    print(f"GTF: {gtf_path}")
    print(f"Long-read SAM: {long_sam}")
    if short_sam:
        print(f"Short-read SAM: {short_sam}")
    
    # Build correct miniQuant command format
    # python main.py quantify -gtf <gtf> -lrsam <long_sam> -o <output> [options]
    cmd = [
        'python', MINIQUANT_MAIN,
        'quantify',  # key: add subcommand
        '-gtf', gtf_path,
        '-lrsam', long_sam,
        '-o', output_dir
    ]
    
    # Add short-read file
    if short_sam:
        cmd.extend(['-srsam', short_sam])
    
    # Add other parameters
    alpha = str(kwargs.get('alpha', '1.0'))
    em_choice = kwargs.get('EM_choice', 'LR')
    threads = str(kwargs.get('threads', 4))
    
    cmd.extend(['--alpha', alpha])
    cmd.extend(['--EM_choice', em_choice])
    cmd.extend(['--threads', threads])
    
    # Handle pre-trained model path
    if alpha.lower() == 'adaptive':
        pretrained_model = kwargs.get('pretrained_model_path', 'cDNA-ONT')
        cmd.extend(['--pretrained_model_path', pretrained_model])
    
    # Add other optional parameters
    if 'filtering' in kwargs:
        cmd.extend(['--filtering', str(kwargs['filtering'])])
    if 'multi_mapping_filtering' in kwargs:
        cmd.extend(['--multi_mapping_filtering', str(kwargs['multi_mapping_filtering'])])
    
    print(f"\nExecuting command: {' '.join(cmd)}")
    
    try:
        print("Starting miniQuant...")
        print("Note: For large datasets, this may take several hours. Please be patient.")
        print("The program will run until completion, no timeout limit.")
        print("To cancel, press Ctrl+C\n")
        
        # Remove timeout parameter to allow unlimited runtime
        result = subprocess.run(cmd, cwd=output_dir, capture_output=True, text=True)
        
        print("\n--- miniQuant execution completed ---")
        if result.stdout:
            print("STDOUT (last 2000 chars):")
            print(result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout)
        
        if result.stderr:
            print("STDERR (last 2000 chars):")
            print(result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr)
        print("--- End of output ---")
        
        if result.returncode != 0:
            print(f"miniQuant execution failed (return code: {result.returncode})")
            return False
        else:
            print("miniQuant executed successfully")
            
            # Check output files
            print("\nChecking output files...")
            if os.path.exists(output_dir):
                files = [f for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))]
                if files:
                    print(f"Generated {len(files)} files:")
                    for f in sorted(files)[:10]:
                        fp = os.path.join(output_dir, f)
                        size = os.path.getsize(fp)
                        print(f"  {f} ({size:,} bytes)")
                    if len(files) > 10:
                        print(f"  ... and {len(files)-10} more files")
                else:
                    print("Output directory is empty")
            
            return True
            
    except KeyboardInterrupt:
        print("\nUser interrupted execution")
        return False
    except Exception as e:
        print(f"Error while running miniQuant: {e}")
        return False

def classify_transcript(all_exon_overlap_length, transcript_length, tss_200_overlap, tes_200_overlap, 
                       total_te_proportion, te_overlap_threshold=10, te_ratio_threshold=80.0, 
                       te_feature_threshold=50):
    """
    Classify transcript
    
    Parameters:
    - all_exon_overlap_length: total overlap length between entire transcript and TE (bp)
    - transcript_length: full transcript length (bp)
    - tss_200_overlap: overlap length between TSS 200bp window and TE (bp)
    - tes_200_overlap: overlap length between TES 200bp window and TE (bp)
    - total_te_proportion: TE proportion in transcript (%)
    - te_overlap_threshold: overlap length threshold (bp), default 10bp
    - te_ratio_threshold: TE proportion in transcript threshold (%), default 80%
    - te_feature_threshold: TSS/TES window overlap length threshold (bp), default 50bp
    
    Returns:
    - classification: classification result (str)
    - sub_classifications: list of sub-classifications (list)
    """
    
    # 1. Determine if it is Gene-alone
    if all_exon_overlap_length <= te_overlap_threshold:
        return 'Gene-alone', []
    
    # 2. Determine if it is TE-alone
    te_ratio = (all_exon_overlap_length / transcript_length * 100) if transcript_length > 0 else 0
    if te_ratio > te_ratio_threshold and tss_200_overlap > te_feature_threshold:
        return 'TE-alone', []
    
    # 3. TE-Gene classification (includes three subclasses)
    sub_classifications = []
    
    if tss_200_overlap > te_feature_threshold:
        sub_classifications.append('TSS')
    
    if tes_200_overlap > te_feature_threshold:
        sub_classifications.append('TE-terminating')
    
    if tss_200_overlap <= te_feature_threshold and tes_200_overlap <= te_feature_threshold:
        sub_classifications.append('TE-exonization')
    
    return 'TE-Gene', sub_classifications

def generate_enhanced_integrated_output(transcript_output, first_exon_overlaps, all_exon_overlaps, 
                                       output_path, transcript_info, first_exon_threshold=50.0, 
                                       total_threshold=50.0, skip_quantification=False, 
                                       te_overlap_threshold=10, te_ratio_threshold=80.0,
                                       te_feature_threshold=50):
    """Generate enhanced integrated output
    
    New parameters:
    - te_overlap_threshold: overlap length threshold between transcript and TE, used to distinguish Gene-alone from others (bp)
    - te_ratio_threshold: TE proportion in full transcript threshold, used to identify TE-alone (%)
    - te_feature_threshold: overlap length threshold between TSS/TES 200bp windows and TE (bp)
    """
    if skip_quantification:
        print("Generating enhanced integrated output (fast mode - no TPM values)...")
    else:
        print("Generating enhanced integrated output...")
    
    # Read transcript quantification results (or generate dummy data)
    transcript_abundance = None
    
    if skip_quantification:
        print("Using dummy TPM values for analysis...")
        # Create dummy transcript abundance data
        transcript_ids = list(transcript_info.keys())
        transcript_abundance = pd.DataFrame({
            'Isoform': transcript_ids,
            'TPM': [0.0] * len(transcript_ids)  # all transcripts use same dummy TPM
        })
        print(f"Created dummy quantification data: {len(transcript_abundance)} transcripts")
    else:
        # Original reading logic
        result_files = ["Isoform_abundance.out", "LR_EM_expression.out", "EM_expression.out", 
                       "LR_EM_expression.out", "SR_EM_expression.out"]
        
        found_file = None
        for filename in result_files:
            filepath = f"{transcript_output}/{filename}"
            if os.path.exists(filepath):
                try:
                    transcript_abundance = pd.read_csv(filepath, sep='\t')
                    found_file = filename
                    print(f"Read: {filename} ({len(transcript_abundance)} rows)")
                    break
                except Exception as e:
                    print(f"Failed to read {filename}: {e}")
        
        if transcript_abundance is None:
            print("\nNo transcript quantification result file found")
            print(f"Searched files: {result_files}")
            print(f"In directory: {transcript_output}")
            return
    
    # Determine column names
    id_col = None
    for col in ['Isoform', 'transcript_id', 'isoform_id']:
        if col in transcript_abundance.columns:
            id_col = col
            break
    
    if not id_col:
        print(f"No transcript ID column found")
        print(f"Available columns: {list(transcript_abundance.columns)}")
        return
    
    tpm_col = None
    for col in ['TPM', 'tpm', 'Tpm']:
        if col in transcript_abundance.columns:
            tpm_col = col
            break
    
    # Generate results
    all_results = []
    te_initiated_results = []
    
    for _, row in transcript_abundance.iterrows():
        transcript_id = row[id_col]
        tpm = row[tpm_col] if tpm_col else 0
        gene_id = transcript_info.get(transcript_id, {}).get('gene_id', '')
        transcript_length = transcript_info.get(transcript_id, {}).get('length', 0)
        
        result_row = {
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'transcript_TPM': tpm,
            'transcript_length': transcript_length
        }
        
        first_exon_te_proportion = 0
        total_te_proportion = 0
        num_tes = 0
        te_classes_str = ''
        
        # First calculate first exon proportion
        if transcript_id in first_exon_overlaps:
            fe_data = first_exon_overlaps[transcript_id]
            # Use corrected, deduplicated total overlap length
            total_unique_fe_overlap = fe_data['total_unique_overlap']
            first_exon_te_proportion = (
                total_unique_fe_overlap / fe_data['first_exon_length'] * 100
            ) if fe_data['first_exon_length'] > 0 else 0
        
        # Then calculate whole transcript proportion + collect detailed information for each TE
        loci = []
        subfamilies = []
        families = []
        classes = []
        lengths = []
        overlap_lengths = []  # new: overlap length per TE
        
        if transcript_id in all_exon_overlaps:
            ae_data = all_exon_overlaps[transcript_id]
            total_te_proportion = (
                ae_data['total_overlap_length'] / ae_data['total_exon_length'] * 100
            ) if ae_data['total_exon_length'] > 0 else 0
            num_tes = ae_data['num_tes']
            te_classes_str = ';'.join(ae_data['te_classes'])
            
            # Extract information of all TEs from te_overlap_details
            te_details = ae_data.get('te_overlap_details', {})
            for te_key, te_data in te_details.items():
                # te_key = (repeat_name, class)
                _, te_class = te_key
                loci.append(te_data.get('locus', 'Unknown'))
                subfamilies.append(te_data.get('subfamily', 'Unknown'))
                families.append(te_data.get('family', 'Unknown'))
                classes.append(te_class)
                lengths.append(str(te_data.get('te_length', 0)))
                overlap_lengths.append(str(te_data.get('overlap_length', 0)))  # new: overlap length
        
        # Get TSS/TES 200bp window overlap with TE
        tss_200_overlap = 0
        tes_200_overlap = 0
        if transcript_id in first_exon_overlaps:
            tss_200_overlap = first_exon_overlaps[transcript_id]['tss_200_overlap']
            tes_200_overlap = first_exon_overlaps[transcript_id]['tes_200_overlap']
        
        # Get total overlap length between whole transcript and TE, for classification
        all_exon_overlap_length = 0
        if transcript_id in all_exon_overlaps:
            all_exon_overlap_length = all_exon_overlaps[transcript_id]['total_overlap_length']
        
        # Classify transcript
        transcript_classification, sub_classifications = classify_transcript(
            all_exon_overlap_length=all_exon_overlap_length,
            transcript_length=transcript_length,
            tss_200_overlap=tss_200_overlap,
            tes_200_overlap=tes_200_overlap,
            total_te_proportion=total_te_proportion,
            te_overlap_threshold=te_overlap_threshold,
            te_ratio_threshold=te_ratio_threshold,
            te_feature_threshold=te_feature_threshold
        )
        
        result_row.update({
            'num_TEs': num_tes,
            'TE_class': te_classes_str,
            'first_exon_TE_proportion': first_exon_te_proportion,
            'total_transcript_TE_proportion': total_te_proportion,
            'tss_200_overlap': tss_200_overlap,
            'tes_200_overlap': tes_200_overlap,
            'transcript_classification': transcript_classification,
            'transcript_sub_classification': ','.join(sub_classifications) if sub_classifications else '',
            # new 5 columns: multiple TEs separated by semicolons
            'locus': ';'.join(loci) if loci else '',
            'subfamily': ';'.join(subfamilies) if subfamilies else '',
            'family': ';'.join(families) if families else '',
            'class': ';'.join(classes) if classes else '',
            'length': ';'.join(lengths) if lengths else '',
            'overlap_length': ';'.join(overlap_lengths) if overlap_lengths else ''  # new: overlap length per TE
        })
        
        all_results.append(result_row)

        # TE-initiated transcripts
        meets_first_exon = (transcript_id in first_exon_overlaps and first_exon_te_proportion >= first_exon_threshold)
        meets_total = (transcript_id in all_exon_overlaps and total_te_proportion >= total_threshold)
        
        if meets_first_exon or meets_total:
            # Get dominant TE information
            dominant_te_all = find_dominant_te_in_transcript(transcript_id, all_exon_overlaps)
            dominant_te_first_exon = find_dominant_te_in_first_exon(transcript_id, first_exon_overlaps)
            
            # Construct TE-initiated transcript row
            te_initiated_row = {
                'transcript_id': transcript_id,
                'gene_id': gene_id,
                'transcript_TPM': tpm,
                'transcript_length': transcript_length,
                'first_exon_TE_proportion': first_exon_te_proportion,
                'total_transcript_TE_proportion': total_te_proportion,
                'tss_200_overlap': tss_200_overlap,
                'tes_200_overlap': tes_200_overlap,
                'transcript_classification': transcript_classification,
                'transcript_sub_classification': ','.join(sub_classifications) if sub_classifications else '',
                'meets_first_exon_threshold': meets_first_exon,
                'meets_total_threshold': meets_total,
                
                # Original initiating TE information (TE with longest overlap in first exon)
                'TE': dominant_te_first_exon['repeat_name'] if dominant_te_first_exon else 'Unknown',
                'TE_subfamily': dominant_te_first_exon['subfamily'] if dominant_te_first_exon else 'Unknown',
                'TE_family': dominant_te_first_exon['family'] if dominant_te_first_exon else 'Unknown',
                'TE_class': dominant_te_first_exon['class'] if dominant_te_first_exon else 'Unknown',
                
                # New: dominant TE in whole transcript
                'dominant_TE_all_repeat_name': dominant_te_all['repeat_name'] if dominant_te_all else 'Unknown',
                'dominant_TE_all_class': dominant_te_all['class'] if dominant_te_all else 'Unknown',
                'dominant_TE_all_family': dominant_te_all['family'] if dominant_te_all else 'Unknown',
                'dominant_TE_all_subfamily': dominant_te_all['subfamily'] if dominant_te_all else 'Unknown',
                'dominant_TE_all_locus': dominant_te_all['locus'] if dominant_te_all else 'Unknown',
                'dominant_TE_all_length': dominant_te_all['te_length'] if dominant_te_all else 0,
                'dominant_TE_all_overlap_length': dominant_te_all['overlap_length'] if dominant_te_all else 0,
                'dominant_TE_all_transcript_proportion': (dominant_te_all['overlap_length'] / transcript_length * 100) if dominant_te_all and transcript_length > 0 else 0,
                
                # New: dominant TE in first exon detailed information
                'dominant_TE_first_exon_repeat_name': dominant_te_first_exon['repeat_name'] if dominant_te_first_exon else 'Unknown',
                'dominant_TE_first_exon_class': dominant_te_first_exon['class'] if dominant_te_first_exon else 'Unknown',
                'dominant_TE_first_exon_family': dominant_te_first_exon['family'] if dominant_te_first_exon else 'Unknown',
                'dominant_TE_first_exon_subfamily': dominant_te_first_exon['subfamily'] if dominant_te_first_exon else 'Unknown',
                'dominant_TE_first_exon_locus': dominant_te_first_exon['locus'] if dominant_te_first_exon else 'Unknown',
                'dominant_TE_first_exon_te_length': dominant_te_first_exon['te_length'] if dominant_te_first_exon else 0,
                'dominant_TE_first_exon_overlap_length': dominant_te_first_exon['overlap_length'] if dominant_te_first_exon else 0,
                'dominant_TE_first_exon_te_overlap_proportion': dominant_te_first_exon.get('te_overlap_proportion', 0) if dominant_te_first_exon else 0,  # new: overlap proportion relative to TE total length
            }
            
            # Add first exon length and proportion information, including TSS/TES 200bp window overlaps
            if transcript_id in first_exon_overlaps:
                first_exon_length = first_exon_overlaps[transcript_id]['first_exon_length']
                te_initiated_row['first_exon_length'] = first_exon_length
                te_initiated_row['dominant_TE_first_exon_proportion'] = (dominant_te_first_exon['overlap_length'] / first_exon_length * 100) if dominant_te_first_exon and first_exon_length > 0 else 0
                # New: TSS/TES 200bp window overlaps
                #te_initiated_row['tss_200_overlap'] = first_exon_overlaps[transcript_id]['tss_200_overlap']
                #te_initiated_row['tes_200_overlap'] = first_exon_overlaps[transcript_id]['tes_200_overlap']
            else:
                te_initiated_row['first_exon_length'] = 0
                te_initiated_row['dominant_TE_first_exon_proportion'] = 0
                #te_initiated_row['tss_200_overlap'] = 0
                #te_initiated_row['tes_200_overlap'] = 0
            
            te_initiated_results.append(te_initiated_row)
    
    # Save results
    all_df = pd.DataFrame(all_results)
    output_file1 = f"{output_path}/transcript_quantification_with_TE.tsv"
    all_df.to_csv(output_file1, sep='\t', index=False)
    print(f"Saved all transcripts: {output_file1} ({len(all_results)} rows)")
    
    # Output classification information
    print("\nTranscript classification parameters:")
    print(f"  - te_overlap_threshold: {te_overlap_threshold} bp (Gene-alone threshold)")
    print(f"  - te_ratio_threshold: {te_ratio_threshold}% (TE-alone threshold)")
    print(f"  - te_feature_threshold: {te_feature_threshold} bp (TSS/TES window threshold)")
    print("\nTranscript classification results:")
    if 'transcript_classification' in all_df.columns:
        class_counts = all_df['transcript_classification'].value_counts()
        for class_name, count in class_counts.items():
            print(f"  - {class_name}: {count}")
    
    if te_initiated_results:
        te_df = pd.DataFrame(te_initiated_results)
        output_file2 = f"{output_path}/TE_derived_transcripts.tsv"
        te_df.to_csv(output_file2, sep='\t', index=False)
        print(f"Saved enhanced TE-derived transcripts: {output_file2} ({len(te_initiated_results)} rows)")
        
        # Output column descriptions
        print("\nEnhanced TE_derived_transcripts.tsv contains the following new columns:")
        print("Dominant TE information (whole transcript):")
        print("  - dominant_TE_all_repeat_name: TE repeat name")
        print("  - dominant_TE_all_class: TE class")
        print("  - dominant_TE_all_family: TE family")
        print("  - dominant_TE_all_subfamily: TE subfamily")
        print("  - dominant_TE_all_locus: TE locus position")
        print("  - dominant_TE_all_length: TE length")
        print("  - dominant_TE_all_overlap_length: overlap length")
        print("  - dominant_TE_all_transcript_proportion: proportion in full transcript")
        print("Dominant TE information (first exon):")
        print("  - dominant_TE_first_exon_*: detailed information of the dominant TE in first exon")
        print("  - dominant_TE_first_exon_proportion: proportion in first exon")
        
    else:
        print("No TE-initiated transcripts found meeting thresholds")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Long-read Transposable Element Quantification')
    parser.add_argument('-g', '--gtf', required=True, help='Transcript GTF annotation file')
    parser.add_argument('-t', '--te-gtf', required=True, help='TE GTF annotation file')
    parser.add_argument('-l', '--long-sam', required=True, help='Long-read SAM file')
    parser.add_argument('-s', '--short-sam', help='Short-read SAM file (optional)')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    
    # miniQuant parameters - support string parameters like adaptive
    parser.add_argument('--alpha', default='1.0', help='Mixing parameter (supports adaptive)')
    parser.add_argument('--EM_choice', default='LR', help='EM algorithm choice')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    
    # TE threshold parameters
    parser.add_argument('--first-exon-threshold', type=float, default=0.0, 
                       help='First exon TE proportion threshold (%)')
    parser.add_argument('--total-threshold', type=float, default=0.0, 
                       help='Full transcript TE proportion threshold (%)')
    
    # Transcript classification parameters
    parser.add_argument('--te-overlap-threshold', type=int, default=10,
                       help='Overlap length threshold between transcript and TE, for Gene-alone distinction (bp)')
    parser.add_argument('--te-ratio-threshold', type=float, default=80.0,
                       help='TE proportion in full transcript threshold for TE-alone identification (%)')
    parser.add_argument('--te-feature-threshold', type=int, default=50,
                       help='TSS/TES 200bp window overlap with TE threshold (bp)')
    
    # New: skip quantification switch
    parser.add_argument('--skip-quantification', action='store_true',
                       help='Skip transcript quantification, only perform TE overlap analysis (fast mode)')
    
    # Other miniQuant parameters
    parser.add_argument('--pretrained_model_path', default='cDNA-ONT')
    parser.add_argument('--filtering', type=bool, default=False)
    parser.add_argument('--multi_mapping_filtering', type=bool, default=False)
    
    args = parser.parse_args()
    
    print("=" * 60)
    if args.skip_quantification:
        print("Integrated transcript quantification and TE detection - fast mode (overlap analysis only)")
    else:
        print("Integrated transcript quantification and TE detection - full mode")
    print("=" * 60)
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    transcript_output = output_dir / "transcript_quantification_temp"
    
    try:
        # 1. Parse annotation files
        print("\nStep 1: Parsing transcript annotation...")
        transcript_exons, transcript_info = parse_transcript_annotation(args.gtf)
        
        print("\nStep 2: Parsing TE annotation...")
        te_intervals = parse_te_annotation(args.te_gtf)
        
        # 2. Run miniQuant (optional)
        if args.skip_quantification:
            print("\nStep 3: Skipping transcript quantification (fast mode)")
            print("Will use dummy TPM values for analysis...")
            transcript_quantification_success = True
        else:
            print("\nStep 3: Running transcript quantification...")
            miniquant_kwargs = {
                'alpha': args.alpha,
                'EM_choice': args.EM_choice,
                'threads': args.threads,
                'pretrained_model_path': args.pretrained_model_path,
                'filtering': args.filtering,
                'multi_mapping_filtering': args.multi_mapping_filtering
            }
            
            transcript_quantification_success = run_miniquant(
                args.gtf, args.long_sam, args.short_sam, 
                str(transcript_output), **miniquant_kwargs
            )
            
            if not transcript_quantification_success:
                print("Transcript quantification failed")
                return
        
        # 3. TE overlap detection
        print("\nStep 4: Detecting TE overlaps...")
        first_exon_overlaps = find_first_exon_TE_overlaps(transcript_exons, te_intervals)
        
        print("\nStep 5: Calculating enhanced TE overlaps...")
        all_exon_overlaps = calculate_enhanced_all_exon_TE_overlaps(transcript_exons, te_intervals)
        
        # 4. Generate integrated output
        print("\nStep 6: Generating enhanced integrated output...")
        generate_enhanced_integrated_output(
            str(transcript_output), first_exon_overlaps, all_exon_overlaps,
            str(output_dir), transcript_info, 
            args.first_exon_threshold, args.total_threshold,
            skip_quantification=args.skip_quantification,
            te_overlap_threshold=args.te_overlap_threshold,
            te_ratio_threshold=args.te_ratio_threshold,
            te_feature_threshold=args.te_feature_threshold
        )
        
        print("\n" + "=" * 60)
        if args.skip_quantification:
            print("Fast analysis completed!")
            print("=" * 60)
            print(f"Results saved in: {output_dir}")
            print("Main output files:")
            print("  1. transcript_quantification_with_TE.tsv - TE analysis for all transcripts (dummy TPM)")
            print("  2. TE_initiated_transcripts.tsv - enhanced TE-initiated transcript details (dummy TPM)")
            print("Note: Fast mode uses dummy TPM values. For real quantification, run full mode.")
        else:
            print("Analysis completed!")
            print("=" * 60)
            print(f"Results saved in: {output_dir}")
            print("Main output files:")
            print("  1. transcript_quantification_with_TE.tsv - quantification results for all transcripts")
            print("  2. TE_initiated_transcripts.tsv - enhanced TE-initiated transcript details")
        
    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()