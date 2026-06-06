import pandas as pd
import re
import argparse
import numpy as np
import os

def cal_TE_exp_add_thres(transcript_TE_annotation, transcript_quantification_with_TE, output_dir, percent_threshold=0.5):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read transcript_quantification_with_TE.tsv to get transcript TPM and TE mapping
    quant_df = pd.read_csv(transcript_quantification_with_TE, sep='\t')
    
    # Extract unique transcript TPM
    transcript_tpm = {}
    TE_locus_to_transcripts = {}
    for _, row in quant_df.iterrows():
        transcript_id = row['transcript_id']
        tpm = float(row['transcript_TPM'])
        if transcript_id not in transcript_tpm:
            transcript_tpm[transcript_id] = tpm
        
        locus_str = row.get('locus', '')
        if pd.isna(locus_str):
            locus_str = ''
        if locus_str:
            locus_list = locus_str.split(';')
            length_str = row.get('length', '')
            if pd.isna(length_str):
                length_str = ''
            length_list = length_str.split(';') if length_str else []
            overlap_str = row.get('overlap_length', '')
            if pd.isna(overlap_str):
                overlap_str = ''
            overlap_list = overlap_str.split(';') if overlap_str else []
            for i, locus in enumerate(locus_list):
                if i < len(overlap_list) and i < len(length_list):
                    try:
                        overlap = float(overlap_list[i])
                        length = float(length_list[i])
                        prop = overlap / length if length > 0 else 0
                    except ValueError:
                        prop = 0
                else:
                    prop = 0
                if locus not in TE_locus_to_transcripts:
                    TE_locus_to_transcripts[locus] = []
                TE_locus_to_transcripts[locus].append((transcript_id, prop))

    # Read TE annotation file (GTF format)
    TE_info = {}
    with open(transcript_TE_annotation, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            attributes = parts[8]
            attr_dict = {}
            for attr in attributes.split(';'):
                if attr.strip():
                    key_value = attr.strip().split(' ', 1)
                    if len(key_value) == 2:
                        key = key_value[0]
                        value = key_value[1].strip('"')
                        attr_dict[key] = value
            locus = attr_dict.get('locus', '')
            subfamily = attr_dict.get('gene_id', '')  # Assuming family_id is subfamily
            family = attr_dict.get('family_id', '')
            te_class = attr_dict.get('class_id', '')
            length = int(attr_dict.get('consensus_length', 0))
            transcript_id_te = attr_dict.get('transcript_id', '')
            TE_info[locus] = {
                'subfamily': subfamily,
                'family': family,
                'class': te_class,
                'length': length,
                'transcript_id': transcript_id_te
            }

    # Calculate TE loci expression
    TE_exp = []
    for locus, info in TE_info.items():
        te_exp_row = {
            'TE_loci': info['transcript_id'],
            'TE_subfamily': info['subfamily'],
            'TE_family': info['family'],
            'TE_class': info['class'],
            'TE_length': info['length']
        }
        transcripts = TE_locus_to_transcripts.get(locus, [])
        total_exp = 0.0
        for transcript, prop in transcripts:
            if transcript in transcript_tpm and prop >= percent_threshold:
                total_exp += transcript_tpm[transcript]
        te_exp_row['TE_expression'] = total_exp
        TE_exp.append(te_exp_row)

    # Write TE loci exp
    te_loci_df = pd.DataFrame(TE_exp)
    te_loci_df.to_csv(os.path.join(output_dir, 'TE_loci_exp_TPM.tsv'), sep='\t', index=False)

    # Aggregate to subfamily
    subfamily_info = {}
    subfamily_transcripts = {}
    for locus, info in TE_info.items():
        subfam = info['subfamily']
        if subfam not in subfamily_info:
            subfamily_info[subfam] = {'TE_family': info['family'], 'TE_class': info['class']}
        transcripts = TE_locus_to_transcripts.get(locus, [])
        for transcript, prop in transcripts:
            if transcript in transcript_tpm and prop >= percent_threshold:
                subfamily_transcripts.setdefault(subfam, set()).add(transcript)

    subfamily_exp = {}
    for subfam, transcripts in subfamily_transcripts.items():
        fam_cls = subfamily_info.get(subfam, {})
        total_exp = sum(transcript_tpm[transcript] for transcript in transcripts)
        subfamily_exp[subfam] = {
            'TE_family': fam_cls.get('TE_family', ''),
            'TE_class': fam_cls.get('TE_class', ''),
            'TE_expression': total_exp
        }

    subfam_df = pd.DataFrame(list(subfamily_exp.values()))
    subfam_df.insert(0, 'TE_subfamily', list(subfamily_exp.keys()))
    subfam_df.to_csv(os.path.join(output_dir, 'TE_subfamily_exp_TPM.tsv'), sep='\t', index=False)

    # Aggregate to family
    family_exp = {}
    for subfam, data in subfamily_exp.items():
        fam = data['TE_family']
        cls = data['TE_class']
        if fam not in family_exp:
            family_exp[fam] = {'TE_class': cls, 'TE_expression': 0.0}
        family_exp[fam]['TE_expression'] += data['TE_expression']

    fam_df = pd.DataFrame(list(family_exp.values()))
    fam_df.insert(0, 'TE_family', list(family_exp.keys()))
    fam_df.to_csv(os.path.join(output_dir, 'TE_family_exp_TPM.tsv'), sep='\t', index=False)

    # Aggregate to class
    class_exp = {}
    for fam, data in family_exp.items():
        cls = data['TE_class']
        if cls not in class_exp:
            class_exp[cls] = {'TE_expression': 0.0}
        class_exp[cls]['TE_expression'] += data['TE_expression']

    class_df = pd.DataFrame(list(class_exp.values()))
    class_df.insert(0, 'TE_class', list(class_exp.keys()))
    class_df.to_csv(os.path.join(output_dir, 'TE_class_exp_TPM.tsv'), sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate TE expression from transcript TPM data.")
    parser.add_argument("-t", "--te_anno", help="Path to TE annotation file (GTF format)")
    parser.add_argument("-i", "--quant", help="Path to transcript quantification with TE mapping (TSV format)")
    parser.add_argument("-o", "--out_dir", help="Output directory for TE expression files")
    parser.add_argument("-p", "--thresh", type=float, help="Percentage threshold for overlap proportion (e.g., 0.01 for 1%)")

    args = parser.parse_args()

    cal_TE_exp_add_thres(
        args.te_anno,
        args.quant,
        args.out_dir,
        args.thresh
    )
