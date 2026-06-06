import argparse
from TrEESR import TrEESR
from TransELS import TransELS
from EM import EM,EM_SR,EM_hybrid
import config
import os
import sys
# import os
# os.system("taskset -p 0xfffff %d" % os.getpid())
# affinity_mask = os.sched_getaffinity(0)
# os.sched_setaffinity(0, affinity_mask)


# ---------------------------------------------------------------------------
# cal_TE dispatcher (called before bulk EM config is set up)
# ---------------------------------------------------------------------------

def _run_cal_te(args):
    """
    Dispatch handler for the  cal_TE  subcommand.
    Called early (before the bulk EM config block) so that cal_TE does not
    require any bulk-EM arguments (READ_JUNC_MIN_MAP_LEN, normalize_sr_A, …).
    """
    # Make sure TE_analyse is importable from this directory
    _this_dir = os.path.dirname(os.path.abspath(__file__))
    _te_dir = os.path.join(_this_dir, "TE_analyse")
    for d in [_this_dir, _te_dir]:
        if d not in sys.path:
            sys.path.insert(0, d)

    from TE_analyse.cal_te_main import run_te_annotation

    te_table_path = run_te_annotation(
        gtf=args.gtf_annotation_path,
        te_gtf=args.te_gtf_path,
        output_dir=args.output_path,
        skip_quantification=args.skip_quantification,
        threads=args.threads,
        first_exon_threshold=args.first_exon_threshold,
        total_threshold=args.total_threshold,
        te_overlap_threshold=args.te_overlap_threshold,
        te_ratio_threshold=args.te_ratio_threshold,
        te_feature_threshold=args.te_feature_threshold,
    )

    # --- Optional spatial TE analysis ---
    if args.sp_quant is not None:
        from TE_analyse.te_quantification import run_sc_te_analysis

        st_te_output_dir = os.path.join(args.output_path, "ST_TE_output")
        run_sc_te_analysis(
            te_table_path=te_table_path,
            quant_dir=args.sp_quant,
            output_dir=st_te_output_dir,
            percent_threshold=args.percent_threshold,
            output_loci=args.output_loci,
            barcode_label="spot_barcode",
        )
        print(f"\n[cal_TE] Spatial TE output: {st_te_output_dir}", flush=True)

    # --- Optional single-cell TE analysis ---
    if args.sc_quant is not None:
        from TE_analyse.te_quantification import run_sc_te_analysis

        sc_te_output_dir = os.path.join(args.output_path, "SC_TE_output")
        run_sc_te_analysis(
            te_table_path=te_table_path,
            quant_dir=args.sc_quant,
            output_dir=sc_te_output_dir,
            percent_threshold=args.percent_threshold,
            output_loci=args.output_loci,
            barcode_label="cell_barcode",
        )
        print(f"\n[cal_TE] Single-cell TE output: {sc_te_output_dir}", flush=True)

    # --- Optional bulk TE quantification ---
    if args.bulk_quant is not None:
        import pandas as pd
        from TE_analyse.calculate_te_expression import cal_TE_exp_add_thres

        print("\n[cal_TE] Bulk TE quantification...", flush=True)

        # Merge real TPM values into the TE annotation table
        te_df = pd.read_csv(te_table_path, sep='\t', low_memory=False)
        quant_df = pd.read_csv(args.bulk_quant, sep='\t')

        def _detect_col(df, candidates):
            for c in candidates:
                if c in df.columns:
                    return c
            lower2col = {c.lower(): c for c in df.columns}
            for c in candidates:
                if c.lower() in lower2col:
                    return lower2col[c.lower()]
            tpm_like = [c for c in df.columns if 'tpm' in c.lower()]
            return tpm_like[0] if tpm_like else None

        quant_id_col  = _detect_col(quant_df, ['transcript_id', 'Isoform', 'isoform_id', 'target_id', 'Name'])
        quant_tpm_col = _detect_col(quant_df, ['TPM', 'tpm', 'Tpm', 'IsoTPM', 'transcript_TPM'])

        if quant_id_col is None or quant_tpm_col is None:
            print("[cal_TE] WARNING: Cannot auto-detect transcript ID or TPM column in "
                  f"--bulk_quant file. Skipping bulk TE quantification.", flush=True)
        else:
            print(f"[cal_TE] Bulk quant — ID col: {quant_id_col}, TPM col: {quant_tpm_col}", flush=True)
            quant_map = quant_df.set_index(quant_id_col)[quant_tpm_col]
            new_tpm   = te_df['transcript_id'].map(quant_map)
            mask      = new_tpm.notna()
            te_df.loc[mask, 'transcript_TPM'] = new_tpm[mask].astype(float)
            print(f"[cal_TE] TPM updated for {mask.sum()}/{len(te_df)} transcripts", flush=True)

            updated_te_table = os.path.join(args.output_path, 'transcript_quantification_with_TE_TPM.tsv')
            te_df.to_csv(updated_te_table, sep='\t', index=False)

            bulk_te_output = os.path.join(args.output_path, 'Bulk_TE_output')
            cal_TE_exp_add_thres(
                transcript_TE_annotation=args.te_gtf_path,
                transcript_quantification_with_TE=updated_te_table,
                output_dir=bulk_te_output,
                percent_threshold=args.percent_threshold,
            )
            print(f"\n[cal_TE] Bulk TE output: {bulk_te_output}", flush=True)

    print("\n[cal_TE] All done.", flush=True)


def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="Isoform quantification tools",add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help',dest="subparser_name")

    # -----------------------------------------------------------------------
    # cal_TE subcommand
    # -----------------------------------------------------------------------
    parser_te = subparsers.add_parser(
        'cal_TE',
        help='Calculate TE association for each transcript and (optionally) per-spot TE metrics'
    )
    req_te = parser_te.add_argument_group('required arguments')
    req_te.add_argument('-gtf', '--gtf_annotation_path', type=str, required=True,
                        help='Path to transcript GTF annotation file')
    req_te.add_argument('-te_gtf', '--te_gtf_path', type=str, required=True,
                        help='Path to TE GTF annotation file (from create_simple_TE_gtf.py)')
    req_te.add_argument('-o', '--output_path', type=str, required=True,
                        help='Output directory')

    opt_te = parser_te.add_argument_group('transcript TE annotation options')
    opt_te.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of parallel worker processes for the TE overlap '
                             'computation.  Transcripts are grouped by chromosome and '
                             'distributed across workers.  threads=1 (default) runs '
                             'sequentially with no subprocess overhead.  '
                             'Recommended: use the number of chromosomes (e.g. 22 for '
                             'human autosomes) or the number of available CPU cores, '
                             'whichever is smaller.')
    opt_te.add_argument('--skip_quantification', action='store_true', default=True,
                        help='Skip transcript quantification and use dummy TPM=0 '
                             '(default: True for cal_TE; real TPM can be merged later '
                             'with replace_tpm_from_quant.py)')
    opt_te.add_argument('--first_exon_threshold', type=float, default=50.0,
                        help='First-exon TE proportion threshold (%%) for TE-derived '
                             'transcript detection [default: 50.0]')
    opt_te.add_argument('--total_threshold', type=float, default=50.0,
                        help='Full-transcript TE proportion threshold (%%) for TE-derived '
                             'transcript detection [default: 50.0]')
    opt_te.add_argument('--te_overlap_threshold', type=int, default=10,
                        help='Minimum TE overlap length (bp) to distinguish Gene-alone '
                             'from TE-containing transcripts [default: 10]')
    opt_te.add_argument('--te_ratio_threshold', type=float, default=80.0,
                        help='TE proportion in full transcript (%%) for TE-alone '
                             'identification [default: 80.0]')
    opt_te.add_argument('--te_feature_threshold', type=int, default=50,
                        help='TE overlap length (bp) within TSS/TES 200 bp window '
                             'for sub-classification [default: 50]')

    opt_te_cal = parser_te.add_argument_group(
        'TE quantification options (bulk / per-spot / per-cell)'
    )
    opt_te_cal.add_argument('--percent_threshold', type=float, default=0.5,
                            help='Minimum TE-overlap proportion (overlap_length / TE_length) '
                                 'to associate a transcript with a TE. '
                                 'Applied in both bulk and per-spot/cell modes [default: 0.5]')
    opt_te_cal.add_argument('--bulk_quant', type=str, default=None,
                            help='Path to bulk transcript quantification TSV file '
                                 '(must contain transcript_id and TPM columns, e.g. output '
                                 'from miniQuant quantify or isoQuant). '
                                 'If provided, TE expression is computed at locus / subfamily / '
                                 'family / class level and written to <output_path>/Bulk_TE_output/.')
    opt_te_cal.add_argument('--sp_quant', type=str, default=None,
                            help='Path to directory containing the isoform MEX output '
                                 'from a prior miniQuant quantify --sp_mode run '
                                 '(e.g. /path/to/SP_output). '
                                 'If provided, per-spot TE metrics are computed '
                                 'and written to <output_path>/ST_TE_output/.')
    opt_te_cal.add_argument('--sc_quant', type=str, default=None,
                            help='Path to directory containing the isoform MEX output '
                                 'from a prior miniQuant quantify --sc_mode run '
                                 '(e.g. /path/to/SC_output). '
                                 'If provided, per-cell TE metrics are computed '
                                 'and written to <output_path>/SC_TE_output/.')
    opt_te_cal.add_argument('--output_loci', action='store_true', default=True,
                            help='Also output spot_te_loci_counts.tsv '
                                 '(spot/cell × individual TE locus matrix; may be very large)')

    # -----------------------------------------------------------------------
    parser_TrEESR = subparsers.add_parser('cal_K_value', aliases=['TrEESR'],help='Calculate K values')
    # parser_TransELS = subparsers.add_parser('quantify', aliases=['TransELS'],help='Isoform quantification')
    
    requiredNamed_TrEESR = parser_TrEESR.add_argument_group('required named arguments for calculation of K value')
    requiredNamed_TrEESR.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TrEESR.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    optional_TrEESR = parser_TrEESR.add_argument_group('optional arguments')
    optional_TrEESR.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",required=False)
    optional_TrEESR.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=False)
    optional_TrEESR.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    optional_TrEESR.add_argument('--sr_region_selection',type=str, default='read_length',help="SR region selection methods [default:read_length][read_length,num_exons,real_data]")
    optional_TrEESR.add_argument('--singular_values_tol',type=float,default=0,help="Singular value tolerence")
    optional_TrEESR.add_argument('--filtering',type=str,default='False', help="Whether the very short long reads will be filtered[default:True][True,False]")
    optional_TrEESR.add_argument('--READ_JUNC_MIN_MAP_LEN',type=int, default=1,help="minimum mapped read length to consider a junction")
    optional_TrEESR.add_argument('--same_struc_isoform_handling',type=str, default='merge',help="How to handle isoforms with same structures within a gene[default:merge][merge,keep]")
    optional_TrEESR.add_argument('--multi_exon_region_weight',type=str, default='regular',help="The weight in matrix A for multi_exon_region[default:regular][regular,minus_inner_region]")
    optional_TrEESR.add_argument('--output_matrix_info',type=str, default='False',help="Whether output matrix info [default:False] [True,False]")
    optional_TrEESR.add_argument('--normalize_sr_A',type=str, default='True',help="Whether normalize sr A [default:True] [True,False]")
    optional_TrEESR.add_argument('--keep_sr_exon_region',type=str, default='nonfullrank',help="Keep exon region for SR if using real data to filter region nonfullrank: only keep zero count exon region in non fulll rank gene [default:nonfullrank][nonfullrank,all,none]")
    optional_TrEESR.add_argument('--use_weight_matrix',type=str, default='False',help="Whether use weight matrix[default:True][True,False]")
    optional_TrEESR.add_argument('--normalize_lr_A',type=str, default='True',help="Whether normalize lr A [default:True] [True,False]")
    optional_TrEESR.add_argument('--add_full_length_region',type=str, default='all',help="Whether add full length region[default:all] [all,nonfullrank,none]")
    optional_TrEESR.add_argument('--sr_design_matrix',type=str, default='weight',help="How to calculate design matrix [default:weight][weight,binary]")
    weight_path = os.path.dirname(os.path.realpath(__file__))+'/weights/nanosim_weight_dict.pkl'
    # assert os.path.exists(weight_path)
    optional_TrEESR.add_argument('--region_weight_path',type=str, default=None,help="Mili LR region weight path")

   
    parser_EM = subparsers.add_parser('quantify', aliases=['EM'],help='Isoform quantification by EM algorithm')
    requiredNamed_EM = parser_EM.add_argument_group('required named arguments for isoform quantification')
    requiredNamed_EM.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_EM.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    optional_EM = parser_EM.add_argument_group('optional arguments')
    optional_EM.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=False,default=None)
    optional_EM.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",default=None)
    optional_EM.add_argument('-srfastq','--short_read_fastq', type=str, help="The path of short read fastq file",default=None)
    optional_EM.add_argument('-sr_m1','--short_read_mate1_fastq', type=str, help="The path of short read mate 1 fastq file",default=None)
    optional_EM.add_argument('-sr_m2','--short_read_mate2_fastq', type=str, help="The path of short read mate 2 fastq file",default=None)

    optional_EM.add_argument('-ref_genome','--reference_genome', type=str, help="The path of reference genome file",default=None)
    optional_EM.add_argument('--SR_quantification_option', type=str, help="SR quantification option[Options: Mili, kallisto,Salmon, RSEM] [default:kallisto]",default='kallisto')
    # optional_EM.add_argument('--kallisto_index', type=str, help="kallisto index",default='/fs/project/PCON0009/Yunhao/Project/Mili/Annotation/kallistoIndex/gencode.v39.transcripts.clean.dedup.m')
    optional_EM.add_argument('--alpha',type=str,default='adaptive', help="Alpha[default:adaptive]: SR and LR balance parameter")
    optional_EM.add_argument('--beta',type=str, default='1e-6',help="Beta[default:1e-6]: L2 regularization parameter")
    optional_EM.add_argument('--filtering',type=str,default='False', help="Whether the very short long reads will be filtered[default:False][True,False]")
    optional_EM.add_argument('--multi_mapping_filtering',type=str,default='best', help="How to filter multi-mapping reads[default:best][unique_only,best]")
    optional_EM.add_argument('--training',type=str,default='False', help="Generate training dict")
    optional_EM.add_argument('--DL_model',type=str,default=None,help='DL model to use')
    optional_EM.add_argument('--assign_unique_mapping_option',type=str,default='linear_model',help='How to assign unique mapping reads [Options:linear_model,manual_assign] [default:linear_model]')
    optional_EM.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    optional_EM.add_argument('--READ_JUNC_MIN_MAP_LEN',type=int, default=1,help="minimum mapped read length to consider a junction")
    optional_EM.add_argument('--use_weight_matrix',type=str, default='False',help="Whether use weight matrix[default:True][True,False]")
    optional_EM.add_argument('--normalize_lr_A',type=str, default='True',help="Whether normalize lr A [default:True] [True,False]")
    # optional_EM.add_argument('--same_struc_isoform_handling',type=str, default='keep',help="How to handle isoforms with same structures within a gene[default:merge][merge,keep]")
    optional_EM.add_argument('--add_full_length_region',type=str, default='all',help="Whether add full length region[default:all] [all,nonfullrank,none]")
    optional_EM.add_argument('--multi_exon_region_weight',type=str, default='regular',help="The weight in matrix A for multi_exon_region[default:regular][regular,minus_inner_region]")
    optional_EM.add_argument('--sr_design_matrix',type=str, default='weight',help="How to calculate design matrix [default:weight][weight,binary]")
    optional_EM.add_argument('--output_matrix_info',type=str, default='False',help="Whether output matrix info [default:False] [True,False]")
    optional_EM.add_argument('--normalize_sr_A',type=str, default='True',help="Whether normalize sr A [default:False] [True,False]")
    optional_EM.add_argument('--sr_region_selection',type=str, default='read_length',help="SR region selection methods [default:real_data][read_length,num_exons,real_data]")
    optional_EM.add_argument('--keep_sr_exon_region',type=str, default='nonfullrank',help="Keep exon region for SR if using real data to filter region nonfullrank: only keep zero count exon region in non fulll rank gene [default:nonfullrank][nonfullrank,all,none]")
    optional_EM.add_argument('--region_weight_path',type=str, default=None,help="Mili LR region weight path")
    optional_EM.add_argument('--EM_choice',type=str, default='LR',help="EM_choice[SR,LR,hybrid]")
    optional_EM.add_argument('--iter_theta',type=str, default='False',help="Whether use updated theta to re-calculate conditional prob [True,False]")
    optional_EM.add_argument('--kde_path',type=str, default='/fs/project/PCON0009/Au-scratch2/haoran/_projects/long_reads_rna_seq_simulator/models/kde_H1-hESC_dRNA',help="KDE model path")
    optional_EM.add_argument('--eff_len_option',type=str, default='kallisto',help="Calculation of effective length option [kallisto,RSEM]")
    optional_EM.add_argument('--EM_SR_num_iters',type=int, default=200,help="Number of EM SR iterations")
    optional_EM.add_argument('--EM_output_frequency',type=int, default=200,help="Frequency(in itertations) of outputting EM results")
    optional_EM.add_argument('--pretrained_model_path',type=str, default='cDNA-ONT',help="The pretrained model path to identify the alpha")
    optional_EM.add_argument('--alpha_df_path',type=str, default=None,help="Alpha df path")
    optional_EM.add_argument('--inital_theta','--initial_theta',type=str, default='uniform',help="initial_theta [LR,SR,LR_unique,SR_unique,uniform,hybrid,hybrid_unique,random]")
    optional_EM.add_argument('--inital_theta_eps','--initial_theta_eps',type=float, default=0.0,help="initial_theta eps [float]")
    optional_EM.add_argument('--eps_strategy',type=str, default='add_eps_small',help="how to add initial_theta eps [add_eps_all,add_eps_small]. (add_eps_small: add isoform with theta < eps with eps. add_eps: add eps to all isoforms)")
    optional_EM.add_argument('--isoform_start_end_site_tolerance',type=int, default=20,help="Isoform Start and end site tolerance for mapping long reads")
    optional_EM.add_argument('--junction_site_tolerance',type=int, default=5,help="Junction site tolerance for mapping long reads")
    optional_EM.add_argument('--read_len_dist_sm_dict_path',type=str, default=None,help="The path of read length distribution for long reads")
    optional_EM.add_argument('--LR_cond_prob_calc',type=str, default='form_2',help="How to calculate LR length distribution [form_1,form_2]")
    optional_EM.add_argument('--singular_values_tol',type=float,default=0,help="Singular value tolerence")

    # Bulk mode arguments
    bulk_EM = parser_EM.add_argument_group('bulk mode arguments')
    bulk_EM.add_argument('--bulk_mode', action='store_true', default=False,
                         help="Enable bulk mode: standard transcript-level quantification without "
                              "barcode/UMI extraction (default mode when no other mode flag is given)")

    # Single-cell mode arguments
    sc_EM = parser_EM.add_argument_group('single-cell mode arguments')
    sc_EM.add_argument('--sc_mode',action='store_true',default=False,
                       help="Enable single-cell mode: extract cell barcodes/UMIs and output a cell × isoform count matrix")
    sc_EM.add_argument('--barcode_in_readname',action='store_true',default=True,
                       help="Barcode/UMI embedded in read name by flexiplex (default: True). Use --no-barcode_in_readname for CB/UB SAM tags")
    sc_EM.add_argument('--no_barcode_in_readname',dest='barcode_in_readname',action='store_false',
                       help="Read barcode from CB:Z and UMI from UB:Z SAM tags instead of read name")
    sc_EM.add_argument('--barcode_separator',type=str,default='_',
                       help="Separator used by flexiplex between readname/barcode/UMI in read name (default: '_')")
    sc_EM.add_argument('--cb_tag',type=str,default='CB',
                       help="SAM tag for cell barcode when not in read name (default: CB)")
    sc_EM.add_argument('--umi_tag',type=str,default='UB',
                       help="SAM tag for UMI when not in read name (default: UB)")
    sc_EM.add_argument('--umi_dedup_hamming',type=int,default=0,
                       help="Hamming distance for UMI deduplication (0=exact match, 1=1-mismatch; default: 0)")
    sp_EM = parser_EM.add_argument_group('spatial transcriptomics mode arguments')
    sp_EM.add_argument('--sp_mode', action='store_true', default=False,
                       help="Enable spatial transcriptomics mode: extract spot barcodes/UMIs and output "
                            "a spot × isoform count matrix. Use --barcode_in_readname / --no_barcode_in_readname "
                            "to control barcode extraction (same as single-cell mode).")
    sp_EM.add_argument('--barcode_whitelist', type=str, default=None,
                       help="Path to spatial barcode whitelist file (one barcode per line). "
                            "Only barcodes present in this list are retained in the output. "
                            "For Visium, use the tissue barcode list; for Visium HD use the filtered "
                            "tissue positions barcodes. Barcodes not in the whitelist are discarded "
                            "before quantification.")
    sp_EM.add_argument('--tissue_positions', type=str, default=None,
                       help="Path to tissue_positions.csv (Visium format: columns barcode, in_tissue, "
                            "array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres; with or without "
                            "header). When provided, spatial coordinates are appended as leading columns to "
                            "spot_isoform_counts.tsv, enabling direct import into Squidpy / SpatialDE.")

    args = parser.parse_args()

    # -------------------------------------------------------------------
    # cal_TE: early dispatch — does not need bulk-EM config setup
    # -------------------------------------------------------------------
    if getattr(args, 'subparser_name', None) == 'cal_TE':
        _run_cal_te(args)
        return
    # -------------------------------------------------------------------

    if args.filtering == 'True':
        args.filtering = True
    else:
        args.filtering = False
    # config.same_struc_isoform_handling = args.same_struc_isoform_handling
    config.output_path = args.output_path
    config.threads = args.threads
    config.same_struc_isoform_handling = 'keep'
    config.READ_JUNC_MIN_MAP_LEN = args.READ_JUNC_MIN_MAP_LEN
    config.multi_exon_region_weight = args.multi_exon_region_weight
    config.sr_region_selection = args.sr_region_selection
    config.region_weight_path = args.region_weight_path
    config.sr_design_matrix = args.sr_design_matrix
    if args.output_matrix_info == 'True':
        config.output_matrix_info = True
    else:
        config.output_matrix_info = False
    config.keep_sr_exon_region = args.keep_sr_exon_region
    if args.normalize_sr_A == 'True':
        config.normalize_sr_A = True
    else:
        config.normalize_sr_A = False
    if args.normalize_lr_A == 'True':
        config.normalize_lr_A = True
    else:
        config.normalize_lr_A = False
    if args.use_weight_matrix == 'True':
        config.use_weight_matrix = True
    else:
        config.use_weight_matrix = False
    config.add_full_length_region = args.add_full_length_region
    config.singular_values_tol = args.singular_values_tol
    # config.kallisto_index = args.kallisto_index
    # print('\n'.join(f'{k}={v}' for k, v in vars(args).items()))
    if args.subparser_name in ['cal_K_value','TrEESR']:
        print('[INFO] Calculate K values')
        TrEESR(args.gtf_annotation_path,args.output_path,args.short_read_sam_path,args.long_read_sam_path,args.sr_region_selection,args.filtering,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN)
    elif args.subparser_name in ['quantify','EM']:
        config.kde_path = args.kde_path
        if args.training == 'True':
            args.training = True
        else:
            args.training = False
        print('[INFO] Isoform quantification by miniQuant',flush=True)
        if (args.short_read_sam_path is None) or (args.alpha == 1.0):
            args.alpha = 1.0
            args.SR_quantification_option = 'Mini'
        if (args.alpha == 'adaptive'):
            alpha = 'adaptive'
        else:
            try:
                alpha = float(args.alpha)
            except:
                raise Exception('Alpha given is not numeric')
        if (args.beta == 'adaptive'):
            beta = 'adaptive'
        else:
            try:
                beta = float(args.beta)
            except:
                raise Exception('Beta given is not numeric')
        # if args.SR_quantification_option not in ['Mili','kallisto','Salmon','RSEM']:
        #     raise Exception('SR_quantification_option is not valid.Options: [Mili, kallisto,Salmon, RSEM]')
        if (args.multi_mapping_filtering is None) or (not args.multi_mapping_filtering in ['unique_only','best']):
            args.multi_mapping_filtering = 'no_filtering'
        SR_fastq_list = []
        if args.short_read_fastq is not None:
            SR_fastq_list = [args.short_read_fastq]
        elif args.short_read_mate1_fastq is not None:
            SR_fastq_list = [args.short_read_mate1_fastq,args.short_read_mate2_fastq]
        if args.DL_model is None:
            args.DL_model = args.SR_quantification_option + '.pt'
        config.EM_SR_num_iters = args.EM_SR_num_iters
        config.inital_theta_eps = args.inital_theta_eps
        config.EM_output_frequency = args.EM_output_frequency
        config.isoform_start_end_site_tolerance = args.isoform_start_end_site_tolerance
        config.junction_site_tolerance = args.junction_site_tolerance
        config.eps_strategy = args.eps_strategy
        config.read_len_dist_sm_dict_path = args.read_len_dist_sm_dict_path
        config.LR_cond_prob_calc = args.LR_cond_prob_calc
        # Mode config
        config.bulk_mode = args.bulk_mode
        # SC mode config
        config.sc_mode = args.sc_mode
        config.barcode_in_readname = args.barcode_in_readname
        config.barcode_separator = args.barcode_separator
        config.cb_tag = args.cb_tag
        config.umi_tag = args.umi_tag
        config.umi_dedup_hamming = args.umi_dedup_hamming
        # ST mode config
        config.sp_mode = args.sp_mode
        config.barcode_whitelist_path = args.barcode_whitelist
        config.tissue_positions_path = args.tissue_positions
        active_modes = sum([config.bulk_mode, config.sc_mode, config.sp_mode])
        if active_modes > 1:
            raise ValueError('--bulk_mode, --sc_mode and --sp_mode are mutually exclusive.')
        if config.sp_mode:
            # Spatial mode reuses the barcode extraction path of SC mode
            config.sc_mode = True
        if args.pretrained_model_path in ['cDNA-ONT','dRNA-ONT','cDNA-PacBio']:
            args.pretrained_model_path = os.path.dirname(os.path.realpath(__file__))+'/pretrained_models/' + args.pretrained_model_path +'/'
        config.pretrained_model_path = args.pretrained_model_path
        if args.EM_choice == 'SR':
            config.eff_len_option = args.eff_len_option
            args.long_read_sam_path = None
            args.alpha = 0
            args.inital_theta = 'SR'
            config.alpha = args.alpha
            config.alpha_df_path = args.alpha_df_path
            config.inital_theta = args.inital_theta
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)
        elif args.EM_choice == 'hybrid':
            # args.alpha = 0.5
            config.alpha = args.alpha
            config.alpha_df_path = args.alpha_df_path
            if args.alpha_df_path is None:
                config.alpha_df_path = args.output_path +'/hybrid_alpha.tsv'
            config.inital_theta = args.inital_theta
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)
        else:
            if args.EM_choice == 'LR':
                args.EM_choice = 'LIQA_modified'
            args.short_read_sam_path = None
            args.alpha = 1
            args.inital_theta = 'LR'
            config.alpha = args.alpha
            config.alpha_df_path = args.alpha_df_path
            config.inital_theta = args.inital_theta
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)    
    else:
        parser.print_help()
if __name__ == "__main__":
    parse_arguments()
