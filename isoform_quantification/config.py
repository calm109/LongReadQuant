use_weight_matrix = True
same_struc_isoform_handling = 'keep'
add_full_length_region = 'all'
READ_JUNC_MIN_MAP_LEN = 1
multi_exon_region_weight = 'regular'
output_matrix_info = False
normalize_sr_A = False
sr_region_selection = 'real_data'
keep_sr_exon_region = 'nonfullrank'
normalize_lr_A = True
kallisto_index = ''
kde_path = ''
region_weight_path = ''
eff_len_option = 'kallisto'
EM_SR_num_iters = 1000
alpha_df_path = None
alpha = 0.5
EM_output_frequency = 50
isoform_start_end_site_tolerance = 20
junction_site_tolerance = 5
eps_strategy = 'add_eps_small'
pseudo_count_SR = None
pseudo_count_LR = 1
std_f_len = None
mean_f_len = None
read_len_dist_sm_dict_path = None
LR_cond_prob_calc = 'form_1'
sr_design_matrix = 'regular'
output_path = None
pretrained_model_path = None
threads = None
singular_values_tol = 0

# Bulk mode
bulk_mode = False

# Single-cell mode
sc_mode = False
barcode_in_readname = True   # True: flexiplex embeds barcode in read name
barcode_separator = '_'      # separator between readname / barcode / UMI
cb_tag = 'CB'                # SAM tag for cell barcode (used when barcode_in_readname=False)
umi_tag = 'UB'               # SAM tag for UMI (used when barcode_in_readname=False)
umi_dedup_hamming = 0        # hamming distance for UMI deduplication (0 = exact match)

# Spatial transcriptomics mode
sp_mode = False
barcode_whitelist_path = None  # path to spatial barcode whitelist file (one barcode per line)
tissue_positions_path = None   # path to tissue_positions.csv (Visium format)