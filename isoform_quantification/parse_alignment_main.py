from parse_alignment import map_read,parse_read_line
from patch_mp import patch_mp_connection_bpo_17560
from parse_annotation_main import check_valid_region
from collections import defaultdict
import traceback
from operator import itemgetter, attrgetter
from functools import partial
import numpy as np
import time
import random
import multiprocessing as mp
import os
from util import check_region_type
import config

# from memory_profiler import profile
# def parse_alignment_iteration(alignment_file_path,gene_points_dict,gene_interval_tree_dict, filtered_gene_regions_dict,
#                     start_pos_list, start_gname_list, end_pos_list, end_gname_list,
#                     READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST,map_f,line_nums):
def debuginfoStr(info):
    print(info,flush=True)
    with open('/proc/self/status') as f:
        memusage = f.read().split('VmRSS:')[1].split('\n')[0][:-3]
    mem = int(memusage.strip())/1024
    print('Mem consumption: '+str(mem),flush=True)
def parse_alignment_iteration(alignment_file_path, READ_JUNC_MIN_MAP_LEN,map_f,temp_queue,long_read,aln_line_marker):
    os.nice(10)
    start_file_pos,num_lines = aln_line_marker
    with open(alignment_file_path, 'r') as aln_file:
        local_gene_regions_read_count = {}
        local_gene_regions_read_length = {}
        local_gene_regions_read_pos = {}
        aln_file.seek(start_file_pos)
        line_num_ct = 0
        max_buffer_size = 1e1
        buffer_size = 0
        for line in aln_file:
            line_num_ct += 1
            if line_num_ct > num_lines:
                break
            try:
                if line[0] == '@':
                    continue
                fields = line.split('\t')
                if (fields[2] == '*'):
                    continue
                aln_line = parse_read_line(line)
                mapping = map_f(points_dict,interval_tree_dict, filtered_gene_regions_dict,
                    start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                    READ_JUNC_MIN_MAP_LEN, CHR_LIST,aln_line)
                if (mapping['read_mapped']):
                    random.seed(mapping['read_name'])
                    for mapping_area in [random.choice(mapping['mapping_area'])]:
                        rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
                        if rname not in local_gene_regions_read_count:
                            local_gene_regions_read_count[rname],local_gene_regions_read_length[rname] = {},{}
                            local_gene_regions_read_pos[rname] = {}
                        if gname not in local_gene_regions_read_count[rname]:
                            local_gene_regions_read_count[rname][gname],local_gene_regions_read_length[rname][gname] = {},{}
                            local_gene_regions_read_pos[rname][gname] = {}
                        if region_name not in local_gene_regions_read_count[rname][gname]:
                            local_gene_regions_read_count[rname][gname][region_name],local_gene_regions_read_length[rname][gname][region_name] = 0,[]
                            local_gene_regions_read_pos[rname][gname][region_name] = []
                        local_gene_regions_read_count[rname][gname][region_name] += 1 
                        # if long_read:
                        local_gene_regions_read_length[rname][gname][region_name].append(mapping['read_length'])
                        local_gene_regions_read_pos[rname][gname][region_name].append(mapping)
                    buffer_size += 1
            except Exception as e:
                tb = traceback.format_exc()
                print(Exception('Failed to on ' + line, tb))
                continue
            if buffer_size > max_buffer_size:
                temp_queue.put((local_gene_regions_read_count,local_gene_regions_read_length,local_gene_regions_read_pos))
                local_gene_regions_read_count,local_gene_regions_read_length = {},{}
                local_gene_regions_read_pos = {}
                buffer_size = 0
        if buffer_size > 0:
            temp_queue.put((local_gene_regions_read_count,local_gene_regions_read_length,local_gene_regions_read_pos))
    return 
def mapping_listener(temp_queue,gene_regions_read_count,gene_regions_read_length,gene_regions_read_pos):
    num_mapped_lines = 0
    num_lines = 0
    while True:
        msg = temp_queue.get()
        if msg == 'kill':
            break
        else:
            local_gene_regions_read_count,local_gene_regions_read_length,local_gene_regions_read_pos = msg
            for rname in local_gene_regions_read_count:
                for gname in local_gene_regions_read_count[rname]:
                    for region in local_gene_regions_read_count[rname][gname]:
                        num_mapped_lines += local_gene_regions_read_count[rname][gname][region]
                        gene_regions_read_count[rname][gname][region] += local_gene_regions_read_count[rname][gname][region]
                        gene_regions_read_length[rname][gname][region] += local_gene_regions_read_length[rname][gname][region]
                        gene_regions_read_pos[rname][gname][region] += local_gene_regions_read_pos[rname][gname][region]

            # for mapping in local_all_mappings:
            #     num_lines += 1
            #     if len(mapping['gene_candidates'])>0:
            #         num_mapped_to_gene += 1
            #     if (mapping['read_mapped']):
            #         num_mapped_lines += 1
            #         for mapping_area in [random.choice(mapping['mapping_area'])]:
            #             rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
            #             if region_name in gene_regions_read_count[rname][gname]:
            #                 gene_regions_read_count[rname][gname][region_name] += 1 
            #                 gene_regions_read_length[rname][gname][region_name].append(mapping['read_length'])
            #         read_lens.append(mapping['read_length'])
            #         read_names.update(local_read_names)
    return gene_regions_read_count,gene_regions_read_length,num_mapped_lines,gene_regions_read_pos

# @profile
def get_aln_line_marker(alignment_file_path,threads):
    with open(alignment_file_path, 'r') as aln_file:
        line_offset = []
        offset = 0
        for line in aln_file:
            if line[0] != '@':
                line_offset.append(offset)
            offset += len(line)
    num_aln_lines = len(line_offset)
    chunksize, extra = divmod(num_aln_lines, threads)
    if extra:
        chunksize += 1
    aln_line_marker = []
    for i in range(threads):
        aln_line_marker.append((line_offset[i*chunksize],chunksize))
    return aln_line_marker
