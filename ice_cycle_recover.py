import re
import sys
import time
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO
import os
import argparse

#
# This script examines the outputs from each split cluster run to find clusters which may have been split by the length cutoffs.
# It then generates a new fasta file containing the fasta sequences for the flnc sequences that were involved in these split cluster scenarios.
# Run the output fasta through cluster again and then merge.
#


ap = argparse.ArgumentParser(description='This script examines the outputs from each split cluster run to find clusters which \
may have been split by the length cutoffs. It then generates a new fasta file containing the fasta sequences for the flnc sequences \
that were involved in these split cluster scenarios. Run the output fasta through cluster again and then merge.')

ap.add_argument('-r', type=str, nargs=1, help='File with list of cluster report csv files')
ap.add_argument('-f', type=str, nargs=1, help='File with list of split fasta files used for clustering')
ap.add_argument('-o', type=str, nargs=1, help='Output fasta file name')

opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0
if not opts.r:
    print("Cluster report file list missing")
    missing_arg_flag = 1
if not opts.f:
    print("Fasta file list missing")
    missing_arg_flag = 1
if not opts.o:
    print("Output file name missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")

filelist_file = opts.r[0]
fastalist_file = opts.f[0]
outfile_name = opts.o[0]

########################################
filelist_file_contents = open(filelist_file).read().rstrip("\n").split("\n")
fastalist_file_contents = open(fastalist_file).read().rstrip("\n").split("\n")
outfile = open(outfile_name,"w")

fasta_dict = {} # fasta_dict[fasta id] = fasta seq

for filepath in fastalist_file_contents:
    print("going through fasta")
    for seq_record in SeqIO.parse(filepath, "fasta"):
        seq_name = str(seq_record.id)
        
        seq_string = str(seq_record.seq)
        seq_string = seq_string.upper()
        seq_length = len(seq_string)
        
        fasta_dict[seq_name] = seq_string

cluster_dict = {} # cluster_dict[file num][cluster id][flnc id] = 1

flnc_cluster_dict = {} # flnc_cluster_dict[flnc id][file cluster id] = 1

multi_cluster_flnc_dict = {} # multi_cluster_flnc_dict[flnc id] = number of occurrences
file_count = 0

flnc_list = []

for filepath in filelist_file_contents:
    
    report_file_content = open(filepath).read().rstrip("\n").split("\n")
    file_count += 1
    file_id = str(file_count)
    cluster_dict[file_id] = {}
    
    for line in report_file_content:
        if line.startswith("cluster_id"):
            continue
        
        line_split = line.split(",")
        
        cluster_id = line_split[0]
        flnc_id = line_split[1]
        
        if cluster_id not in cluster_dict[file_id]:
            cluster_dict[file_id][cluster_id] = {}
        
        cluster_dict[file_id][cluster_id][flnc_id] = 1
        
        if flnc_id not in flnc_cluster_dict:
            flnc_cluster_dict[flnc_id] = {}
        
        file_cluster_id = file_id + "_" + cluster_id
        flnc_cluster_dict[flnc_id][file_cluster_id] = 1
        
        if flnc_id not in multi_cluster_flnc_dict:
            multi_cluster_flnc_dict[flnc_id] = 0
            flnc_list.append(flnc_id)
        
        multi_cluster_flnc_dict[flnc_id] += 1

uniq_flnc_dict = {} # uniq_flnc_dict[flnc] = 1

for flnc_id in flnc_list:
    
    if multi_cluster_flnc_dict[flnc_id] > 1:
        file_cluster_id_list  = flnc_cluster_dict[flnc_id].keys()
        
        #collect all clusters that flnc belongs to
        for file_cluster_id in file_cluster_id_list:
            file_id = file_cluster_id.split("_")[0]
            cluster_id = file_cluster_id.split("_")[1]
            
            #collect all flnc in cluster
            for this_flnc_id in cluster_dict[file_id][cluster_id]:
                uniq_flnc_dict[this_flnc_id] = 1

uniq_flnc_list = uniq_flnc_dict.keys()
uniq_flnc_list.sort()

for flnc_id in uniq_flnc_list:
    flnc_seq = fasta_dict[flnc_id]
    
    outline = ">" + flnc_id
    outfile.write(outline)
    outfile.write("\n")
    
    outline = flnc_seq
    outfile.write(outline)
    outfile.write("\n")


