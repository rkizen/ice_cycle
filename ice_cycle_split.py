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
# This script partitions FLNC fasta files by transcript length for running with ICE.
# Allows for window overlap and remerging for second ICE run to resolve split clusters.
#


ap = argparse.ArgumentParser(description='This script partitions FLNC fasta files by transcript length for running with ICE. \
Allows for window overlap and remerging for second ICE run to resolve split clusters.')
ap.add_argument('-f', type=str, nargs=1, help='Fasta file containing transcript sequences')
ap.add_argument('-w', type=int, nargs=1, help='Window size for length of sequences')
ap.add_argument('-o', type=int, nargs=1, help='Overlap size for length of sequences')
ap.add_argument('-p', type=str, nargs=1, help='Prefix for output split fasta files')

opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0
if not opts.f:
    print("fasta file missing")
    missing_arg_flag = 1
if not opts.w:
    print("window size missing")
    missing_arg_flag = 1
if not opts.o:
    print("overlap size missing")
    missing_arg_flag = 1
if not opts.p:
    print("prefix missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")

fasta_file = opts.f[0]
window_size = int(opts.w[0])
overlap_size = int(opts.o[0])
outfile_prefix = opts.p[0]

#check that window size is an integer
try:
    window_size += 1
    window_size = window_size -1
except TypeError:
    print("Window size provided is not an integer")

#check that overlap size is an integer
try:
    overlap_size += 1
    overlap_size = overlap_size -1
except TypeError:
    print("Overlap size provided is not an integer")

#check that window size is larger than overlap size
if window_size <= overlap_size:
    print("Window size must be larger than overlap size")
    sys.exit()

num_windows = 100000 / window_size
overlap_half = overlap_size/2
window_low_list = []
window_high_list = []

window_dict = {} # window_dict[window range][fasta header] = fasta sequence 


for i in xrange(num_windows):
    if i == 0:
        window_low = 0
        window_low_list.append(window_low)
        window_high = window_size + overlap_size
        window_high_list.append(window_high)
        window_dict[window_low] = {}
        continue
    window_bin = window_size * i
    window_low = (window_size * i) 
    window_high = (window_size * (i+1)) + overlap_size
    window_low_list.append(window_low)
    window_high_list.append(window_high)
    
    window_dict[window_low] = {}
    

for seq_record in SeqIO.parse(fasta_file, "fasta"):
    seq_name = str(seq_record.id)
    seq_string = str(seq_record.seq)
    seq_string = seq_string.upper()
    seq_length = len(seq_string)
    
    window_bin = (seq_length / window_size) * window_size
    window_over = seq_length % window_size
    
    window_dict[window_bin][seq_name] = seq_string
    
    if window_over < overlap_size:
        next_window_bin = window_bin + window_size
        window_dict[next_window_bin][seq_name] = seq_string
 
 
for window_bin in window_low_list:
    
    seq_list = window_dict[window_bin].keys()
    if len(seq_list) == 0:
        continue
    
    seq_list.sort()
    new_outfile_name = outfile_prefix + "_" + str(window_bin) + "_" + str(window_bin + window_size + overlap_size) + ".fa"
    outfile = open(new_outfile_name,"w")
    
    
    for fasta_header in seq_list:
        fasta_seq = window_dict[window_bin][fasta_header]
    
        header_line = ">" + fasta_header
        outfile.write(header_line)
        outfile.write("\n")
        outfile.write(fasta_seq)
        outfile.write("\n")


