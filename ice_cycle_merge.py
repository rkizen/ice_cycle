import re
import sys
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO
import os
import argparse

#
# This script merges the outputs from the recover cluster run and the split cluster runs
#

ap = argparse.ArgumentParser(description='This script merges the outputs from the recover cluster run and the split cluster runs')

ap.add_argument('-rr', type=str, nargs=1, help='Recover cluster report csv file')
ap.add_argument('-rf', type=str, nargs=1, help='Recover cluster fasta file')

ap.add_argument('-sf', type=str, nargs=1, help='File list of fasta files for split run inputs')

ap.add_argument('-cr', type=str, nargs=1, help='File list of cluster report csv files for split runs')
ap.add_argument('-cf', type=str, nargs=1, help='File list of cluster fasta files for split run outputs')

ap.add_argument('-p', type=str, nargs=1, help='Prefix for output files')

opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0
if not opts.rr:
    print("Recover cluster report csv file missing")
    missing_arg_flag = 1
if not opts.rf:
    print("Recover cluster fasta file missing")
    missing_arg_flag = 1
if not opts.sf:
    print("File list of fasta files for split run inputs missing")
    missing_arg_flag = 1
if not opts.cr:
    print("File list of cluster report csv files for split runs missing")
    missing_arg_flag = 1
if not opts.cf:
    print("File list of cluster fasta files for split run outputs missing")
    missing_arg_flag = 1
if not opts.p:
    print("Prefix for output files missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")


recovercluster_file = opts.rr[0]
recover_fasta_file = opts.rf[0]
flnc_fastalist_file = opts.sf[0]
filelist_file = opts.cr[0]
cluster_fastalist_file = opts.cf[0]
outfile_prefix = opts.p[0]



########################################

recovercluster_file_contents = open(recovercluster_file).read().rstrip("\n").split("\n")
filelist_file_contents = open(filelist_file).read().rstrip("\n").split("\n")
flnc_fastalist_file_contents = open(flnc_fastalist_file).read().rstrip("\n").split("\n")
cluster_fastalist_file_contents = open(cluster_fastalist_file).read().rstrip("\n").split("\n")


outfile_fasta_name = outfile_prefix + ".fa"
outfile_report_name = outfile_prefix + "_cluster_report.FL.csv"
outfile_split_name = outfile_prefix + "_split_clusters.txt"

outfile_fasta = open(outfile_fasta_name,"w")
outfile_report = open(outfile_report_name,"w")
outfile_split = open(outfile_split_name,"w")

outfile_report.write("cluster_id,read_id,read_type")
outfile_report.write("\n")

#print header for splut report file
outline = "\t".join(["rec_cluster_id","pre_cluster_line","num_clusters","diff_cluster_flag","cluster_support_line"])
outfile_split.write(outline)
outfile_split.write("\n") 

#####################################################
#Use this to make sure no sequences were lost
flnc_fasta_dict = {} # flnc_fasta_dict[fasta id] = 1

for filepath in flnc_fastalist_file_contents:
    #print("going through flnc fasta")
    for seq_record in SeqIO.parse(filepath, "fasta"):
        seq_name = str(seq_record.id)
        
        seq_string = str(seq_record.seq)
        seq_string = seq_string.upper()
        seq_length = len(seq_string)
        
        flnc_fasta_dict[seq_name] = 1

#####################################################
#collect cluster fasta sequences
cluster_fasta_dict = {} # cluster_fasta_dict[cluster id] = seq
cluster_id_dict = {} # cluster_id_dict[cluster id] = fasta cluster id
cluster_header_dict = {} # cluster_header_dict[cluster id] = fasta header 

file_count = 0

for filepath in cluster_fastalist_file_contents:
    
    file_count += 1
    file_id = str(file_count)
    
    #print("going through cluster fasta")
    for seq_record in SeqIO.parse(filepath, "fasta"):
        seq_name = str(seq_record.id)
        seq_desc = str(seq_record.description)

        
        seq_string = str(seq_record.seq)
        seq_string = seq_string.upper()
        seq_length = len(seq_string)
        cluster_fasta_id = seq_name.split()[0]
        cluster_id = cluster_fasta_id.split("/")[0]
        new_cluster_id = file_id + "_" + cluster_id
        
        new_cluster_fasta_id = cluster_fasta_id.replace(cluster_id,new_cluster_id)
        new_seq_name = seq_desc.replace(cluster_id,new_cluster_id)
        
        cluster_fasta_dict[new_cluster_id] = seq_string
        cluster_id_dict[new_cluster_id] = new_cluster_fasta_id
        cluster_header_dict[new_cluster_id] = new_seq_name

#####################################################  
#collect recover cluster fasta sequences
recover_cluster_fasta_dict = {} # recover_cluster_fasta_dict[fasta id] = seq
recover_cluster_id_dict = {} # recover_cluster_id_dict[cluster id] =  
recover_cluster_header_dict = {} # recover_cluster_header_dict[cluster id] = fasta header

#print("going through recover cluster fasta")
for seq_record in SeqIO.parse(recover_fasta_file, "fasta"):
    seq_name = str(seq_record.id)
        
    seq_string = str(seq_record.seq)
    seq_string = seq_string.upper()
    seq_length = len(seq_string)
    
    seq_desc = str(seq_record.description)

    cluster_fasta_id = seq_name.split()[0]
    cluster_id = cluster_fasta_id.split("/")[0]
    new_cluster_id = "rec_" + cluster_id
        
    new_cluster_fasta_id = cluster_fasta_id.replace(cluster_id,new_cluster_id)
    new_seq_name = seq_desc.replace(cluster_id,new_cluster_id)

    recover_cluster_fasta_dict[new_cluster_id] = seq_string
    recover_cluster_id_dict[new_cluster_id] = new_cluster_fasta_id
    recover_cluster_header_dict[new_cluster_id] = new_seq_name
    
#####################################################

#keep track of all flnc used in recover clustering
rec_clust_flnc_dict = {} # rec_clust_flnc_dict[flnc id] = 1
#use this to make sure you dont repeat clusters
uniq_rec_cluster_dict = {} # uniq_rec_cluster_dict[cluster id] = 1

for line in recovercluster_file_contents:
    if line.startswith("cluster_id"):
        continue
        
    line_split = line.split(",")
        
    cluster_id = line_split[0]
    flnc_id = line_split[1]
    new_cluster_id = "rec_" + cluster_id    
    
    rec_clust_flnc_dict[flnc_id] = new_cluster_id
    

    
    report_line = ",".join([new_cluster_id,flnc_id,"FL"])
    outfile_report.write(report_line)
    outfile_report.write("\n")
        
    # only print out if this is uniq fasta line
    if new_cluster_id not in uniq_rec_cluster_dict:
        uniq_rec_cluster_dict[new_cluster_id] = 1
    
        cluster_seq = recover_cluster_fasta_dict[new_cluster_id]
    
        cluster_fasta_header = recover_cluster_header_dict[new_cluster_id]

        cluster_fasta_header = ">" + cluster_fasta_header
        
        outfile_fasta.write(cluster_fasta_header)
        outfile_fasta.write("\n")
        outfile_fasta.write(cluster_seq)
        outfile_fasta.write("\n")
        
        
    
#use this to keep track of clusters with flnc that are included in recovered clusters
rec_redundant_cluster_dict = {} # rec_redundant_cluster_dict[cluster id] = 1
#use this to make sure you dont repeat clusters
uniq_cluster_dict = {} #uniq_cluster_dict[cluster id] = 1

cluster_flnc_dict = {} # cluster_flnc_dict[normal cluster id][flnc id] = 1

split_cluster_dict = {} # split_cluster_dict[recovery cluster id][normal cluster id] = 1


file_count = 0
for filepath in filelist_file_contents:
    
    report_file_content = open(filepath).read().rstrip("\n").split("\n")
    file_count += 1
    file_id = str(file_count)

    
    # collect clusters with flnc that are in recovery clusters
    for line in report_file_content:
        if line.startswith("cluster_id"):
            continue
        
        line_split = line.split(",")
        
        cluster_id = line_split[0]
        flnc_id = line_split[1]
        
        new_cluster_id = file_id + "_" + cluster_id
        
        if flnc_id in rec_clust_flnc_dict:
            rec_redundant_cluster_dict[new_cluster_id] = 1
            rec_cluster_id = rec_clust_flnc_dict[flnc_id]
            
            # this is for figuring out which rec clusters actually merged different clusters
            if rec_cluster_id not in split_cluster_dict:
                split_cluster_dict[rec_cluster_id] = {}
            split_cluster_dict[rec_cluster_id][new_cluster_id] = 1
            

        
    # print out fasta and report for non-recovery clusters
    for line in report_file_content:
        if line.startswith("cluster_id"):
            continue
        
        line_split = line.split(",")
        
        cluster_id = line_split[0]
        flnc_id = line_split[1]
        new_cluster_id = file_id + "_" + cluster_id
        
        #store info about cluster to flnc match
        if new_cluster_id not in cluster_flnc_dict:
            cluster_flnc_dict[new_cluster_id] = {}
        cluster_flnc_dict[new_cluster_id][flnc_id] = 1
        
        
        if new_cluster_id in rec_redundant_cluster_dict and flnc_id not in rec_clust_flnc_dict:
            print("Error with ambiguously claimed flnc sequence")
            print(file_id + " " + new_cluster_id + " " + flnc_id )
            sys.exit()
            
        # skip this if it is already represented in the recovery clusters
        if new_cluster_id in rec_redundant_cluster_dict:
            continue
        
        report_line = ",".join([new_cluster_id,flnc_id,"FL"])
        outfile_report.write(report_line)
        outfile_report.write("\n")
        
        # only print uniq cluster fasta
        if new_cluster_id not in uniq_cluster_dict:
            uniq_cluster_dict[new_cluster_id] = 1
    
            cluster_seq = cluster_fasta_dict[new_cluster_id]
            cluster_fasta_header = cluster_header_dict[new_cluster_id]
      
            cluster_fasta_header = ">" + cluster_fasta_header
            
            outfile_fasta.write(cluster_fasta_header)
            outfile_fasta.write("\n")
            outfile_fasta.write(cluster_seq)
            outfile_fasta.write("\n")
            
            

rec_cluster_list = split_cluster_dict.keys()        
rec_cluster_list.sort()



for rec_cluster_id in rec_cluster_list:
    
    if rec_cluster_id not in split_cluster_dict:
        continue
    pre_cluster_list = split_cluster_dict[rec_cluster_id].keys()
    pre_cluster_list.sort()
    
    #collect cluster support info
    cluster_support_list = []
    for pre_cluster_id in pre_cluster_list:
        pre_flnc_list = cluster_flnc_dict[pre_cluster_id].keys()
        pre_flnc_list.sort()
        pre_flnc_line = ",".join(pre_flnc_list)
        cluster_support_list.append(pre_flnc_line)
    
    cluster_support_line = ";".join(cluster_support_list)
    
    diff_cluster_flag = "same_cluster"
    
    #check if merged clusters contain same reads or if they have some different reads
    if len(pre_cluster_list) > 1:
        pre_cluster1 = pre_cluster_list[0]
        pre_flnc1_list = cluster_flnc_dict[pre_cluster1].keys()
        pre_flnc1_list.sort()
        
        for i in xrange(len(pre_cluster_list)):
            pre_cluster2 = pre_cluster_list[i]
            pre_flnc2_list = cluster_flnc_dict[pre_cluster2].keys()
            pre_flnc2_list.sort()
            
            if pre_flnc1_list != pre_flnc2_list:
                diff_cluster_flag = "diff_cluster"
        
        num_clusters = len(pre_cluster_list)

        pre_cluster_line  = ",".join(pre_cluster_list)
        outline = "\t".join([rec_cluster_id,pre_cluster_line,str(num_clusters),diff_cluster_flag,cluster_support_line])
        outfile_split.write(outline)
        outfile_split.write("\n") 
    
    
    

            


####################################################################################
