import os
import numpy as np

def construct_bins(bin_size, ref, cov_matrix, rna_loc_list):
    print("Convert the genes to genome bins with the size of %skb..."%bin_size)
    bin_size *= 1000
    f_in = open("genome_size_%s.txt"%ref)
    bin_dic = {}
    bin_list = [] 
    for line in f_in:
        chr_num = line.split("\t")[0]
        length = int(line.split("\t")[1])
        start = 1
        while (length-start) > bin_size-1:
            bin_ = chr_num + ":" + str(start) + "-" + str(start+bin_size-1)
            bin_list.append(bin_)
            bin_dic[bin_] = []
            start = start + bin_size
        bin_ = chr_num + ":" + str(start) + "-" + str(length)
        bin_list.append(bin_)
        bin_ = chr_num + ":" + str(start) + "-" + str(start+bin_size-1)
        bin_dic[bin_] = []
    f_in.close()
    
    for i in range(len(rna_loc_list)):
        rna_loc = rna_loc_list[i]
        chrom = rna_loc.split(":")[0]
        start = float(rna_loc.split(":")[1].split("-")[0])
        stop = float(rna_loc.split(":")[1].split("-")[1])
        rna_mean_loc = (start + stop) / 2
        bin_ = chrom + ":" + str(int(rna_mean_loc/bin_size)*bin_size+1) + "-" + str(int(rna_mean_loc/bin_size)*bin_size+bin_size)
        bin_dic[bin_].append(cov_matrix[i])
    
    bin_list_new = []
    cn_matrix = []
    count = 0
    for bin_ in bin_dic:
        cn_data = bin_dic[bin_]
        if len(cn_data) != 0:
            bin_list_new.append(bin_list[count])
            cn_matrix.append(np.mean(np.array(cn_data), axis=0))
        count += 1
    cn_matrix = np.array(cn_matrix)

    return cn_matrix, bin_list_new
        
