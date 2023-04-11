import argparse
import logging
import os
import numpy as np
import pandas as pd
from multiprocessing import Pool
from functools import partial

from convert2loc import construct_locs
from convert2bin import construct_bins
from preprocess import filter_cell_bin, normalize
from baseline_correct import Kmeans_clu_method, Ward_clu_method, GMM_method, synthetic_method
from segmentation import chr_segment
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def main():
    parser = argparse.ArgumentParser()
    # required parameters
    parser.add_argument("--input_file", default=None, type=str, required=True, help="The count data from scRNA-seq (gene x cell).")
    parser.add_argument("--output_dir", default=None, type=str, required=True, help="The output dictionary.")

    # other paramters
    parser.add_argument("--ensembl", default=86, type=int, help="The used Ensembl release. (default: 86)")
    parser.add_argument("--bin_size", default=220, type=float, help="The size of genomic bins (kb) (default: 220).")
    parser.add_argument("--baseline_estimation", default="Ward_clu_method", type=str, help="The baseline estimation method (default: Ward_clu_method).")
    parser.add_argument("--pca_components", default=3, type=int, help="The n_components of PCA for cell clustering (default: 3).")
    parser.add_argument("--n_clusters", default=6, type=int, help="The number of clusters for baseline detection (default: 6)")
    parser.add_argument("--K", default=5, type=float, help="The topK distances used for congruent map construction (default: 5).")
    parser.add_argument("--sigma", default=0.5, type=float, help="The sigma of Gaussian function for congruent map construction (default: 0.5).")
    parser.add_argument("--n_threads", default=20, type=int, help="The number of threads for segmentation (default: 20).")
    parser.add_argument("--ref", default="hg38", type=str, help="The reference used for bin construction (default: hg38).")
    parser.add_argument("--seed", default=111, type=int, help="random seed for GMM initialization (default: 111).") 
    
    parser.add_argument("--gmm_dp", action="store_true", help="Use dynamic programming for GMM.")
    parser.add_argument("--KR_norm", action="store_true", help="Use KR normlization for depth congruent map.")
    parser.add_argument("--overwrite_output_dir", action="store_true", help="Overwrite the content of the output directory.")
    parser.add_argument("--save_figures", action="store_true", help="Save the figures for analysis.")
    
    args = parser.parse_args()
    
    if (os.path.exists(args.output_dir) and os.listdir(args.output_dir) and not args.overwrite_output_dir):
        raise ValueError("Output directory ({}) already exists and is not empty.".format(args.output_dir))
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)


    logger = logging.getLogger(__name__)

    print("***** Begin running SeCNV *****")
    
    print("Step1: Sort the genes according to the genomic coordinates with Ensembl release %s."%args.ensembl) 
    construct_locs(args.input_file, os.path.join(args.output_dir, "RNA_expression_loc.csv"), args.ensembl)
    
    print("Step2: Filtering cells and genes.")
    # read data
    data = pd.read_csv(os.path.join(args.output_dir, "RNA_expression_loc.csv"), index_col=0)
    rna_loc_list = data.columns.values.tolist()
    sample_list = data._stat_axis.values.tolist()
    chr_name = [each.split(":")[0] for each in rna_loc_list]
    cov_matrix = data.to_numpy()
    chr_name = np.array(chr_name)
    rna_loc_list = np.array(rna_loc_list)
    sample_list = np.array(sample_list)
    cov_matrix, chr_name, rna_loc_list, sample_list = filter_cell_bin(cov_matrix.T, chr_name, rna_loc_list, sample_list)
    
    print("Step3: Data normalization.")
    cov_matrix = normalize(cov_matrix)
    
    print("Step4: Baseline profile estimation.") 
    cov_matrix, normal_sample_list = eval(args.baseline_estimation)(cov_matrix, sample_list, args.pca_components, args.n_clusters, args.output_dir, args.save_figures, args.seed, args.gmm_dp)
    if args.save_figures:
        plt.figure(figsize=(6,5))
        plt.clf()
        sns.set()
        sns.heatmap(cov_matrix.T, center=0, cmap=sns.diverging_palette(220, 20, as_cmap=True))
        plt.savefig(os.path.join(args.output_dir, "cnv_profile_beforeseg.png"))
        plt.close()
    pd.DataFrame(normal_sample_list).to_csv(os.path.join(args.output_dir, "normal_list.csv"))
'''
    print("Step5: Converting genes to genomic bins with %s."%args.ref)
    bin_cov_matrix, bin_list = construct_bins(args.bin_size, args.ref, cov_matrix, rna_loc_list)
    
    print("Step6: Breakpoints detection using %s thread(s)."%args.n_threads)
    chosen_chr = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    chr_name = [each.split(":")[0] for each in bin_list]
    chr_name = np.array(chr_name)
    partial_chr_segment = partial(chr_segment, chr_name=chr_name, cov_matrix=bin_cov_matrix, K=args.K, sigma=args.sigma, KR_norm=args.KR_norm, output_dir=args.output_dir, save_figures=args.save_figures)
    pool = Pool(args.n_threads)
    new_bin_cov_matrix = np.vstack(pool.map(partial_chr_segment, chosen_chr))
    pool.close()
    pool.join()
    if args.save_figures:
        plt.figure(figsize=(6,5))
        plt.clf()
        sns.set()
        sns.heatmap(new_bin_cov_matrix.T, center=0, cmap=sns.diverging_palette(220, 20, as_cmap=True))
        plt.savefig(os.path.join(args.output_dir, "cnv_profile_afterseg.png"))
        plt.close()
    
    pd.DataFrame(new_bin_cov_matrix, columns=sample_list, index=bin_list).to_csv(os.path.join(args.output_dir, "relative_results_bin_by_cell.csv"))
    if args.save_figures:
        chrom_loc = {}
        chrom_loc["chr0"] = 0
        for i in range(len(bin_list)):
            bin_ = bin_list[i]
            chrom = bin_.split(":")[0]
            chrom_loc[chrom] = i
        
        plt.figure(figsize=(6,5))
        plt.clf()
        sns.set()
        sns.heatmap(new_bin_cov_matrix.T, center=0, cmap=sns.diverging_palette(220, 20, as_cmap=True))
        for loc in chrom_loc:
            plt.vlines(chrom_loc[loc], 0 , len(sample_list), colors="grey", linestyles="--")
        plt.savefig(os.path.join(args.output_dir, "cnv_profile_inbin.png"))
        plt.close()

        normal_index = [sample_list[i] in normal_sample_list for i in range(len(sample_list))]
        abnormal_index = [sample_list[i] not in normal_sample_list for i in range(len(sample_list))]
        new_bin_cov_matrix_sort = np.hstack((new_bin_cov_matrix[:, normal_index], new_bin_cov_matrix[:, abnormal_index]))
        plt.figure(figsize=(6,5))
        plt.clf()
        sns.set()
        sns.heatmap(new_bin_cov_matrix_sort.T, center=0, cmap=sns.diverging_palette(220, 20, as_cmap=True))
        for loc in chrom_loc:
            plt.vlines(chrom_loc[loc], 0 , len(sample_list), colors="grey", linestyles="--")
        plt.savefig(os.path.join(args.output_dir, "cnv_profile_inbin_sort.png"))
        plt.close()
'''

if __name__ == "__main__":
    main()
