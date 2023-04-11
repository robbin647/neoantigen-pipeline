import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def filter_cell_bin(data, chr_name, rna_loc_list, sample_list):
    logger.info("There are %s cells and %s genes in the dataset."%(data.shape[1], data.shape[0]))
    
    idx = np.argwhere(np.count_nonzero(data, axis=0) < 200)
    data = np.delete(data, idx, axis=1)
    sample_list = np.delete(sample_list, idx, axis=0)
    
    logger.info("Remove cells with fewer than 200 genes, remaining %s cells."%data.shape[1])
    
    idx = np.argwhere(np.count_nonzero(data, axis=1) < (data.shape[1]*0.05))
    data = np.delete(data, idx, axis=0)
    chr_name = np.delete(chr_name, idx, axis=0)
    rna_loc_list = np.delete(rna_loc_list, idx, axis=0)
    
    logger.info("Remove genes detected in less than 0.05 of the cells, remaining %s genes."%data.shape[0])
    
    chosen_chr = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    for chrom in chosen_chr:
        chrom_data = data[np.where(chr_name == chrom)]
        idx = np.argwhere(np.count_nonzero(chrom_data, axis=0) < 5)
        data = np.delete(data, idx, axis=1)
        sample_list = np.delete(sample_list, idx, axis=0)
    
    logger.info("Remove cells with less than 5 genes in chromosomes, remaining %s cells."%data.shape[1])

    return data, chr_name, rna_loc_list, sample_list


def normalize(data):
    
    data = np.log(np.sqrt(data) + np.sqrt(data+1))
    data = data - np.mean(data, axis=0)
    import rpy2
    from rpy2 import robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r = robjects.r
    r['library']('dlm')
    model = r['dlmModPoly'](order=1, dV=0.16, dW=0.001)
    for i in range(data.shape[1]):
        data[:, i] = r['dlmSmooth'](data[:, i], model)[0][1:]
    data = data - np.mean(data, axis=0)

    return data
