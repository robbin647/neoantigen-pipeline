import numpy as np
import pandas as pd
from pyensembl import EnsemblRelease

def construct_locs(rna_file, loc_file, resembl):
    # transfer the gene to genome locations
    RNA_data = pd.read_csv(rna_file, index_col=0, sep="\t")

    cell_list = RNA_data.columns.values.tolist()
    gene_list = RNA_data._stat_axis.values.tolist()
    expression_data = RNA_data.to_numpy().T

    # Use EnsemblRelease
    data = EnsemblRelease(resembl)
    all_gene_list = data.gene_names()

    loc_list = [] # gene location
    choose_list = [] # chosed gene
    sort_list = [] # gene location for sorting
    for i in range(len(gene_list)):
        gene = gene_list[i]
        if gene in all_gene_list:
            contig = data.genes_by_name(gene)[0].contig
            start = data.genes_by_name(gene)[0].start
            end = data.genes_by_name(gene)[0].end
            if contig in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]:
                loc_list.append("chr" + contig + ":" + str(start) + "-" + str(end))
                choose_list.append(i)
                sort_list.append((int(contig), (start+end)/2))

    # sort loc_list and choose_list according to sort_list
    index_sort = [i[0] for i in sorted(enumerate(sort_list), key=lambda x:x[1])]
    loc_list_sort = []
    choose_list_sort = []
    for idx in index_sort:
        loc_list_sort.append(loc_list[idx])
        choose_list_sort.append(choose_list[idx])

    expression_data = expression_data[:, choose_list_sort]
    pd.DataFrame(expression_data, index=cell_list, columns=loc_list_sort).to_csv(loc_file)
