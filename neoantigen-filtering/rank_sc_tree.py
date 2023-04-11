
'''
A tree construction script to output newick trees


INPUT:
 a 2D peptide vs. cell type matrix that 
  1) comes from matrix correlation of the following two:
        peptide type vs. gene symbol matrix (== pVACseq table)
        cell type vs. gene symbol matrix (single-cell coverage matrux)
        
  2) is in csv format

OUTPUT:
  1) XXX.newick: the constructed tree from the matrix in newick format
  2) XXX.top_n_neoantigens.csv: the top n neoantigens in the tree
'''
import argparse
import os
# import umap
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances, cosine_similarity
from sklearn.preprocessing import StandardScaler
from collections import deque
from skbio import DistanceMatrix
from skbio.tree import nj

############################ Global Variable ###################################

myParser = argparse.ArgumentParser("Build newick tree from peptide type vs. gene symbol matrix")
myParser.add_argument("--in", required=True, metavar="IN_MATRIX", dest="file_in", type=str, help="The input peptide vs. gene symbol matrix file")
myParser.add_argument("--sample-id", required=True, metavar="SAMPLE_ID", dest="sample_id", type=str, help="The sample ID. To be used to name the output files")
myParser.add_argument("--dim", required=False, metavar="DIMENSIONS", dest="dimensions", type=int, help="Number of features to use when calculating distances. Default is to use all features")
myParser.add_argument("--max-sim", required=False, default="1e6", metavar="MAX_SIMILARITY", dest="max_sim", type=str, help="the maximum value allowed in similarity matrix. Default is 10^6")
myParser.add_argument("--top-n", required=False, default=10, metavar="NO_OF_TOP_NEOANTIGENS", dest="top_n_neoantigen", type=int, help="Number of top neoantigens to output. Default is 10")
myParser.add_argument("--top-n-out",required=False, default=None, metavar="TOP_N_NEOANTIGEN_OUT_FILE", dest="top_n_neoantigen_out", type=str, help="The output file for the top n neoantigens. Default is <sample ID>.top_n_neoantigens.txt")
myParser.add_argument("--newick-out", required=False, default=None, metavar="NEWICK_OUT_FILE", dest="newick_out", type=str, help="The output file for the newick tree. Default is <sample ID>.newick")
args = myParser.parse_args()

#SAMPLE_INPUT = r'D:\My_code\fyp\sample_tmp\OV002-T_R.average.filtered.csv'
SAMPLE_INPUT = args.file_in

if args.dimensions is not None:
    DIMENSIONS = int(args.dimensions) 
else:
    DIMENSIONS = None

if float(args.max_sim) != 1e6:
    MAX_SIMILARITY = float(args.max_simi)
else:
    MAX_SIMILARITY = 1e6 # the  maximum value allowed in similarity matrix


TOP_N_NEOANTIGEN = int(args.top_n_neoantigen) if args.top_n_neoantigen != 10 else 10

if args.top_n_neoantigen_out is not None:
    # create parent directory if not exist
    if not os.path.exists(os.path.dirname(args.top_n_neoantigen_out)):
        os.makedirs(os.path.dirname(args.top_n_neoantigen_out))
    TOP_N_NEOANTIGEN_OUT = os.path.abspath(args.top_n_neoantigen_out)
else:
    TOP_N_NEOANTIGEN_OUT = os.path.abspath(args.sample_id + ".top_n_neoantigens.txt")

if args.newick_out is not None:
    # create parent directory if not exist
    if not os.path.exists(os.path.dirname(args.newick_out)):
        os.makedirs(os.path.dirname(args.newick_out))
    NEWICK_OUT = os.path.abspath(args.newick_out)
else:
    NEWICK_OUT = os.path.abspath(args.sample_id + ".newick")
#############################################################################

def main():
    '''
    Step 1: take in the matrix
    '''
    ## transpose needed so that cell types are rows
    #M = pd.read_csv(tsv, keep_default_na=False).drop("GeneSymbol", axis=1).T
    M = pd.read_csv(SAMPLE_INPUT, keep_default_na=False).set_index('Peptide')
    feature_nums = len(M.columns)
    peptides = M.index 
    

    ## to center data and convert to unit variance (may take a long time)
    M = StandardScaler(with_mean=False).fit_transform(M)
    ## other dimensionality reduction
    # if DIMENSIONS is not None and  DIMENSIONS < feature_nums:
    #     M = umap.UMAP(n_components=DIMENSIONS).fit_transform(M)

    '''
    Step 2: Find an algorithm that
    2-1: takes in whole matrix
    2-2: do clustering on each row
    2-3: return a tree-like structure
    '''

    # compute Euclidean distance between each rows
    distance_matrix = euclidean_distances(M)
    # distance_matrix = cosine_similarity(M)
    distance_matrix[distance_matrix == 0] = 1 / MAX_SIMILARITY ## replace all zeros in distance_matrix
    distance_matrix = 1 / distance_matrix  ## get the inverse of euclidean distance (=similarity)
    diagonal_indices = np.diag_indices(M.shape[0]) ## ---\
    distance_matrix[diagonal_indices] = 0          ## ---/ set all diagonal elements to zero (to avoid DistanceMatrixError)
    distance_matrix = np.round(distance_matrix, 9) ## to solve "not symmetric" error

    dm = DistanceMatrix(distance_matrix, peptides)
    tree = nj(dm)
    
    '''
    Step 3:
        do a level order traversal for "tree" and print out the first TOP_N_NEOANTIGEN nodes
    '''
    with open(TOP_N_NEOANTIGEN_OUT, "wt") as neoantigen_out:
        print("Peptide\tLength\tParent", file=neoantigen_out)
        count = 0
        for node in tree.levelorder(include_self=True):
            if node.name is not None:
                print(f"{node.name}\t{node.length}\t{node.parent}", file=neoantigen_out)
                count += 1
            if count == TOP_N_NEOANTIGEN:
                break
    # visited_nodes = []
    # expanded_nodes = deque()
    # expanded_nodes.append(tree)
    # while len(visited_nodes) < NUM_NODES:
    #     node = expanded_nodes.popleft()
    #     if len(node.children) == 0: # leaf node
    #         visited_nodes.append(node)
    #     else:
    #         if node.name is not None: # named internal node
    #             visited_nodes.append(node)
    #         else: # unnamed internal node
    #             for child in node.children:
    #                 expanded_nodes.append(child)
    # for node in visited_nodes:
    #     print(f"MyName:{node.name}; MyLength:{node.length}; MyParent:{node.parent};")

    '''
    Step 4: Convert the tree to newick format
    '''
    tree.write(NEWICK_OUT)
    

if __name__ == '__main__':
    main()
