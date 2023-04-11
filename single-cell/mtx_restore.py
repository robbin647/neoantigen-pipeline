import csv
import gzip
import os
import scipy.io

######### DEBUG BEGIN Yang Tianxia 20220826 ##########
import argparse

myParser = argparse.ArgumentParser("Create Single-Cell Matrix based on market matrix file")
myParser.add_argument('sample_id', metavar='SAMPLE_ID', type=str, help='the name of the sample')
myParser.add_argument('working_dir', metavar='WORKING_DIR', type=str, default=os.getcwd(), help='the directory of input and output files (Default is current working directory)')

# define sample id
args = myParser.parse_args()

sample_id = args.sample_id
######################################################

# define MEX directory
matrix_dir = os.path.abspath(args.working_dir)

#############DEBUG##############
#print(f"Your sample id: {sample_id}; your working dir: {matrix_dir}")
################################

# read in MEX format matrix as table
mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))
 
# list of transcript ids, e.g. 'ENSG00000243485'
features_path = os.path.join(matrix_dir, "features.tsv.gz")
feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of gene names, e.g. 'MIR1302-2HG'
gene_names = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of feature_types, e.g. 'Gene Expression'
feature_types = [row[2] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]

import pandas as pd
 
# transform table to pandas dataframe and label rows and columns
matrix = pd.DataFrame.sparse.from_spmatrix(mat)
matrix.columns = barcodes
matrix.insert(loc=0, column="GeneSymbol", value=gene_names)
 
# display matrix
print(matrix)
# save the table as a CSV (note the CSV will be a very large file)
matrix.to_csv(os.path.join(matrix_dir,f"{sample_id}.tsv"), index=False, sep='\t')
