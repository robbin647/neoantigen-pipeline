USAGE = """
Replace the 'GeneSymbol' column of the Single Cell matrix 
    with the 'Best Peptide' column of the pVACseq table of predicted peptides 
e.g. In a Single Cell matrix : 
   GeneSymbol  Barcode1  Barcode2  Barcode3 
   Gene1       1         2         3 
   Gene2       4         5         6 
   Gene3       7         8         9 
And if the pVACseq table of predicted peptides is:
   Gene  ... Best Peptide
   Gene2 ... AAAAA
   Gene3 ... BBBBB
Then the output matrix will be:
   Peptide  Barcode1  Barcode2  Barcode3
   AAAAA    4         5         6
   BBBBB    7         8         9
   (Note that 'Gene1' does not exist in the pVACseq table, so it is not included in the output matrix)
"""

import pandas as pd
import os
from argparse import ArgumentParser

############### Get command line arguments ################
myParser = ArgumentParser(description="Inner join Single Cell matrix with pVACseq gene-peptides table on genes")
myParser.add_argument('--peptide', action='store', type=str, required=True, help='pVACseq gene-peptides table')
myParser.add_argument('--sc', action='store', type=str, required=True, help='Single Cell gene-by-cell matrix')
myParser.add_argument('--out', action='store', type=str, required=True, help='Output CSV of the peptide-cell matrix')
myParser.add_argument('--filtered-aggr', action='store', type=str, required=False, default="", help='Optional. If this option is set, will write the pVACseq gene-peptides table filtered to contain only the neoantigens that are in the Single Cell matrix')

args = myParser.parse_args()
if not os.path.exists(args.peptide):
    raise FileNotFoundError('Cannot find the pVACseq gene-peptides table:', args.peptide)
if not os.path.exists(args.sc):
    raise FileNotFoundError('Cannot find the Single Cell gene-by-cell matrix:', args.sc)
AGGREGATED_TSV = os.path.abspath(args.peptide)
SC_FILTERED = os.path.abspath(args.sc)
SC_PEP_CELL = os.path.abspath(args.out)
if args.filtered_aggr != "":
    OUT_FILTERED_AGGR_TABLE = os.path.abspath(args.filtered_aggr)
###############################

def matrix_correlate(AGGREGATED_TSV, SC_FILTERED, SC_PEP_CELL):
    '''Input
         AGGREGATED_TSV: pVACseq gene-peptides table (contains candidate neoantigens)
         SC_FILTERED: Single Cell gene-by-cell coverage matrix
       Output
         SC_PEP_CELL: is the single cell gene-by-cell coverage matrix where "GeneSymbol" column values are replaced by the neoantigen peptides
    '''
    aggregated_table = pd.read_table(AGGREGATED_TSV, sep='\t')
    single_cell = pd.read_csv(SC_FILTERED)
    genes = aggregated_table.loc[:, 'Gene'].astype('str').to_list()
    genes_and_bp = aggregated_table.loc[:, ['Gene', 'Best Peptide']].astype('str')
    sc_genes = single_cell.loc[:, 'GeneSymbol'].astype('str').to_list()
    # find the intersection of genes and sc_genes
    common_genes = list(set(genes).intersection(set(sc_genes)))
    # find the rows of single_cell that contain the common genes
    common_rows = []
    for i, r in single_cell.iterrows():
        if r['GeneSymbol'] in common_genes:
            tmp = single_cell.loc[i, :]
            tmp.loc['GeneSymbol'] = genes_and_bp.loc[genes_and_bp['Gene'] == r['GeneSymbol'], 'Best Peptide'].to_list()[0]
            common_rows.append(tmp)
    sc_pep_cell = pd.DataFrame(common_rows).rename(columns={'GeneSymbol': 'Peptide'})
    sc_pep_cell.to_csv(SC_PEP_CELL, sep=',', index=False)

def extract_peptide(AGGREGATED_TSV, SC_PEP_CELL, OUT_FILTERED_AGGR_TABLE):
    '''Input:
         AGGREGATED_TSV: pVACseq gene-peptides table (contains candidate neoantigens)
         SC_PEP_CELL: is the single cell gene-by-cell coverage matrix where "GeneSymbol" column values are replaced by the neoantigen peptides  
        Output:
        OUT_FILTERED_AGGR_TABLE: is the pVACseq gene-peptides table filtered to contain only the neoantigens that are in the Single Cell matrix

       Map the filetered neoantigens in the peptide-cell matrix back to the pVACseq gene-peptides table.
       And write out each row in the pVACseq gene-peptides table whose 'Best Peptide' field matches the filtered neoantigen.
    '''
    aggregated_table = pd.read_table(AGGREGATED_TSV, sep='\t')
    single_cell = pd.read_csv(SC_PEP_CELL)
    peptides = single_cell.loc[:, 'Peptide'].astype('str').to_list()
    # find the rows of aggregated_table that contain the peptides
    common_rows = []
    for i, r in aggregated_table.iterrows():
        if r['Best Peptide'] in peptides:
            common_rows.append(aggregated_table.loc[i, :])
    filtered_aggregated_table = pd.DataFrame(common_rows)
    # rename all columns in filtered_aggregated_table to remove spaces (required for visualization)
    filtered_aggregated_table.columns = [c.replace(' ', '_') for c in filtered_aggregated_table.columns]
    filtered_aggregated_table.to_csv(OUT_FILTERED_AGGR_TABLE, sep='\t', index=False)

if __name__ == '__main__':
    matrix_correlate(AGGREGATED_TSV, SC_FILTERED, SC_PEP_CELL)
    if OUT_FILTERED_AGGR_TABLE!= "":
        extract_peptide(AGGREGATED_TSV, SC_PEP_CELL, OUT_FILTERED_AGGR_TABLE)
