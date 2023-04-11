#!/usr/bin/python3

## Input: a tsv file 
# 
# # Output: another tsv file 
# 
#  

from genericpath import isfile
import tarfile
import numpy as np
import pandas as pd
import os
import gzip 
import argparse

################ Global Variables ################################
myParser = argparse.ArgumentParser(description="Doing average on Gene-Celltype matrix")
myParser.add_argument('--id', required=True, metavar='SAMPLE_ID', dest='sample_id', type=str, help='the name of the sample')
myParser.add_argument('--in', required=True, metavar='INPUT', dest='input_path',type=str, help='path to input file (can be relative or absolute path)')
myParser.add_argument('--tmp', required=True, metavar='TMP_DIR', dest='temp_dir', type=str, help='temporary directory for intermediary files')
myParser.add_argument('--out', required=True, metavar='OUTPUT', dest='output_path', type=str, help='path to output file (can be relative or absolute path)')

args = myParser.parse_args()

sample_id = args.sample_id
input_path = os.path.abspath(args.input_path)
output_path = os.path.abspath(args.output_path)
temp_dir = os.path.abspath(args.temp_dir)
tmp_files = []

################### Utility Functions ########################
# Determine if two gene names are considered the "same".
# Two gene names are the same if they both contain a dot '.'
# and all parts before the last dot is exactly the same 
#
# e.g. 'A1CF.1' and 'A1CF.2' are the same (common part: A1CF, 
# different after the first dot)
#       'A1CF.1' and 'A1CC.1' are NOT the same (no even have
# a common part before the dot )
#
def is_same_symbol(s1, s2):
    '''Input: s1: str
              s2: str
       Return bool
    '''
    if s1 == s2:
        return True
    if s1.rfind('.') > -1 and s2.rfind('.') > -1: # find the last occurrence of dot
        #if s1.split('.')[0] == s2.split('.')[0]: # if everything before the last dot is the same
        if s1[:s1.rfind('.')] == s2[:s2.rfind('.')]:
            return True
    return False

# Given a dataframe ``df``, the column on which to search ``col_name``
# and the index from where to start checking `start_idx``,
# return a tuple of positional indices (start_idx, end_idx)
# ``start_idx``, ``end_idx``: positional index to determine the start and ending 
# indices of consecutive duplicate records in the matrix (including the beginning and ending indices)
def get_consec_dup_index(df, col_name,start_idx):
    '''Input: df: pandas.DataFrame
              col_name: str
              start_idx: int
        Output: 2-Tuple of (start_idx, end_idx)
    '''
    end_idx = start_idx + 1
    col_to_search = df[col_name]
    while end_idx < len(col_to_search) and is_same_symbol(col_to_search.iloc[start_idx], col_to_search.iloc[end_idx]):
        end_idx += 1
    return start_idx, end_idx-1

# Given a GeneSymbol (a string), just trim everything after the last dot 
def cleanse_genesymbol(name):
    '''
    Input: 
       name : str the GeneSymbol 
    Output: str
    '''
    if name.rfind('.') > 0:
        return name[:name.rfind('.')]
    return name
  
#################### Main Function ###############################

if os.path.isdir(os.path.abspath(temp_dir)): # change working directory
  os.chdir(os.path.abspath(temp_dir))

if args.input_path.endswith('.tar.gz'):
  ### unzip the tar file and place the file temporarily under the tmp directory
  tar = tarfile.open(input_path, mode='r:gz')
  tar.extractall()
  tar.close()
  ### change the input_path to the decompressed file
  tmp_files.append(os.path.join(temp_dir, os.path.basename(args.input_path)[:-7]))
  input_path = os.path.join(temp_dir, os.path.basename(args.input_path)[:-7])


# read the first line (all columns)
with open(input_path, 'rt') as fin:
    header = fin.readline()
    col_list = header.split('\t')
    dtype_list = [(col_list[0], 'object')] + [(col_list[i], 'int') for i in range(1, len(col_list))]

mtx = pd.read_csv(input_path, sep='\t', header=0, dtype=dtype_list)

print("Source matrix", mtx)

# Must sort the matrix by the index col (GeneSymbol)
mtx.sort_values(axis='index', by='GeneSymbol', inplace=True)

# Pandas dtype option accepts a dictionary
# So I need to construct a dictionary 
dtype_option = {'GeneSymbol': "object"}
for col in mtx.columns:
    if col == 'GeneSymbol': 
        continue
    dtype_option[col] = 'int32'
        
mtx_result = pd.DataFrame(columns=mtx.columns)
mtx_result = mtx_result.astype(dtype_option, copy=False) # set column data types

cursor = 0
# all arithmetic operation will be done over all columns except 'GeneSymbol'
operate_on_columns = col_list.copy()
operate_on_columns.remove('GeneSymbol') 

# iterate along the 'GeneSymbol' column
# stop and do average when "same" GeneSymbols are encountered

while cursor < mtx.shape[0]:
    (start, end) = get_consec_dup_index(mtx, 'GeneSymbol', cursor)
    if end > start: ### we meet a consecutive duplicate region
        # Apply the column-wise sum to these rows
        tmp = mtx.iloc[start].drop('GeneSymbol', inplace=False) ## type=Series
                
        for i in np.arange(start+1, end+1):
            tmp = tmp + (mtx.iloc[i].drop('GeneSymbol', inplace=False)) ## type=Series
        tmp = tmp / (end - start + 1) ## type=Series
        # cleanse the 'GeneSymbol' name
        new_genesymbol_name = cleanse_genesymbol((mtx.iloc[start])['GeneSymbol']) ## type=str
        # Construct a new row with stripped GeneSymbol, and the average values
        tmp['GeneSymbol'] = new_genesymbol_name ## type=Series
        a_new_row = pd.DataFrame(tmp.to_dict(), columns=mtx.columns, index=[0,]) ## type=DataFrame
        a_new_row = a_new_row.astype(dtype_option, copy=False) ## type=DataFrame
        
        # Append new row to the result matrix
        mtx_result = pd.concat([mtx_result, a_new_row], ignore_index=True)
    else:
        row_to_append = mtx.take([cursor,], axis=0)
        # even if no duplicate, still need to cleanse the 'GeneSymbol' name 
        row_to_append.iat[0, 0] = cleanse_genesymbol(row_to_append.iat[0, 0]) ## assume [0,0] is the GeneSymbol of the first row
        mtx_result = pd.concat([mtx_result,row_to_append ], ignore_index=True)
    cursor = end + 1
    
print('The averaged matrix has shape of', mtx_result.shape)

# Save output to tsv file
mtx_result.to_csv(output_path, sep='\t', index=False)
# clean temp dir
if len(tmp_files) > 0:
    for f in tmp_files:
        os.unlink(f)
        print(f"Cleaned up {f}")
