#!/usr/bin/python3

import numpy as np
import pandas as pd
import argparse 

import os

######################### Global Variables ###################
myParser = argparse.ArgumentParser(description="Remove normal cell columns from matrix file")
myParser.add_argument("--list", required=True, metavar='LIST_FILE', dest='list_file', type=str, help='path to the csv file which is a list of normal cells')
myParser.add_argument("--in", required=True, metavar='INPUT_MATRIX_FILE', dest='input_matrix', type=str, help='a tsv file that contains input matrix')
myParser.add_argument("--out", required=True, metavar='OUTPUT_MATRIX_FILE', dest='output_matrix', type=str, help='a csv file that contains input matrix with columns deleted')

args = myParser.parse_args()

NORMAL_LIST = os.path.abspath(args.list_file)
INPUT_MATRIX = os.path.abspath(args.input_matrix)
OUTPUT_MATRIX = os.path.abspath(args.output_matrix)

print("Using list: " + NORMAL_LIST)
print("Using matrix file: " + INPUT_MATRIX)
print("Using output file: ", OUTPUT_MATRIX)

normal_list = pd.read_csv(NORMAL_LIST, header=0, sep=',')


######################### Main Function  ###################
# get to knoe the column header first
with open(INPUT_MATRIX, 'rt') as fin:
    a_line = fin.readline()
    col_list = a_line.split('\t')
    dtype_list = [('GeneSymbol', 'str')]
    for item in col_list:
        if item == 'GeneSymbol':
            continue
        dtype_list.append((item.strip(), 'int32'))

mtx = pd.read_csv(INPUT_MATRIX, sep='\t', header=0, dtype=dtype_list)

mtx.drop(axis=1, labels=normal_list.loc[:, '0'], inplace=True, errors="raise")
print('The shape of matrix after deleting normal cells: ', mtx.shape)
mtx.to_csv(OUTPUT_MATRIX, sep=',', index=False)
print(f"Successfully saved result to {OUTPUT_MATRIX}")