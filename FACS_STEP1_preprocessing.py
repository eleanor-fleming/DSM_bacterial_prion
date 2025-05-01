## PURPOSE

# This script generates processes FACS read counts data to make it ammenable to mean fluorescence intensity estimation using the
#fitdiscen R package.


## dependencies

# standard packages
import os
import warnings
import numpy as np
import pandas as pd
import functools as ft


## general set-up
warnings.simplefilter("ignore") # supress warnings
outdir = "./FACS_analysis_outputs/"
os.makedirs(outdir, exist_ok=True) #create directory for outputs


## user-specified variables
list_of_sample_names = [] # names (as strings) of FACS samples, eg ['bin1', 'bin2']
list_of_counts_paths = [] # paths (as strings) to csv files containing FACS seq data, one for each bin
list_of_cell_counts = [] # list of integers cooresponding to total cells sorted into each bin. in same order as sample_names
path_to_CVT = '' # path to csv file containing barcode-variant mapping data
list_of_cell_counts = [] # list of int values of cells sorted into each FACS bin
final_filename = '' # string specifying how to name final csv file (eg './050125_rep1_FACS_norm_pseudocounts1.csv')


## import FACS data: barcode count csv for each FACS bin sample
FACS_df = import_FACS_barcode_counts(list_of_sample_names, list_of_counts_paths)

## merge FACS dataframe with variant identities dataframe (CVT)
data = merge_with_CVT(path_to_CVT, FACS_df)

## normalize read counts across all FACS bins
data = normalize_read_counts(data,list_of_sample_names,list_of_cell_counts )

## add 1 pseudocount to each counts column to avoid barcode loss (see function description)
data = add_pseudocount(data, list_of_sample_names, 1)

# drop unneeded counts columns
data=drop_extra_counts_col(data,list_of_sample_names)

# add variant type column
data = add_type(data)

# save as csv
data.to_csv(final_filename)
