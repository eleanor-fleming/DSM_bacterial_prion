## PURPOSE

# This script uses variant-barcode read counts data from the PIN-independence screen to calculate a true positive population enrichment score for each data entry (either individual barcode-variants or variant-barcodes grouped by codon substitution).


# dependencies
import warnings
import pandas as pd
import functools as ft


## general set-up
warnings.simplefilter("ignore") # supress warnings


# user-specificed variables
list_of_sample_names = [] # names (as strings) of screening output populations (eg 'true_positive', 'false_positive', 'pale'

list_of_counts_paths = [] # paths (as strings) to csv files containing read count data, one for screening population

path_to_CVT = '' # path (as string) to csv file containing barcode-variant mapping data

path_to_meanF_scores = '' # path (as string) to csv file containing z-test results for each barcode-variant mean fluorescence

list_of_cell_counts = [] # list of int values of colonies collected in each screening output population

output_filename = '' # string specifying the name to give the enrichment score dataframe when saving as a csv file


# import barcode count csvs for three screen samples
counts_df = import_barcode_counts(list_of_sample_names, list_of_counts_paths)

# merge count data with variant identities dataframe
data = merge_with_CVT(path_to_CVT, counts_df)

# normalize read counts across all samples
data = normalize_read_counts(data,list_of_sample_names, cutoff=0)

# estimate cell count per barcode
data = estimate_cell_count(data, list_of_sample_names, list_of_cell_counts)

# remove compromised barcods
data = remove_compromised_barcodes(data, path_to_meanF_scores)

# group barcodes by codon substitution
grouped = group_by_variant(data, variant_type_filter=None)

# calculate enrichment score, no weight on false positives, cut-off of 1 true positive CFU
final = enrichment_score(grouped, weight, cutoff=0)

# save as csv
final.to_csv(output_filename)