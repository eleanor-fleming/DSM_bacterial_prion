## PURPOSE

# This script performs Z-tests to compare all variant-barcode mean YFP fluorescence estimates (meanF) to that of the WT population (all barcoded WT variants). 


## dependencies
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import norm
import statsmodels.stats.multitest as multi


## user-specified variables
input_data_path = '' # path to csv file containing barcode-variant mean fluorescence data
alpha = 0.05 # alpha value to use for Z-tests
output_filename = '' # string specifying how to name final csv file (eg 'date_replicate_FACS_Ztest.csv')


# import mean flu data
data = pd.read_csv(input_data_path)

# filter for WT only
WT_data = get_variant_type(data, 'WT')

# calculate WT population mean flu
WT_data, pop_mu, var_w = calc_weighted_ave_meanF(WT_data)

# perform Z-test
data = Z_test(data, pop_mu, var_w, alpha, filter_WT='no')

# save Z-test results as csv
data.to_csv(output_filename)

