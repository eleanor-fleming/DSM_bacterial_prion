## Functions for steps following MLE of mean fluorescence in R markdown


## These functions collectively require the following dependencies:

# standard packages: numpy, pandas, scipy stats, scipy.stats norm, statsmodels.stats.multitest, 


def get_variant_type(df, variant):
    
    """
    PURPOSE:
    This function filters a pd df of data to only include an individual variant type (eg single mutants, SSB1-1).
    
    INPUTS:
    1) "df" -- pd dataframe. Must include a column titled "variant_type."
    2) "variant" -- type of variant you want data for, as a string. Options include: "WT", "SSB1-1", "singles", or "mutliptles".
        The "WT" option will also return any synonymous mutants in the dataset.
        
    OUTPUT:
    1) "df_out" -- a filtered copy of the original dataframe, only containing rows pertaining to the variant type
        of interest.
        
    """
    
    if variant == 'WT':
        want = ['WT', 'synonymous']
    else:
        want = [variant]
        
    df_out = df.loc[df['variant_type'].isin(want)]
    df_out = df_out.reset_index(drop=True)
    
    return df_out


def calc_weighted_ave_meanF(df):
    
    """
       
    PURPOSE:
    This function the weighted average of meanF for a dataframe. Weighting is done according to cell count
    from the FACS sorting. It also calculates the standard error of the weighted mean.
    
    INPUTS:
    1) "df" -- pd dataframe. Must include columns titled "total_count" (some of cell counts across all FACS bins),
        "ML_meanF" -- (maximum likelihood estimate of mean YFP fluorescence as determined by fitdistcen package), and
        "var_ML_mean" -- (variance of MLE mean, also from fitdistcen package).

        
    OUTPUT:
    1) "df_out" -- a  copy of the original dataframe, with two additional columns; "weight" indicating relative weight
        each barcode receives in mean calculation, and 'var_weighted_meanF' indicating weighted variance of each
        barcode's mean flu value.
    2) "pop_mu" -- weighted average mean YFP fluorescence across population (all barcodes in input dataframe).
    3) "var_w" -- variance of "pop_mu".
       
    """
        
    # make a copy of df
    df_out = df.copy(deep=True)
    
    # calculate weights (cell count / total WT cell count across all barcodes)
    df_out['weight'] = df_out['total_count'] /(df_out['total_count'].sum())
    
    # calculate pop mu (weighted ave of all individual sample mus)
    pop_mu = sum(df_out['ML_meanF']*df_out['weight'])
    print('Weighted mean YFP fluorescence of population: ' + str(round(pop_mu,3)))
    
    # caclulate standard error of weighted mean
    df_out['var_weighted_meanF'] = (df_out['weight']**2)*df_out['var_ML_meanF']
    var_w = round(df_out['var_weighted_meanF'].sum(),6)
    print('Variance of population mean: ' + str(var_w))
    
    return df_out, pop_mu, var_w


def Z_test(df, pop_mu, var_w, alpha, filter_WT='no'):
    
    """
    
    PURPOSE:
    This function performs a Z-test (with correction for multiple testing) to determine variant-barcodes with 
    mean fluorescence values statistically different from the population mean calculated for all WT-barcodes. 
    
    INPUTS:
    1) "df" -- pd dataframe. Must include columns titled "ML_meanF" (maximum likelihood estimate of mean 
        YFP fluorescence as determined by fitdistcen package), and "var_ML_mean" (variance of MLE mean, 
        also from fitdistcen package).
    2)  "pop_mu" -- weighted average mean YFP fluorescence across WT-barcode population.
    3) "var_w" -- variance of "pop_mu".
    4) "filter_WT" -- option to exclude WT data from Z-test. Default is no.

    OUTPUT:
    1) "df_out" -- a  copy of the original dataframe, with the following additional columns; "Z_score" indicating 
        the Z_score for each variant-barcode, "p-value" containing the unadjusted p-value from each variant-barcode
        Z-test against the WT population mean, and "p-value_corrected" indicating the p-values adjusted for multiple
        testing using the Benjamini Hochberg method.
    
    """
    
    # make copy of df
    df_out = df.copy(deep=True)
    
    # get mean flu and var of mean values
    means = df['ML_meanF'].values 
    varss = df['var_ML_meanF'].values
    
    # calculate Z-score
    nums = means - pop_mu
    denoms = np.sqrt((varss + var_w))
    Z_score = nums/denoms
    df_out['Z_score'] = Z_score
    
    # calc p-value for two-tailed test
    p_value = stats.norm.sf(abs(Z_score))*2
    # correct for multiple testing using BH method which will limit false negatives (more concerning to me)
    p_corrected = multi.multipletests(p_value, 
                                      alpha=alpha, 
                                      method='fdr_bh', 
                                      maxiter=1, 
                                      is_sorted=False, 
                                      returnsorted=False) 
    
    # add p-value as col
    df_out['p-value'] = p_value
    df_out['p-value_corrected'] = p_corrected[1]
    df_out['reject_boolean'] = p_corrected[0]
    print('Number of barcodes with normal (False) and abberent (True) mean YFP fluorescence: ' + str(df_out.reject_boolean.value_counts()))
    
    # sort df by z-score before returning
    df_out = df_out.sort_values('p-value_corrected')
    
    return df_out 
    
