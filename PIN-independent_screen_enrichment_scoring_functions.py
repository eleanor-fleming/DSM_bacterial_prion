## Functions for data processing and enrichment score calculation for PIN-independent screening data


## These functions collectively require the following dependencies:

# standard packages: warnings,pandas, functools


# functions

def import_barcode_counts(list_of_sample_names, list_of_counts_paths):
    
    """
    
    PURPOSE:
    This function imports the barcode counts dataframes (from a list of paths), modifies column names
    (to match CVT dataframe, needed later), and merges all counts dataframes into a single df.
    
    INPUTS:
    1) "list_of_sample_names" -- list of strings. Each string specifies the name of one screening sample.
        Will be used to rename the counts columns. 
       (Eg 'blue' --> col name will be 'blue_count')
    2) "list_of_counts_paths" --  list of strings, each specifying the path to the directory in which 
        the barcode counts csv is stored (eg ["./counts/replicate1/blue", "./counts/replicate1/pale"]).
        Must be in same order as list_of_sample_names.
    
    OUTPUTS:
    1) df -- a dataframe containing read counts for each barcode in each sample. 
       Compatible with CVT (can be merged).
    
    """
    
    to_merge = []
    
    for i in range(0, len(list_of_sample_names)):
        df = pd.read_csv(list_of_counts_paths[i]) # read in data file as pd df
        new_name = list_of_sample_names[i] + '_count'
        df = df.rename(columns={'Barcode': 'barcode', 'Count': new_name}) # make lowercase because that's what CVT requires
        df = df.drop('Unnamed: 0', axis=1) # drop unnamed col
        to_merge.append(df)
    
    merged = ft.reduce(lambda left, right: pd.merge(left, right, on='barcode'), to_merge)
    
    print('Imported counts data for ' + str(merged.shape[0]) + ' variant-barcodes.')
    
    return merged


def merge_with_CVT(path_to_CVT, df):
    
    """
    
    PURPOSE:
    This function imports the variant-barcode mapping dataframe ("CVT"), and merges it with all FACS sample dataframes.
    
    INPUTS:
    1) "path_to_CVT" --  a string specifying the path to the directory in which variant-barcode mapping dataframe
    ("CVT") csv is stored (eg "./replicate1/CVT").
    2)  "df" -- df with  all of the samples' data (output from 'import_barcode_counts' function).
    
    OUTPUTS:
    1) merged -- a dataframe containing barcode sequences, read counts for each barcode in each sample, 
    and variant mapped to each barcode. 
    
    """
    
    # import CVT, save only necessary columns
    CVT = pd.read_csv(path_to_CVT, index_col=0)
    
    CVT = CVT.loc[:, ['barcode', 
                     'codon_substitutions', 
                     'n_codon_substitutions', 
                     'aa_substitutions', 
                     'n_aa_substitutions']]
    
    
    merged = ft.reduce(lambda left, right: pd.merge(left, right, on='barcode'), [CVT, df])
    
    return merged


def normalize_read_counts(df, list_of_sample_names, cutoff=0):
    
    """
    
    PURPOSE: 
    Normalize to account for read depth and barcode composition difference across samples.
    
    INPUTS:
    1) "df" -- pandas dataframe containing read counts from all samples. This is 
        likely output from the 'merge_with_CVT' function.
    2) "list_of_sample_names" -- list of strings. Each string specifies the name of one sample. 
        Need to be the same names as used in the 'import_barcode_counts' function.
    3) "cutoff" -- int, specifying minimum number of reads an entry must have in a screening population to be
        considered for analysis. Below this value, read count will be reduced to 0. This is intended to limit the
        effects of sequencing noise.
    
    OUTPUTS:
    1) df_out -- same dataframe with new columns containing normalized read counts (one per original read col).
    
    """
    
    names = []
    for n in list_of_sample_names:
        names.append(n + '_count') # get name of column with read counts
        
    # make deep copy of data_df
    df_out = df.copy(deep=True)
    
    # ignore read values below 10 -> likely noisey
    for name in names:
        df_out[name] = df_out[name].mask(df_out[name] < cutoff, 0)
    
    # perform normalization
    for name in names:
        new_name = 'norm_' + name
        df_out[new_name] = (df_out[name]/df_out[name].sum())
        
    return df_out


def estimate_cell_count(df, list_of_sample_names, list_of_cell_counts):
    
    """
    
    PURPOSE: 
    Estimate number of cells per barcode based on normalized read counts.
    
    INPUTS:
    1) "df" -- pandas dataframe containing read counts from all samples. This is 
        likely output from the 'normalize_read_counts' function.
    2) "list_of_sample_names" -- list of strings. Each string specifies the name of one sample. 
        Need to be the same names as used in the 'import_barcode_counts' function.
    3) "list_of_cell_counts" -- list of integer values of cell counts (sorted cells) for each sample.
        Same order as list_of_sample_names.
        
    OUTPUTS:
    1) df_out -- same dataframe with new columns with estimated cell counts (one per original read col).  
    
    """
    
    names = []
    
    # make deep copy of data_df
    df_out = df.copy(deep=True)
    
    # get names of columns with read counts
    for n in list_of_sample_names:
        names.append('norm_' + n + '_count')
        
    # estimate cell count
    for i in range(0, len(names)):
        name = names[i]
        new_name = ('est_cells_' + list_of_sample_names[i])
        df_out[new_name] = df_out[name]*list_of_cell_counts[i]
    
    return df_out


def remove_compromised_barcodes(df, path_to_meanF_scores):
    
    """
    
    PURPOSE:
    This function removes variant-barcodes that have statistically different mean YFP fluorescence compared to
    the mean of the WT-barcode population.
    
    INPUTS:
    1) "df" -- pandas dataframe for which you would like to remove the compromised barcodes. Must contain 
        'barcode' column.
    2) "path_to_meanF_scores" -- string, specifying the path to the .csv file containing Z-test results for
        barcode-variants.
        
    OUTPUTS:
    1) df_out -- filtered dataframe containing only barcode-variants with "WT-like" mean YFP fluorescence.
    
    """
    
    # copy df
    temp = df.copy(deep=True)
    #print(str(temp.shape[0]) + ' barcodes before filtering.')
    
    # import Z-score dataframe
    scores = pd.read_csv(path_to_meanF_scores, index_col=0)
    
    # merge with data
    merge_on = ['barcode', 'n_codon_substitutions', 'codon_substitutions', 
                'aa_substitutions', 'n_aa_substitutions']
    
    merged = pd.merge(temp, scores, on=merge_on, how='outer')
    
    # filter out bad barcodes based on Z-test boolean
    df_out = merged[merged['reject_boolean'] != True]
    
    print('Removed compromised variant-barcodes. Number of variant-barcodes remaining: ' + str(df_out.shape[0]))
    
    return df_out
    
    
def add_type(df):
    
    """"
    
    PURPOSE:
    This function adds a new column to your dataframe, specifying the type of variant. 
    Options are: 'WT', 'synonymous', 'single', 'multiple', or 'SSB1-1'. SSB1-1 is our positive control variant.
    
    INPUTS:
    1) "df" -- pandas dataframe to which you would like to add Type col. Must contain 'n_condon_substitutions' 
        and 'n_aa_substitutions' columns.
    
    OUTPUTS:
    1) df_out -- same dataframe with added 'variant_type' column.
    
    """
    
    c = df['n_codon_substitutions'].tolist()
    a = df['n_aa_substitutions'].tolist()
    
    types = []
    for i in range(0, len(c)):
        
        if a[i] == 0:
            if c[i] == 0:
                types.append('WT')
            else:
                types.append('synonymous')
             
        elif a[i] == 1:
            types.append('single')
           
        elif a[i] >1:
            if a[i] >4:
                types.append('SSB1-1')
            else:
                types.append('multiple')
           
    df_out = df.copy(deep=True)
        
    df_out['variant_type'] = types
        
    return df_out
    
    
def group_by_variant(df, variant_type_filter=None):
    
    """
    
    PURPOSE:
    This function groups barcode-variants by variant identity (codon substitution).
    
    INPUTS: 
    1) "df" -- dataframe on which to perform grouping
    2) "variant_type_filter" -- string, option to filter by a variant type. Options include, "WT", "SSB1-1",
        "single", "multiple", "SSB1-1", or None. Note: "WT" will include synonymous mutants, and "multiple" 
        will NOT include SSB1-1. In case of 'all' (the default), no filtering will occur.
        
    OUTPUTS:
    1) "df_out" -- dataframe with grouped data (cell count, read count, aa_substitution, condon_substitution
        only).
    
    """
    
    # perform filtering
    if variant_type_filter == None:
        temp = df.copy(deep=True)
    elif variant_type_filter == 'WT':
        want = ['WT', 'synonymous']
        temp  = df[df['variant_type'].isin(want)]
    else:
        want = [variant_type_filter]
        temp  = df[df['variant_type'].isin(want)]
        
    # drop unneeded columns
    to_keep = ['codon_substitutions', 'aa_substitutions', 
               'est_cells_true_positive', 'est_cells_false_positive',
               'est_cells_pale' ]
    temp = temp[to_keep]
    
    # group by codon change
    df_out = temp.groupby(['codon_substitutions', 'aa_substitutions']).sum('numeric_only')
    
    # clean up output df
    df_out = df_out.reset_index()
    
    print('Grouped variant-barcodes by codon substitution. Number of grouped data entries: ' + str(df_out.shape[0]))

    return df_out


def drop_extra_columns(df):
    
    """
    
    PURPOSE:
    This function cleans up your dataframe by removing all columns extraneous to enrichment scoring function.
    This is specificically for use when you are NOT performing barocde-variant grouping.
    
    INPUT:
    1) "df" -- dataframe to clean. Assumes working with part 1 screening data.
    
    OUTPUT:
    2) "df_out" -- dataframe with only minimally required columns needed for enrichment scoring function:
        'barcode', 'codon_substitutions', 
        'aa_substitutions', 'est_cells_true_positive', 'est_cells_false_positive', 'est_cells_pale'
        
    """
    
    df_out = df.copy(deep=True)
    
    to_keep = ['barcode','codon_substitutions', 'aa_substitutions', 'est_cells_true_positive', 
               'est_cells_false_positive','est_cells_pale' ]
    
    df_out = df_out[to_keep]
    
    return df_out

def calculate_enrichment_score(df, weight=0, cutoff=0):
    
    """
    
    PURPOSE:
    This function calculates the enrichment score for each entry in your dataframe. The enrichment score is the
    fraction of true positive CFUs divided by the fraction of pale CFUs belonging to each entry in the dataframe.
    
    There is the option of calculating a "combination enrichment score" which also takes into account each entry's
    representation in the false positive group. This is designed to capture two types of events. First, is 
    counter-screening failures. If your counter screening technique is prone to false negatives, then consider 
    using a weight value >0 to account for reads in the false positive group. We see evidence of this happening 
    with our positive control. The second event is heterogeneity. Although our positive control presents with 100%
    blue (prion-containing) colonies, it is possible there are variants that exhibit different dynamics of 
    conversion, propagation and resolution, cumulatively resulting in <100% blue colonies. Such variants would be
    assumed to accumulate more false positive reads due to loss of the prion state during the extened counter-
    screening procedure.
    
    INPUTS:
    1) "df" -- pd dataframe with screening data, including cell counts for true postive and pale populations.
    2) "weight" -- float, indicating amount of weight to be given to the false positive sample when calculating
        the additive combination enrichment score. Default is 0 in which case the combination score is identical 
        to the original score.
    3) "cutoff" -- int, minimum number of true positive colonies required, to filter final output. Default is 0.
    
    OUTPUTS:
    1) "df_out" -- pd dataframe with an enrichment score and a combination enrichment score calculated for each 
        entry in the original dataframe. Additional new columns include the intermediate calculations such as
        'true_fraction', 'pale_fraction', and 'false_fraction'.
    
    """
    
    df_copy = df.copy(deep=True)
    
    # calc total cells in each screening population
    true_total = sum(df_copy['est_cells_true_positive'])
    false_total = sum(df_copy['est_cells_false_positive'])
    pale_total = sum(df_copy['est_cells_pale'])
    
    # calculate fraction representations
    df_copy['true_fraction'] = (df_copy['est_cells_true_positive']/true_total)
    df_copy['false_fraction'] = (df_copy['est_cells_false_positive']/false_total)
    df_copy['pale_fraction'] = (df_copy['est_cells_pale']/pale_total)
    
    # calculate enrichment scores
    df_copy['True_frac_enrichment'] = df_copy['true_fraction']/df_copy['pale_fraction']
    df_copy['False_frac_enrichment'] = df_copy['false_fraction']/df_copy['pale_fraction']

    # calculate combo enrichment score
    df_copy['Combo_enrichment'] = df_copy['True_frac_enrichment']+(weight*df_copy['False_frac_enrichment'])
    
    # filter by cutoff
    df_out = df_copy.loc[df_copy['est_cells_true_positive'] >= cutoff]
    
    print('Number of variants with colony count > cut-off: ' + str(df_out.shape[0]))
    
    # clean up output
    df_out = df_out.sort_values('Combo_enrichment', ascending=False)
    df_out = df_out.reset_index()
    
    return df_out
