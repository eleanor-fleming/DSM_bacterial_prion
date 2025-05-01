## Functions for steps proceeding MLE of mean fluorescence in R markdown


## These functions collectively require the following dependencies:

# standard packages: os, warnings, numpy, pandas, functools

def import_bccounts(sample_name, counts_path):
    
    """
    
    PURPOSE:
    This function imports the barcode counts dataframe and modifies column names to match CVT dataframes.
    
    INPUTS:
    1) "sample_name" -- string specifying the name of the sample. Will be used to rename the counts col. 
       (Eg 'bin1' --> col name will be 'bin1_counts')
    2) "counts_path" --  a string specifying the path to the directory in which the barcode counts csv is stored 
       (eg "./counts/replicate1/bin2")
    
    
    OUTPUTS:
    1) df -- a dataframe containing read counts for each barcode. Compatible with CVT (can be merged).
    
    
    """
  
    df = pd.read_csv(counts_path) # read in data file as pd df
    new_name = sample_name + '_count'
    df = df.rename(columns={'Barcode': 'barcode', 'Count': new_name}) # make lowercase because that's what CVT requires
    df = df.drop('Unnamed: 0', axis=1) # drop unnamed col
    
    return df


def merge_all_FACs(path_to_CVT, list_of_counts_dfs):
    
    """
    
    PURPOSE:
    This function imports the variant-barcode mapping dataframe ("CVT"), and merges it with all FACS sample dataframes.
    
    INPUTS:
    1) "path_to_CVT" --  a string specifying the path to the directory in which variant-barcode mapping dataframe ("CVT") csv is 
    stored (eg "./replicate1/CVT")
    2)  "list_of_counts_df" -- a list specifying  all of the FACS sample dataframes (eg [bin1_df, bin2_df]
    
    
    OUTPUTS:
    1) df -- a dataframe containing barcode sequences, read counts for each barcode in each FACS sample, and variant mapped to 
    each barcode. 
    
    """
    
    # import CVT, save only necessary columns
    CVT = pd.read_csv(path_to_CVT, index_col=0)
    
    df = CVT.loc[:, ['barcode', 
                     'codon_substitutions', 
                     'n_codon_substitutions', 
                     'aa_substitutions', 
                     'n_aa_substitutions']]
    
    list_of_counts_dfs.append(df)
    merged = ft.reduce(lambda left, right: pd.merge(left, right, on='barcode'), list_of_counts_dfs)
    
    
    return merged

# function for adding variant type
def add_type(df):
    
    """"
    
    PURPOSE:
    This function adds a new column to your dataframe, specifying the type of variant. 
    Options are: 'WT', 'synonymous', 'single', 'multiple', or 'SSB1-1'. SSB1-1 is our positive control variant.
    
    INPUTS:
    1) "df" -- pandas dataframe to which you would like to add Type col. Must contain 'n_condon_substitutions' 
        and 'n_aa_substitutions' columns.
    
    OUTPUTS:
    1) new_df -- same dataframe with added 'variant_type' column.
    
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
           
    df_c = df.copy()
        
    df_c['variant_type'] = types
        
    return df_c

def normalize_read_counts(data_df, list_of_sample_names, list_of_cell_counts):
    
    """
    
    PURPOSE: 
    Normalization done such that read:CFU equivalency is the same across all samples. 
    
    INPUTS:
    1) "data_df" -- pandas dataframe containing read counts from all FACS bins. This is 
    likely output from the 'merge_all_FACS_data' function.
    2) "list_of_sample_names" -- list of strings. Each string specifies the name of one FACS sample. 
        Need to be the same names as used in the 'import_bccounts' function.
    3) "list_of_cell_counts" -- list of integer values of cell counts (sorted cells) for each FACS bin.
        Same order as list_of_sample_names
    
    OUTPUTS:
    1) df -- same dataframe with new columns with normalized read counts (one per original read col).
    
    """
    
    # set up variables
    names = []
    per_CFUs = []
    norms = {}
    
    # make deep copy of data_df
    df = data_df.copy()
    
    # get names of columns with read counts
    for n in list_of_sample_names:
        names.append(n + '_count')
    
    # calc expected value of reads per cfu for each sample
    for i in range(0, len(names)):
        name = names[i]
        per_CFUs.append(sum(df[name])/list_of_cell_counts[i])
        
    # perform read count normalization and clean bc count dfs
    low = min(per_CFUs) # get lowest reads per cfu
    print('lowest value of reads per CFU: ' + str(low))
    
    # calculate each sample's normalization factor by dividng lowest by expected value of sample                   
    for i in range(0, len(names)):
        name = names[i]
        norms[name] = low/per_CFUs[i] 
    print('Normalization value: ' + str(norms))
    
    # perform normalization
    for name in names:
        new_name = ('norm_' + name)
        df[new_name] = round(df[name].apply(lambda x: x*norms[name]),1) # add normalized read count column
      
    return df

def add_pseudocount(norm_df, list_of_sample_names, p):
    
    """
    
    PURPOSE: 
    Add pseudocount(s) to all normalized read columns to avoid data being tossed from 0 values. 
    The MLE implenetation requires at least 1 read in 2 bins. By adding psuedocounts, we avoid data loss.
    
    INPUTS:
    1) "norm_df" -- pandas dataframe containing normalized read counts from all FACS bins. This is 
    likely output from the 'normalize_read_counts' function.
    2) "list_of_sample_names" -- list of strings. Each string specifies the name of one FACS sample. 
        Need to be the same names as used in the 'import_bccounts' function.
    3) "p" -- int, psuedocount(s) to be added
    
    OUTPUTS:
    1) df -- same dataframe with new columns with normalized read counts + pseudocount (one per original read col).
    
    """
    
    #set up variables
    names = []
    
    
    # deep copy df
    df = norm_df.copy()
    
    # get names of columns with normalized read counts
    for n in list_of_sample_names:
        names.append('norm_' + n + '_count')
        
    
    for i in range(0, len(names)):
        name = names[i]
        new_name = ('pseudo_'+ name +'_count')
        df[new_name] = df[name] + p
        df[new_name] = df[new_name].astype(int) # convert to ints
        
    return df
    
    

