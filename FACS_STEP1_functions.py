## Functions for steps proceeding MLE of mean fluorescence in R markdown


## These functions collectively require the following dependencies:

# standard packages: os, warnings, numpy, pandas, functools



def import_FACS_barcode_counts(list_of_sample_names, list_of_counts_paths):
    
    """
    
    PURPOSE:
    This function imports the barcode counts dataframes (from a list of paths), modifies column names
    (to match CVT dataframe, needed later), and merges all counts dataframes into a single df.
    
    INPUTS:
    1) "list_of_sample_names" -- list of strings. Each string specifies the name of one FACS sample.
        Will be used to rename the counts columns. 
       (Eg 'bin1' --> col name will be 'bin1_counts')
    2) "list_of_counts_paths" --  list of strings, each specifying the path to the directory in which 
        the barcode counts csv is stored (eg ["./counts/replicate1/bin1", "./counts/replicate1/bin2"]).
        Must be in same order as list_of_sample_names.
    
    
    OUTPUTS:
    1) df -- a dataframe containing read counts for each barcode in each FACS sample. 
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
    
    return merged


def merge_with_CVT(path_to_CVT, df):
    
    """
    
    PURPOSE:
    This function imports the variant-barcode mapping dataframe ("CVT"), and merges it with all FACS sample dataframes.
    
    INPUTS:
    1) "path_to_CVT" --  a string specifying the path to the directory in which variant-barcode mapping dataframe
    ("CVT") csv is stored (eg "./replicate1/CVT").
    2)  "df" -- df with  all of the FACS sample dataframes (output from 'import_FACS_barcode_counts' function).
    
    
    OUTPUTS:
    1) merged -- a dataframe containing barcode sequences, read counts for each barcode in each FACS sample, 
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


def normalize_read_counts(df, list_of_sample_names):
    
    """
    
    PURPOSE: 
    Normalize to account for read depth and barcode composition difference across samples.
    
    INPUTS:
    1) "df" -- pandas dataframe containing read counts from all FACS bins. This is 
        likely output from the 'merge_all_FACS_data' function.
    2) "list_of_sample_names" -- list of strings. Each string specifies the name of one FACS sample. 
        Need to be the same names as used in the 'import_FACS_barcode_counts' function.
    
    OUTPUTS:
    1) df_out -- same dataframe with new columns with normalized read counts (one per original read col).
    
    """
    
    # make deep copy of data_df
    df_out = df.copy()
    
    # perform normalization
    for n in list_of_sample_names:
        name = n + '_count' # get name of column with read counts
        new_name = 'norm_' + name
        df_out[new_name] = (df_out[name]/df_out[name].sum())
        
    return df_out


def estimate_cell_count(df, list_of_sample_names, list_of_cell_counts):
    
    """
    
    PURPOSE: 
    Estimate number of cells per barcode based on normalized read counts.
    
    INPUTS:
    1) "df" -- pandas dataframe containing read counts from all FACS bins. This is 
        likely output from the 'normalize_read_counts' function.
    2) "list_of_sample_names" -- list of strings. Each string specifies the name of one FACS sample. 
        Need to be the same names as used in the 'import_FACS_barcode_counts' function.
    3) "list_of_cell_counts" -- list of integer values of cell counts (sorted cells) for each FACS bin.
        Same order as list_of_sample_names.
        
    OUTPUTS:
    1) df_out -- same dataframe with new columns with estimated cell counts (one per original read col).  
    
    """
    
    names = []
    
    # make deep copy of data_df
    df_out = df.copy()
    
    # get names of columns with read counts
    for n in list_of_sample_names:
        names.append('norm_' + n + '_count')
        
    # estimate cell count
    for i in range(0, len(names)):
        name = names[i]
        new_name = ('est_cells_' + list_of_sample_names[i])
        df_out[new_name] = df_out[name]*list_of_cell_counts[i]
        
    return df_out
    


def add_pseudocount(df, list_of_sample_names, p):
    
    """
    
    PURPOSE: 
    Add pseudocount(s) to all normalized read columns to avoid data being tossed from 0 values. 
    The MLE implenetation requires at least 1 read in 2 bins. By adding psuedocounts, we avoid data loss.
    
    INPUTS:
    1) "df" -- pandas dataframe containing normalized read counts from all FACS bins. This is 
        likely output from the 'estimate_cell_counts' function.
    2) "list_of_sample_names" -- list of strings. Each string specifies the name of one FACS sample. 
        Need to be the same names as used in the 'import_FACS_barcode_counts' function.
    3) "p" -- int, psuedocount(s) to be added
    
    OUTPUTS:
    1) df_out -- same dataframe with new columns with normalized read counts + pseudocount (one per original read col).
    
    """
    
    #set up variables
    names = []
    
    # deep copy df
    df_out = df.copy()
    
    # get names of columns with normalized read counts
    for n in list_of_sample_names:
        names.append('est_cells_' + n)
        
    # add pseudocount
    for i in range(0, len(names)):
        name = names[i]
        new_name = 'pseudo_'+ name
        df_out[new_name] = df_out[name] + p
        df_out[new_name] = df_out[new_name].astype(int) # convert to ints
        
    return df_out



def drop_extra_counts_col(df,list_of_sample_names):
    
    """
    
    PURPOSE: 
    Remove unneeded FACS barcode counts columns (eg original counts col, and normalized counts col).
    This function assumes you have gone through the normalization and pseudocounting steps already.
    
    INPUTS:
    1) "df" -- pandas dataframe containing normalized and pseudocounted read counts from all FACS bins. This is 
        likely output from the 'add_pseudocount' function.
    2) "list_of_sample_names" -- list of strings. Each string specifies the name of one FACS sample. 
        Need to be the same names as used in the 'import_FACS_barcode_counts' function.
    
    
    OUTPUTS:
    1) df_out -- same dataframe with only one counts column per FACS sample. This is the normalized+peudocount added
        column. It will be renamed by sample name, eg 'bin1_cell_count'.
    
    """
    #set up variables
    final_names = {}
    to_remove = []
    
    # deep copy df
    df_out = df.copy()
    
    # get names of columns to remove
    for n in list_of_sample_names:
        to_remove.append(n + '_count') # original data
        to_remove.append('norm_' + n + '_count') # normalized only
        to_remove.append('est_cells_' + n) # normalized only
        k = 'pseudo_est_cells_' + n
        final_names[k] = n + '_cell_count' # dictionary for renaming pseudocounts columns
        
    df_out = df_out.drop(columns=to_remove, errors='ignore')
    df_out = df_out.rename(columns=final_names)
    
    return df_out


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
           
    df_out = df.copy()
        
    df_out['variant_type'] = types
        
    return df_out