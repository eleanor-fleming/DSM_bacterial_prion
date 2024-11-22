## Functions for steps proceeding and through consensus building


## These functions collectively require the following dependencies:

# standard packages: os, warnings, numpy, pandas, plotnine, matplotlib.pyplot, 
# alignparse specific: alignparse.consensus, alignparse.minimap2, alignparse.constants, alignparse.targets
# dms_variants specific: dms_variants.plotnine_themes, dms_variants.utils, dms_variants.codonvarianttable


def align(target_names, alignment_parameters, run_names, library_name, lib_dir, align_outdir):
    
    """
    PURPOSE:
    This function preforms all the steps up to and including read alignment.
    
    INPUTS:
    1) "target_names" -- a list of target file names as strings with the .gb extension (eg "WT_coding.gb")
    2)  "alignment_parameters" -- a string file name for the .yaml file that specifies the alignment 
        parameters specific to your target (eg "R1_R2_final.yaml")
    3) "run_names": a list of run file names as strings each WITHOUT the .fasta extension (eg "lib1_R1") 
    4) "library_name": a string specifying which library the read files pertain to (eg "rep1")
    5) "library_directory": a string specifying the directory in which the read files are stored (eg "./input/rep1/")
    6) "align_outdir": a string specifying the directory in which to save alignment output files (eg "./lib1_alignemnt")
    
    OUTPUTS:
    1) readstats -- a dataframe summarizing counts of aligned, filtered, and unmapped reads
    2) aligned -- a dictionary in which keys are the target names and values are dataframes specifying information
        on each read successfully aligned to that target
    3) filtered -- a dictionary organized identically to the aligned variable, but pertaining to reads that did not
        align (eg reads that were filtered due to failure to meet the alignment parameters specified in the 
        "alignmentparameters_file")
    
    NOTES:
    Make sure you specify 3 variables when calling the function. 
    
    """
            
    # STEP 1: create target object
    targets = alignparse.targets.Targets(
    seqsfile=target_names,
    feature_parse_specs= alignment_parameters)
    print('Target created.')

    # STEP 2: create mapper object, using optimal parameters for amplicon alignment
    mapper = alignparse.minimap2.Mapper(alignparse.minimap2.OPTIONS_CODON_DMS)
    print(f"Using `minimap2` {mapper.version} with these options:\n" + " ".join(mapper.options))
    
    # STEP 3: create read files dataframe (required by align_and_parse function)
    libraries = [library_name] * len(run_names)
    df_runs = pd.DataFrame(
    {
        "name": run_names,
        "library": libraries,
        "fastq": [f"{library_directory}/{name}.fastq" for name in run_names],
    })

    # STEP 4: perform alignment
    %time readstats, aligned, filtered = targets.align_and_parse(
        df=df_runs,
        mapper=mapper,
        outdir=align_outdir,
        name_col="name",
        queryfile_col="fastq", 
        overwrite=True)
    
    return readstats, aligned, filtered


def pair_aligned_reads(aligned):
    
    """
    
    PURPOSE:
    This function re-pairs aligned reads based on read names. 
    
    INPUTS:
    1) 'aligned' -- a dictionary in which keys are the target names and values are dataframes specifying
        information on each read successfully aligned to that target. (This the 2nd output from the align function.)
    
    OUTPUTS:
    1) dataframe with columns specifying read names and additional target features (eg barcode sequence, barcode
        accuracy, gene sequence, gene accuracy), as specified in the target.gb files provided to the 
        align_and_parse function.

    """
    
    # get dfs out of aligned dictionary
    coding = aligned['WT_amplicon_coding']
    noncoding = aligned['WT_amplicon_noncoding']

    # merge on query name, only keeping pairs (name found in both dfs)
    paired = pd.merge(coding, noncoding, on=['query_name'])
    print('Number of reads aligned to coding strand: ' + str(coding.shape[0]))
    print('Number of reads aligned to noncoding strand: ' + str(noncoding.shape[0]))
    print('Number of paired aligned reads: ' + str(aligned_df.shape[0]))

    # clean up paired_df
    paired.drop(['name_x', 'name_y'], inplace=True, axis=1)
    paired = paired.rename(
        columns={'query_clip5_x':'R1_q5clip',
                 'query_clip3_x': 'R1_q3clip',
                 'query_clip5_y': 'R2_q5clip',
                 'query_clip3_y': 'R2_q3clip',
                 'barcode_sequence':'barcode'})
    
    return paired


def QC_filter(paired, error_cutoff, num_bins=25):
    
    """
    
    PURPOSE:
    This function performs quality filtering to identify reads likely to contain sequencing errors in
    features of interest in order to avoid false-positive variant calls in the consensus building step. As coded, 
    this function assumes the two features of interest are the gene and the barcode.
    
    INPUTS:
    1) "paired" -- the dataframe containing paried aligned reads and relevant feature information (output from
        the pair_aligned_reads function).
    2) "error_cutoff" -- numerical value indicating QC value to use as the filtering threshold (eg to retain
        reads with scores of QC30 and above, the error_cutoff would be 0.001, cooresponding to an allowance of
        1 mutation per 1000 nucleotides.
    3) "num_bins" -- integer specifying the number of bins (for binning error rates) to use for histogram. Default
        value is 25.
        
    VISUALIZATIONS:  
    1) histogram of read counts at binned error rates for each feature of interest.
    2) table displaying retained vs not retained read counts (based on error_cutoff).
    
    OUTPUTS:
    1)  modified input dataframe containing a column detailing "retained" status (to be used as a filter in 
        consensus building set that follows).
    
    """
    
    # calculate error rates (1-accuracy) for barcodes and genes
    paried = paired.assign(
        barcode_error=lambda x: np.clip(1 - x["barcode_accuracy"], 1e-7, None), # avoid 0 to accomodate log plot
        gene_error=lambda x: np.clip(1 - x["gene_accuracy"], 1e-7, None),
    )
    
    # mark reads that pass threshold ("error_cutoff") as "retained" (new col)
    paired = aligned_df.assign(
        retained=lambda x: ((x["gene_error"] <= error_cutoff) & (x["barcode_error"] <= error_cutoff)
        )
    )
    
    # Visualization 1: histogram of read counts at binned error rates for each feature of interest
    p = (
        ggplot(
            aligned_df.melt(
                value_vars=["barcode_error", "gene_error"],
                var_name="feature_type",
                value_name="error rate",
            ),
            aes("error rate"),
        )
        + geom_histogram(bins=num_bins)
        + geom_vline(xintercept=error_cutoff, linetype="dashed", color=CBPALETTE[1], size=2)
        + facet_grid(" ~ feature_type")
        + theme(figure_size=(4.5, 2))
        + ylab("number of reads")
        + scale_x_log10()
    )

    p
    
    # Visualization 2: 
    (
        paired.groupby(["retained"])
        .size()
        .rename("number of reads")
        .reset_index()
    )
    
    return paired
    
    
    
def filter_indels(paired):
    
    """
    
    PURPOSE:
    This function identifies reads containing indel(s) in any features of interest (in this case just the gene),
    such that indel-containing reads can be excluded during the consensus building step. 
    
    INPUTS:
    1) 
    
    The function returns
    a the paired dataframe with two new columns; 'n_indel_col' specifying the number of indels in the gene sequence,
    and 'has_indel', boolean whether each read contains at least one indel (True) or none (False).
    
    This function produces 1 visualization:
    
    1) table specifying reads counts with and without indels
    
    """
    
    # mark reads that have at least 1 indel (new col)
    paired = alignparse.consensus.add_mut_info_cols(
    paired, mutation_col="gene_mutations", n_indel_col="n_indels"
    )
    paired = paired.assign(has_indel=lambda x: x["n_indels"] > 0)
        
    # Visualization 1:table specifying reads counts with and without indels
    (
    paired.query("retained")
        .groupby(["has_indel"])
        .size()
        .rename("number_reads")
        .reset_index()
    )
    
    return paired


 def generate_consensus(paired, max_sub_diffs, max_indel_diffs, max_minor_sub_frac, max_minor_indel_frac, 
                        min_support, exclude_indels, library_name, WT_gene_sequence, num_bins=30):

    """
    
    PURPOSE:
    This function generates a consensus gene sequence for each barcode (when possible) from the paired reads 
    dataframe. 
    
    INPUTS:
    1) 'max_sub_diffs' -- specifying the maximum number of substitutions to allow in a gene sequence. Paired reads 
        (barcode, gene) exceeding this cutoff will be discared.
    2) "max_indel_diffs" -- similar to max_sub_diffs, the parameter specifies the maximum number of indels to allow.
    3) "max_minor_sub_frac" -- specifies the maximum fraction of reads that can disagree (different set of gene 
        mutations) with the majority of reads sharing the same barcode. Any barcode with a "minority fraction" 
        exceeding this cutoff will be discarded.
    4) "max_minor_indel_frac" -- similar to the max_minor_sub_frac, but pertaining to indels.
    5) "min_support" -- minimum number of reads each barcode must have before being considered.
    6) "exclude_indels" -- boolean specifying if you want to exclude barcodes pertaining to reads containing
        indels from dataset (following consensus generation).
    6) "library_name" -- string of name of library (eg "rep1") to which the reads belong. This value is used to
        generate custom file names for saving outputs.
    7) "WT_gene_sequence" -- string (all caps) specifying the DNA sequence of target gene.
    8) "num_bins" -- integer specifying the number of bins (for binning variant call support) to use for 
        histogram 1 (see below). Default value is 30.
    
    VISUALIZATIONS:
    1) histogram of barcode counts per variant call support bin (num reads used to build consensus).
    2) frequency of mutantions across gene length but mutation type (synonymous, misesnse, nonsense).
    3) histogram of number of variant possessing number of CODON muts (bins, eg 1 mutation ,2 mutation)
    
    OUTPUTS:
    1) "consensus" -- dataframe with the following columns; barcode sequence, variant call support (num of reads used 
        to build consensus), and various other attributes of the cooresponding gene sequence (
        eg codon substitution identity and and number of indels.) This dataframe is saved as a csv and is not 
        retured by the function.
    2) "barcode_variant" -- dataframe expands on the consensus dataframe, adding amino acid mutation identities (
        translations of the codon substitutions), number of codon changes, number of amino acid changes, and drops
        the number of indels column. This dataframe is saved as a csv and also returned by the function.
    
    """
    
    # build generate consensus gene sequence for each barcode, stored as df
    consensus, dropped = alignparse.consensus.simple_mutconsensus(
    paired.query("retained"),
    group_cols=("library", "barcode"),
    mutation_col="gene_mutations",
    max_sub_diffs=max_sub_diffs,
    max_indel_diffs=max_indel_diffs,
    max_minor_sub_frac=max_minor_sub_frac,
    max_minor_indel_frac=max_minor_indel_frac,
    min_support=min_support,
    )
    
    # add mutation information column to consensus df
    consensus = alignparse.consensus.add_mut_info_cols(
    consensus,
    mutation_col="gene_mutations",
    sub_str_col="substitutions",
    n_indel_col="number_of_indels",
    overwrite_cols=True,
    )
    
    # OPTIONAL filter indel-containing gene/barcodes from consensus df
    if exclude_indels == True:
        consensus = consensus.loc[consensus['number_of_indels'] == 0]
        
    # save consensus as csv (to use dmsvariants function)
    csv_filename = './' + library_name + '_consensus.csv'
    consensus.to_csv(csv_filename) 

    # translate codon mutations to amino acid mutations and clean up, using codonvarianttable function
    CVT = dms_variants.codonvarianttable.CodonVariantTable(
        barcode_variant_file=csv_filename, 
        geneseq=WT_gene_sequence, 
        substitutions_col='substitutions'
    )
    
    # create and save barcode variant dataframe
    barcode_variant_df = CVT.barcode_variant_df
    variant_df_filename = './' + library_name + '_bcvariants.csv' 
    barcode_variant_df.to_csv(variant_df_filename) # save barcode_variant_df as csv
    
    # Visualization 1: histogram barcode counts per variant call support bins (num reads used to build consensus)
    %matplotlib inline
    x = np.log(consensus['variant_call_support'].tolist())
    plt.hist(x, density=False, bins=num_bins)  # density=False would make counts
    plt.ylabel('No. barcode-variant pairs')
    plt.xlabel('No. reads used to build consensus (log)')
    
    # Visualization 2: frequency of mutantions across gene length but mutation type (synonymous, misesnse, nonsense)
    p = CVT.plotMutFreqs("single", "codon", samples=None)
    p
    
    # Visualization 3: histogram of number of variant possessing number of CODON muts 
    # (bins are preset, eg 1 mutation ,2 mutations, 3 mutations, etc.)
    p2 = CVT.plotNumMutsHistogram("codon", samples=None)
    p2 = p2 + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
    p2

    
    return consensus, barcode_variant



def distinguish_rare_codons(barcode_variant, num_bins=15):
    
    """
    
    PURPOSE:
    This function creates alternate amino acid symbols for any rare codons used in the library. As coded, this
    function assumes the only rare codons used are CTA for leucine and AGG for arginine. 
    
    INPUTS:
    1) "barcode_variant" -- dataframe with columns specifying barcode sequences, gene codon mutations, gene amino
        acid mutations. (This is 2nd output from the generate_consensus function).
    2) "num_bins" -- integer specifying the number of bins (for binning # barcodes per variant) to use for 
        histogram (see below). Default value is 15.
    
    VISUALIZATIONS:
    1) histogram of barcodes per single amino acid substitution variant.
    
    OUTPUTS:
    1) singles -- dataframe including only barcode-variants with single amino acid substitutions. This dataframe
        retains all the columns in the input dataframe, plus an additional column called 'new_aa_substitutions' 
        which distinguishes rare leucines with a lowercase l and rare arginines with a lowercase r.
    
    """
    
    # filter to only include single amino acid mutants
    singles = barcode_variant_df[barcode_variant_df.n_aa_substitutions == 1]

    # add rare codon information in new mutations column
    new_aa_substitutions = [] # to make new col from
    cta_count = 0 # rare codon counts
    agg_count = 0
    for index, row in singles.iterrows():
        codon = row['codon_substitutions']
    
        if 'CTA' in codon: # if it has a rare leu codon...
            cta_count +=1
            mut = row['aa_substitutions']
            new_aa_substitutions.append(mut[0:-1] + 'l') # symbol for rare Leu: lower case "l"
        
        elif 'AGG' in codon: # if has rare arg mutations
            agg_count +=1
            mut = row['aa_substitutions']
            new_aa_substitutions.append(mut[0:-1] + 'r') # symbol for rare Arg: lower case "r"
        
        else:
            new_aa_substitutions.append(row['aa_substitutions']) # if not rare, keep existing notation
    
    # add as new column to singles dataframe
    singles['new_aa_substitutions'] = new_aa_substitutions
    
    # Visualization 1: histogram of barcodes per variant
    bins = list(range(1,num_bins+1)) # set bins
    singles_bc_counts = []
    for mut in set(new_aa_substitutions): 
        singles_bc_counts.append(new_aa_substitutions.count(mut)) # get unique singles
    plt.hist(singles_bc_counts, bins=bins, align='left')
    plt.title('Barcodes per single mutant variant')
    plt.xlabel('Number of barcodes')
    plt.ylabel('Number of variants with this barcode count')
    plt.xticks(range(num_bins+1))
    plt.show()
    
    return singles
