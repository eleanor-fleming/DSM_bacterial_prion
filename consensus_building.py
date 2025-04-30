## PURPOSE

# This script generates the mappings of barcodes to gene sequences. It requires reads of sufficient length to cover both the barcode and gene sequences. The script aligns reads, performs quality filtering, then generates consensus gene sequences. Optional additional steps include filtering the dataset to only include single amino acid substitution mutants and signifying  rare codons (leucine and arginine) with distinct amino acid symbols.


## dependencies

# standard packages
import os
import warnings
import numpy as np
import pandas as pd
from plotnine import *
import matplotlib.pyplot as plt

# alignparse specific
import alignparse.consensus
import alignparse.minimap2
from alignparse.constants import CBPALETTE
import alignparse.targets

# dms_variants specific
import dms_variants.plotnine_themes
import dms_variants.utils
import dms_variants.codonvarianttable

# custom functions
from consensus_building_functions import align
from consensus_building_functions import pair_aligned_reads
from consensus_building_functions import QC_filter
from consensus_building_functions import filter_indels
from consensus_building_functions import generate_consensus
from consensus_building_functions import distinguish_rare_codons


## general set-up
warnings.simplefilter("ignore") # supress warnings
theme_set(dms_variants.plotnine_themes.theme_graygrid()) # set plotting theme
outdir = "./consensus_building/"
os.makedirs(outdir, exist_ok=True) #create directory for outputs


## user-specified variables

# alignment variables
target_names = [] # names (as strings) of targets (maps) to use in alignemnt, in list format
alignment_parameters = '' # name of .yaml file that specifies alignment parameters, string
run_names = [] # names (as strings) of sequencing runs, in list format
library_name = '' # name for library the reads come from (eg "replicate 1"), string
library_directory = '' # directory where fasta read files are stored, string
align_outdir = '' # directory in which to save alignment result files, string

# consensus generation variables
max_sub_diff = None             # do not cap the number of mutations allowed
max_indel_diffs = None          # do not cap the number of indels allowed
max_minor_sub_frac = 0.1        # see function notes for description
max_minor_indel_frac = 0.1      # see function notes for description
min_support = 10                # minimum number of reads required to build a consensus sequence
exclude_indels = True           # exclude indels from dataset (not of interest to us)
WT_gene_sequence = 'CAGAACTCCCAAAATGGCGGTAACTTCGGCAATCAGAACTCCAACTTCGGTGGGGGCAACTTCAACTCGCAGAATAACTTTAACGGTTACGGCTCCAACCAGGGTTACAATGGTCACAACTCGAACAATAACCAGAATTTCAACCAAAACCGCAACTCTGGCTCCAACTTCGGCTACTCTAACCAGAATAACAATAAGACCAACTTCTCTAACAACACTCCGTAA'

                    
## perform read alignment
readstats, aligned, filtered = align(target_names, alignment_parameters, run_names, library_name, library_directory, align_outdir)

## pair aligned reads
paired = pair_aligned_reads(aligned)

# perform quality filtering
pairedQC = QC_filter(paired, 1e-3) # suitable quality threshold for Illumina

# perform indel filtering
paired_final = filter_indels(pairedQC)

# calculate empirical accuracy of gene mutation calls on retained reads with no indels
paired_final['library']=library_name
alignparse.consensus.empirical_accuracy(
    paired_final.query("retained & not has_indel"), mutation_col="gene_mutations")

# generate consensus sequences from high-quality reads, saves resulting dataframe as csv file
barcode_variant = generate_consensus(paired_final, 
                   max_sub_diffs, 
                   max_indel_diffs, 
                   max_minor_sub_frac, 
                   max_minor_indel_frac, 
                   min_support, 
                   exclude_indels, 
                   library_name, 
                   WT_gene_sequence)

# mark rare codons with separate amino acid symbols (applies to single mutants only)
singles = distinguish_rare_codons(barcode_variant, num_bins=15, library_name)