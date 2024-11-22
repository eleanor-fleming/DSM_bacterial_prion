# DMS_bacterial_prion

## DEEP MUTATIONAL SCANNING OF A BACTERIAL PRION PROTEIN

Analysis of deep mutational scanning of barcoded codon variants of Campylobacter hominis single-stranded DNA-binding (SSB) protein prion domain (PrD)

Study and analysis by Eleanor Fleming

## PURPOSE
The purpose of this project is to analyze data collected from a screen or selection of a barcoded deep mutational scanning library. In this case our plasmid library consists of all possible
single amino acid substitutions in a gene encoding a bacterial prion domain, each with at least one unique N16 barcode. The prion domain is expressed with an N-terminal monomeric yellow fluorescent protein (mYFP) fusion to enable estimation of protein stability via fluorescent activated cell sorting (FACS, see below). We screened our library in two related screens for two aspects of prion behavior; conversion (initial spontaneous appearance of amyloid conformation) and propagation (continued inheritance of amyloid form via self-templating reaction). 

## SUMMARY OF WORKFLOW
1) generate the barcode-gene variant mappings which includes generating high quailty consensus gene sequences. [consensus_building_functions.py and consensus_building.py]
2) estimate variant protein stability from barcode counts from FACS data [script in progress].
3) analyze barcode counts in conversion screen outputs to identify variants with enhanced conversion dynamics [script in progress].
4) analyze barcode counts in propagation screen outputs to identify variants with diminished and/or enhanced propagation dynamics [script in progress].

## RUNNING THE ANALYSIS
These scripts make use of functions from the dms_variants and alignparse Python packages written by the Bloom lab. Instructions for how to build the required computing environment are found here: https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS/blob/master/README.md 
