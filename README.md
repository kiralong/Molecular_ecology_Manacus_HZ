# README - Molecular Ecology of a Manacus hybrid zone
Code repository for the article "Telomeres, hybridization, and genome-wide heterozygosity in a manakin hybrid zone"

Authors: Ben J. Vernasco, Kira M. Long, Michael J. Braun, and Jeffrey D. Brawn

Primary Contacts: benvernasco@gmail.com and kiralong778@gmail.com

## Description of each script

# calculate_idnv-obs_het.py
A custom python script that calculates both the SNP only and all sites heterozygosity ratio for each indiviual.

# filter_sumstats_to_whitelist.py
A custom  python script that allows the user to filter snps for generating a whitelist of snps to use in the `stacks` software.

# format_structurefile_for_gghybrid.sh
A Bash script that formats a standard structure file to the desired input for the R package gghybrid.

# gghybrid.HI_230526.R
An R script that runs the R package `gghybrid`.

# popmap_telomere_samples_ONLY.tsv
An example of a popmap used in `stacks`.

# run_populations.sh
A Bash script for running the `stacks` module `populations` on a computing cluster

# run_vcftools_het.sh
A Bash scriipt for running `vcftools` to calculate heterozygosity on a computeing cluster
