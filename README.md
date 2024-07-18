# Molecular Ecology of a _Manacus_ hybrid zone Manuscript Repository

>Ben J. Vernasco*, Kira M. Long*, Michael J. Braun, Jeffrey D. Brawn. (2024) **Genetic and telomeric variability: Insights from a tropical avian hybrid zone**. In Press. _Molecular Ecology_. (*Co-first Authors)

Code repository describing bioinformatics for the assesment of associations between telomere variation, hybridization and genetic variation in _Manacus_ birds.

Authors: Ben J. Vernasco, Kira M. Long, Michael J. Braun, and Jeffrey D. Brawn

Primary Contacts: benvernasco@gmail.com and kiralong778@gmail.com

# Description of each script

## calculate_idnv-obs_het.py
A custom python script that calculates both the SNP only and all sites heterozygosity ratio for each individual. Requires an input vcf with all sites, such as a `populations.all.vcf` from the `populations` module of *Stacks* [(Rochette et al. 2019)](https://catchenlab.life.illinois.edu/stacks/) using the `--vcf-all` flag. This script calculates the observed individual heterozygosity, i.e. the proportion of heterozygous sites across all genotyped sites.

`calculate_idnv-obs_het.py` example usage:
```
calculate_idnv-obs_het.py --vcf /path/to/vcf-all/populations.all.vcf --outdir /path/to/output/directory
```

## filter_sumstats_to_whitelist.py
A custom  python script that allows the user to filter snps for generating a whitelist of snps to use in the *Stacks* [(Rochette et al. 2019)](https://catchenlab.life.illinois.edu/stacks/) software.

## format_structurefile_for_gghybrid.sh
A Bash script that formats a standard structure file to the desired input for the R package gghybrid [(Bailey 2024)](https://doi.org/10.1111/1755-0998.13910).

## gghybrid.HI_230526.R
An R script that runs the R package [gghybrid](https://github.com/ribailey/gghybrid?tab=readme-ov-file).

## popmap_telomere_samples_ONLY.tsv
An example of a popmap used in *Stacks*.

## run_populations.sh
A Bash script for running the *Stacks* module `populations` on a computing cluster using SLURM.

## run_vcftools_het.sh
A Bash script for running `vcftools` to calculate heterozygosity on a computing cluster using SLURM.

## Data availability
Raw RAD-seq reads can be found on NCBI under BioProject [PRJNA893627](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA893627).
All other data used in this publication are on [dryad](https://doi.org/10.5061/dryad.qnk98sfnh).
