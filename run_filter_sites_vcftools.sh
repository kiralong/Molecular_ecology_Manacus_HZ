#!/bin/bash
#SBATCH -p secondary
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -J vcftools_filter_sites

module load gcc

work=/projects/aces/kira2/telomeres_collab/telomere_only_samples/all_sites
invcf=$work/populations.all.vcf.gz
outvcf=$work/populations.all.R80.vcf.gz
cd $work

vcftools --gzvcf $invcf --stdout --max-missing 0.8 --recode | gzip > $outvcf
