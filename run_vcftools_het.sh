#!/bin/bash
#SBATCH -p aces
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 98:00:00
#SBATCH -J vcftools_telomeres_het

module load gcc

work=/projects/aces/kira2/telomeres_collab/telomere_only_samples/snps_only
invcf=$work/populations.snps.vcf.gz
cd $work

vcftools --gzvcf $invcf --het --out output_het
