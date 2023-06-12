#!/bin/bash
#SBATCH -p aces
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -J populations_telomeres_mac3_hwe_structure_wl
#SBATCH -t 68:00:00

module load gcc/7.2.0
thr=${SLURM_NTASKS}

stacks=/projects/catchenlab/local/bin
popmap_path=/projects/aces/kira2/telomeres_collab/popmap/popmap_telomere_samples_ONLY.tsv
gstacks_out=/projects/aces/kira2/chapter_1_HZ_movement/gstacks_runs/220308_gstacks_rm-pcr-dups_ALL
populations_output_path=/projects/aces/kira2/telomeres_collab/populations_runs
whitelist_path=/projects/aces/kira2/telomeres_collab/whitelists/telomeres_p3_r0.8_maf0.03_36261_whitelist.tsv

mac=3
r=0.8
p=3

populations_output=$populations_output_path/populations_telomeres_2910_p${p}_r${r}_maf${maf}_hwe_structure_36k_wl
mkdir -p $populations_output

cmd=(
    $stacks/populations
    --in-path $gstacks_out
    --out-path $populations_output
    --popmap $popmap_path
    --threads $thr
    --min-samples-per-pop $r
    --min-mac $mac
    --min-population $p
    --hwe
    --whitelist $whitelist_path
    --structure
#    --genepop
#    --plink
#    --vcf
#    --ordered-export
#    --vcf-all
)

"${cmd[@]}"
