#!/bin/bash

working_dir=/projects/aces/kira2/telomeres_collab/populations_runs/populations_telomeres_2910_p3_r0.8_maf0.03_hwe_structure_36k_wl
structure_file=$working_dir/populations.structure
output_file=$working_dir/populations.p3.r80.maf03.36k.structure.csv

cat $structure_file | grep -v "^#" | tr '\t' ',' | sed -E 's/,0/,NA/g' |  sed -E 's/,,/INDLABEL,POPID,/' | \
	sed 's/,pr_090/,090PR/' | \
	sed 's/,cg_100/,100CG/' | \
	sed 's/,ss_020/,020SS/' > $output_file
