#!/bin/bash

chr=$1
input_file=chromosome_${chr}_summary_stats.txt
module load plink/1.90

plink --bfile g1000_eur --a1-allele $input_file 4 2 "#" --r yes-really square gz --out LD_chr${chr}_BD 
