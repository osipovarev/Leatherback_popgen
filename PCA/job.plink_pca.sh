#!/bin/bash
#SBATCH --job-name=plink_pca
#SBATCH --output=./logs/plink_pca.out 
#SBATCH --error=./logs/plink_pca.err 
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=24
#SBATCH --partition=cpu

module load uri/main
module load PLINK/2.00a3.7-gfbf-2023a  # loads plink v1.9

PREFIX=ld_pruned_0.2.$1

plink --bfile $PREFIX \
	--allow-extra-chr \
	--chr-set 95 no-xy \
	--threads 24 \
	--pca \
	--out pca.$PREFIX
