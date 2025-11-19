#!/bin/bash
#SBATCH --job-name=plink_prune_ld
#SBATCH --output=./logs/plink_prune_ld.out 
#SBATCH --error=./logs/plink_prune_ld.err 
#SBATCH --time=2:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=12
#SBATCH --partition=cpu

module load uri/main
module load PLINK/2.00a3.7-gfbf-2023a  # loads plink v1.9

PREFIX=$1

# LD pruning with r2 > 0.2
plink --bfile $PREFIX \
      --indep-pairwise 50 10 0.2 \
	--chr-set 95 no-xy \
	--allow-extra-chr \
	--make-bed \
	--threads 12 \
	--out ld_pruned_0.2.$PREFIX
