#!/bin/bash
#SBATCH --job-name=plink_make_bed
#SBATCH --output=./logs/plink_make_bed.out 
#SBATCH --error=./logs/plink_make_bed.err 
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=24
#SBATCH --partition=cpu

module load uri/main
module load PLINK/2.00a3.7-gfbf-2023a  # loads plink v1.9

WDIR=/nese/meclab/Katya/PCA/
VCF=$1
PREFIX=$2

plink --vcf $VCF \
	--make-bed \
	--double-id \
	--allow-extra-chr \
	--chr-set 95 no-xy \
	--threads 24 \
	--out $PREFIX



