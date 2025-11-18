#!/bin/bash
#SBATCH -J subset
#SBATCH -o ./logs/subset/subset_%a_%A.log
#SBATCH -e ./logs/subset/subset_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 12:00:00
#SBATCH -c 12 # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=64000  # Requested Memory

module load bcftools/1.19

VCF=/nese/meclab/Katya/new_Leatherback_samples/snps.DerCor_combined_filtered.vcf.gz
BED=$1
OUT=$2

bcftools view -R $BED $VCF -Oz -o $OUT
