#!/bin/bash
#SBATCH -J bcftools
#SBATCH -o ./logs/bcftools/bcftools_%j_%A.log
#SBATCH -e ./logs/bcftools/bcftools_%j_%A.err
#SBATCH -p cpu
#SBATCH -t 2:00:00
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=4000  # Requested Memory
#SBATCH --mail-type=END

module load bcftools/1.19

WDIR=/nese/meclab/Katya/snpEff_222/

IMPACT=$1
VCF=$WDIR/$IMPACT.snps.DerCor_combined_filtered.snpEff.vcf.gz
SAMPLES=$2

bcftools view -S $SAMPLES $VCF | grep -v "AC=0" | bcftools view -O z -o split_by_group/$IMPACT.${SAMPLES%.samples.lst}.snpEff.vcf.gz

