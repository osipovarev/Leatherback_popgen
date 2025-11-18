#!/bin/bash
#SBATCH -J bcftools
#SBATCH -o ./logs/bcftools/bcftools_%a_%A.log
#SBATCH -e ./logs/bcftools/bcftools_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 2:00:00
#SBATCH -c 8 
#SBATCH --threads-per-core=1
#SBATCH --mem=2000
#SBATCH --array=1-222
#SBATCH --mail-type=END

module load bcftools/1.19

WDIR=/nese/meclab/Katya/ANNOVAR_222/

IMPACT=$1
VCF=$WDIR/$IMPACT.derCor_annovar.vcf.gz
SAMPLES=/nese/meclab/Katya/new_Leatherback_samples/new.dc_222_samples.lst

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLES)

bcftools view -s $SAMPLE $VCF | bcftools view -i 'GT="het"' -O z -o $WDIR/split_by_sample/$IMPACT.het.$SAMPLE.vcf.gz
bcftools view -s $SAMPLE $VCF | bcftools view -i 'GT="hom"' -O z -o $WDIR/split_by_sample/$IMPACT.hom.$SAMPLE.vcf.gz

