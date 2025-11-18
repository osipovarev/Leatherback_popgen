#!/bin/bash
#SBATCH -J bcftools
#SBATCH -o ./logs/bcftools/bcftools_%a_%A.log
#SBATCH -e ./logs/bcftools/bcftools_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 2:00:00
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=2000  # Requested Memory
#SBATCH --array=1-222

module load bcftools/1.19

VCF=$1

SAMPLES=/nese/meclab/Katya/new_Leatherback_samples/new.dc_222_samples.lst
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLES)

i=snps
hom=$(bcftools view -s $SAMPLE $VCF |  bcftools view -v ${i} -i 'GT!="0/0" && GT="hom"' | wc -l)
het=$(bcftools view -s $SAMPLE $VCF |  bcftools view -v ${i} -i 'GT="het"' | wc -l )
echo -e "$SAMPLE\t$i\t$hom\thom\t$VCF"; \
echo -e "$SAMPLE\t$i\t$het\thet\t$VCF"; \

