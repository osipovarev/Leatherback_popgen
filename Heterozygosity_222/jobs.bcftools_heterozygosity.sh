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
#SBATCH --mail-type=END

module load bcftools/1.19

WDIR=/nese/meclab/Katya/Heterozygosity_222/

VCF=$1
SAMPLES=/nese/meclab/Katya/new_Leatherback_samples/new.dc_222_samples.lst

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLES)

hets=$(bcftools view -s $SAMPLE $VCF | bcftools view -i 'GT="het"' | grep -v "AC=0" | grep -v ^# | wc -l)
total=$(bcftools view -s $SAMPLE $VCF | grep -v "AC=0" | grep -v ^# | wc -l)
echo -e "$SAMPLE\t$hets\t$total"
