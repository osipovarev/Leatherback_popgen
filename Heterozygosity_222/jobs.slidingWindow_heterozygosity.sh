#!/bin/bash
#SBATCH -J slidingWindow
#SBATCH -o ./logs/slidingWindow/slidingWindow_%a_%A.log
#SBATCH -e ./logs/slidingWindow/slidingWindow_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 1:00:00
#SBATCH -c 4  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=1000  # Requested Memory
#SBATCH --array=1-222
#SBATCH --mail-type=END

module load conda/latest
conda activate /work/pi_lkomoroske_umass_edu/.conda/envs/hetwindow

WDIR=/nese/meclab/Katya/new_Leatherback_samples/
CHROMS=$WDIR/chrom.sizes

VCF=$1
NAME=$(echo -e $VCF | awk -F"/" '{print $NF}' | sed 's/.vcf.gz//')

SAMPLES=/nese/meclab/Katya/new_Leatherback_samples/dc_222_samples.lst
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLES)

S_VCF=$WDIR/split_by_sample/$SAMPLE.$NAME.vcf.gz
SlidingWindowHet.py $S_VCF 1000000 100000 $CHROMS

