#!/bin/bash
#SBATCH -J bcftools
#SBATCH -o ./logs/bcftools/bcftools_%a_%A.log
#SBATCH -e ./logs/bcftools/bcftools_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 2:00:00
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=2000  # Requested Memory
#SBATCH --array=1-153
#SBATCH --mail-type=END

module load bcftools/1.15

WDIR=/nese/meclab/Katya/ROH_turtles/
SCRATCH=/scratch/workspace/lkomoroske_umass_edu-Dc_WGR/
SAMPLES=$SCRATCH/Dc_SNPs_May24_samplelist.txt

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLES)

bcftools roh --AF-dflt 0.4 $SCRATCH/snpEff/ind_vcf_forsnpEff/$SAMPLE.vcf > $WDIR/ROH_per_individual/$SAMPLE.roh.txt

