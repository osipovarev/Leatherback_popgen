#!/bin/bash
#SBATCH -J randomLines
#SBATCH -o ./logs/randomLines/randomLines_%a_%A.log
#SBATCH -e ./logs/randomLines/randomLines_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 2:00:00
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=64000  # Requested Memory

module load bcftools/1.19

VCF=$1
OUT=$2

bcftools sort <(cat <(zgrep "^#" $VCF) <(zgrep -v "^#" $VCF | randomLines stdin 100000 stdout)) -Oz -o $OUT

