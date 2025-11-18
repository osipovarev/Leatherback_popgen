#!/bin/bash
#SBATCH -J WGR_snpEff
#SBATCH -o ./logs/snpEff/WGR_snpEff_%a_%A.log
#SBATCH -e ./logs/snpEff/WGR_snpEff_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 10:00:00
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=32000  # Requested Memory
#SBATCH --mail-type=END

module load snpeff/2017-11-24
module load uri/main
module load tabixpp/1.1.2-GCC-12.3.0

WDIR=/nese/meclab/Katya/snpEff_222/
CONFIG=$WDIR/snpEff.config
VCF=/nese/meclab/Katya/new_Leatherback_samples/snps.DerCor_combined_filtered.vcf.gz
GENINPUT=rDerCor1.pri.cur.20210524

snpEff eff -nodownload -c $CONFIG $GENINPUT $VCF > $WDIR/snps.DerCor_combined_filtered.snpEff.vcf
bgzip $WDIR/snps.DerCor_combined_filtered.snpEff.vcf
tabix $WDIR/snps.DerCor_combined_filtered.snpEff.vcf.gz
