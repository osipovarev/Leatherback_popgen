#!/bin/bash
#SBATCH -J WGR_ANNOVAR
#SBATCH -o ./logs/ANNOVAR/WGR_ANNOVAR_%a_%A.log
#SBATCH -e ./logs/ANNOVAR/WGR_ANNOVAR_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 24:00:00
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=96000  # Requested Memory
#SBATCH --mail-type=END


WDIR=/nese/meclab/Katya/ANNOVAR_222/
DB=/nese/meclab/Katya/Scripts/annovar/turtledb/
VCF=/nese/meclab/Katya/new_Leatherback_samples/snps.DerCor_combined_filtered.vcf.gz
BUILD=derCor
PROTOCOL=genes

table_annovar.pl $VCF $DB -out derCor_annovar -buildver $BUILD -protocol $PROTOCOL -operation g -remove -polish -vcfinput -nastring .

