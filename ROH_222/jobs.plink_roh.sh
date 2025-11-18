#!/bin/bash
#SBATCH -J ROH_PLINK_filt_DerCor
#SBATCH -o ./logs/ROH_PLINK_filt_DerCor.log
#SBATCH -e ./logs/ROH_PLINK_filt_DerCor.err
#SBATCH -p cpu
#SBATCH -t 04:00:00
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=96000  # Requested Memory
#SBATCH --mail-type=END


module load uri/main
module load PLINK/2.00a3.7-gfbf-2023a


q=20
plink --homozyg --bfile DerCor_combined_filtered \ 
	--allow-extra-chr \
	--chr-set 95 no-xy \
	--threads 32 \
	--out DerCor_combined_filtered_ROH_${q} \
	--homozyg-snp 20 \
	--homozyg-kb 50 \
	--homozyg-window-snp $q \
	--homozyg-window-het 1 \
	--homozyg-window-missing 5 \
	--homozyg-window-threshold 0.01; \


