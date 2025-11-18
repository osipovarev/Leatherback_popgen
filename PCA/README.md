# PCA analysis of new derCor samples


## 1. Convert VCF to PLINK 
```
INDIR=/nese/meclab/Katya/new_Leatherback_samples/

VCF=$INDIR/derCor_42_snpArcher_raw.vcf.gz
sbatch job.plink_make_bed.sh $VCF derCor_42

VCF=$INDIR/derCor_28_snpArcher_raw.vcf.gz
batch job.plink_make_bed.sh $VCF derCor_28

VCF=/nese/meclab/Lisa/Dc_WGR_analyses_docsMay2024/Final_filtered_combvcf/DerCor_combined_filtered.vcf
sbatch job.plink_make_bed.sh $VCF derCor_153

VCF=$INDIR/DerCor_combined_filtered.vcf.gz
sbatch job.plink_make_bed.sh $VCF derCor_222
```

derCor_42:  15837364 variants and 41 samples pass filters and QC.
derCor_28:  4811884 variants and 28 samples pass filters and QC.
derCor_153: 10787173 variants and 153 samples pass filters and QC.
derCor_222: 17691434 variants and 222 samples pass filters and QC.
derCor_222: 13590825

## 2. Linkage Disequilibrium (LD) Pruning
```
sbatch job.plink_prune_ld.sh derCor_42
sbatch job.plink_prune_ld.sh derCor_28
sbatch job.plink_prune_ld.sh derCor_153
sbatch job.plink_prune_ld.sh derCor_222
```
### r=0.2:
ld_pruned_0.2.derCor_42.log:  11369768 of 15837364 variants removed.
ld_pruned_0.2.derCor_28.log:   4184835 of  4811884 variants removed.
ld_pruned_0.2.derCor_153.log:  6426681 of 10787173 variants removed.
ld_pruned_0.2.derCor_222.log:  9630303 of 17691434 variants removed.
ld_pruned_0.2.derCor_222.log:  7245510 of 13590825 variants removed

### r=0.1:
ld_pruned_0.1.derCor_153.log:  7351951 of 10787173 variants removed.
ld_pruned_0.1.derCor_222.log: 11555824 of 17691434 variants removed.


## 3. PCA
```
sbatch job.plink_pca.sh derCor_42
sbatch job.plink_pca.sh derCor_28
sbatch job.plink_pca.sh derCor_153
sbatch job.plink_pca.sh derCor_222
```

