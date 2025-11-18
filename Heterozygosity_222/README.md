## Calculate genome-wide heterozygosity

### for all variant sites: 
```
VCF=/nese/meclab/Katya/new_Leatherback_samples/DerCor_combined_filtered.vcf.gz

sbatch --array=1-222 jobs.bcftools_heterozygosity.sh $VCF

cat logs/bcftools/bcftools_*.log > all_samples.hets.variants.tsv
rm -r logs/
```


### for snps only: 
```
VCF=/nese/meclab/Katya/new_Leatherback_samples/snps.DerCor_combined_filtered.vcf.gz

sbatch --array=1-222 jobs.bcftools_heterozygosity.sh $VCF

cat logs/bcftools/bcftools_*.log > all_samples.hets.snps.tsv
rm -r logs/
```

### combine into one table
```
cat <(awk '{print $0"\tvariants"}' all_samples.hets.variants.tsv) <(awk '{print $0"\tsnps"}' all_samples.hets.snps.tsv) > all_samples.hets.all.tsv
```



## Calculate window-based heterozigosity

### Split by sample
```
VCF=../new_Leatherback_samples/snps.DerCor_combined_filtered.vcf.gz
sbatch --array=1-222 jobs.bcftools_split_samples.sh $VCF

VCF=../new_Leatherback_samples/DerCor_combined_filtered.vcf.gz
sbatch --array=1-222 jobs.bcftools_split_samples.sh $VCF
``` 


### Calculate heterozygosity
```
# test
conda activate /work/pi_lkomoroske_umass_edu/.conda/envs/hetwindow
VCF=/nese/meclab/Katya/new_Leatherback_samples/snps.DerCor_combined_filtered.vcf.gz

sbatch --array=1-222 jobs.slidingWindow_heterozygosity.sh $VCF
```

### Combine into one table
```
for S in $(cat $SAMPLES); do f=$(echo $WDIR/split_by_sample/$S.snps.DerCor_combined_filtered.vcf.gz_het_1000000win_100000step.txt); cat $f | awk -v S=$S '{print $0"\t"S}'; done > all_sliding_window_hets.tsv
```

