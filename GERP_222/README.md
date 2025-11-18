## Subset the original VCF for GERP high & low scoring regions (conding/noncoding)
```
## GERP >= 1.0
sbatch job.subset_vcf.sh coding.derCor.gerp_1.0.bed coding.gerp_1.0.derCor_222.vcf.gz
sbatch job.subset_vcf.sh noncoding.derCor.gerp_1.0.bed noncoding.gerp_1.0.derCor_222.vcf.gz

## GERP < 1.0
sbatch job.subset_vcf.sh coding.derCor.gerp_less_1.0.bed coding.gerp_less_1.0.derCor_222.vcf.gz
sbatch job.subset_vcf.sh noncoding.derCor.gerp_less_1.0.bed noncoding.gerp_less_1.0.derCor_222.vcf.gz
```

## Get a random subset of 100k variants from the original VCF for GERP lower scoring regions: conding/noncoding (aka controls)
```
for i in coding.gerp_less noncoding.gerp noncoding.gerp_less; \
do \
    sbatch job.subset_vcf_random.sh ${i}_1.0.derCor_222.vcf.gz 100k_${i}_1.0.derCor_222.vcf.gz; \
done
```


## Make summary table
gerp_less_1 = LOW impact
gerp_1 = HIGH impact
```
VCF=coding.gerp_1.0.derCor_222.vcf.gz
VCF=100k_coding.gerp_less_1.0.derCor_222.vcf.gz
VCF=100k_noncoding.gerp_1.0.derCor_222.vcf.gz
VCF=100k_noncoding.gerp_less_1.0.derCor_222.vcf.gz

sbatch --array=1-222 jobs.bcftools_by_state.sh $VCF

OUT_SUMMARY=all_samples.gerp.derCor_222.by_state_var.tsv

cat logs/bcftools/bcftools*.log | sed 's/.derCor_222.vcf.gz//' | sed 's/100k_//' | sed 's/.gerp_less_1.0/\tLOW/' | sed 's/.gerp_1.0/\tHIGH'  > $OUT_SUMMARY

rm -r logs/
```

