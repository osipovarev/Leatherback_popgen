## 1. Annotate variants with snpEff
```
sbatch job.snpEff_combined.sh
```

### 2. Split combined annotated VCF by impact; then split by sample
```
sbatch job.by_impact_combined.sh

sbatch --array=1-222 jobs.bcftools_split_samples.sh HIGH
sbatch --array=1-222 jobs.bcftools_split_samples.sh MODERATE
sbatch --array=1-222 jobs.bcftools_split_samples.sh LOW
sbatch --array=1-222 jobs.bcftools_split_samples.sh MODIFIER
```



### 3. Count number of snps and indels of homozygotes and heterozygotes by impact
```
bash count_var_by_impact.sh
```


### 4. Make combined table with total number of homozygotes
```
paste <(grep -v -w LOF all_samples.hom_het_by_impact.tsv ) <(cut -f4 snpEff_combined/all_samples.hom_het_by_impact.tsv)  > all_samples.hom_het_by_impact_with_ref.tsv
```


