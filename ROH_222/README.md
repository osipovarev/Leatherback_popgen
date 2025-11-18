# Running ROH analysis iwth PLINK

## 1. Make bed file for PLINK
```
sbatch job.plink_make_bed.sh /nese/meclab/Katya/new_Leatherback_samples/DerCor_combined_filtered.vcf.gz DerCor_combined_filtered
```
Total genotyping rate is 0.995047.
13590825 variants and 222 samples pass filters and QC.


## 2. Get ROHs with PLINK
q=20 - this paragmeter is coming from optimization
```
sbatch jobs.plink_roh.sh 
```


## 3. Count ROHs by length class 
```
bash get_binned_ROH.sh ./ > all_samples_binned_ROH.tsv
```
