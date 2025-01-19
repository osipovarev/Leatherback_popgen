### Count total number of variants by impact
```
scratch
cd snpEff/snpEff_out
WDIR=/nese/meclab/Katya/snpEff_turtles/

for f in $(ls *.g.genes.txt); \
do i=${f%.g.genes.txt}; \
	high=$(cut -f5 $f | tail -n +3 | sum_stdin.py); \
	low=$(cut -f6 $f | tail -n +3 | sum_stdin.py); \
	mod=$(cut -f7 $f | tail -n +3 | sum_stdin.py); \
	modif=$(cut -f5,6,7,8 $f | awk '$0=$4-$1-$2-$3{print}' | tail -n +3 | sum_stdin.py); \
	 echo -e "$i\t$high\t$low\t$mod\t$modif"; \
done > $WDIR/all_samples.total_variant_by_impact.tsv

```

### Count number of homozygotes and heterozygotes by impact
```
WDIR=/nese/meclab/Katya/snpEff_turtles/snpEff_vcfs_by_impact/
SCRATCH=/scratch/workspace/lkomoroske_umass_edu-Leatherback_WGRv2/
SAMPLES=$SCRATCH/Dc_SNPs_May24_samplelist.txt

### wrote jobs.bcftools_by_impact.sh
sbatch --array=1-153 jobs.bcftools_by_impact.sh


for SAMPLE in $(cat $SAMPLES); \
do \
	for STATE in hom het; \
	do \
		for IMPACT in HIGH MODERATE LOW; \
		do \
			NUMBER=$(grep -v ^# $WDIR/$IMPACT.$STATE.$SAMPLE.snpEff.vcf | wc -l); \
			echo -e "$SAMPLE\t$IMPACT\t$STATE\t$NUMBER"; \
		done; \
	done; \
done > all_samples.hom_het_by_impact.tsv

### add LOF (subset of HIGH) variants separately
IMPACT=HIGH

for SAMPLE in $(cat $SAMPLES); \
do \
	for STATE in hom het; \
	do \
		NUMBER=$(grep -v ^# $WDIR/$IMPACT.$STATE.$SAMPLE.snpEff.vcf | grep -w LOF | wc -l); \
	echo -e "$SAMPLE\tLOF\t$STATE\t$NUMBER"; \
	done; \
done >> all_samples.hom_het_by_impact.tsv

```

### Count number of snps and indels of homozygotes and heterozygotes by impact
```
WDIR=/nese/meclab/Katya/snpEff_turtles/snpEff_combined/snpEff_out/
module load bcftools/1.19

bash count_var_by_impact.sh
```




### Split combined annotated VCF by impact; then split by sample
```
sbatch job.by_impact_conbined.sh

sbatch --array=1-153 jobs.bcftools_split_samples.sh HIGH
sbatch --array=1-153 jobs.bcftools_split_samples.sh MODERATE
sbatch --array=1-153 jobs.bcftools_split_samples.sh LOW
sbatch --array=1-153 jobs.bcftools_split_samples.sh MODIFIER
```



### Count number of snps and indels of homozygotes and heterozygotes by impact
```
WDIR=/nese/meclab/Katya/snpEff_turtles/snpEff_combined/snpEff_out/
module load bcftools/1.19

bash count_var_by_impact.sh
```


### Make combined table with total number of homozygotes
```
paste <(grep -v -w LOF all_samples.hom_het_by_impact.tsv ) <(cut -f4 snpEff_combined/all_samples.hom_het_by_impact.tsv)  > all_samples.hom_het_by_impact_with_ref.tsv
```


