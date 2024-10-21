
# 1. Genetic load analysis

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
	modif=$(cut -f8 $f | tail -n +3 | sum_stdin.py); \
	 echo -e "$i\t$high\t$low\t$mod\t$modif"; \
done > $WDIR/all_samples.total_variant_by_effect.tsv

```

### Count number of homozygotes and heterozygotes by impact
```
WDIR=/nese/meclab/Katya/snpEff_turtles/snpEff_vcfs_by_impact/
SCRATCH=/scratch/workspace/lkomoroske_umass_edu-Dc_WGR/
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


# 2. Inbreeding and genetic load

### Calculation of F_ROH for each individual of leatherback turtle

genome size = 2,164,762,090 bp

```
WDIR=/nese/meclab/Katya/ROH_turtles/
SCRATCH=/scratch/workspace/lkomoroske_umass_edu-Dc_WGR/
SAMPLES=$SCRATCH/Dc_SNPs_May24_samplelist.txt

sbatch --array=1-153 jobs.bcftools_roh.sh


for SAMPLE in $(cat $SAMPLES); do L_ROH=$(grep ^RG ROH_per_individual/$SAMPLE.roh.txt | awk '{sum += $6} END {print sum}'); echo -e "$SAMPLE\t$L_ROH"; done > all_samples.L_ROH.tsv

```

