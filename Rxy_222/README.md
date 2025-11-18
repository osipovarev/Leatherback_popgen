# Rxy analysis between larger_declining and small_stable populations


## 1. Prepare sample sets by demographic trend group
```
WDIR=/nese/meclab/Katya/new_Leatherback_samples
cut -f2,15,16,17,30,34 dc_222_samples_info.tsv > simplified_dc_222_samples_info.tsv
INFO=simplified_dc_222_samples_info.tsv

for i in $(cut -f5 $INFO | tail -n +2 | s -u); do g -v exclude $INFO | g $i | cut -f1 | awk '{if($0 ~ /^dc_/) print $0 "_bb_noq"; else print $0}' > $i.samples.lst; done

INFO=dc_222_samples_info.tsv
for i in $(cut -f32 $INFO | tail -n +2 | s -u); do g -v exclude $INFO | g -w $i | cut -f2 | awk '{print $0 "_bb_noq.g"}' > $i.samples.lst; done
```


## 2. Split annotated snp VCF into sample groups
```
WDIR=/nese/meclab/Katya/Rxy_222/

mv ../new_Leatherback_samples/*.samples.lst .

md split_by_group/

for IMPACT in HIGH MODERATE LOW MODIFIER; \
do \
    for SAMPLES in larger_declining small_stable; \
    do  \
        sbatch job.bcftools_split_samples.sh $IMPACT $SAMPLES.samples.lst; \
    done; \
done
```


## 3. Get genotypes; count derived alleles; split into blocks for jackknife
```
md split_by_block/

for IMPACT in HIGH MODERATE LOW MODIFIER; \
do \
    for SAMPLES in larger_declining small_stable; \
    do  \

        ## get genotypes
        bcftools query -f '%CHROM\t%POS[\t%GT]\n' split_by_group/$IMPACT.$SAMPLES.snpEff.vcf.gz > split_by_group/$IMPACT.$SAMPLES.gts; \

        ## count derived alleles
        bash count_derived_alleles.awk split_by_group/$IMPACT.$SAMPLES.gts > split_by_group/$IMPACT.$SAMPLES.counts; \

        ## split into 100 blocks
        total=$(wc -l < split_by_group/$IMPACT.$SAMPLES.counts); \
        blocksize=$((total / 100)); \
        split -l $blocksize split_by_group/$IMPACT.$SAMPLES.counts split_by_block/$IMPACT.$SAMPLES.block_ ; \

    done; \
done
```


## 4. Run block jackknife approach
block jackknife = divided the SNP data sets into 100 contiguous blocks and then recomputed the statistic on all of the data except for the data from that block.
```
for IMPACT in HIGH MODERATE LOW MODIFIER; \
do \
    for SAMPLES in larger_declining small_stable; \
    do  \
        for b in split_by_block/$IMPACT.$SAMPLES.block_*; \
        do  \
            grep -v -F -f $b split_by_group/$IMPACT.$SAMPLES.counts > tmp_counts; \
            cut -f3 tmp_counts | sum_stdin.py >> $IMPACT.$SAMPLES.jackknife; \
        done; \
    done; \
done    
```


## 5. Compute Rxy between larger_declining and small_stable populations; normalize
```
S1=small_stable
S2=larger_declining

s1=$(wc -l < $S1.samples.lst)
s2=$(wc -l < $S2.samples.lst)

for IMPACT in HIGH MODERATE LOW MODIFIER; \
do \
    paste $IMPACT.$S1.jackknife $IMPACT.$S2.jackknife | awk -v s1=$s1 -v s2=$s2 '{print $1/$2/s1*s2}' > $IMPACT.$S1.$S2.rxy; \
done

## normalize by Rxy in neutral regions
for IMPACT in HIGH MODERATE LOW; do paste $IMPACT.$S1.$S2.rxy MODIFIER.$S1.$S2.rxy | awk '{print $1/$2}' > norm.$IMPACT.$S1.$S2.rxy; done

## make final table: Rxy comparison between populations
for IMPACT in HIGH MODERATE LOW; do cat norm.$IMPACT.$S1.$S2.rxy | awk -v i=$IMPACT '{print $0"\t"i}'; done > norm.$S1.$S2.tsv
```

## 6. Make summary table for all population comparisons
```
for f in $(ls norm*tsv); do c=$(echo -e "$f" | sed 's/norm.//' | sed 's/.tsv//'); awk -v c=$c '{print $0"\t"c}' $f; done > summary_rxy.all_populaitons.tsv
```








