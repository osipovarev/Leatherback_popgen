# Admixture analysis

## 1. Copy ld-pruned bed file from PCA_222 analysis
```
BED=ld_pruned_0.2.derCor_222
cp ../PCA_222/$BED.bed ./
cp ../PCA_222/$BED.bim ./
cp ../PCA_222/$BED.fam ./

## fix non-integer chr names
sed -i 's/SUPER_//' $BED.bim

## remove missing genotypes
plink --bfile $BED --geno 0.99 --maf 0.001 --chr-set 95 no-xy --make-bed --out filt.$BED
```


## 2. Run Admixture all 222 samples
```
sbatch --array=2-8 jobs.admixture.sh $BED.bed
```


## 3. Subset Western / Eastern Pacific samples
```
WDIR=/nese/meclab/Katya/new_Leatherback_samples/
INFO=$WDIR/dc_222_samples_info.tsv

grep 'Western Pacific' $INFO | cut -f2 | awk '{print $0"_bb_noq.g\t"$0"_bb_noq.g"}' > WP.samples.txt
grep 'Eastern Pacific' $INFO | cut -f2 | awk '{print $0"_bb_noq.g\t"$0"_bb_noq.g"}' > EP.samples.txt

module load uri/main
module load PLINK/2.00a3.7-gfbf-2023a # loads plink v1.9


BED=ld_pruned_0.2.derCor_222
for i in WP EP; \
do \
    plink --bfile $BED --keep $i.samples.txt -make-bed --double-id --allow-extra-chr --chr-set 95 no-xy --threads 24 --out $i.$BED; \
    plink --bfile $i.$BED --geno 0.99 --maf 0.001 --chr-set 95 no-xy --make-bed --out filt.$i.$BED \
done
```


## 4. Run Admixture on Western /Eastern Pacific samples
```
BED=filt.WP.ld_pruned_0.2.derCor_222.bed
sbatch --array=2-8 jobs.admixture.sh $BED

BED=filt.EP.ld_pruned_0.2.derCor_222.bed
sbatch --array=2-8 jobs.admixture.sh $BED
```

## CV errors
### WP
error (K=2): 0.30173
error (K=3): 0.30886
error (K=4): 0.32499
error (K=5): 0.33848
error (K=6): 0.35652
error (K=7): 0.36792
error (K=8): 0.38532

### EP
(K=3): 0.43859
(K=4): 0.50333
(K=5): 0.52752
(K=6): 0.57116
(K=7): 0.58702
(K=8): 0.65475

### all
(K=2): 0.19556
(K=3): 0.13749
(K=4): 0.18372
(K=5): 0.14005
(K=6): 0.14260
(K=7): 0.14460
(K=8): 0.14634

### Prepare files to plot Admixture results
```
for g in EP WP; \
do \
    for i in {2..8}; \
    do \
        paste <(cut -f1 filt.$g.ld_pruned_0.2.derCor_222.nosex | sed 's/_bb_noq.g//') filt.$g.ld_pruned_0.2.derCor_222.$i.Q | sed 's/ /\t/g' > named.filt.$g.ld_pruned_0.2.derCor_222.$i.Q; \
    done; \
done
```
