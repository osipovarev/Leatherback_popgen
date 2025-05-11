## 3. Run GERP
```
mamba install bioconda::gerp

sbatch --array=1-40 jobs.filter_maf.sh

sbatch --array=1-40 jobs.gerp.sh

```


### 3.1. Get GERP summary stats in histograms
```
# distribution of scores
cat mafs/*rates		  | cut -f2 | python3 hist_stdin.py 20 gerp.rates.hist.pdf

# distribution of sum scores per conserved element
cat mafs/*maf.rates.elems | cut -f4 | python3 hist_stdin.py 100 elem.rates.hist.pdf

# distribution of lenghts of conserved elements
cat mafs/*maf.rates.elems | cut -f6 | python3 hist_stdin.py 50 elem.length.hist.pdf
```


### 3.2. Separate GERP results into coding and non-coding regions
```
sbatch --array=1-40 jobs.bedtools.sh

cat mafs/coding*bed | cut -f4 | python3 hist_stdin.py 50 coding.gerp.rates.hist.pdf

cat mafs/noncoding*bed | cut -f4 | python3 hist_stdin.py 50 noncoding.gerp.rates.hist.pdf
```


### 3.3. Get high score GERP. Overlap with annotation
```
ANNODIR=/n/holylfs05/LABS/informatics/Lab/project-eosipova/Leatherback/Annotation/
ANNO=$ANNODIR/renamed.derCor.ncbi.bed
EXONS=$ANNODIR/exons.renamed.derCor.ncbi.bed

## Split annotation into coding exons
bed12ToBed6 -i $ANNO > $EXONS


## Extract positions with GERP score = 1.0
### coding
bedtools merge -i <(cat coding.*.maf.rates.bed | awk '$4>=1{print }' | sort -k1,1 -k2,2n) > ../coding.derCor.gerp_1.0.bed

### noncoding
bedtools merge -i <(cat noncoding.*.maf.rates.bed | awk '$4>=1{print }' | sort -k1,1 -k2,2n) > ../noncoding.derCor.gerp_1.0.bed


## Overlap coding GERP>=1 and CDS annotation
cat $EXONS | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
279574900

cat coding.derCor.gerp_1.0.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
27468972

bedtools intersect -a coding.derCor.gerp_1.0.bed -b $EXONS | uniq -u | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
10053840

## Total noncoding bases
cat noncoding.derCor.gerp_1.0.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

```
