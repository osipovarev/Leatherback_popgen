# ANNOVAR: running gene-based variant annotation of Leatherback samples


### 1. Prepare custom database from derCor gene annotation
```
WDIR=/nese/meclab/Katya/Scripts/annovar/
md turtledb/
cd turtledb/

# link to gene annotation
ln -s /nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105/Annotation/GCF_009764565.3_rDerCor1.pri.v4_genomic_rename.gtf genes.gtf

# link to genome
ln -s /nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105/rDerCor1.pri.cur.20210524.fasta genome.fa

gtfToGenePred -genePredExt -ignoreGroupsWithoutExons genes.gtf stdout | awk '{print "585\t"$0}' > derCor_genes.txt

retrieve_seq_from_fasta.pl --format refGene --seqfile genome.fa derCor_genes.txt --out derCor_genesMrna.fa

```


### 2. Run ANNOVAR annotation
```
sbatch job.annovar_combined.sh 

mv derCor_annovar.derCor_multianno.vcf derCor_annovar.vcf
VCF=derCor_annovar.vcf

module load uri/main
module load tabixpp/1.1.2-GCC-12.3.0

bgzip $VCF
tabix $VCF.gz
```


### 3. Assign impact to classes of variants; split VCF by impact
```
frameshift_insertion    HIGH
frameshift_deletion HIGH
frameshift_block_substitution   HIGH
stopgain    HIGH
stoploss    HIGH
nonframeshift_insertion MODERATE
nonframeshift_deletion  MODERATE
nonframeshift_block_substitution    MODERATE
nonsynonymous_SNV   MODERATE
synonymous_SNV  LOW
unknown MODIFIER

## split by impact
sbatch job.by_impact_annovar.sh

## fix header
module load uri/main
module load bcftools/1.19
module load tabixpp/1.1.2-GCC-12.3.0

VCF=derCor_annovar.vcf.gz

echo '##INFO=<ID=.,Number=0,Type=Flag,Description="Placeholder for undefined INFO">' > header.txt


for impact in HIGH MODERATE LOW MODIFIER; \
do \
    bcftools annotate -h header.txt $impact.$VCF -O z -o fixed.vcf.gz; \
    mv fixed.vcf.gz $impact.$VCF; \
done
```


### 4. Split by sample and by state (hom/het)
```
for impact in HIGH MODERATE LOW MODIFIER; \
do \
    sbatch --array=1-222 jobs.bcftools_split_samples.sh $impact; \
done
```


### 5. Count number of snps and indels of homozygotes and heterozygotes by impact
```
bash count_var_by_impact.sh
```


