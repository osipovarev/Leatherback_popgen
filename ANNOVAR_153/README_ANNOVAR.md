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
WDIR=/nese/meclab/Katya/ANNOVAR_turtles/
DB=/nese/meclab/Katya/Scripts/annovar/turtledb/
ANNO=/nese/meclab/Katya/Scripts/annovar/turtledb/derCor_genes.txt
VCF=/nese/meclab/Lisa/Dc_WGR_analyses_docsMay2024/Final_filtered_combvcf/DerCor_combined_filtered.vcf

sbatch job.annovar_conbined.sh
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

VCF=derCor_annovar.derCor_multianno.vcf
bgzip $VCF
tabix $VCF.gz

sbatch job.by_impact_annovar.sh

## fix header
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
sbatch --array=1-153 jobs.bcftools_split_samples.sh HIGH
sbatch --array=1-153 jobs.bcftools_split_samples.sh MODERATE
sbatch --array=1-153 jobs.bcftools_split_samples.sh LOW
sbatch --array=1-153 jobs.bcftools_split_samples.sh MODIFIER
```

### 5. Count number of snps and indels of homozygotes and heterozygotes by impact
```
bash count_var_by_impact.sh
```


