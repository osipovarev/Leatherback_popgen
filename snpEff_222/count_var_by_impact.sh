WDIR=/nese/meclab/Katya/snpEff_222
SAMPLES=/nese/meclab/Katya/new_Leatherback_samples/new.dc_222_samples.lst
OUTFILE=all_samples.snp_indel_hom_het_by_impact.tsv

module load bcftools/1.19
i=snps

for SAMPLE in $(cat $SAMPLES); 
do 
    for STATE in hom het; 
    do 
        for IMPACT in HIGH MODERATE LOW MODIFIER; 
        do 
		if [ $STATE == "hom" ]; 
		then 
			NUMBER=$(bcftools view -v ${i} -i 'GT=="0/0"' $WDIR/split_by_sample/$IMPACT.$STATE.$SAMPLE.snpEff.vcf.gz | grep -v ^# |  wc -l);
			echo -e "$SAMPLE\t$IMPACT\t${STATE}_ref\t$NUMBER\t$i";
			NUMBER=$(bcftools view -v ${i} -i 'GT!="0/0"' $WDIR/split_by_sample/$IMPACT.$STATE.$SAMPLE.snpEff.vcf.gz | grep -v ^# |  wc -l);
		else
                	NUMBER=$(bcftools view -v ${i} $WDIR/split_by_sample/$IMPACT.$STATE.$SAMPLE.snpEff.vcf.gz | grep -v ^# |  wc -l); 
		fi;
                echo -e "$SAMPLE\t$IMPACT\t$STATE\t$NUMBER\t$i"; \
        done; 
    done; 
done > $OUTFILE
