#!/usr/bin/bash

OUTDIR=$1
OUTFILE=$OUTDIR/DerCor_combined_filtered_ROH*.hom
SAMPLES=/nese/meclab/Katya/new_Leatherback_samples/new.dc_222_samples.lst

for s in $(cat $SAMPLES); \
do \
	bin1=$(grep $s $OUTFILE | awk '($9<100){print $9}' | sum_stdin.py); \
	bin2=$(grep $s $OUTFILE | awk '($9>=100 && $9<500){print $9}' | sum_stdin.py); \
	bin3=$(grep $s $OUTFILE | awk '($9>=500){print $9}' | sum_stdin.py); \

	nbin1=$(grep $s $OUTFILE | awk '($9<100){print $9}' | wc -l ); \
	nbin2=$(grep $s $OUTFILE | awk '($9>=100 && $9<500){print $9}' | wc -l ); \
	nbin3=$(grep $s $OUTFILE | awk '($9>=500){print $9}' | wc -l ); \

	echo -e "$s\t$OUTDIR\t$nbin1\t$nbin2\t$nbin3\t$bin1\t$bin2\t$bin3"; \
done

