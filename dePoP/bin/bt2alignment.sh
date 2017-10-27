#!/bin/sh
bam=$1
sam=$bam'.sam'
report=$bam'_report.txt'
uns_bam=$bam'_unsorted.bam'

if [ -f $bam ]; then rm -f $bam; fi
if [ -f $bam.bai ]; then rm -f $bam.bai; fi
if [ -f $report ]; then rm -f $report; fi

bowtie2-align-s -x $2 -U $3 --rg-id SM --rg 'SM:'$4 --threads $5 -S $sam 2>> $report 

samtools view -bS $sam > $uns_bam
samtools 'sort' $uns_bam -o $bam
samtools index $bam

rm $sam $uns_bam
#chmod -w $report $bam $bam.bai
