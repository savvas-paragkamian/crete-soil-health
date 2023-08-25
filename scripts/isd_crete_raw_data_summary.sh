#!/bin/bash

###############################################################################
# script name: isd_crete_raw_data_summary.sh
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to use basic command line tools to summarise the reads
# from fastq.gz files of ISD Crete 2016
###############################################################################
# usage:./isd_crete_raw_data_summary.sh

echo "This script runs some commands to summarise reads."

cd ../ena_data

echo "Calculate number of reads"
zcat *.fastq.gz | gawk 'BEGIN{RS="@" ; FS="\n"}{reads[$1]}END{print "total reads = " length(reads)}'

echo "Find the reads with primers, both forward and reverse"

zcat *1.fastq.gz  | gawk 'BEGIN{RS="@" ; FS="\n"} {match($2,/ACTCCTACGGGAGGCAGCAG/) ; a[RSTART]++}END{for (i in a){print i "\t" a[i]}}' > ../isd_crete_reads_1_primer_summary.tsv

zcat *2.fastq.gz  | gawk 'BEGIN{RS="@" ; FS="\n"}{match($2,/GGACTACHVGGGTWTCTAAT/) ; a[RSTART]++}END{for (i in a){print i "\t" a[i]}}' > ../isd_crete_reads_2_primer_summary.tsv

echo "Find how many reads have Ns and count them per sample"

zcat *.fastq.gz |  gawk 'BEGIN{RS="@" ; FS="\n"}($2 ~ /N/){split($1,sample,".") ; count[sample[1]][gsub(/N/,"")]++ }END{print "sample" "\t" "N_of_Ns" "\t" "total_reads"; for (s in count){for (c in count[s]){print s "\t" c "\t" count[s][c]}}}' > ../isd_crete_reads_Ns-summary.tsv


