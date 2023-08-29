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
# usage:./isd_crete_reads_summary.sh -i ~/isd_crete/ena_data/ -o ~/isd_crete/output/

#################################### USAGE ####################################
## Usage of the script
usage="Use the parameter -i for the path of fastq.gz \\
and -o for the path of the output files.\n \\
Example: ./isd_crete_raw_data_summary.sh -i ~/isd_crete/ena_data/ -o ~/isd_crete/output/ \n"

## User input PDF file
while getopts "i:o:" option
do
   case "$option" in
      i)   files="${OPTARG}";;
      o)   output="${OPTARG}";;
      ?|:)   echo -e "$usage" ; exit 1;;
      *)   echo -e "option ${OPTARG} unknown. \n$usage" ; exit 1 ;;
   esac
done

## Detect if no options were passed
if [ $OPTIND -eq 1 ];
    then echo -e "No options were passed. \n$usage "; exit 1;
fi

## Detect if a single option is missing
if [ -z "${files}" ]; then
    echo -e "Option -f empty. $usage"; exit 1;
fi

if [ -z "${output}" ]; then
    echo -e "Option -d empty. $usage"; exit 1;
fi

## Successful call of the script parameters

## remove any leading and trailing spaces
OUTPUT="$(echo -e "${output}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

echo -e "data: $files \noutput directory: $OUTPUT \n"

########################## BEGIN SCRIPT #####################################
echo "This script runs some commands to summarise reads."

cd $files

## The files output
PRIMER_1=${output}/isd_crete_reads_1_primer_summary.tsv
PRIMER_2=${output}/isd_crete_reads_2_primer_summary.tsv
Ns=${output}/isd_crete_reads_Ns-summary.tsv

zcat *.fastq.gz | gawk -v pr1=$PRIMER_1 -v pr2=$PRIMER_2 -v ns=$Ns 'BEGIN{RS="@" ; FS="\n"}{
reads[$1]=1;
match($2,/ACTCCTACGGGAGGCAGCAG/) ; a[RSTART]++ ; 
match($2,/GGACTACHVGGGTWTCTAAT/) ; b[RSTART]++ ;
if ($2 ~ /N/){split($1,sample,".") ; count[sample[1]][gsub(/N/,"")]++ }
}END{print "total reads = " length(reads);
print "Find the reads with primers, both forward and reverse. \
The start position of the primer is the first column and the number of reads the second";
print "start_position" "\t" "total_reads" > pr1; for (i in a){print i "\t" a[i] >> pr1 };
print "start_position" "\t" "total_reads" > pr2; for (j in b){print j "\t" b[j] >> pr2 };
print "how many reads have Ns and count them per sample";
print "sample" "\t" "N_of_Ns" "\t" "total_reads" > ns; for (s in count){for (c in count[s]){print s "\t" c "\t" count[s][c] >> ns}}}'
########################## END ###############################################
# The following oneliners are 2X slower because the same files are read 
# 3 times
#echo "Calculate number of reads"
#zcat *.fastq.gz | gawk 'BEGIN{RS="@" ; FS="\n"}{reads[$1]=1}END{print "total reads = " length(reads)}'
#echo "Find the reads with primers, both forward and reverse"
#zcat *.fastq.gz  | gawk -v pr1=$PRIMER_1 'BEGIN{RS="@" ; FS="\n"}{match($2,/ACTCCTACGGGAGGCAGCAG/) ; a[RSTART]++}END{for (i in a){print i "\t" a[i] > pr1 }}'
#zcat *.fastq.gz  | gawk 'BEGIN{RS="@" ; FS="\n"}{match($2,/GGACTACHVGGGTWTCTAAT/) ; a[RSTART]++}END{for (i in a){print i "\t" a[i]}}' > $PRIMER_2
#echo "Find how many reads have Ns and count them per sample"
#zcat *.fastq.gz |  gawk 'BEGIN{RS="@" ; FS="\n"}($2 ~ /N/){split($1,sample,".") ; count[sample[1]][gsub(/N/,"")]++ }END{print "sample" "\t" "N_of_Ns" "\t" "total_reads"; for (s in count){for (c in count[s]){print s "\t" c "\t" count[s][c]}}}' > $Ns


