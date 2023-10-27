#!/bin/bash -l

#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --mem=400G
#SBATCH --job-name="isd-species"
#SBATCH --mail-user=s.paragkamian@hcmr.gr
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --requeue

start=`date +%s`

data_path="/home1/s.paragkamian/isd-crete/ena_data"
output_path="/home1/s.paragkamian/isd-crete/dada2_output"
module purge # unloads all previous loads

module load  R/4.1.1 #loads  R/4.1.1
#/home1/s.paragkamian/isd-crete/scripts/isd_crete_dada2_taxonomy.R $data_path $output_path


dada2=`date +%s`
runtime_dada2=$((dada2-start))
echo $((runtime_dada2/60)) " minutes for the DADA2 pipeline" 
echo "start the helper oneliners"

cd $output_path
cd taxonomy

gawk -F"\t" 'BEGIN{print "file" "\t" "asv_id" "\t" "abundance"} \ 
    (NR==1){split($0,asv,"\t")}(NR>1){split($0, data, "\t"); \
    for (i=2; i<=length(data); ++i){print data[1] "\t" "asv"i-1 "\t" data[i]}}' seqtab_nochim.tsv > seqtab_nochim_long.tsv

gawk -F"\t" 'BEGIN{print "asv_id" "\t" "asv"}(NR==1){split($0,asv,"\t"); \
    for (i in asv){print "asv"i "\t" asv[i]}}' seqtab_nochim.tsv > asv_fasta_ids.tsv

cd ../
/home1/s.paragkamian/isd-crete/scripts/isd_crete_reads_summary.sh -i filtered -o $output_path

module purge

end=`date +%s`
runtime=$((end-start))
echo "Job ID: " $SLURM_JOB_ID
echo "Job name: " $SLURM_JOB_NAME
echo $runtime "in seconds" 
echo $((runtime/60)) "in minutes" 
