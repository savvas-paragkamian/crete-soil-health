#!/bin/bash -l

#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --mem=400G
#SBATCH --job-name="isd-crete workflow"
#SBATCH --mail-user=s.paragkamian@hcmr.gr
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --requeue

start=`date +%s`

###################### ENVIRONMENT #########################
module purge # unloads all previous loads

# activate the environment
module load anaconda3/default

conda activate R-4.3.1

# paths
repository_path="/home1/s.paragkamian/isd-crete/"
data_path="/home1/s.paragkamian/isd-crete/ena_data"
output_path="/home1/s.paragkamian/isd-crete/dada2_output"

############################ get ISD Crete data ##############################

cd $repository_path

echo "retrieve the filereport of ENA for the ISD Crete project id PRJEB21776"

curl -o ena_metadata/filereport_read_run_PRJEB21776_tsv.txt \
    'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB21776&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp,bam_ftp&format=tsv'

############## Metadata - Attributes ##############

echo "retrive metadata xml and transform to tsv"
cd scripts/

./get_isd_crete_2016_attributes.py

cd $repository_path
./scripts/ena_xml_to_csv.py

############ sequences - fastq files ############
cd $data_path
./get_isd_crete_2016_fastq.sh

####################### Amplicon to Taxonomy pipeline #########################

############ dada2 pipeline ################
echo "start the DADA2 pipeline for ASV inference"

cd $repository_path
./scripts/isd_crete_dada2_taxonomy.R $data_path $output_path

dada2=`date +%s`
runtime_dada2=$((dada2-start))
echo $((runtime_dada2/60)) " minutes for the DADA2 pipeline" 



############# helper oneliners ##############
echo "start the helper oneliners"

cd $output_path
cd taxonomy

gawk -F"\t" 'BEGIN{print "file" "\t" "asv_id" "\t" "abundance"} \ 
    (NR==1){split($0,asv,"\t")}(NR>1){split($0, data, "\t"); \
    for (i=2; i<=length(data); ++i){print data[1] "\t" "asv"i-1 "\t" data[i]}}' seqtab_nochim.tsv > seqtab_nochim_long.tsv

gawk -F"\t" 'BEGIN{print "asv_id" "\t" "asv"}(NR==1){split($0,asv,"\t"); \
    for (i in asv){print "asv"i "\t" asv[i]}}' seqtab_nochim.tsv > asv_fasta_ids.tsv

cd $repository_path
./scripts/isd_crete_reads_summary.sh -i filtered -o $output_path

####################### Microbial Ecology analysis #########################
# enter the working directory
cd $repository_path

echo "metadata enrichment with spatial and remote-sensing data script"

./scripts/isd_crete_spatial.R

echo "executing biodiversity script"

./scripts/isd_crete_biodiversity.R

echo "compositionality and differential abundance"
#####./scripts/isd_crete_compositionality.R

echo "executing numerical ecology script"

./scripts/isd_crete_numerical_ecology.R

echo "UMAP dimention reduction scripts"

./scripts/isd_crete_umap.py

echo "Network inferrence and network analysis"
### need to fix the Julia environment
./scripts/isd_crete_network.jl

./scripts/isd_crete_network_analysis.R
############################ Visualisation ################################
echo "create all the figures"

./script/figures.R

end=`date +%s`
runtime=$((end-start))
echo "Job ID: " $SLURM_JOB_ID
echo "Job name: " $SLURM_JOB_NAME
echo $runtime "in seconds" 
echo $((runtime/60)) "in minutes" 
