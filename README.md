# ISD - Island Sampling Day Crete 2016

This repository contains the scripts for the analysis of the ISD Crete 2016.
ISD Crete 2016 is a soil microbiome sampling from 72 distinct sites around 
Crete island, Greece. 

During the 18th Genome Standards Consortium workshop the participants along
with members from Hellenic Centre for Marine Research organised the Island 
Sampling Day Crete 2016.

The goal of this sampling was to put metagenomic standards into play.

For more info regarding the methodology of the sampling 
visit the [website](https://lab42open-team.github.io/isd-crete-website/).

## Contents

* [Scripts](#scripts)
* [Data retrieval](#data-retrieval)
* [Sequences](#sequences)
* [Inference and Taxonomy](#inference-and-taxonomy)
* [Metadata](#metadata)
* [Spatial data](#spatial-data)
* [Analysis](#analysis)
* [Software](#software)
* [Hardware](#hardware)
* [Citation](#citation)
* [Licence](#licence)

## Scripts
The scripts of the analysis are in the `scripts` folder and cover the following tasks:

#### task : data retrieval
1. get sequences [here](scripts/get_isd_crete_2016_fastq.sh)
2. get metadata [here](scripts/get_isd_crete_2016_attributes.py)
3. transform metadata [script](https://github.com/savvas-paragkamian/isd-crete/blob/main/scripts/ena_xml_to_csv.py)

#### task : Inference and Taxonomy
1. dada2 [hpc job script](https://github.com/savvas-paragkamian/isd-crete/blob/main/scripts/isd_crete_hpc_job_dada2.sh)
2. dada2 analysis [script](https://github.com/savvas-paragkamian/isd-crete/blob/main/scripts/isd_crete_dada2_taxonomy.R)
3. pema [hpc job script](https://github.com/savvas-paragkamian/isd-crete/blob/main/scripts/isd_crete_pema_asv.sh)

output: 

#### task : Analysis
1. biodiversity [script](https://github.com/savvas-paragkamian/isd-crete/blob/main/scripts/isd_crete_biodiversity.R)
2. spatial data [script](https://github.com/savvas-paragkamian/isd-crete/blob/main/scripts/isd_crete_spatial.R)
3. numerical ecology [script](https://github.com/savvas-paragkamian/isd-crete/blob/main/scripts/isd_crete_numerical_ecology.R)
4. figures [script](https://github.com/savvas-paragkamian/isd-crete/blob/main/scripts/figures.R)

## Data retrieval

The sequences of the ISD Crete 2016 are available in ENA project [PRJEB21776](https://www.ebi.ac.uk/ena/browser/view/PRJEB21776)
and were downloaded using the ENA API. The script is available [here](scripts/get_isd_crete_2016_fastq.sh) for the
raw data and [here](scripts/get_isd_crete_2016_attributes.py) for the metadata.
The metadata are in `xml` format and using this [script](https://github.com/savvas-paragkamian/isd-crete/blob/main/scripts/ena_xml_to_csv.py) we transform them to
tabular format.

From 4 samples DNA quality was low and there aren't any sequences nor metadata.

To find which samples don't appear in metadata run the following command:
```
gawk -F"\t" '($2 ~ /source material identifiers/){split($3,sample,"_"); a[sample[2]][sample[4]]++}END{for (i in a){for (j in a[i]){print i "\t" j "\t" a[i][j]}}}' ena_isd_2016_attributes.tsv
```
The samples DNA wasn't sequenced are: 
* isd_7_site-3_loc_2
* isd_10_site_1_loc_1
* isd_10_site_1_loc_1
* isd_10_site_2_loc_2

## Sequences

In this [script](scripts/isd_crete_reads_summary.sh) there are some oneliners
to provide some basic information regarding the reads.

Total reads (forward and reverse) = 121232490 in 140 samples

Primers used are FWD: 5'-ACTCCTACGGGAGGCAGCAG-3' REV: 5'-GGACTACHVGGGTWTCTAAT-3'

Numbers of Ns in reads : there are many reads with Ns and the `dada` function
doesn't accept so after the filtering they are all removed.

## Inference and Taxonomy
We used PEMA and DADA2 for the clustering of OTUs and ASVs, respectively.

### PEMA

PEMA incorporates state of the art tools for each step of the analysis.

The filtering step is slower than DADA2 because of Trimmomatic.

### DADA2
DADA2 for our dataset required a total of 1121 minutes (18 hours, 41 minutes)
to run on a FAT node of the ZORBAS HPC - IMBBC - HCMR.

FAT node Specs : Intel(R) Xeon(R) Gold 6230 CPU @ 2.10GHz, 40 CPUs, 503gb RAM
OS = Debian 4.19.146-1

### Extra handling of DADA2 output

To transform the matrix of the DADA2 output to a long format use the following
command. This one-liner creates an asv id for each asv sequence, which is the
header of the file. Use this if the DADA2 matrix/array object was saved as
table instead of RDS. This creates a matrix that has N colunm names but N+1
columns.


```
gawk -F"\t" 'BEGIN{print "file" "\t" "asv_id" "\t" "abundance"}(NR==1){split($0,asv,"\t")}(NR>1){split($0, data, "\t"); for (i=2; i<=length(data); ++i){print data[1] "\t" "asv"i-1 "\t" data[i]}}' seqtab_nochim.tsv > seqtab_nochim_long.tsv
```

A companion one-liner that contains the sequences of asvs and their created ids.
```
gawk -F"\t" 'BEGIN{print "asv_id" "\t" "asv"}(NR==1){split($0,asv,"\t"); for (i in asv){print "asv"i "\t" asv[i]}}' seqtab_nochim.tsv > asv_fasta_ids.tsv
```

## Metadata

Each `xml` file, i.e. each sample, has a total of 43 attributes which are 

```
source material identifiers
organism
total nitrogen method
ENA-SUBMISSION
geographic location (region and locality)
environment (biome)
total organic C method
ENA-SUBMITTED-FILES
ENA-EXPERIMENT
ENA-FASTQ-FILES
target subfragment
common name
total organic carbon
DNA concentration
vegetation zone
total nitrogen
soil environmental package
environment (material)
amount or size of sample collected
geographic location (depth)
environment (feature)
ENA-LAST-UPDATE
sample collection device or method
place name
ENA-FIRST-PUBLIC
project name
collection date
target gene
ENA-CHECKLIST
geographic location (country and/or sea)
ENA-STUDY
pcr primers
sequencing method
water content method
storage conditions (fresh/frozen/other)
water content
investigation type
geographic location (elevation)
sample volume or weight for DNA extraction
geographic location (latitude)
geographic location (longitude)
current land use (emp 500 soil)
ENA-RUN
```

## Spatial data

Spatial data from Copernicus system are also used. These are the 
[CORINE Land Cover](https://land.copernicus.eu/pan-european/corine-land-cover/clc2018?tab=download)
and the [Digital Elevation Models](https://www.eea.europa.eu/data-and-maps/data/copernicus-land-monitoring-service-eu-dem).
CORINE Land Cover data are used to identify human pressures on the Natura2000
regions and on the hotspots of the arthropod endemic taxa.

The protected areas of [Natura2000 SCI](https://www.eea.europa.eu/data-and-maps/data/natura-14)
(habitats directive) and [Wildlife refugees](https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA)
are used in this analysis. 

## Analysis



## Software

PEMA
DADA2
vegan
tidyverse

ElementTree

GAWK
R 
Python
Conda/Bioconda
## Hardware

Most computations were performed in the Zorbas HPC facility of [IMBBC-HCMR](https://hpc.hcmr.gr),
see here for more [info](https://doi.org/10.1093/gigascience/giab053).

## Citation

## Licence

GNU GPLv3 license (for 3rd party scripts separate licenses apply).
