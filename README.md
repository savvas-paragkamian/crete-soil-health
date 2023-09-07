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

* [Data retrieval](#data-retrieval)
* [Sequences](#sequences)
* [Taxonomic assignement](#taxonomic-assignement)
* [Metadata](#metadata)
* [Spatial data](#spatial-data)
* [Software](#software)
* [Hardware](#hardware)
* [Citation](#citation)
* [Licence](#licence)

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

## Sequences

In this [script](scripts/isd_crete_reads_summary.sh) there are some oneliners
to provide some basic information regarding the reads.

Total reads (forward and reverse) = 121232490 in 140 samples

Primers used are FWD: 5'-ACTCCTACGGGAGGCAGCAG-3' REV: 5'-GGACTACHVGGGTWTCTAAT-3'

Numbers of Ns in reads : there are many reads with Ns and the `dada` function
doesn't accept so after the filtering they are all removed.

## Taxonomic assignement
We used PEMA and DADA2 for the clustering of OTUs and ASVs, respectively.

### PEMA

PEMA incorporates state of the art tools for each step of the analysis.

The filtering step is slower than DADA2 because of Trimmomatic.

### DADA2
DADA2 for our dataset required a total of 1121 minutes (18 hours, 41 minutes)
to run on a FAT node of the ZORBAS HPC - IMBBC - HCMR.

FAT node Specs : Intel(R) Xeon(R) Gold 6230 CPU @ 2.10GHz, 40 CPUs, 503gb RAM
OS = Debian 4.19.146-1

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
