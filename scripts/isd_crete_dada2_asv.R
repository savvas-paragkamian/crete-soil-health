#!/usr/bin/env Rscript

###############################################################################
# script name: isd_crete_dada2_asv.R
# developed by: Savvas Paragkamian, Johanna Holms
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to use the DADA2 pipeline to infer ASVs 
# from fastq.gz files of ISD Crete 2016
###############################################################################
# usage:./isd_crete_dada2_asv.R
###############################################################################

# packages
library(dada2,lib.loc="/home1/s.paragkamian/software/R/4.1.1")

# Working Environment
# path of the sequences as retrieved from ENA
path <- "/home1/s.paragkamian/isd-crete/ena_data"
output_path <- "/home1/s.paragkamian/isd-crete/dada2_output"

setwd(output_path)

# to view the files in the directory
#list.files(path)

# Forward and reverse fastq filenames have the ENA format: 
# sampleENAid_1.fastq.gz and sampleENAid_2.fastq.gz
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# 1. quality control

if (!dir.exists("quality_plots")){
    dir.create("quality_plots")
    print("quality_plots directory created")
} else{
    print("dir exists")
}
# plots of the Phred score
## function to create individual plots for each file
quality_plots <- function(all_fn) {

    names(all_fn) <- sapply(strsplit(basename(all_fn), "\\."), `[`,1)

    for (i in seq_along(all_fn)){
        
        quality_plot <- plotQualityProfile(all_fn[i])
        ggplot2::ggsave(plot=quality_plot,
                        path="quality_plots/",
                        filename = paste0(names(all_fn[i]), ".png",sep=""),
                        device="png") 
    }
}

all_fn <- c(fnFs,fnRs)
quality_plots(all_fn)
stop("Manual break inserted here")


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# filter and trim sequences
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
head(out)

# 2.
