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
# usage:./isd_crete_dada2_taxonomy.R
###############################################################################

# packages
library(dada2, lib.loc="/home1/s.paragkamian/software/R/4.1.1")

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

################# REMOVE ################################
#fnFs <- head(fnFs)
#fnRs <- head(fnRs)
################# REMOVE ################################

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

print(paste0("Samples processing : ", length(sample.names)))

# 1. quality control
create_dir <- function(dir_name){

    if (!dir.exists(dir_name)){
    dir.create(dir_name)
    print(paste0(dir_name, " directory created", sep=""))
    } else{
        print(paste0(dir_name, " directory exists", sep=""))
    }
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

create_dir("quality_plots")
all_fn <- c(fnFs,fnRs)
quality_plots(all_fn)
#stop("Manual break inserted here")

# 2. filter and trim sequences
print("start the filter and trim")

# Place filtered files in filtered/ subdirectory
create_dir("filtered")
filtFs <- file.path(output_path,
                    "filtered",
                    paste0(sample.names, "_1_filt.fastq.gz"))
filtRs <- file.path(output_path,
                    "filtered",
                    paste0(sample.names, "_2_filt.fastq.gz"))

filtered <- filterAndTrim(fwd=fnFs, filt=filtFs,
                     rev=fnRs, filt.rev=filtRs,
                     truncLen=c(280,230),
                     maxN=0,
                     maxEE=2,
                     truncQ=2,
                     minLen=100,
#                     rm.phix=TRUE,
                     compress=TRUE,
                     multithread=TRUE,
                     verbose=T)

#write.table(filtered,
#            paste0(output_path,"/filtered_summary.tsv", sep=""),
#            sep="\t",
#            col.names = TRUE)

# 3. Learn errors

print("learning errors")

create_dir("errors")

set.seed(100)
## Errors forward
err_F <- learnErrors(filtFs,
                      multithread=TRUE,
                      randomize=TRUE,
                      MAX_CONSIST=20,
                      verbose=1)

saveRDS(err_F, paste0(output_path,"/errors/err_F.rds", sep=""))

g_errors_F <- plotErrors(err_F, nominalQ=TRUE)
ggplot2::ggsave(plot=g_errors_F,
                path="errors/",
                filename = "g_errors_F.png",
                device="png") 

## Errors reverse
err_R <- learnErrors(filtRs,
                      multithread=TRUE,
                      randomize=TRUE,
                      MAX_CONSIST=20,
                      verbose=1)

saveRDS(err_R, paste0(output_path,"/errors/err_R.rds", sep=""))
g_errors_R <- plotErrors(err_R, nominalQ=TRUE)
ggplot2::ggsave(plot=g_errors_R,
                path="errors/",
                filename = "g_errors_R.png",
                device="png") 

# 4. Sample Inferrence
print("sample inference")
create_dir("taxonomy")

dadaFs <- dada(filtFs,
               err=err_F,
               multithread=TRUE)

print("check convergence Fs")
dada2:::checkConvergence(dadaFs[[1]])

dadaRs <- dada(filtRs,
               err=err_R,
               multithread=TRUE)

print("check convergence Rs")
dada2:::checkConvergence(dadaRs[[1]])
# Merge pairs
print("merge pairs")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                      verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

write.table(seqtab,
            paste0(output_path,"/taxonomy/seqtab.tsv", sep=""),
            sep="\t",
            col.names = TRUE,
            row.names=TRUE)

saveRDS(seqtab, paste0(output_path,"/taxonomy/seqtab.RDS", sep=""))
# # Remove chimeras
print("remove chimeras")
seqtab.nochim <- removeBimeraDenovo(seqtab,
                             method="consensus",
                             multithread=TRUE)

write.table(seqtab.nochim,
            paste0(output_path,"/taxonomy/seqtab_nochim.tsv", sep=""),
            sep="\t",
            col.names = TRUE,
            row.names=TRUE)

saveRDS(seqtab.nochim, paste0(output_path,"/taxonomy/seqtab_nochim.RDS", sep=""))
## Summary

getN <- function(x) sum(getUniques(x))

track <- cbind(filtered, 
               sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: 
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,
            paste0(output_path,"/isd_crete_read_track_summary.tsv", sep=""),
            sep="\t",
            col.names = TRUE)

# # Assign Taxonomy
# Add your path to SILVA

print("assign taxonomy")
taxa <- assignTaxonomy(seqtab.nochim,
                       "/home1/s.paragkamian/databases/SILVA_138_SSU/silva_nr99_v138.1_wSpecies_train_set.fa",
                       multithread=20,
                       tryRC = TRUE, 
                       verbose = TRUE)

saveRDS(taxa, paste0(output_path,"/taxonomy/dada2_taxonomy.RDS", sep=""))
print("assignTaxonomy done.")

# Add your path to SILVA
# This function requires more than 250 gb of memory (RAM)!!
taxa <- addSpecies(taxa,
                   "/home1/s.paragkamian/databases/SILVA_138_SSU/silva_species_assignment_v138.1.fa",
                   tryRC = TRUE)

print("addSpecies done.")

print("begin saving data.")
saveRDS(taxa, paste0(output_path,"/taxonomy/dada2_taxa_species.RDS", sep=""))
print("data saved.")
#stop("Manual break inserted here")
