#!/usr/bin/env Rscript

###############################################################################
# script name: isd_crete_numerical_ecology.R
# developed by: Savvas Paragkamian, Johanna Holms
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to use Crete biodiversity data of ASVs, taxonomy and 
# the sample metadata to perform ecological analyses on biodiversity, ordination
# and multivariate comparison.
#
###############################################################################
# OUTPUT:
#
###############################################################################
# usage:./isd_crete_numerical_ecology.R
###############################################################################
#library(mia)
library(phyloseq)
library(vegan)
library(dplyr)
library(tibble)
library(readr)
library(magrittr)
library(tidyr)
library(ggplot2)

################################## Load data ##################################
crete_biodiversity <- read_delim("results/crete_biodiversity_asv.tsv",delim="\t")

crete_biodiversity_matrix <- readRDS("results/crete_biodiversity_matrix.RDS")
tax_tab <- readRDS("results/tax_tab.RDS")

# Metadata

master_metadata_old <- read.delim("Crete/Composite_MetaData_from_master.csv", sep=",")

metadata <- read_delim("results/metadata_all.tsv", delim="\t")

## location pairs of each site
metadata$sites <- gsub("_loc_*.","", metadata$source_material_identifiers)
metadata$location <- gsub(".*(loc_.*)","\\1", metadata$source_material_identifiers)

samples_locations <- metadata %>%
    dplyr::select(ENA_RUN,sites,location) %>%
    pivot_wider(id_cols=sites,
                names_from=location,
                values_from=ENA_RUN)

# differences of old and new data

master_metadata_old$team_site_location_id[which(!(master_metadata_old$team_site_location_id %in% metadata$source_material_identifiers))]

################################## functions ##################################
dist_long <- function(x,method){
    method <- method
    df <- as.data.frame(as.matrix(x)) %>%
    rownames_to_column() %>%
    pivot_longer(-rowname,
                 values_to=method,
                 names_to="colname")

    return(df)
}


####################### Destriptors Statistics ###############################
pw <- pairwise.wilcox.test(metadata_all$Observed, metadata_all$vegetation_zone, p.adjust.method="BH")
pw_e <- pairwise.wilcox.test(metadata_all$Observed, metadata_all$elevation_bin, p.adjust.method="BH")
pw_s <- pairwise.wilcox.test(metadata_all$Shannon, metadata_all$LABEL1, p.adjust.method="BH")
print(pw)


box_shannon <- ggplot(data=metadata_all, mapping=aes(x=LABEL2, y=Shannon))+
    geom_boxplot()+
    geom_jitter(width = 0.2)+
    theme_bw()+
    theme(axis.text.x = element_text(face="bold",
                                     size = 10,
                                     angle = 45,
                                     vjust = 1,
                                     hjust=1))


ggsave("figures/box_shannon.png", 
       plot=box_shannon, 
       device="png", 
       height = 23, 
       width = 23, 
       units="cm")


# correlations of diversity and other numerical metadata
with(metadata_all, cor(Shannon, water_content))


metadata_n <- metadata_all
rownames(metadata_n) <- metadata_all$ENA_RUN
nums <- unlist(lapply(metadata_n, is.numeric), use.names = FALSE)
metadata_n <- metadata_n[,c(nums)]
cc <- cor(metadata_n)

cc_sp <- cor(metadata_n, method="spearman")

write.table(cc_sp,
            "results/metadata_sprearman.tsv",
            sep="\t",
            row.names=T,
            col.names=NA)

############################## Community analysis ###########################
###################### Co-occurrence of samples and ASVs ####################

biodiversity_m <- biodiversity_srs
biodiversity_m[biodiversity_m > 0 ] <- 1
biodiversity_m <- as.matrix(biodiversity_m)

## matrix multiplication takes up a lot of memory and CPU, I had an error
## Error: vector memory exhausted (limit reached?)
## cd ~ ; touch .Renviron 
## echo R_MAX_VSIZE=100Gb >> .Renviron

asv_cooccur <- biodiversity_m %*% t(biodiversity_m)
sample_cooccur <- t(biodiversity_m) %*% biodiversity_m
isSymmetric(sample_cooccur) # is true so we can remove the lower triangle
sample_cooccur[lower.tri(sample_cooccur)] <- NA

sample_cooccur_l <- dist_long(sample_cooccur,"cooccurrence") %>%
    filter(rowname!=colname) %>%
    na.omit()


############################## Dissimilarity ###########################
# use the vegan package, the matrix must be transposed
biodiversity_srs_t <- t(biodiversity_srs)

bray <- vegdist(biodiversity_srs_t,
                method="bray")

jaccard <- vegdist(biodiversity_srs_t,
                method="jaccard",
                binary=TRUE)

aitchison <- vegdist(biodiversity_srs_t,
                method="robust.aitchison")


bray_l <- dist_long(bray, "bray")
jaccard_l <- dist_long(jaccard, "jaccard")
aitchison_l <- dist_long(aitchison, "robust.aitchison")



######################## Site locations comparison #################
dissi_loc <- samples_locations %>%
    left_join(sample_cooccur_l,
              by=c("loc_1"="rowname", "loc_2"="colname")) %>%
    left_join(bray_l,
              by=c("loc_1"="rowname", "loc_2"="colname")) %>%
    left_join(jaccard_l,
              by=c("loc_1"="rowname", "loc_2"="colname")) %>%
    left_join(aitchison_l,
              by=c("loc_1"="rowname", "loc_2"="colname"))


summary(dissi_loc)

#####################################################################
z <- betadiver(biodiversity_srs_t, "z")
mod <- with(metadata_all, betadisper(z, LABEL1))
#sac <- specaccum(biodiversity_srs_t)

# Ordination
nmds <- vegan::metaMDS(biodiversity_srs_t,
                       k=2,
                       distance = "bray",
                       trymax=100)
stressplot(nmds)
ordiplot(nmds,display="sites", cex=1.25)
ordisurf(nmds,metadata_all$dem,main="",col="forestgreen")
ordihull(nmds,display="sites",label=T,  groups=metadata_all$LABEL1, cex=1.25)
ordihull(nmds,display="sites",label=T,  groups=metadata_all$elevation, cex=1.25)


tse <- runUMAP(tse,
               name = "UMAP",
               assay.type = "counts",
               ncomponents = 3)
plotReducedDim(tse, "UMAP",
               colour_by = "Group",
               ncomponents = c(1:3))
#### tests 

nmds <- vegan::metaMDS(t(crete_biodiversity_matrix),
                       k=2,
                       distance = "bray",
                       trymax=100)

stressplot(nmds)
ordiplot(nmds,display="sites", cex=1.25)

# ancombc analysis

assays = SimpleList(counts = crete_biodiversity_matrix)
smd = DataFrame(metadata1)
tax_tab = DataFrame(tax_tab)

## create TSE
tse = TreeSummarizedExperiment(assays = assays,
                               colData = smd,
                               rowData = tax_tab)
## convert TSE to phyloseq
pseq = makePhyloseqFromTreeSummarizedExperiment(tse)

total = median(sample_sums(pseq))
standf = function(x, t=total) round(t * (x / sum(x)))

pseq_std  = transform_sample_counts(pseq,standf )
gpsf = filter_taxa(pseq_std, function(x) sd(x)/mean(x) > 3.0, TRUE)

gpsfb = subset_taxa(gpsf, Order=="Frankiales")

jpeg(file="saving_plot1.jpeg")
plot_bar(gpsfb, "LABEL1", "Abundance", "Family", title="test")
dev.off()


GP = pseq_std
#wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
#GP1 = prune_taxa(wh0, GP)

#GP1 = transform_sample_counts(GP, function(x) 1E6 * x/sum(x))

GP.ord <- ordinate(pseq, "NMDS", "bray")

jpeg(file="saving_plot1.jpeg")
p1 = plot_ordination(GP, GP.ord, type="samples", color="LABEL1", title="taxa")
p1
dev.off()


###
set.seed(123)
output = ancombc2(data = tse, assay_name = "counts", tax_level = "Family",
                  fix_formula = NULL, rand_formula = NULL,
                  p_adj_method = "holm", 
                  prv_cut = 0.10, lib_cut = 100, s0_perc = 0.05,
                  group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)

data(dietswap, package = "microbiome")
    tse = mia::makeTreeSummarizedExperimentFromPhyloseq(dietswap)
    colData(tse)$bmi_group = factor(colData(tse)$bmi_group,
                                    levels = c("obese",
                                               "overweight",
                                               "lean"))
set.seed(123)
# Note that setting max_iter = 1 and B = 1 is only for the sake of speed
# Use default or larger values for max_iter and B for better performance
out = ancombc2(data = tse, assay_name = "counts", tax_level = "Phylum",
               fix_formula = "nationality + timepoint + bmi_group",
               rand_formula = "(timepoint | subject)",
               p_adj_method = "holm",
               prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
               group = "bmi_group", struc_zero = TRUE, neg_lb = TRUE,
               alpha = 0.05, n_cl = 1, verbose = TRUE, global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
               iter_control = list(tol = 1e-2, max_iter = 1, verbose = TRUE),
               em_control = list(tol = 1e-5, max_iter = 1),
               lme_control = lme4::lmerControl(),
               mdfdr_control = list(fwer_ctrl_method = "holm", B = 1),
               trend_control = list(contrast =list(matrix(c(1, 0, -1, 1),nrow = 2, byrow = TRUE)),node = list(2), solver = "ECOS", B = 1))
res_prim = out$res
res_global = out$res_global
res_pair = out$res_pair
res_dunn = out$res_dunn
res_trend = out$res_trend





### highest species biodiversity sample
print("sample with the highest microbial species diversity")
crete_biodiversity_s[which(crete_biodiversity_s$Species==max(crete_biodiversity_s$Species)),]

highest_crete_biodiversity_s <- crete_biodiversity_s %>% 
    left_join(metadata_wide,
              by=c("ENA-RUN"="ENA-RUN")) %>% 
    select(place_name,source_material_identifiers, Species, `ENA-RUN`)

write_delim(highest_crete_biodiversity_s,
            "results/highest_crete_biodiversity_sample_taxonomy.tsv",
            delim="\t")

################################# Taxonomy #####################################
## Phyla distribution

phyla_dist <- crete_biodiversity %>%
    group_by(sampleID,Phylum) %>%
    summarise(sum_abundance=sum(abundance), .groups="keep") %>%
    na.omit(Phylum)

total_phyla_dist <- phyla_dist %>% 
    group_by(Phylum) %>%
    summarise(sum_abundance=sum(sum_abundance),
              n_samples=n()) %>%
    mutate(sampleID="All samples") %>%
    arrange(desc(n_samples))


## create bar plots for each sample at family level, class level, Phylum etc
## shannon per sample
## bray curtis per sample
## 
