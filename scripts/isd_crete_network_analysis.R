#!/usr/bin/env Rscript

###############################################################################
# script name: isd_crete_network_analysis.R
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to use analyse the networks of the soil microbiome of
# Crete.
###############################################################################
# OUTPUT:
#
###############################################################################
# usage:./isd_crete_network_analysis.R
###############################################################################

source("scripts/functions.R")
library(igraph)
library(tibble)
library(readr)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(pheatmap)

############################### load data #####################################
community_matrix <- read_delim("results/network_community_matrix.tsv", delim="\t")
network_taxa_metadata <- read_delim("results/network_taxa_metadata.tsv", delim="\t")# %>% dplyr::select(-reads_srs_sd)
graph <- igraph::read_graph("results/network_output.gml", format="gml")

#for (attrib in colnames(network_taxa_metadata)){
#    print(attrib)
#    graph <- graph |> set_vertex_attr(as.character(attrib),
#                    value=network_taxa_metadata[,attrib])
#
#}

E(graph)$sign <- ifelse(E(graph)$weight > 0, "positive", "negative")
E(graph)$color <- ifelse(E(graph)$sign == "positive", "green", "red")

taxa_cooccur <- cor(community_matrix, method="spearman")
isSymmetric(taxa_cooccur) # is true so we can remove the lower triangle
#taxa_cooccur[lower.tri(taxa_cooccur)] <- NA

taxa_cooccur_l <- dist_long(taxa_cooccur,"cooccurrence") %>%
    filter(rowname!=colname) %>%
    na.omit()

png("figures/network_heatmap_spearman_taxa.png",
    res=300,
    width=70,
    height=30,
    unit="cm")

pheatmap(taxa_cooccur,
         color=colorRampPalette(c("white","skyblue", "palevioletred3"))(20))

dev.off()
############################### graph metrics #################################
is.connected(graph)
igraph::is.directed(graph)

vcount(graph)

is.weighted(graph)

edge_density(graph,loops = F)
reciprocity(graph)

###### only for positive values of weights
diameter(graph_positive)

average.path.length(graph_positive)
############################### subgraph metrics ####################### 
graph_motifs <- motifs(graph, 3)


########################### Centralities ###############################

V(graph)$degree <- degree(graph)
V(graph)$strength <- strength(graph)
V(graph)$betweenness <- igraph::betweenness(graph,weights = NA,directed = TRUE)
V(graph)$closeness_all <- igraph::closeness(graph,weights = NA, mode = "all")
V(graph)$transitivity_l <- igraph::transitivity(graph,type = "local",weights = NA)
V(graph)$transitivity_g <- igraph::transitivity(graph,type = "global",weights = NA)

graph_tbl <- graph |>
    as_tbl_graph() |>
    activate(nodes) |>
    left_join(network_taxa_metadata, by=c("label"="asv_id")) |>
    activate(edges) |>
    mutate(weight_original = weight) |>
    mutate(weight = abs(weight))

####################### clustering ###########################
graph_info <- cluster_infomap(graph_tbl)
graph_louvain <- cluster_louvain(graph_tbl)

V(graph_tbl)$louvain <- membership(graph_louvain)
############################## plot #####################################

png(file="figures/network_asv_sign.png",
    width = 50,
    height = 30,
    res=300,
    units = "cm",
    bg="white")

plot(graph_tbl, layout=layout_in_circle,
     vertex.label=NA,
     vertex.size=V(graph_tbl)$degree/10,
     edge.color=E(graph_tbl)$color,
     edge.size=0.001,
     vertex.frame.color=V(graph_tbl)$louvain,
     vertex.color=V(graph_tbl)$louvain,
     main = "Crete soil microbial interactome")
dev.off()


gg <- ggraph(graph_tbl, layout = 'fr', weights=abs(weight)) + 
    geom_edge_link(aes(color=color)) + 
    geom_node_point(mapping=aes(colour=Phylum, size=degree)) + 
    scale_edge_color_manual(values=c("red"="palevioletred3", "green"="darkolivegreen4"))+
    theme(legend.position = 'bottom') +
    theme_bw() +
    theme(
#        panel.background = element_rect(fill='transparent'), #transparent panel bg
        panel.border = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
#        legend.background = element_rect(fill='transparent'), #transparent legend bg
#        legend.box.background = element_rect(fill='transparent'), #transparent legend panel
        line = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom")


ggsave("figures/network_fr.png",
       plot=gg,
       device="png",
       height = 50,
       width = 53,
       dpi = 300,
       units="cm")

gg_circle <-  ggraph(graph_tbl, circular=T) + 
  geom_edge_link(aes(color=color)) + 
  geom_node_point(mapping=aes(colour=Phylum)) + 
  scale_edge_color_manual(values=c("red"="palevioletred3", "green"="darkolivegreen4"))+
  theme(legend.position = 'bottom') +
  theme_bw()

ggsave("figures/network_circular.png",
       plot=gg_circle,
       device="png",
       height = 50,
       width = 53,
       units="cm")

gg_matrix <- ggraph(graph, 'matrix') + 
  geom_edge_point(aes(colour = color, size=abs(weight)), mirror = TRUE) + 
  scale_edge_color_manual(values=c("red"="red", "green"="green"))+
  theme(legend.position = 'bottom')

ggsave("figures/network_matrix.png",
       plot=gg_matrix,
       device="png",
       height = 50,
       width = 53,
       units="cm")
