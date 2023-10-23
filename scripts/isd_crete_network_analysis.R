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

library(igraph)
library(tibble)
library(readr)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidygraph)
library(ggraph)

############################### load data #####################################
genera_phyla_stats <- read_delim("results/genera_phyla_stats.tsv", delim="\t") %>% dplyr::select(-reads_srs_sd)
graph <- read_graph("results/network_output.gml", format="gml")

genera_phyla_stats_g <- genera_phyla_stats %>% dplyr::filter(Genus %in% V(graph)$label)

for (attrib in colnames(genera_phyla_stats)){
    print(attrib)
    set_vertex_attr(graph,
                    var,
                    value=genera_phyla_stats[,attrib])

}

E(graph)$sign <- ifelse(E(graph)$weight > 0, "positive", "negative")


graph_positive <- graph
E(graph_positive)$weight <- abs(E(graph)$weight)


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



##############as_tbl_graph#################### Centralities ###############################

V(graph)$degree <- degree(graph)
V(graph)$strength <- strength(graph)

graph_tbl <- graph %>%
    as_tbl_graph() %>%
    activate(nodes) %>% 
    left_join(genera_phyla_stats, by=c("label"="Genus"))

############################## plot #####################################
E(graph)$color <- ifelse(E(graph)$sign == "positive", "green", "red")

png(file="figures/network_genera_sign.png",
    width = 50,
    height = 30,
    res=300,
    units = "cm",
    bg="white")

plot(graph, layout=layout_in_circle,
     vertex.label=NA,
     vertex.size=V(graph)$degree/10,
     edge.color=E(graph)$color,
     edge.size=abs(E(graph)$weight)/100,
     vertex.frame.color="coral2",
     vertex.color="coral2",
     main = "Crete soil microbial interactome")
dev.off()


ggraph(graph_tbl, layout = 'fr') + 
  geom_edge_link(aes(color=color)) + 
  geom_node_point(aes(size = centrality_pagerank())) + 
  theme(legend.position = 'bottom')

