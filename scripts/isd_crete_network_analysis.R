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



########################### Centralities ###############################

V(graph)$degree <- degree(graph)
V(graph)$strength <- strength(graph)

graph_tbl <- graph |>
    as_tbl_graph() |>
    activate(nodes) |>
    left_join(network_taxa_metadata, by=c("label"="asv_id")) 

############################## plot #####################################

png(file="figures/network_asv_sign.png",
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

gg_circle <-  ggraph(graph_tbl, circular=T, weights=abs(weight)) + 
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
