using Graphs
using FlashWeave

#data_path = "results/community_matrix_genera.tsv"
data_path = string("results/network_genera_community_matrix.tsv")

#meta_data_path = "results/sample_metadata.tsv"
netw_results = learn_network(data_path, sensitive=true, heterogeneous=false)

save_network("results/network_output.gml", netw_results, detailed=true)

