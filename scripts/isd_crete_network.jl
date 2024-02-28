using Graphs
using FlashWeave

data_path = string("results/network_community_matrix.tsv")
metadata_path = string("results/sample_metadata.tsv")

#meta_data_path = "results/sample_metadata.tsv"
netw_results = learn_network(data_path, sensitive=true, heterogeneous=false, verbose=true, max_k=3, alpha=0.05, FDR=true)

#netw_results_heterogenous = learn_network(data_path, metadata_path, sensitive=false, heterogeneous=true, verbose=true, max_k=6, alpha=0.01, FDR=true)

save_network("results/network_output.gml", netw_results, detailed=true)

