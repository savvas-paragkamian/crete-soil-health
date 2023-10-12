#!/usr/bin/env python3

###############################################################################
# script name: isd_crete_umap.py
# developed by: Savvas Paragkamian
# framework: ISD Crete 2016
###############################################################################
# GOAL:
# Aim of this script is to apply the umap dimention reduction method 
# to the metagenomic data.
###############################################################################
# usage:./isd_crete_umap.py
###############################################################################
import pandas as pd
import numpy as np
import os, sys
import umap

community_matrix_path = "results/genera_samples_matrix.tsv"
samples_umap_file = "results/umap_samples.tsv"
taxa_umap_file = "results/umap_taxa.tsv"

community_matrix = pd.read_csv(community_matrix_path, sep="\t")
community_array = community_matrix.drop(columns=['ENA_RUN'])
# taxa, the colnames to a new dataframe
cols_list = community_array.columns.tolist()
df_taxa = pd.DataFrame(cols_list)
#print(df_taxa)
# samples
df_samples = community_matrix[['ENA_RUN']]
#print(df_samples)
community_array = community_array.to_numpy()
#print(community_array.shape)

# for taxa UMAP the array needs to be transposed 
community_array_t = np.transpose(community_array)
#print(community_array_t.shape)

def umap_function(array, df_ids):
    
    np.random.seed(123)
    # umap function
    embedding = umap.UMAP(n_neighbors=80,
                          min_dist=0.3,
                          n_components=3,
                          metric='braycurtis').fit_transform(array)
    
    #add NumPy matrix as new columns in DataFrame
    df = pd.concat([df_ids, pd.DataFrame(embedding)], axis=1)
    df.columns =['id', 'UMAP1', 'UMAP2', 'UMAP3']

    return(df)

# samples
print("UMAP samples is ongoing")
samples_umap = umap_function(community_array,df_samples)
# save to file
samples_umap.to_csv(samples_umap_file, header=True, index=None, sep='\t')

# taxa
print("UMAP taxa is ongoing")
taxa_umap = umap_function(community_array_t,df_taxa)
# save to file
taxa_umap.to_csv(taxa_umap_file, header=True, index=None, sep='\t')

