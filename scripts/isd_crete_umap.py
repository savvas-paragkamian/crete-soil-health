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
from sklearn.datasets import load_digits

community_matrix_path = "results/genera_samples_matrix.tsv"
community_matrix = pd.read_csv(community_matrix_path, sep="\t")
community_array = community_matrix.drop(columns=['ENA_RUN'])
community_array = community_array.to_numpy()
print(type(community_matrix))
print(type(community_array))
digits = load_digits()
print(type(digits.data))

embedding = umap.UMAP(n_neighbors=80,
                      min_dist=0.3,
                      n_components=3,
                      metric='braycurtis').fit_transform(community_array)
print(embedding)
