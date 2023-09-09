#!/usr/bin/Rscript

###############################################################################
# script name: isd_crete_spatial.R
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to use spatial data to enrich the metadata of the 
# ISD Crete project
#
###############################################################################
# OUTPUT:
#
###############################################################################
# usage:./isd_crete_spatial.R
###############################################################################
library(sf)
library(dplyr)
library(tibble)
library(readr)
library(magrittr)
library(tidyr)
