###############################################################
# Genomic Prediction Pipeline
# GBLUP Model - To Predict GEBVs using SNP Markers
# Author: Dr. Kavyashree N M
###############################################################

# --- 1. Load Libraries ---
# rrBLUP is the industry standard for starting with Genomic Selection
if(!require(rrBLUP)) install.packages("rrBLUP")
library(rrBLUP)
library(tidyverse)
