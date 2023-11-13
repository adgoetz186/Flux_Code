import os
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf

# Takes in a list of symbol gene names and creates file formatted to be used by HGNC checker

# Navigates to main folder
os.chdir("../../..")

# Loads scRNAseq data

print(lf.RNA_Seq_Load("Data/Input/scRNA-seq-data/A549_sciCAR_data.mat"))