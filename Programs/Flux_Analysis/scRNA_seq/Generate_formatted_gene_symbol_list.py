import os
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf

# Takes in a list of symbol gene names and creates file formatted to be used by HGNC checker

# Navigates to main folder
os.chdir("../../..")

# Loads scRNAseq data
scRNAseq_gene_list = list(lf.RNA_Seq_Load("Data/Input/scRNA-seq-data/A549_sciCAR_data.mat")["RNA"]["Features"])

with open("Data/Input/nomenclature_files/symbol_list/A549_gene_names",'w') as symbol_file:
	for i in scRNAseq_gene_list:
		symbol_file.write(f"{i},\n")

