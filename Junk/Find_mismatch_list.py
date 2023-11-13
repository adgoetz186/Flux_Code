import json
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf
import os
from collections import Counter
import random
import numpy as np
import mygene
import gzip
import re



mg = mygene.MyGeneInfo()
os.chdir("..")


scRNAseq_gene_list = lf.RNA_Seq_Load("Data/Input/Symbol_Handling/scRNA_seq-data/A549_sciCAR_data.mat")["RNA"]["Features"]

#mygene_search_results = mg.querymany(hgnc_genes_list,scopes='hgnc',species = "human")
#symbol_genes_list = [i["symbol"] for i in mygene_search_results]
with open("Data/Input/Symbol_Handling/Model_Gene_Symbols/recon2_2","r") as gene_file:
	model_gene_list = eval(gene_file.readlines()[0])
	
# ignore case
# model_gene_list = [i.upper() for i in model_gene_list]
# scRNAseq_gene_list = [i.upper() for i in scRNAseq_gene_list]

model_genes_with_no_scRNAseq_entry = []

for i in model_gene_list:
	if not i in scRNAseq_gene_list:
		model_genes_with_no_scRNAseq_entry.append(i)

with open("Data/Input/Symbol_Handling/Model_scRNA_mismatch_symbols/recon2_2","w") as mismatch_file:
	mismatch_file.write(str(model_genes_with_no_scRNAseq_entry))