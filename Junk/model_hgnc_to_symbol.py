import json
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf
import os
import random
import numpy as np
import mygene
import gzip
import re



mg = mygene.MyGeneInfo()
os.chdir("..")


#rna = lf.RNA_Seq_Load("Data/Input/scRNA_seq-data/A549_sciCAR_data.mat")["RNA"]["Features"]
mat = lf.loadmat("Data/Input/MAT_Models/recon2_2.mat")
hgnc_genes_list = list(set([int(i.split(":")[-1]) for i in mat['MODEL1603150001']["genes"]]))
mygene_search_results = mg.querymany(hgnc_genes_list,scopes='hgnc',species = "human")
symbol_genes_list = [i["symbol"] for i in mygene_search_results]
with open("Data/Input/Symbol_Handling/Model_Gene_Symbols/recon2_2","w") as gene_file:
  gene_file.write(str(symbol_genes_list))