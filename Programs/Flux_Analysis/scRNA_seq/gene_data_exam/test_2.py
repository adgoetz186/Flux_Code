import json
import gzip
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import numpy as np
import json
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf
import os
import random
import numpy as np
import mygene
import gzip
import re



mg = mygene.MyGeneInfo()





#test = pd.read_csv("matrix.mtx.gz",compression='gzip')
test_2f = pd.read_csv("barcodesf.tsv.gz", compression='gzip', delimiter="\t", header=None)
test_3f = pd.read_csv("featuresf.tsv.gz", compression='gzip', delimiter="\t", header=None)

test_2 = pd.read_csv("barcodes.tsv.gz", compression='gzip', delimiter="\t", header=None)
test_3 = pd.read_csv("features.tsv.gz", compression='gzip', delimiter="\t", header=None)
#print(test)
print(test_2f)
print(test_2)


#rna = lf.RNA_Seq_Load("Data/Input/scRNA_seq-data/A549_sciCAR_data.mat")["RNA"]["Features"]

mygene_search_results = mg.querymany(test_3[0].tolist()[:10])
print(mygene_search_results)
input()
symbol_genes_list = [i["symbol"] for i in mygene_search_results]
with open("Data/Input/Symbol_Handling/Model_Gene_Symbols/recon2_2","w") as gene_file:
  gene_file.write(str(symbol_genes_list))

print(test_3f)
print(test_3)
print(np.alltrue((test_3==test_3f).to_numpy()))
scRNAseq_gene_list = test_3[1].tolist()
list_of_dup = []
count = 0
for i in scRNAseq_gene_list:
	count+=1
	if count % 100 == 0:
		print(count/len(scRNAseq_gene_list))
	for j in scRNAseq_gene_list:
		if i.lower() == j.lower() and i != j:
			if [i,j] not in list_of_dup and [j,i] not in list_of_dup:
				print(count, len(list_of_dup)+1, len(scRNAseq_gene_list))
				list_of_dup.append([i,j])
print(list_of_dup)
input()


print(np.size(np.argwhere(cell_names=="A549")))
print(test_2.columns)
full_data = test_2.to_numpy()
print(int(full_data[0][0].split()[0]),int(full_data[0][0].split()[1]))
scRNA_seq = np.zeros((int(full_data[0][0].split()[0]),int(full_data[0][0].split()[1])))
entry_values = test_2.to_numpy()[1:]
for i in entry_values:
	val = i[0].split()
	scRNA_seq[int(val[0])-1,int(val[1])-1] = int(val[2])
print(np.sum(scRNA_seq, axis=0))
print(np.size(np.argwhere(np.sum(scRNA_seq, axis=0))))
print(np.sum(scRNA_seq,axis = 1))
