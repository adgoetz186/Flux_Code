import json
import gzip
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import numpy as np


test = pd.read_csv("GSM3271040_RNA_sciCAR_A549_cell.txt.gz",compression='gzip')
test_2 = pd.read_csv("/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Input/scRNA-seq-data/GSM3271042_RNA_only_A549_gene_count.txt.gz",compression='gzip',delimiter=",")
test_3 = pd.read_csv("/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Input/scRNA-seq-data/GSM3271042_RNA_only_A549_gene.txt.gz",compression='gzip')
print(test)
cell_names = test["cell_name"].to_numpy()
print(np.size(np.argwhere(cell_names=="A549")))
print(test["experiment"].tolist())
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
