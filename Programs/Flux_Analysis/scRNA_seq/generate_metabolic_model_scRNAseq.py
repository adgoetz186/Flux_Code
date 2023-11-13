import json
import gzip
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import numpy as np

# Creates scRNAseq matrix for metabolic genes using HGNC nomenclature
#

RNA_Seq_Location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Input/scRNA-seq-data/A549_sciCAR_data.mat"
Conversion_json_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Input/nomenclature_files/to_hgnc_dictionary/to_hgnc_A549_genes.json"
Minimal_Model_Location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_A549/min_model"
scRNA_Seq_Location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Intermediate/scRNAseq_for_model/A549_Recon_2_2_hgnc.csv"

test = pd.read_csv("/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Input/scRNA-seq-data/GSM3271042_RNA_only_A549_cell.txt.gz",compression='gzip')
test_2 = pd.read_csv("/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Input/scRNA-seq-data/GSM3271042_RNA_only_A549_gene_count.txt.gz",compression='gzip')
test_3 = pd.read_csv("/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Input/scRNA-seq-data/GSM3271042_RNA_only_A549_gene.txt.gz",compression='gzip')
print(test["treatment_time"])
print(test_3)
print(set(test["cell_name"].tolist()))
print(np.size(np.argwhere(test["treatment_time"].to_numpy())==0))
print(test_3)
gene_list = test_3["gene_short_name"].tolist()
input()

list_of_dup = []
count = 0
for i in gene_list:
	count+=1
	for j in gene_list:
		if i.lower() == j.lower() and i != j:
			if [i,j] not in list_of_dup and [j,i] not in list_of_dup:
				print(count, len(list_of_dup)+1, len(gene_list))
				print(i,j)
				list_of_dup.append([i,j])
print("Done")
input()
input()
# Loads scRNAseq data
scRNAseq_gene_list = list(lf.RNA_Seq_Load(RNA_Seq_Location)["RNA"]["Features"])
print(len(scRNAseq_gene_list))
input()
print(lf.RNA_Seq_Load(RNA_Seq_Location)["RNA"])
print("start")
list_of_dup = []
count = 0
#for i in scRNAseq_gene_list:
#	count+=1
#	for j in scRNAseq_gene_list:
#		if i.lower() == j.lower() and i != j:
#			if [i,j] not in list_of_dup and [j,i] not in list_of_dup:
#				print(count, len(list_of_dup)+1, len(scRNAseq_gene_list))
#				list_of_dup.append([i,j])
print(list_of_dup)
print(len(list_of_dup))
input()
rna_mat = lf.RNA_Seq_Load(RNA_Seq_Location)["RNA"]["data"].todense()

print(np.min(np.sum(rna_mat,axis = 0)))
print(np.max(np.sum(rna_mat,axis = 0)))
print(np.min(np.sum(rna_mat,axis = 1)))
print(np.min(np.shape(rna_mat)))
print(scRNAseq_gene_list)
input()
ind = 0
for i in scRNAseq_gene_list:
	print(f"{i} count is {np.sum(rna_mat[ind])}")
	ind+=1
input()

# If data is not in HGNC nomenclature the below should be true and the relevant "to_hgnc" json file should be
# called.
use_conversion_list = True

# Converts to HGNC using the json file
if use_conversion_list:
	with open(Conversion_json_location, "r") as conversion_file:
		conversion_dict = json.load(conversion_file)
	print(conversion_dict)
	for i in range(len(scRNAseq_gene_list)):
		try:
			scRNAseq_gene_list[i] = conversion_dict[scRNAseq_gene_list[i]]
		except KeyError:
			print(f"Could not find match in search for {scRNAseq_gene_list[i]}\n"
			      f"This seems to result from HGNC checker not liking the name")
			scRNAseq_gene_list[i] = ''

# Access model and obtain model genes
recon_flux_model = Flux_Balance_Model()

recon_flux_model.load_pkl_model(Minimal_Model_Location)
genes_in_model = recon_flux_model.get_used_genes()



# Gets scRNAseq gene list and data matrix
scRNAseq_mat = lf.RNA_Seq_Load(RNA_Seq_Location)["RNA"]["data"].todense()

# Takes the data for t = 0 hours, this is often not needed as scRNAseq data typically has no temporal dependence.
scRNAseq_mat_0h = scRNAseq_mat[:,:702]

# This shows genes with no hits were dropped
#print(np.min(np.sum(scRNAseq_mat,axis=1)))

# Creates matrix with each row corresponding to specific gene found in metabolic model
# For genes which exist in model but not in scRNAseq data it is assumed that gene had no hits
metabolic_scRNAseq_mat_0h = np.zeros((len(genes_in_model),np.shape(scRNAseq_mat_0h)[1]))
list_of_genes_in_model_not_in_data = []
for i in range(len(genes_in_model)):
	try:
		metabolic_scRNAseq_mat_0h[i] = scRNAseq_mat_0h[scRNAseq_gene_list.index(genes_in_model[i])]
	except ValueError:
		list_of_genes_in_model_not_in_data.append(genes_in_model[i])
		continue
print(f"There were {len(list_of_genes_in_model_not_in_data)} genes in the model which were not in the scRNAseq data\nThose genes are:")
print(list_of_genes_in_model_not_in_data)
print(metabolic_scRNAseq_mat_0h)

# Saves array with gene labels in HGNC format for each row
#name_array = np.reshape(np.array(genes_in_model),(-1,1))
#array_to_save = np.hstack((name_array,metabolic_scRNAseq_mat_0h))
rna_seq_df = pd.DataFrame(metabolic_scRNAseq_mat_0h,index = genes_in_model,columns = [f"cell_{i}" for i in range(np.shape(metabolic_scRNAseq_mat_0h)[1])])
rna_seq_df.to_csv(scRNA_Seq_Location)


