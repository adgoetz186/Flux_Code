import json
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf
import os
import mygene

os.chdir("..")
mg = mygene.MyGeneInfo()
print(os.getcwd())
scRNAseq_gene_list = lf.RNA_Seq_Load("Data/Input/scRNA-seq-data/A549_sciCAR_data.mat")["RNA"]["Features"]

with open("Data/Input/Symbol_Handling/Model_scRNA_mismatch_symbols/recon2_2", "r") as gene_file:
	mismatch_gene_list = eval(gene_file.readlines()[0])
	
with open("Data/Input/Symbol_Handling/Model_Gene_Symbols/recon2_2", "r") as gene_file:
	model_gene_list = eval(gene_file.readlines()[0])

with open("Data/Input/Symbol_Handling/Model_scRNA_mismatch_database_dictionary/recon2_2.json","r") as database:
	data_dict = json.load(database)

list_of_gene_info = data_dict["DocumentSummarySet"]["DocumentSummary"]

relevant_alias_converter = {}
missing_gene_list = []
count = 0
for i in mismatch_gene_list:
	gene_match = False
	for j in range(len(list_of_gene_info)):
		if i == list_of_gene_info[j]["Name"]:
			alias_list = list_of_gene_info[j]["OtherAliases"].replace(" ", "").split(",")
			print(alias_list)
			for k in alias_list:
				if k in scRNAseq_gene_list:
					gene_match = True
					relevant_alias_converter[i] = k
	if not gene_match:
		missing_gene_list.append(i)
print(relevant_alias_converter)
print(len(relevant_alias_converter))
print(relevant_alias_converter["B4GAT1"])
print(missing_gene_list)
input()

with open("Data/Input/Symbol_Handling/Model_scRNA_convertion_dict/recon2_2.json","w") as relevant_alias_converter_file:
	json.dump(relevant_alias_converter,relevant_alias_converter_file)

for i in range(len(missing_gene_list)):
	print(i+1, missing_gene_list[i])

print(len(missing_gene_list))

matching_gene_list_RNAseq = []
for i in model_gene_list:
	if i in scRNAseq_gene_list:
		matching_gene_list_RNAseq.append(i)
	else:
		try:
			if relevant_alias_converter[i] in scRNAseq_gene_list:
				matching_gene_list_RNAseq.append(relevant_alias_converter[i])
		except KeyError:
			print(f"{i} returned no matches")
print(matching_gene_list_RNAseq)
print(len(matching_gene_list_RNAseq))
input()
matching_gene_list_model = []
for i in model_gene_list:
	if i in scRNAseq_gene_list:
		matching_gene_list_model.append(i)
	else:
		try:
			if relevant_alias_converter[i] in scRNAseq_gene_list:
				matching_gene_list_model.append(i)
		except KeyError:
			print(f"{i} returned no matches")
#print(matching_gene_list_model)
print(len(matching_gene_list_model))
			