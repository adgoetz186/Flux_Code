import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf
import json
import os

os.chdir("..")

scRNAseq_gene_list = lf.RNA_Seq_Load("Data/Input/scRNA_seq-data/A549_sciCAR_data.mat")["RNA"]["Features"]

with open("Data/Input/Symbol_Handling/Model_scRNA_convertion_dict/recon2_2.json","r") as converter_file:
	converter = json.load(converter_file)
	
with open("Data/Input/Symbol_Handling/Model_Gene_Symbols/recon2_2", "r") as gene_file:
	model_gene_list = eval(gene_file.readlines()[0])
	
model_gene_list_matching = []
model_gene_list_nonmatching = []
for i in model_gene_list:
	if i in scRNAseq_gene_list:
		model_gene_list_matching.append(i)
	else:
		try:
			if converter[i] in scRNAseq_gene_list:
				model_gene_list_matching.append(i)
		except KeyError:
			model_gene_list_nonmatching.append(i)
RNA_gene_list_matching = []
for i in model_gene_list:
	if i in scRNAseq_gene_list:
		RNA_gene_list_matching.append(i)
	else:
		try:
			if converter[i] in scRNAseq_gene_list:
				RNA_gene_list_matching.append(converter[i])
		except KeyError:
			continue
			
with open("Data/Input/Symbol_Handling/Matching_Gene_List/recon2_2", "w") as write_file:
	write_file.write(str(model_gene_list_matching) + "\n")
	write_file.write(str(RNA_gene_list_matching) + "\n")
	write_file.write(str(model_gene_list_nonmatching))
print(model_gene_list_matching)
print(len(model_gene_list_matching))
print(model_gene_list_nonmatching)
print(len(model_gene_list_nonmatching))
print(RNA_gene_list_matching)
print(len(RNA_gene_list_matching))
