import copy
import os
from functools import partial
import re
import pandas as pd
import numpy as np
from pathlib import Path
import scipy.io as sio
import time
import json
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model

# _____ Setting the CWD to be Flux_Code BEGIN _____
# Flux_Code path here
path_to_FC = ""
if path_to_FC == "":
	try:
		# Obtains the location of the Cell_signaling_information folder if it is in the cwd parents
		path_to_FC = Path.cwd().parents[
			[Path.cwd().parents[i].parts[-1] for i in range(len(Path.cwd().parts) - 1)].index(
				"Flux_Code")]
	except ValueError:
		print("Flux_Code not found in cwd parents, trying sys.path")
		try:
			# Obtains the location of the Cell_signaling_information folder if it is in sys.path
			path_to_CSI = Path(sys.path[[Path(i).parts[-1] for i in sys.path].index("Flux_Code")])
		except ValueError:
			print("Flux_Code not found in sys.path "
			      "consult 'Errors with setting working directory' in README")
else:
	path_to_CSI = Path(path_to_FC)
os.chdir(path_to_FC)
# _____ Setting the CWD to be Flux_Code END _____

def create_entrez_to_express_dict(gene_names,gene_array,entrez_dict):
	mouse_gene_list = list(gene_names[:,0])
	
	list_unmatch = []
	list_m_ens_unmatch = []

	for i in entrez_dict.keys():
		if "mouse_ensemble" in entrez_dict[i].keys():
			if entrez_dict[i]["mouse_ensemble"] in mouse_gene_list:
				match_ind = mouse_gene_list.index(entrez_dict[i]["mouse_ensemble"])
				entrez_dict[i]["expression_val"] = gene_array[match_ind]/entrez_dict[i]["degen"]
				print(i, entrez_dict[i],1)
				#print(avg_exp[mouse_gene_list.index(entrez_dict[i]["mouse_ensemble"])], entrez_dict[i]["degen"],entrez_dict[i]["expression_val"])
			else:
				print(i, entrez_dict[i],2)
				entrez_dict[i]["expression_val"] = 0
				list_m_ens_unmatch.append([i,entrez_dict[i]["mouse_ensemble"]])
		else:
			entrez_dict[i]["expression_val"] = 0
			print(i, entrez_dict[i],3)
			list_unmatch.append(i)
	print(list_unmatch)
	print(list_m_ens_unmatch)
	return entrez_dict

def name_to_express(entrez_dict,matchobj):
	# converts from gene name to RNA expression
	entrez_name = matchobj.group(0)
	print(entrez_name)
	print(entrez_dict[entrez_name])
	print(str(entrez_dict[entrez_name]["expression_val"]))
	return str(entrez_dict[entrez_name]["expression_val"])

def strip_pare(matchobj):
	# strips parenthesis which encompass single value
	return matchobj.group(0)[1:-1]

def convert_and(matchobj):
	# converts all pure and sections to value
	ands = matchobj.group(0)
	ands = ands.replace(" ","")
	ands = np.average(np.array([float(i) for i in ands.split("and")]))
	return str(ands)

def convert_or_internal(matchobj):
	# converts all pure or statements enclosed in parenthesis to value
	ors = matchobj.group(0)[1:-1]
	ors = ors.replace(" ","")
	if ors.count("and") > 0:
		print(ors)
		print("ERROR")
		input()
	print(ors.split("or"))
	ors = np.sum(np.array([float(i) for i in ors.split("or")]))
	return str(ors)

def convert_or_external(matchobj):
	# converts final or statements to value, must be ran at very end
	ors = matchobj.group(0)
	ors = ors.replace(" ","")
	if ors.count("and") > 0:
		print(ors)
		print("ERROR")
		input()
	print(ors.split("or"))
	ors = np.sum(np.array([float(i) for i in ors.split("or")]))
	return str(ors)

def convert_gr_Rule_to_RAS(gr_rule,entrez_dict):
	# converts name to expression
	print(gr_rule)
	name_to_real_express = partial(name_to_express, entrez_dict)
	gr_rule = re.sub("\d+\.\d+", name_to_real_express, gr_rule)
	print(gr_rule)
	int_or_altered = True
	while int_or_altered:
		gr_rule_cp_2 = copy.deepcopy(gr_rule)
		and_altered = True
		while and_altered:
			gr_rule_cp = copy.deepcopy(gr_rule)
			gr_rule = re.sub("\(\d+\.?\d*e?-?\d*\)",strip_pare,gr_rule)
			gr_rule = re.sub("(\d+\.?\d*e?-?\d* ?(and) ?)+\d+\.?\d*e?-?\d*", convert_and, gr_rule)
			if gr_rule_cp == gr_rule:
				and_altered = False
		gr_rule = re.sub("\(\d+\.?\d*e?-?\d*\)", strip_pare, gr_rule)
		gr_rule = re.sub("\((\d+\.?\d*e?-?\d* ?(or) ?)+\d+\.?\d*e?-?\d*\)", convert_or_internal, gr_rule)
		if gr_rule_cp_2 == gr_rule:
			int_or_altered = False
	gr_rule = re.sub("(\d+\.?\d*e?-?\d* ?(or) ?)+\d+\.?\d*e?-?\d*",convert_or_external,gr_rule)
	if gr_rule == "":
		print(np.nan)
		return np.nan
	else:
		print(float(gr_rule))
		return float(gr_rule)
		
		
		
Raw_gene_file_location = Path("Data/scRNA/Raw/B61-D2-GEX-B5")

entrez_dict_location = Path("Data/Models/Gene_Data_For_Models/entrez_dict/RECON1.json")



with open(entrez_dict_location, ) as entrez:
	entrez_dict = json.load(entrez)

features = np.loadtxt(Raw_gene_file_location/"features.tsv.gz",delimiter="\t",dtype = str)
gene_mat = np.array(sio.mmread(Raw_gene_file_location / "matrix.mtx.gz").todense())
print(gene_mat)
print(list(features[:,0]).index("ENSMUSG00000035472"))
print(np.average(gene_mat[list(features[:,0]).index("ENSMUSG00000035472"),:]))
print(np.shape(gene_mat))
print(np.shape(features))
input()

flux_model_dict = {}
for filename in os.listdir(HR_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(HR_model_location / filename)
	flux_model_dict[filename.split("_")[0]] = recon_flux_model
gene_array = np.average(gene_mat,axis = 1)


barcodes = np.loadtxt(Raw_gene_file_location/"barcodes.tsv.gz",delimiter="\t",dtype = str)
features = np.loadtxt(Raw_gene_file_location/"features.tsv.gz",delimiter="\t",dtype = str)


