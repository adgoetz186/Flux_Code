import json
import os
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

# Performs initial alignment of t-cell sizes with experimental data
# Each conditions error is found uniquely here allowing for one or more conditions to be grouped to find more than
# one cell size. At the time of writting this we only allow for one cell size.

# _____ Setting the CWD to be Flux_Code BEGIN _____
# Flux_Code path here
path_to_FC = ""
if path_to_FC == "":
	try:
		# Obtains the location of the Cell_signaling_information folder if it is in the cwd parents
		path_to_FC = Path.cwd().parents[[Path.cwd().parents[i].parts[-1] for i in range(len(Path.cwd().parts) - 1)].index(
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

scRNA_path = Path("Data/scRNA/A549_Raw/A549_sciCAR_data.mat")
bulk_output_path = Path("Data/scRNA/A549_Bulk/matrix.npy")
entrez_symbol_array_path = Path("Data/scRNA/A549_Bulk/entrez_symbol_features.npy")
#entrez_features_output_path = Path("Data/scRNA/A549_Bulk/entrez_features.npy")
human_gene_info_path = Path("Data/scRNA/full_gene_convert_info/hgnc_complete_set.json")

# bigg recon1 is used for its more complete gene descriptions
recon1_bigg_path = Path("Data/Models/json_models/Bigg_format/RECON1.json")

with open(human_gene_info_path,'rb')  as jfile:
	t = json.load(jfile)

with open(recon1_bigg_path,'rb')  as jfile:
	bigg_recon1 = json.load(jfile)







scRNA_mat = sio.loadmat(scRNA_path)
#print(scRNA_mat.keys())
print(scRNA_mat['RNA'])
avg_expression = np.average(scRNA_mat['RNA'][0][0][0].todense(),axis=1)
avg_expression = np.reshape(avg_expression,(1,-1))
avg_expression = np.asarray(avg_expression)[0]

print(len([i[0][0] for i in scRNA_mat['RNA'][0][0][1]]))
print(len([i[0][0] for i in scRNA_mat['RNA'][0][0][1]]))
feature_array = np.array([i[0][0] for i in scRNA_mat['RNA'][0][0][1]])
print(np.shape(scRNA_mat['RNA'][0][0][0].todense()))
print(len(list(scRNA_mat['RNA'][0][0][2][0])))
print(feature_array)
print(len(list(feature_array)))
print(np.min(avg_expression))
print(np.shape(avg_expression))
print('oksa')
input()
upper_case_fa = [i.upper() for i in list(feature_array)]
print(upper_case_fa)


print(t['response'].keys())
symbol_dict = {i['symbol']:i for i in t['response']['docs']}
print(symbol_dict)
key_list = []
for i in t['response']['docs']:
	key_list += list(i.keys())
key_list = list(set(key_list))
key_list.sort()
print(key_list)

use_prev = True

print(bigg_recon1['genes'][0])
# make 3 dicts for each level of search

for i in feature_array:
	if "TPI1P2" in i:
		print(13, i)
		input()
man_dict = {}
man_dict['284208_AT1'] = 'B3GNTL1'
man_dict['2351_AT1'] = 'FOLR2L'
man_dict['57733_AT1'] = 'GBA3'
man_dict['6537_AT1'] = 'SLC6A10P'
man_dict['286297_AT1'] = 'LOC286297'
man_dict['2713_AT1'] = 'GK3'
man_dict['2679_AT1'] = 'GGT3P'
man_dict['2974_AT1'] = 'GUCY1B2'
man_dict['26237_AT1'] = 'ENO1B'
man_dict['7359_AT1'] = 'UGP2'
man_dict['22305_AT1'] = 'Vmn2r37'
man_dict['22305_AT2'] = 'Vmn2r37'
man_dict['2655_AT1'] = 'GGCT'
man_dict['26062_AT1'] = 'HYAL6P'
man_dict['4844_AT1'] = 'NOS2P2'
man_dict['4845_AT1'] = 'NOS2P1'
man_dict['80068_AT1'] = 'MTHFD2L'
man_dict['285216_AT1'] = 'MTHFD2P1'
man_dict['348477_AT1'] = 'MIA3'
man_dict['5240_AT1'] = 'PGP'
man_dict['8781_AT1'] = 'PSPHP1'
man_dict['2544_AT1'] = 'SLC37A4'
man_dict['9954_AT1'] = 'HS3ST3A2'
man_dict['9952_AT1'] = 'HS3ST3B2'
man_dict['286016_AT1'] = 'TPI1P2'

symbol_dict = {}
alias_dict = {}
curate_sym_dict = {}

print('done')
input()
for human_gene_info in bigg_recon1['genes']:

	if 'ncbigene' in list(human_gene_info['annotation'].keys()):

		if human_gene_info['name'] in feature_array:
			symbol_dict[human_gene_info['name']] = human_gene_info['annotation']['ncbigene'][0]
		elif 'refseq_synonym' in human_gene_info['annotation'].keys():
			found_it = False
			for alias in human_gene_info['annotation']['refseq_synonym']:
				if alias in feature_array:
					alias_dict[alias] = human_gene_info['annotation']['ncbigene'][0]
					found_it = True
	else:
		if human_gene_info['id'] in man_dict.keys():
			curate_sym_dict[man_dict[human_gene_info['id']]] = human_gene_info['id'].split("_")[0]
print(curate_sym_dict)
print('done')
input()
#		if 'prev_symbol' in human_gi_keys:
#			for prev in human_gene_info['prev_symbol']:
#				prev_dict[prev] = human_gene_info['entrez_id']
#print(symbol_dict)
#input()
#print(alias_dict)
#print(prev_dict)
#input()
symbols = symbol_dict.keys()
alias = alias_dict.keys()
#prevs = prev_dict.keys()

use_prev = True
# A much faster
entrez_feature_list = []
prev = 0.001
#feature_array = ["PGP","SLC37A4","MTHFD2P1","MTHFD2L"]
for data_name_feature in feature_array:
	#print(len(entrez_feature_list)/len(feature_array))
	if len(entrez_feature_list)/len(feature_array) > prev:
		print(len(entrez_feature_list)/len(feature_array))
		prev+=0.001
	#print(data_name_feature)
	#print(len(entrez_feature_list))
	prev_list_len = len(entrez_feature_list)
	count = 0
	if data_name_feature in symbol_dict:
		entrez_feature_list.append(symbol_dict[data_name_feature])
		print(1, data_name_feature,symbol_dict[data_name_feature])
	elif data_name_feature in alias_dict:
		if list(alias_dict.keys()).count(data_name_feature) == 1:
			entrez_feature_list.append(alias_dict[data_name_feature])
			print(2, data_name_feature, alias_dict[data_name_feature],data_name_feature)
		else:
			print("AMBIGUOUS")
			input()
	elif data_name_feature in curate_sym_dict:
		entrez_feature_list.append(curate_sym_dict[data_name_feature])
		print(entrez_feature_list[-1],data_name_feature)
		print('ok')
		input()

	#elif data_name_feature in prev_dict:
	#	entrez_feature_list.append(prev_dict[data_name_feature])
	else:
		entrez_feature_list.append('NO_ENTREZ')

	if prev_list_len+1 != len(entrez_feature_list):
		print(entrez_feature_list[-1], data_name_feature,count)
		raise Exception(f"Something went wrong, the size of the list has grown more than expected\nlist_0 = {prev_list_len}, list_1 = {len(entrez_feature_list)}")
print('done')

input()
entrez_symbol_array = np.transpose(np.vstack((np.array(entrez_feature_list),feature_array)))
print(entrez_symbol_array)

rxns = ["SLC18A3", "AGXT", "ARG1","CH25H" ,"CDO1" ,"CYP4A11","CYP4A11", "CYP4A11", "ALDOB", "GGCT", "GCTG", "HSD3B1", "HSD3B1", "PDXP", "SLC5A3", "PLA2G3", "MTHFD2L","FLJ13105","LOC285216","MTHFD2P1", "SULT1A3", "CYP11B1", "CYP11B1","CYP11B2", "PDXP", "PCK1", "PGP" , "G6PT2","SLC37A4", "PDXP", "AGXT", "TH", "UGT1A8"]
for rxn in rxns:
	print(1, rxn, rxn in upper_case_fa)
	if rxn in upper_case_fa:
		print(2, rxn, rxn in feature_array)
		if rxn in feature_array:
			print(3, rxn, avg_expression[list(feature_array).index(rxn)])
			print(feature_array[list(feature_array).index(rxn)])
			print(entrez_symbol_array[list(feature_array).index(rxn)])
#mouse
input()

# it looks like the above is working too well, a gene will be mapped to its newest entrez, and sometimes that is different from the old one. Maybe if it can be shown to be rare, just do it manually
np.save(bulk_output_path,avg_expression)
np.save(entrez_symbol_array_path,entrez_symbol_array)
#print([i[0] for i in scRNA_mat['RNA'][0][0][2][0]])