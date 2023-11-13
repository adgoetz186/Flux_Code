import copy

import numpy as np
import matplotlib.pyplot as plt
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import os
import pickle
import pandas as pd
import scipy.io as sio
from pathlib import Path
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

input_path = Path("Data/HR/HR_bi_directional_points/recon_1b_t_cells_gene_Day_2")
HR_model_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/HR_ready_Day_2")
subsystem_convert = Path("Data/subsystem_converter/subsystem_convert.pkl")





# loads flux models

recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(HR_model_location / os.listdir(HR_model_location)[0])
recon_flux_model.convert_model_to_bidirectional_flux()

# loads reaction names
lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
names = rxn_list
names.append("total_cost")
print(names)

start = time.time()
with open(subsystem_convert, "rb") as readfile:
	subsystem_convert = pickle.load(readfile)
print(subsystem_convert)
subsystem_dict = {}
for i in subsystem_convert.keys():
	if i in names:
		pathway = ""
		if not isinstance(subsystem_convert[i], str):
			pathway = "None"
		else:
			pathway = subsystem_convert[i]
		if pathway in subsystem_dict.keys():
			subsystem_dict[pathway].append(i)
		else:
			subsystem_dict[pathway] = [i]
	else:
		print("no")

print(subsystem_dict.keys())
print(subsystem_dict['Pentose Phosphate Pathway'])
print(subsystem_dict['Glycolysis/Gluconeogenesis'])
print(subsystem_dict['Citric Acid Cycle'])
glycolosis = ["GLCt1r","HEX1","PGI","G6PDH2r","TALA",'TKT1','TKT2',"PFK","PGK","FBA","GAPD","ENO","PYK","PGM","LDH_L","PYRt2m","TPI","L_LACt2r","L_LACt4r","EX_glc(e)","EX_lac_L(e)","ME2"]
mthgxl_c_rxns = ["AACTOOR","ALR2","LALDO2","MGSA","MGSA2"]
h2o2_c_rxns = ["H2O2syn","GTHP","PRDX"]
TCA = ['ACITL', 'ACONT', 'ACONTm', 'AKGDm', 'CITL', 'CSm', 'FUM', 'FUMm', 'ICDHxm', 'ICDHy', 'ICDHyrm', 'MDH', 'MDHm', 'SUCD1m', 'SUCOAS1m', 'SUCOASm',"EX_o2(e)","AKGMALtm","ME2","O2t","H2O2syn"]
lipid = ["EX_glc(e)","HEX1","PGI","PFK","FBA","TPI","GLCt1r","G3PD1","G3PD2m","CITtam","CSm","ACONT","ACCOAC","ACACT1r"]
ppp = ['DRBK', 'DRPA', 'G6PDH1rer', 'G6PDH2r', 'G6PDH2rer', 'GND', 'PGL', 'PPM', 'PRPPS', 'RPE', 'RPI', 'TALA', 'TKT1', 'TKT2']



B6_list = []
TC_list = []

dtf_dict = {}
# generates dtf
for file_name in os.listdir(input_path):
	flux_samples = np.load(input_path / file_name)
	print(np.shape(flux_samples))
	if "B6" in file_name:
		B6_list.append(np.average(flux_samples,axis=0))
	elif "TC" in file_name:
		TC_list.append(np.average(flux_samples,axis=0))


B6 = np.vstack(B6_list)
TC = np.vstack(TC_list)

dtf_dict["B6_mice"] = pd.DataFrame(B6,columns=names)
dtf_dict["TC_mice"] = pd.DataFrame(TC,columns=names)
print(dtf_dict["B6_mice"].loc[:,"H2O2syn"])
print(dtf_dict["TC_mice"].loc[:,"H2O2syn"])

print(dtf_dict["B6_mice"].loc[:,"ALR"])
print(dtf_dict["TC_mice"].loc[:,"ALR"])
print(np.shape(B6))
print(np.shape(TC))

# Test statistic is B6 is larger
# p value is 2 tail, use significance of p value then use test stat for sign
B6_adj = copy.deepcopy(B6)[:,:-1]
TC_adj = copy.deepcopy(TC)[:,:-1]

print(B6_adj)
for row in B6_adj:
	row/=np.sum(row)
for row in TC_adj:
	row/=np.sum(row)
t_test = st.ttest_ind(B6_adj,TC_adj,equal_var=False)
print(t_test)
#t_test = st.ttest_ind(B6,TC,equal_var=False)
print(np.argsort(t_test[1]))
print(t_test[1][np.argsort(t_test[1])])
print((np.arange(np.size(t_test[1]))+1)/np.size(t_test[1])*1)
print(t_test[1][np.argsort(t_test[1])] <= (np.arange(np.size(t_test[1]))+1)/np.size(t_test[1])*.3)
hp = np.argwhere(t_test[1][np.argsort(t_test[1])] <= (np.arange(np.size(t_test[1]))+1)/np.size(t_test[1])*.6)[-1]
#print(np.argsort(t_test[1])[:hp[0]+1])

plt.scatter(np.arange(np.size(t_test[1]))+1,t_test[1][np.argsort(t_test[1])])
plt.scatter(np.arange(np.size(np.argsort(t_test[1])[:hp[0]+1]))+1,t_test[1][np.argsort(t_test[1])[:hp[0]+1]])
plt.plot(np.arange(np.size(t_test[1]))+1,(np.arange(np.size(t_test[1]))+1)/np.size(t_test[1])*.6)
plt.show()
test_stat = list(t_test[0][np.argsort(t_test[1])[:hp[0]+1]])
rxn_name = list(np.array(names)[np.argsort(t_test[1])[:hp[0]+1]])
rxn_t_dict_B6 = dict(zip(rxn_name,test_stat))

for rxn in rxn_t_dict_B6.keys():
	if rxn_t_dict_B6[rxn] < 0:
		rxn_t_dict_B6[rxn] = "red"
	else:
		rxn_t_dict_B6[rxn] = "blue"
rxn_t_dict_TC = dict(zip(rxn_name,test_stat))
for rxn in rxn_t_dict_TC.keys():
	if rxn_t_dict_TC[rxn] < 0:
		rxn_t_dict_TC[rxn] = "blue"
	else:
		rxn_t_dict_TC[rxn] = "red"
print(rxn_t_dict_B6)
#print(rxn_t_dict_B6["PYRt2m"])
plt.show()
rxn_list = glycolosis

#print(recon_flux_model.met_dict.keys())
recon_flux_model.metabolite_info("mthgxl[c]",dtf_dict["B6_mice"])
recon_flux_model.metabolite_info("mthgxl[c]",dtf_dict["TC_mice"])
edge = []
weight = []
for rxn in rxn_list:
	ph_edge,ph_weight = recon_flux_model.element_edge_vert_generator(rxn,rxn_dtf = dtf_dict["B6_mice"],element="C",break_element_list=["atp[c]","adp[c]","nad[c]","nadh[c]","h[c]","nadph[c]","nadp[c]","coa[m]"])
	edge += ph_edge
	weight += ph_weight
recon_flux_model.generate_rxn_met_graph(edge,weight,color_edge_of_vert=rxn_t_dict_B6,vert_leak_tol=1e-2)

edge = []
weight = []
for rxn in rxn_list:
	ph_edge,ph_weight = recon_flux_model.element_edge_vert_generator(rxn,rxn_dtf = dtf_dict["TC_mice"],element="C",break_element_list=["atp[c]","adp[c]","nad[c]","nadh[c]","h[c]","nadph[c]","nadp[c]","coa[m]"])
	edge += ph_edge
	weight += ph_weight
recon_flux_model.generate_rxn_met_graph(edge,weight,color_edge_of_vert=rxn_t_dict_TC,vert_leak_tol=1e-2)


