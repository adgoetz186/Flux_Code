import os
import time
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model


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

model_file_location = Path("Data/Models/pkl_models/recon1b_t_cell/exp_aligned/")
output_model_folder = Path("Data/Minimal_Model_Data/recon1b_t_cell/")
for filename in os.listdir(model_file_location):
	model_file_path = model_file_location/filename
	recon1_exp_align = Flux_Balance_Model()
	recon1_exp_align.load_pkl_model(model_file_path)
	for met in recon1_exp_align.rxn_dict["biomass_reaction"]["S"]:
		print(met,recon1_exp_align.get_met_mol_weight(met),recon1_exp_align.rxn_dict["biomass_reaction"]["S"][met])
	input()
	print(recon1_exp_align.rxn_dict["biomass_reaction"])
	print(recon1_exp_align.get_rxn_mass_change("biomass_reaction"))
	print(len(recon1_exp_align.rxn_dict))
	print(len(recon1_exp_align.met_dict))
	for rxn in recon1_exp_align.rxn_dict.keys():
		rxn_chem_change = recon1_exp_align.get_rxn_element_change(rxn)
		#print(recon1_exp_align.get_rxn_element_change("ARTPLM1"))
		if rxn_chem_change["Ra"] != 0 and "EX" not in rxn and "DM" not in rxn:
			print(rxn, rxn_chem_change)
			
	print(recon1_exp_align.rxn_dict["maintenance_ATP"])
	print(recon1_exp_align.rxn_dict["CK"])
	input()
	mm_list = {}
	mm_list["adp[c]"] = 427.201
	mm_list["ala-L[c]"] = 89.094
	mm_list["arg-L[c]"] = 174.204
	mm_list["asn-L[c]"] = 132.119
	mm_list["asp-L[c]"] = 133.103
	mm_list["atp[c]"] = 507.18
	mm_list["chsterol[c]"] = 1493.9
	mm_list["ctp[c]"] = 479.12
	mm_list["cys-L[c]"] = 121.16
	mm_list["datp[n]"] = 487.15
	mm_list["dctp[n]"] = 467.156
	mm_list['dgtp[n]'] = 467.156
	for i in recon1_exp_align.rxn_dict["biomass_reaction"]["S"].keys():
		print(i)
		input()