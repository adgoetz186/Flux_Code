import os
import pandas as pd
import numpy as np
import scipy.stats as st
import sys
from pathlib import Path
import shutil


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

ms_data = pd.read_csv(Path("Data/experimental_alignment_data/recon_1b_t_cells/Data_For_All_Mice/Expt_Met_Flux_Data.csv"))
bm_data = pd.read_csv(Path("Data/experimental_alignment_data/recon_1b_t_cells/Data_For_All_Mice/Biomass_Data.csv"))
cell_size_path = Path("Data/experimental_alignment_data/recon_1b_t_cells/Data_For_All_Mice/cell_size.txt")
Allowed_Exchange_Reactions_path = Path("Data/experimental_alignment_data/recon_1b_t_cells/Data_For_All_Mice/Allowed_Exchange_Reactions")
output_folder_dir = Path("Data/experimental_alignment_data/recon_1b_t_cells/Individual_Mice/")
mouse_number_filename_dict = {1:"B6-1",2:"B6-2",3:"B6-3",4:"B6-4",5:"TC-5",6:"TC-6",7:"TC-7",8:"TC-8"}
exp_metabolite_exchange_rxn_dict = {}
exp_metabolite_exchange_rxn_dict["Lactate"] = "EX_lac_L(e)"
exp_metabolite_exchange_rxn_dict["Glucose"] = "EX_glc(e)"
exp_metabolite_exchange_rxn_dict["Glutamine"] = "EX_gln_L(e)"
exp_metabolite_exchange_rxn_dict["Serine"] = "EX_ser_L(e)"
exp_metabolite_exchange_rxn_dict["Asparagine"] = "EX_asn_L(e)"
exp_metabolite_exchange_rxn_dict["Leucine"] = "EX_leu_L(e)"
exp_metabolite_exchange_rxn_dict["Lysine"] = "EX_lys_L(e)"
exp_metabolite_exchange_rxn_dict["Isoleucine"] = "EX_ile_L(e)"
exp_metabolite_exchange_rxn_dict["Cystine"] = "EX_cys_L(e)"
exp_metabolite_exchange_rxn_dict["Arginine"] = "EX_arg_L(e)"
exp_metabolite_exchange_rxn_dict["Threonine"] = "EX_thr_L(e)"
exp_metabolite_exchange_rxn_dict["Phenylalanine"] = "EX_phe_L(e)"
exp_metabolite_exchange_rxn_dict["Methionine"] = "EX_met_L(e)"
exp_metabolite_exchange_rxn_dict["Tyrosine"] = "EX_tyr_L(e)"
exp_metabolite_exchange_rxn_dict["Histidine"] = "EX_his_L(e)"
exp_metabolite_exchange_rxn_dict["Tryptophan"] = "EX_trp_L(e)"
exp_metabolite_exchange_rxn_dict["Proline"] = "EX_pro_L(e)"
exp_metabolite_exchange_rxn_dict["Aspartic_Acid"] = "EX_asp_L(e)"
exp_metabolite_exchange_rxn_dict["Glycine"] = "EX_gly(e)"
exp_metabolite_exchange_rxn_dict["Alanine"] = "EX_ala_L(e)"
exp_metabolite_exchange_rxn_dict["Glutamic_Acid"] = "EX_glu_L(e)"
#Optional
exp_metabolite_exchange_rxn_dict["Valine"] = "EX_val_L(e)"
exp_metabolite_exchange_rxn_dict["L-ornithine"] = "EX_orn(e)"


mice_numbers = list(set(ms_data["Mouse Number"]))
mice_numbers.sort()

for mouse_number in mice_numbers:
	if not os.path.isdir(output_folder_dir / mouse_number_filename_dict[mouse_number]):
		os.mkdir(output_folder_dir / mouse_number_filename_dict[mouse_number])
	shutil.copy(cell_size_path,output_folder_dir / mouse_number_filename_dict[mouse_number] / "cell_size.txt")
	shutil.copy(Allowed_Exchange_Reactions_path, output_folder_dir / mouse_number_filename_dict[mouse_number] / "Allowed_Exchange_Reactions")

for mouse_number in mice_numbers:
	reaction_list = []
	value_list = []
	unit_list = []
	mouse_data = ms_data.loc[ms_data["Mouse Number"] == mouse_number]
	for metabolite in exp_metabolite_exchange_rxn_dict.keys():
		reaction_list.append(exp_metabolite_exchange_rxn_dict[metabolite])
		value_list.append(mouse_data[metabolite].to_numpy()[0])
		unit_list.append("fmol/(cell*hr)")
	mouse_data_bm = bm_data.loc[bm_data["Mouse Number"] == mouse_number]
	print(mouse_data_bm)
	print(mouse_data_bm["Biomass"])
	print(mouse_data_bm["Biomass"].to_numpy()[0])
	output_dataframe_bm = pd.DataFrame({"reaction": ["biomass_reaction"], "value": [mouse_data_bm["Biomass"].to_numpy()[0]], "units": ["mmol/(gDw hr)"]})
	
	output_dataframe = pd.DataFrame({"reaction":reaction_list,"value":value_list,"units":unit_list})
	output_dataframe.to_csv(output_folder_dir / mouse_number_filename_dict[mouse_number] / "Measured_Rates.csv",
	                        index=False)
	
	output_dataframe_bm.to_csv(output_folder_dir / mouse_number_filename_dict[mouse_number] / "Biomass_Data.csv",
	                        index=False)


	

