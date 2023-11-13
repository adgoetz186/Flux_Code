import os
import pandas as pd
import numpy as np
import scipy.stats as st
import sys
import shutil


ms_data = pd.read_csv("/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_t_cell/Data_For_All_Mice/Met_Flux_Data.csv")
cell_size_path = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_t_cell/Data_For_All_Mice/cell_size.txt"
Allowed_Exchange_Reactions_path = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_t_cell/Data_For_All_Mice/Allowed_Exchange_Reactions"
output_folder_dir = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_t_cell/Individual_Mice/"
mouse_number_filename_dict = {1:"B6-1",2:"B6-2",3:"B6-3",4:"B6-4",5:"TC-5",6:"TC-6",7:"TC-7",8:"TC-8"}
exp_metabolite_exchange_rxn_dict = {}
exp_metabolite_exchange_rxn_dict["Lactate"] = "EX_lac_L_b"
exp_metabolite_exchange_rxn_dict["Glucose"] = "EX_glc_D_b"
exp_metabolite_exchange_rxn_dict["Glutamine"] = "EX_gln_L_b"
exp_metabolite_exchange_rxn_dict["Serine"] = "EX_ser_L_b"
exp_metabolite_exchange_rxn_dict["Asparagine"] = "EX_asn_L_b"
exp_metabolite_exchange_rxn_dict["Leucine"] = "EX_leu_L_b"
exp_metabolite_exchange_rxn_dict["Lysine"] = "EX_lys_L_b"
exp_metabolite_exchange_rxn_dict["Isoleucine"] = "EX_ile_L_b"
exp_metabolite_exchange_rxn_dict["Cystine"] = "EX_cys_L_b"
exp_metabolite_exchange_rxn_dict["Arginine"] = "EX_arg_L_b"
exp_metabolite_exchange_rxn_dict["Threonine"] = "EX_thr_L_b"
exp_metabolite_exchange_rxn_dict["Phenylalanine"] = "EX_phe_L_b"
exp_metabolite_exchange_rxn_dict["Methionine"] = "EX_met_L_b"
exp_metabolite_exchange_rxn_dict["Tyrosine"] = "EX_tyr_L_b"
exp_metabolite_exchange_rxn_dict["Histidine"] = "EX_his_L_b"
exp_metabolite_exchange_rxn_dict["Tryptophan"] = "EX_trp_L_b"
exp_metabolite_exchange_rxn_dict["Proline"] = "EX_pro_L_b"
exp_metabolite_exchange_rxn_dict["Aspartic_Acid"] = "EX_asp_L_b"
exp_metabolite_exchange_rxn_dict["Glycine"] = "EX_gly_b"
exp_metabolite_exchange_rxn_dict["Alanine"] = "EX_ala_L_b"
exp_metabolite_exchange_rxn_dict["Glutamic_Acid"] = "EX_glu_L_b"
#Optional
exp_metabolite_exchange_rxn_dict["Valine"] = "EX_val_L_b"
exp_metabolite_exchange_rxn_dict["L-ornithine"] = "EX_orn_b"


mice_numbers = list(set(ms_data["Mouse Number"]))
mice_numbers.sort()

for mouse_number in mice_numbers:
	if not os.path.isdir(output_folder_dir+mouse_number_filename_dict[mouse_number]):
		os.mkdir(output_folder_dir+mouse_number_filename_dict[mouse_number])
	shutil.copy(cell_size_path,output_folder_dir+mouse_number_filename_dict[mouse_number]+"/cell_size.txt")
	shutil.copy(Allowed_Exchange_Reactions_path, output_folder_dir + mouse_number_filename_dict[mouse_number]+"/Allowed_Exchange_Reactions")

for mouse_number in mice_numbers:
	reaction_list = []
	value_list = []
	unit_list = []
	mouse_data = ms_data.loc[ms_data["Mouse Number"] == mouse_number]
	for metabolite in exp_metabolite_exchange_rxn_dict.keys():
		reaction_list.append(exp_metabolite_exchange_rxn_dict[metabolite])
		value_list.append(mouse_data[metabolite].to_numpy()[0])
		unit_list.append("fmol/(cell*hr)")
	output_dataframe = pd.DataFrame({"reaction":reaction_list,"value":value_list,"units":unit_list})
	output_dataframe.to_csv(output_folder_dir+mouse_number_filename_dict[mouse_number]+"/Measured_Rates.csv",index=False)
