import os
import pandas as pd
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model


model_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_A549/exp_aligned_alpha_given"

dataframe_location_folder = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Minimal_Model_Data/recon2_2_A549/"

minimal_model_output_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_A549/min_model"

dataframe_list = []
for file in os.listdir(dataframe_location_folder):
	dataframe_list.append(pd.read_csv(dataframe_location_folder+file,index_col=0,header=0))

recon2_2_exp_align = Flux_Balance_Model()
recon2_2_exp_align.load_pkl_model(model_file_location)

flux_matrix = recon2_2_exp_align.build_essential_flux_matrix_from_dataframes(dataframe_list)
print(len(recon2_2_exp_align.reaction_names))
recon2_2_exp_align.generate_essential_flux_model_from_essential_flux_matrix(flux_matrix)
recon2_2_exp_align.save_model_as_pkl(minimal_model_output_location)
print(len(recon2_2_exp_align.reaction_names))