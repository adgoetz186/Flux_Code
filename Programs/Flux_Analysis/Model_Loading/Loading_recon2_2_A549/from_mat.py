from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
import numpy as np
import sys

print(os.listdir())
# Loads the recon model from the specified mat file and saves it as a dictionary of pickled objects.
# The new saved model can be loaded with load_pkl_model()

# File locations should end with "/"
# File names should not end with file specifier (.txt, .py, etc)
mat_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/mat_models/recon2_2"


output_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_A549/raw"

# This program assumes certain standard names for various parts of the model:
#   standard names: Stochiometric matrix = "S", lower bounds = "lb", upper bounds = "ub",
#   rhs of Sv = b constraint = "b", dictionary of pinched reactions = "pinched_reactions",
#   list of reaction names = "reaction_names", gene reaction rule list = "grRules",
#   list of gene names = "genes", list of metabolite names = "metabolite_names"
#   list of metabolite compositions = "metabolite_comp"
# If a loaded model uses different names, the loader wont know how to handle the alternative names
# This dictionary serves as a translator, the keys should be the standard names given above and the values should be
# the corresponding keys in the mat file. Running the program with no model_header_dict will provide a list of
# mat file keys.
model_comp={"S":"S", "lb" : "lb", "ub" : "ub", "b" : "b","pinched_reactions" : None, "metabolite_names" : "mets", "reaction_names" : "rxns", "grRules" : "grRules", "genes" : "genes", "metabolite_comp" : "metFormulas"}
model_name_to_save = "Recon2_2"

recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_mat_model(mat_file_location,stored_model_name=model_name_to_save,model_comp=model_comp)

for i in range(len(recon_flux_model.lb)):
    if recon_flux_model.ub[i] == np.inf:
        recon_flux_model.update_reaction_bounds(i,"keep",1000)
    if recon_flux_model.lb[i] == -np.inf:
        recon_flux_model.update_reaction_bounds(i,-1000,"keep")

recon_flux_model.save_model_as_pkl(output_file_location)