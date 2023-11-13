from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
import numpy as np


os.chdir("../../../../")
recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_mat_model("Flux_Code/Data/Input/MAT_Models/recon2_2.mat",stored_model_name="Recon2_2",model_comp={"S":"S", "lb" : "lb", "ub" : "ub", "b" : "b","pinched_reactions" : None, "metabolite_names" : "mets", "reaction_names" : "rxnNames", "grRules" : "grRules", "genes" : "genes", "metabolite_comp" : "metFormulas"})

# checking for metabolites
metabolites = recon_flux_model.metabolite_names
i1 = metabolites.index("gln_L_c")
i2 = metabolites.index("glu_L_c")
i3 = metabolites.index("h2o_c")
i4 = metabolites.index("nh4_c")
i5 = metabolites.index("ser_L_c")
i6 = metabolites.index("ser_L_m")

GLUN_S_col = np.zeros((np.shape(recon_flux_model.S)[0],1))
GLUN_S_col[i1] = -1
GLUN_S_col[i2] = 1
GLUN_S_col[i3] = -1
GLUN_S_col[i4] = 1

SERtm_S_col = np.zeros((np.shape(recon_flux_model.S)[0],1))
SERtm_S_col[i5] = -1
SERtm_S_col[i6] = 1

#Potentially remove HGNC:4331
recon_flux_model.add_reaction(GLUN_S_col,-1000,1000,reaction_name = "GLUN",Gene_Reaction_Rule = "HGNC:29570 or HGNC:4331")
recon_flux_model.add_reaction(SERtm_S_col,-1000,1000,reaction_name = "SERtm",Gene_Reaction_Rule = "HGNC:16085")

recon_flux_model.reaction_info(1)

recon_flux_model.reaction_info("GLUN")
recon_flux_model.reaction_info("SERtm")