from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
import numpy as np

# Adds reactions which are important for model behavior (GLUN and SERtm)

# File locations should end with "/"
# File names should not end with file specifier (.txt, .py, etc)
raw_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_A549/raw"

output_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_A549/rxn_added"


recon_flux_model = Flux_Balance_Model()

recon_flux_model.load_pkl_model(raw_file_location)
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

recon_flux_model.reaction_info("GLUN")
recon_flux_model.reaction_info("SERtm")

tol = 1e-7
recon_flux_model.update_reaction_bounds("DM_atp_c_",2-tol ,2+tol)
recon_flux_model.update_reaction_bounds("biomass_reaction",0.03025-tol , 0.03025+tol )
recon_flux_model.update_reaction_bounds("EX_biomass_c", 0.03025-tol , 0.03025+tol )
recon_flux_model.update_reaction_bounds("biomass_other", 0.0016335-tol , 0.0016335+tol )
recon_flux_model.update_reaction_bounds("CITtbm","keep",0)
recon_flux_model.update_reaction_bounds("CK","keep",0)
recon_flux_model.update_reaction_bounds("CLFORtex",0,0)



recon_flux_model.save_model_as_pkl(output_file_location)