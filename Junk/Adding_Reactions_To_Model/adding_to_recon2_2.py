import numpy as np
import os
import sys
import pickle

os.chdir(os.path.dirname(sys.argv[0]))
os.chdir("../..")
flux_modelb = pickle.load(open("Data/Intermediate/Raw_MDL_Models/recon2_2.mdl", "rb"))

# 'GLUN' comparable reaction report:
# [['gln_L_c', -1.0], ['glu_L_c', 1.0], ['h2o_c', -1.0], ['nh4_c', 1.0]]
# 'SERtm' comparable reaction report:
# [['ser_L_c', -1.0], ['ser_L_m', 1.0]]

metabolites = flux_modelb.species_names
i1 = metabolites.index("gln_L_c")
i2 = metabolites.index("glu_L_c")
i3 = metabolites.index("h2o_c")
i4 = metabolites.index("nh4_c")
i5 = metabolites.index("ser_L_c")
i6 = metabolites.index("ser_L_m")

print(np.shape(flux_modelb.S)[0])

GLUN_S_col = np.zeros((np.shape(flux_modelb.S)[0],1))
GLUN_S_col[i1] = -1
GLUN_S_col[i2] = 1
GLUN_S_col[i3] = -1
GLUN_S_col[i4] = 1

SERtm_S_col = np.zeros((np.shape(flux_modelb.S)[0],1))
SERtm_S_col[i5] = -1
SERtm_S_col[i6] = 1

print(len(flux_modelb.reaction_names))
print(len(flux_modelb.grRules))
input()
print(np.shape(GLUN_S_col))

print(flux_modelb.grRules)
print(isinstance(flux_modelb.grRules,list))
input()

dm = flux_modelb.S.todense()
print(dm)

flux_modelb.reaction_info(1)

#Potentially remove HGNC:4331
flux_modelb.add_reaction(GLUN_S_col,-1000,1000,reaction_name = "GLUN",Gene_Reaction_Rule = "HGNC:29570 or HGNC:4331")
flux_modelb.add_reaction(SERtm_S_col,-1000,1000,reaction_name = "SERtm",Gene_Reaction_Rule = "HGNC:16085")
print(len(flux_modelb.reaction_names))
print(len(flux_modelb.grRules))
flux_modelb.reaction_info(1)

flux_modelb.reaction_info("GLUN")
flux_modelb.reaction_info("SERtm")

print(flux_modelb.reaction_names[-2:])
print(len(flux_modelb.reaction_names))
save_object(flux_modelb, "Data/Intermediate/Raw_MDL_Models/recon2_2_GLUN_and_SERtm_added.mdl")