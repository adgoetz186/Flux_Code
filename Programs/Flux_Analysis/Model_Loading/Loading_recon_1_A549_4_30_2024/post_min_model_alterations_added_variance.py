import os
import sys
import numpy as np
import pandas as pd
import copy as cp
from pathlib import Path

# Tweaks minimal model with biological knowledge

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
from Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
# _____ Setting the CWD to be Flux_Code END _____

# removes reactions that cant carry flux
# removes curated reactions:
# ALR,

min_model_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/pos_min_model.json")
output_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024")
measured_rates_location = Path("Data/experimental_alignment_data/recon_1_A549/Measured_Rates.txt")



recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(min_model_file_location)

# Add check for pinch, just look at reactions


print(recon_flux_model.test_feasibility())
print(recon_flux_model.reaction_info("ALR"))
print(recon_flux_model.reaction_info("PGI_inverted"))
print(recon_flux_model.reaction_info("PGI"))
print(recon_flux_model.reaction_info("SBTR"))


del recon_flux_model.model_dict["rxn_dict"]["ALR"]
del recon_flux_model.model_dict["rxn_dict"]["PGI_inverted"]
del recon_flux_model.model_dict["rxn_dict"]["SBTR"]
print(recon_flux_model.test_feasibility())

input()

list_of_pinched_exchange = pd.read_csv(measured_rates_location)["reaction"].to_list()
print(list_of_pinched_exchange)
print(recon_flux_model.test_feasibility())
input()

recon_model_MK1 = cp.deepcopy(recon_flux_model)
lb, ub, S, b, rxn_list, met_list = recon_model_MK1.dicts_to_mats()

biomass_rxn = (lb[rxn_list.index("biomass_reaction")] + ub[rxn_list.index("biomass_reaction")])/2
print(biomass_rxn,lb[rxn_list.index("biomass_reaction")],ub[rxn_list.index("biomass_reaction")])
print(recon_flux_model.test_feasibility())
input(0)
CV_T = np.sqrt(5.57) / 8.6
# https://www.pnas.org/doi/full/10.1073/pnas.2116260119 (CD8, but seems to have similar average)
# gives inter division time which seems similar to doubling time, time it takes from cell creation to cell doubling
# total variance is 5.57 h^2, mean is 8.6 h
# https://www.spandidos-publications.com/10.3892/etm.2019.7873
CV_H = 6.7/33.7
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4802345/
# For now we will use a CV of 0.2 which lies between these values (though much closer to HeLa cells)
CV = 0.2
print(CV_T,CV_H)
input(0)

pinch_dt = np.log(2) / biomass_rxn
sqrt_dt = pinch_dt * CV
recon_model_MK1.update_reaction_bounds("biomass_reaction", np.log(2) / (pinch_dt + sqrt_dt),np.log(2) / (pinch_dt - sqrt_dt))
recon_model_MK1.reaction_info("biomass_reaction")
print(recon_flux_model.test_feasibility())
input(2)


for i in list_of_pinched_exchange:
    if i in recon_model_MK1.model_dict['rxn_dict'].keys():
        recon_model_MK1.update_reaction_bounds(i, lb[rxn_list.index(i)]*.5, ub[rxn_list.index(i)]*2)
        recon_model_MK1.reaction_info(i)
    elif (i+"_inverted") in recon_model_MK1.model_dict['rxn_dict'].keys():
        recon_model_MK1.update_reaction_bounds((i+"_inverted"),lb[rxn_list.index((i+"_inverted"))] * .5,
                                              ub[rxn_list.index((i+"_inverted"))] * 2)
        recon_model_MK1.reaction_info((i+"_inverted"))
print(recon_flux_model.test_feasibility())
input(3)
rxn_bound_dict = recon_model_MK1.fva()
list_of_rxns = list(recon_model_MK1.model_dict["rxn_dict"].keys())
list_of_rxns_to_remove = []
for rxn in list_of_rxns:
    if abs(recon_model_MK1.model_dict["rxn_dict"][rxn]["ub"] - recon_model_MK1.model_dict["rxn_dict"][rxn]["lb"])<1e-12 and abs(recon_model_MK1.model_dict["rxn_dict"][rxn]["ub"])<1e-12:
        list_of_rxns_to_remove.append(rxn)


#print(list_of_rxns_to_remove)
#for i in list_of_rxns_to_remove:
#    del recon_model_MK1.model_dict["rxn_dict"][i]



for rxn_name in recon_model_MK1.model_dict["rxn_dict"].keys():
    recon_model_MK1.update_reaction_bounds(rxn_name,rxn_bound_dict[rxn_name]['lb'],rxn_bound_dict[rxn_name]['ub'])
print(recon_model_MK1.test_feasibility())
recon_model_MK1.model_dict["model_name"] = "CurVar_pos_min_model"
print("done")


for i in list_of_pinched_exchange:
    if i in recon_model_MK1.model_dict['rxn_dict'].keys():
        if abs(recon_model_MK1.model_dict['rxn_dict'][i]['lb']*4 - recon_model_MK1.model_dict['rxn_dict'][i]['ub']) > abs(recon_model_MK1.model_dict['rxn_dict'][i]['lb'])/1e3:
            recon_model_MK1.reaction_info(i)
            print(abs(recon_model_MK1.model_dict['rxn_dict'][i]['lb']*4 - recon_model_MK1.model_dict['rxn_dict'][i]['ub']), recon_model_MK1.model_dict['rxn_dict'][i]['lb']*4,recon_model_MK1.model_dict['rxn_dict'][i]['lb']*2 , recon_model_MK1.model_dict['rxn_dict'][i]['ub'])
            input('warn')
    elif (i+"_inverted") in recon_model_MK1.model_dict['rxn_dict'].keys():
        if abs(recon_model_MK1.model_dict['rxn_dict'][(i+"_inverted")]['lb'] * 4 - recon_model_MK1.model_dict['rxn_dict'][(i+"_inverted")][
            'ub']) >  abs(recon_model_MK1.model_dict['rxn_dict'][(i+"_inverted")]['lb'])/1e3:
            recon_model_MK1.reaction_info((i+"_inverted"))
            print(abs(recon_model_MK1.model_dict['rxn_dict'][(i+"_inverted")]['lb']*4 - recon_model_MK1.model_dict['rxn_dict'][(i+"_inverted")]['ub']), recon_model_MK1.model_dict['rxn_dict'][(i+"_inverted")]['lb'] * 4,recon_model_MK1.model_dict['rxn_dict'][(i+"_inverted")]['lb']*2, recon_model_MK1.model_dict['rxn_dict'][(i+"_inverted")]['ub'])
            input('warn')
print(recon_model_MK1.test_feasibility())
input()
recon_model_MK1.save_model_as_fast_key_json(output_file_location)
