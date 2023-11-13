import os
import pickle
import gurobipy as gp
import copy
import copy as cp
import scipy.optimize as so
from gurobipy import GRB
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Save_And_Load import save_object
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import numpy as np
from matplotlib import pyplot as plt
os.chdir("../..")


flux_modelb = pickle.load(open("Data/Intermediate/Raw_MDL_Models/recon2_2_GLUN_and_SERtm_added.mdl", "rb"))

flux_modelb.reaction_info("EX_orn_b")

flux_modelb.pinch_restricted_exchange_reactions("Data/Input/Experimental_Data/recon2_2/Allowed_Exchange_Reactions", restore_essential=False)

flux_modelb.reaction_info("EX_orn_b")


tol = 1e-5
flux_modelb.update_reaction_bounds("DM_atp_c_",1.9 ,2.1)
flux_modelb.update_reaction_bounds("biomass_reaction",0.0301 , 0.0303 )
flux_modelb.update_reaction_bounds("EX_biomass_c", 0.0301 , 0.0303 )
#flux_modelb.update_reaction_bounds("biomass_other", 0.0016335-tol , 0.0016335+tol )
flux_modelb.update_reaction_bounds("biomass_other", 0.0015 , 0.0017 )

print(flux_modelb.test_feasibility())


flux_modelb.update_reaction_bounds("CITtbm","keep",0)
flux_modelb.update_reaction_bounds("CK","keep",0)
flux_modelb.update_reaction_bounds("CLFORtex",0,0)
for i in range(len(flux_modelb.lb)):
    if flux_modelb.ub[i] == np.inf:
        flux_modelb.update_reaction_bounds(i,"keep",1000)
    if flux_modelb.lb[i] == -np.inf:
        flux_modelb.update_reaction_bounds(i,-1000,"keep")
print(flux_modelb.test_feasibility())
for i in range(len(flux_modelb.lb)):
    if len(str(flux_modelb.ub[i]).split(".")[-1]) > 4:
        print(i)
input()
input()
"Data/Input/Experimental_Data/Measured_Rates.txt"
# There did seem to be some risk of not allowed fluxes being just very close to 0
# Pinching wont fix this
# Write option to pinch to zero if zero falls between lb and ub
#flux_modelb.pinch_reactions(10**(-5))
alpha_list = 10**np.linspace(-2.5,-2.2,30)
obj_list = []
error_list = []
exchange_reaction_dict = flux_modelb.get_exchanged_reaction_info()
print(flux_modelb.get_exchanged_reaction_info())
exchange_reaction_index_list = [flux_modelb.reaction_names.index(exchange_reaction_dict[i]) for i in list(exchange_reaction_dict.keys())]


# A549 has 244 pg of protein
# the fraction of biomass protein is 0.706
mass_of_A549 = 244*10**(-12)/0.706

alpha = np.log10(10**(-12)*(1/mass_of_A549))


exchange_reactions = np.zeros((len(exchange_reaction_index_list),np.size(alpha_list)))


select_reaction_names = ["ENO","AKGDm","CSm","FUM","FUMm","ICDHxm"]
select_reaction_index = [flux_modelb.reaction_names.index(i) for i in select_reaction_names]

select_reactions = np.zeros((len(select_reaction_index),np.size(alpha_list)))

for alpha_ind in range(np.size(alpha_list)):
#    for i in range(len(flux_modelb.reaction_names)):
#        if "EX_urea_b" in flux_modelb.reaction_names[i]:
#            flux_modelb.update_reaction_bounds(flux_modelb.reaction_names[i], "keep", 1*alpha_list[alpha_ind]*0.1)
    X,Obj,ind,exp = flux_modelb.fit_to_experimental_data("Data/Input/Experimental_Data/recon2_2/Measured_Rates.txt",alpha_list[alpha_ind],1)
    counter = 0
    for i in exchange_reaction_index_list:
        exchange_reactions[counter,alpha_ind] = X[i]
        counter+=1
    counter = 0
    for i in select_reaction_index:
        select_reactions[counter, alpha_ind] = X[i]
        counter += 1
    error = 0
    for i in ind:
        error += abs(X[i]-alpha_list[alpha_ind]*exp[i])/abs(alpha_list[alpha_ind]*exp[i])
    error_list.append(error)
    obj_list.append(Obj)
plt.plot(np.log10(alpha_list),np.log10(obj_list))
plt.xlabel("Alpha (log10)")
plt.ylabel("objective function (log10)")
plt.show()

plt.plot(np.log10(alpha_list),np.log10(error_list))
plt.xlabel("Alpha (log10)")
plt.ylabel("sum_i |alpha * e[i] - model.X[i]|/|alpha * e[i]|")
print(np.log10(alpha_list)[np.argmin(np.log10(error_list))])
plt.show()


counter = 0
for i in select_reaction_index:
    plt.plot(np.log10(alpha_list),select_reactions[counter])
    plt.title(f"{flux_modelb.reaction_names[i]}")
    plt.xlabel("Alpha (log10)")
    plt.ylabel("Model Flux")
    plt.show()
    counter += 1

counter = 0
for i in exchange_reaction_index_list:
    if np.sum(exchange_reactions[counter]) != 0:
        plt.plot(np.log10(alpha_list),exchange_reactions[counter])
        plt.title(f"{flux_modelb.reaction_names[i]}")
        plt.xlabel("Alpha (log10)")
        plt.ylabel("Model Flux")
        plt.show()
    counter += 1
print("done")

save_object(flux_modelb, "Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl")
