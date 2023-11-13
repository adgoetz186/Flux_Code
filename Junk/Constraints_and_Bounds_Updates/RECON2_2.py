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

tol = 1e-7
flux_modelb.update_reaction_bounds("DM_atp_c_",2-tol ,2+tol)
flux_modelb.update_reaction_bounds("biomass_reaction",0.03025-tol , 0.03025+tol )
flux_modelb.update_reaction_bounds("EX_biomass_c", 0.03025-tol , 0.03025+tol )
flux_modelb.update_reaction_bounds("biomass_other", 0.0016335-tol , 0.0016335+tol )

flux_modelb.update_reaction_bounds("CITtbm","keep",0)
flux_modelb.update_reaction_bounds("CK","keep",0)
flux_modelb.update_reaction_bounds("CLFORtex",0,0)
for i in range(len(flux_modelb.lb)):
    if flux_modelb.ub[i] == np.inf:
        flux_modelb.update_reaction_bounds(i,"keep",10000)
    if flux_modelb.lb[i] == -np.inf:
        flux_modelb.update_reaction_bounds(i,-10000,"keep")
"Data/Input/Experimental_Data/Measured_Rates.txt"
print(flux_modelb.test_feasibility())

# There did seem to be some risk of not allowed fluxes being just very close to 0
# Pinching wont fix this
# Write option to pinch to zero if zero falls between lb and ub
#flux_modelb.pinch_reactions(10**(-5))
#alpha_list = 10**np.linspace(20,22,30)
obj_list = []
error_list = []
exchange_reaction_dict = flux_modelb.get_exchanged_reaction_info()
print(flux_modelb.get_exchanged_reaction_info())
exchange_reaction_index_list = [flux_modelb.reaction_names.index(exchange_reaction_dict[i]) for i in list(exchange_reaction_dict.keys())]


# A549 has 244 pg of protein
# the fraction of biomass protein is 0.706
mass_of_A549 = 244*10**(-12)/0.706

alpha = 10**(-12)*(1/mass_of_A549)

print(np.log10(alpha))

exchange_reactions = np.zeros((len(exchange_reaction_index_list),1))


select_reaction_names = ["ENO","AKGDm","CSm","FUM","FUMm","ICDHxm"]
select_reaction_index = [flux_modelb.reaction_names.index(i) for i in select_reaction_names]

select_reactions = np.zeros((len(select_reaction_index),1))



X,Obj,ind,exp = flux_modelb.fit_to_experimental_data("Data/Input/Experimental_Data/recon2_2/Measured_Rates.txt",alpha,1)

for i in ind:
    flux_modelb.update_reaction_bounds(flux_modelb.reaction_names[i], X[i] - tol, X[i] + tol)

counter = 0
for i in exchange_reaction_index_list:
    exchange_reactions[counter,0] = X[i]
    counter+=1
counter = 0
for i in select_reaction_index:
    select_reactions[counter, 0] = X[i]
    counter += 1
error = 0
for i in ind:
    error += abs(X[i]-alpha*exp[i])/abs(alpha*exp[i])
print(error)
error_list.append(error)
obj_list.append(Obj)
labels = [flux_modelb.reaction_names[i] for i in ind]
plot_experiment = [exp[i] for i in ind]
plot_model = [X[i]/alpha for i in ind]

etplot_experiment = [alpha*exp[i] for i in ind]
etplot_model = [X[i] for i in ind]

error_term = 0
for i in range(len(etplot_model)):
    error_term+= (etplot_experiment[i]-etplot_model[i])**2/etplot_experiment[i]**2
print(error_term)
x = np.arange(len(labels))  # the label locations

width = 0.25
fig, ax = plt.subplots()
rects1 = ax.bar(x - width / 2, plot_experiment, width, label='Experiment')
rects2 = ax.bar(x + width / 2, plot_model, width, label='Model')

ax.set_ylabel('Excretion Flux (fmol/(cell * hr))')
ax.set_title('Model and Experimental Agreement')
print(x)
print(labels)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=90)
ax.legend()
# ax.bar_label(rects1, padding=3)
# ax.bar_label(rects2, padding=3)

fig.tight_layout()
plt.hlines(0, -1, 22, color="black")
plt.show()

print(obj_list)
print(10**obj_list[0])
plt.plot(np.log10(alpha),np.log10(obj_list))
plt.xlabel("Alpha (log10)")
plt.ylabel("objective function (log10)")
plt.show()

print(error_list)
plt.plot(np.log10(alpha),np.log10(error_list))
plt.xlabel("Alpha (log10)")
plt.ylabel("sum_i |alpha * e[i] - model.X[i]|/|alpha * e[i]|")
plt.show()

counter = 0
for i in exchange_reaction_index_list:
    if exchange_reactions[counter] != 0:
        print(f"{flux_modelb.reaction_names[i]}", exchange_reactions[counter])
    counter += 1
print("done")
print(flux_modelb.test_feasibility())
input()
save_object(flux_modelb, "Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl")
