import os
import pickle
import gurobipy as gp
import copy
import copy as cp
import scipy.optimize as so
from gurobipy import GRB
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import numpy as np
from matplotlib import pyplot as plt

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

data_file_location = "Data/experimental_alignment_data/recon2_2_A549/"
data_file_name = "align_search_log"

model_file_location = "Data/Models/pkl_models/recon2_2_A549/"
model_file_name = "exp_aligned_alpha_given"

with open(f"{data_file_location}{data_file_name}.pkl", "rb") as readfile:
    model_dict = pickle.load(readfile)

print(model_dict.keys())
X_array = model_dict["model_alignment"]
best_ind = model_dict["best_alpha_index"]
alpha_list = model_dict["alpha_list"]
exp = model_dict["experimental_fluxes"]
Obj_array = model_dict["model_alignment_performances"]
ind = model_dict["experimental_flux_index"]
print(ind)
recon2_2_rxn_added = Flux_Balance_Model()
recon2_2_rxn_added.load_pkl_model(f"{model_file_location}{model_file_name}")

exchange_reaction_dict = recon2_2_rxn_added.get_exchanged_reaction_info()
exchange_reaction_index_list = [recon2_2_rxn_added.reaction_names.index(exchange_reaction_dict[i]) for i in list(exchange_reaction_dict.keys())]

select_reaction_names = ["ARGNm","AGMTm","UREAtm","ENO","AKGDm","CSm","FUM","FUMm","ICDHxm"]
select_reaction_index = [recon2_2_rxn_added.reaction_names.index(i) for i in select_reaction_names]

select_reactions = np.zeros((len(select_reaction_index),np.size(alpha_list)))
exchange_reactions = np.zeros((len(exchange_reaction_index_list),np.size(alpha_list)))


X = X_array[best_ind]
alpha = alpha_list[best_ind]
labels = [recon2_2_rxn_added.reaction_names[i] for i in ind]
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

input()

error_list = []
for alpha_ind in range(np.shape(X_array)[0]):
    counter = 0
    for i in exchange_reaction_index_list:
        exchange_reactions[counter,alpha_ind] = X_array[alpha_ind,i]
        counter+=1
    counter = 0
    for i in select_reaction_index:
        select_reactions[counter, alpha_ind] = X_array[alpha_ind,i]
        counter += 1
    error = 0
    for i in ind:
        error += abs(X_array[alpha_ind,i]-alpha_list[alpha_ind]*exp[i])/abs(alpha_list[alpha_ind]*exp[i])
    error_list.append(error)
plt.plot(np.log10(alpha_list),np.log10(Obj_array))
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
    plt.title(f"{recon2_2_rxn_added.reaction_names[i]}")
    plt.xlabel("Alpha (log10)")
    plt.ylabel("Model Flux")
    plt.show()
    counter += 1

counter = 0
for i in exchange_reaction_index_list:
    if np.sum(exchange_reactions[counter]) != 0:
        plt.plot(np.log10(alpha_list),exchange_reactions[counter])
        plt.title(f"{recon2_2_rxn_added.reaction_names[i]}")
        plt.xlabel("Alpha (log10)")
        plt.ylabel("Model Flux")
        plt.show()
    counter += 1
print("done")