import pickle
import numpy as np
import matplotlib.pyplot as plt


file_location = f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_A549/align_search_log.pkl"
with open(file_location, "rb") as readfile:
	model_dict = pickle.load(readfile)
index_of_measured_fluxes = model_dict["experimental_flux_index"]
experimental_fluxes = model_dict["experimental_fluxes"]
alpha_list = model_dict["alpha_list"]
best_alpha_ind = model_dict["best_alpha_index"]
model_alignment_array = model_dict["model_alignment"]
performances = model_dict["model_alignment_performances"]
plt.plot(alpha_list,performances)
print(best_alpha_ind)
plt.title(f"Best_alpha = {alpha_list[best_alpha_ind]}\nMass = {1/alpha_list[best_alpha_ind]}")
plt.xlabel("alpha")
plt.ylabel("squared error")
plt.scatter(alpha_list[best_alpha_ind],performances[best_alpha_ind])
plt.show()
measured_exp_flux = experimental_fluxes[index_of_measured_fluxes]
model_fit_exp_flux =  model_alignment_array[best_alpha_ind][index_of_measured_fluxes]/alpha_list[best_alpha_ind]
print(alpha_list[best_alpha_ind])
print(measured_exp_flux)
print(model_fit_exp_flux)
print(np.sum((model_fit_exp_flux-measured_exp_flux)**2/(measured_exp_flux)**2))
print(performances[best_alpha_ind])