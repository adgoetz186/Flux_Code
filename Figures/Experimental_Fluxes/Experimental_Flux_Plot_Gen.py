import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import matplotlib.patches as mpatches
import scipy.stats as st

os.chdir("../../..")
flux_dataframe = pd.read_csv("Flux_Code/Data/Input/Experimental_Data/TC_B6/Fluxes.csv",header = 0)
metabolite_name_list = list(flux_dataframe.keys())[2:]
position = 0
c3 = "#0021A5"
c4 = "#FA4616"
c0 = "black"


B6 = np.array([0.021488297,0.019910581,0.020568351,0.020277418])
TC = np.array([0.024513342,0.025816285,0.024089732,0.021298621])
print(st.ranksums(B6,TC))
input()

generated_metabolites = []
for j in metabolite_name_list:
	metabolite_to_plot = ""
	largest_range = 0
	for i in metabolite_name_list:
		try:
			metabolite_range = np.max(flux_dataframe[i][:])-np.min(flux_dataframe[i][:])
			if metabolite_range > largest_range and i not in generated_metabolites:
				metabolite_to_plot = i
				largest_range = metabolite_range
		except TypeError:
			print(i,np.max(flux_dataframe[i][:]),np.min(flux_dataframe[i][:]))
			if metabolite_to_plot == "":
				metabolite_to_plot = i
		print(metabolite_to_plot)
	generated_metabolites.append(metabolite_to_plot)
	try:
		plt.boxplot(flux_dataframe[metabolite_to_plot][:4],positions = [position], whis = (0, 100), showfliers = False, patch_artist = True,
		boxprops = dict(facecolor=c3, color=c3), capprops = dict(color=c3), whiskerprops = dict(color=c3),
		flierprops = dict(color=c3, markeredgecolor=c3), medianprops = dict(color=c0))
		
		plt.boxplot(flux_dataframe[metabolite_to_plot][4:], positions=[position+.25], whis=(0, 100), showfliers=False, patch_artist=True,
		            boxprops=dict(facecolor=c4, color=c4), capprops=dict(color=c4), whiskerprops=dict(color=c4),
		            flierprops=dict(color=c4, markeredgecolor=c4), medianprops=dict(color=c0))
		position += 1
	except TypeError:
		print(metabolite_to_plot)
B6_patch = mpatches.Patch(color=c3, label='B6')
TC_patch = mpatches.Patch(color=c4, label='TC')
plt.legend(handles=[B6_patch,TC_patch])
plt.xticks([.25/2+i for i in range(len(generated_metabolites))],generated_metabolites)
plt.xlabel("Metabolite")
plt.ylabel("Flux (fmol/(cell*hr))")
plt.show()
sns.boxplot(data=flux_dataframe, y="Leucine", x="Strain")
plt.show()