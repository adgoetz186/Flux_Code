import os
import pickle
import sys
from pathlib import Path
from pympler import asizeof
import gurobipy as gp
import copy
import re
import pandas as pd
from gurobipy import GRB
import scipy as sc
import scipy.linalg
import time
import scipy.sparse as sp
import numpy as np
import random
from matplotlib import pyplot as plt
from Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model

def minimal_flux_list_multimodel(flux_model_list, min_model_matrix_location):
	"""
	Takes a list of flux model objects and an essential flux matrix and finds a core model that applies to all models.
	Models must have the same reaction name order
	"""
	
	matrix_list = []
	header_list = []
	for filename in os.listdir(min_model_matrix_location):
		if "header" in filename:
			header_list.append(np.load(min_model_matrix_location / filename))
			print(len(header_list[-1]))
		else:
			matrix_list.append(np.load(min_model_matrix_location / filename))
	# Checks to make sure all headers are the same
	try:
		all_header_same = np.alltrue(np.array([np.alltrue(i == header_list[0]) for i in header_list]))
		if not all_header_same:
			raise ValueError
	except ValueError:
		print("All minimal model headers are not the same, this is unexpected and will likely lead to unexpected and undesired results")

	all_mat = np.average(np.vstack(matrix_list), axis=0)
	print(np.shape(all_mat))
	
	essentiality_scores = {}
	for rxn in header_list[0]:
		print(np.where(header_list[0] == rxn))
		essentiality_scores[rxn] = all_mat[np.where(header_list[0] == rxn)][0]
	
	essential_values = list(set(essentiality_scores.values()))
	essential_values.sort(reverse=True)
	list_or_rxns_to_test = []
	for i in essential_values:
		for rxn in header_list[0]:
			if essentiality_scores[rxn] == i:
				list_or_rxns_to_test.append(rxn)
	saved_ub_list = [{} for i in flux_model_list]
	saved_lb_list = [{} for i in flux_model_list]

	# remove all non constrained fluxes from all models
	for flux_model_ind in range(len(flux_model_list)):
		for rxn in list_or_rxns_to_test:
			if (flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]["ub"] < 0 or flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]["lb"] > 0):
				continue
			else:
				saved_ub_list[flux_model_ind][rxn] = copy.deepcopy(flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]['ub'])
				saved_lb_list[flux_model_ind][rxn] = copy.deepcopy(flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]['lb'])
				flux_model_list[flux_model_ind].update_reaction_bounds(rxn, 0, 0)
	essential_flux_main = []
	count = 0
	# The current code structure slowly adds reactions back but does not stop as soon as one final reaction is found
	# Rather all other similarly ranked reactions are also tested and those which make the model feasible are
	# accepted
	accepted = 0
	current_ind = 0
	essential_reaction_scores = []
	for rxn in list_or_rxns_to_test:
		# restores fluxes of reactions
		for flux_model_ind in range(len(flux_model_list)):
			if (flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]["ub"] < 0 or flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]["lb"] > 0):
				continue
			else:
				flux_model_list[flux_model_ind].update_reaction_bounds(rxn, saved_lb_list[flux_model_ind][rxn], saved_ub_list[flux_model_ind][rxn])
		# finds which models are feasible
		
		# This prevents the code from testing each model if it hits one thats infeasible
		if flux_model_list[current_ind].test_feasibility():
			current_ind += 1
			if current_ind == (len(flux_model_list)):
				essential_flux_main.append(rxn)
				break
			while flux_model_list[current_ind].test_feasibility():
				current_ind +=1
				if current_ind == (len(flux_model_list)):
					essential_flux_main.append(rxn)
					break
		print(accepted,current_ind)
		essential_flux_main.append(rxn)
		accepted+=1
		
		# if there are still infeasible models the tested reactions are all essential and the algorithm goes to the
		# next set of reactions (less essential reactions)
		
		
	for flux_model_ind in range(len(flux_model_list)):
		model_names = list(flux_model_list[flux_model_ind].model_dict["rxn_dict"].keys())
		for rxn in model_names:
			if rxn not in essential_flux_main:
				del flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]
	return flux_model_list


def minimal_flux_list_multimodel_model_differences(flux_model_list, min_model_matrix_location):
	"""
	Takes a list of flux model objects and an essential flux matrix and finds a core model that applies to all models.
	Models must have the same reaction name order
	"""
	
	load_dict = {}
	for filename in os.listdir(min_model_matrix_location):
		load_dict[filename.split("_")[0]] = {}
	
	for filename in os.listdir(min_model_matrix_location):
		if "header" in filename:
			load_dict[filename.split("_")[0]]["header"] = list(np.load(min_model_matrix_location / filename))
		else:
			load_dict[filename.split("_")[0]]["matrix_"+filename.split("_")[-1].replace(".npy","")] = np.load(min_model_matrix_location / filename)
	rxn_names = []
	for key in load_dict.keys():
		for rxn_name in load_dict[key]["header"]:
			if rxn_name not in rxn_names:
				rxn_names.append(rxn_name)
	print(rxn_names)
	print(load_dict)

	
	essentiality_scores = {}
	for rxn in rxn_names:
		div = 0
		essential_count = 0
		for model in load_dict.keys():
			for model_key in load_dict[model].keys():
				if "matrix" in model_key:
					div+=np.shape(load_dict[model][model_key])[0]
					if rxn in load_dict[model]["header"]:
						essential_count += np.sum(load_dict[model][model_key][:,load_dict[model]["header"].index(rxn)])
		
		essentiality_scores[rxn] = essential_count/div

	essential_values = list(set(essentiality_scores.values()))
	essential_values.sort(reverse=True)

	list_of_rxns_to_test = []
	for i in essential_values:
		for rxn in rxn_names:
			if essentiality_scores[rxn] == i:
				list_of_rxns_to_test.append(rxn)
	saved_ub_list = [{} for i in flux_model_list]
	saved_lb_list = [{} for i in flux_model_list]
	
	# remove all non constrained fluxes from all models
	for flux_model_ind in range(len(flux_model_list)):
		for rxn in list_of_rxns_to_test:
			if rxn in flux_model_list[flux_model_ind].model_dict["rxn_dict"].keys():
				if (flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]["ub"] < 0 or
						flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]["lb"] > 0):
					continue
				else:
					saved_ub_list[flux_model_ind][rxn] = copy.deepcopy(
						flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]['ub'])
					saved_lb_list[flux_model_ind][rxn] = copy.deepcopy(
						flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]['lb'])
					flux_model_list[flux_model_ind].update_reaction_bounds(rxn, 0, 0)
	essential_flux_main = []
	count = 0
	# The current code structure slowly adds reactions back but does not stop as soon as one final reaction is found
	# Rather all other similarly ranked reactions are also tested and those which make the model feasible are
	# accepted
	accepted = 0
	current_ind = 0
	essential_reaction_scores = []
	for rxn in list_of_rxns_to_test:
		# restores fluxes of reactions
		for flux_model_ind in range(len(flux_model_list)):
			if rxn in flux_model_list[flux_model_ind].model_dict["rxn_dict"].keys():
				if (flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]["ub"] < 0 or
						flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]["lb"] > 0):
					continue
				else:
					flux_model_list[flux_model_ind].update_reaction_bounds(rxn, saved_lb_list[flux_model_ind][rxn],
					                                                       saved_ub_list[flux_model_ind][rxn])
		# finds which models are feasible
		
		if current_ind == (len(flux_model_list)):
			# This handles the case where the final model became feasible inside of a while loop. Python does not allow double breaks
			break
		elif flux_model_list[current_ind].test_feasibility():
			# This prevents the code from testing each model if it hits one thats infeasible
			print(rxn)
			current_ind += 1
			if current_ind == (len(flux_model_list)):
				print(1)
				#print(flux_model_list[current_ind].test_feasibility())
				essential_flux_main.append(rxn)
				break
			while flux_model_list[current_ind].test_feasibility():
				print(2)
				current_ind += 1
				print(current_ind,len(flux_model_list))
				if current_ind == (len(flux_model_list)):
					print(3)
					essential_flux_main.append(rxn)
					break
		print(accepted, current_ind)
		essential_flux_main.append(rxn)
		accepted += 1
	
	# if there are still infeasible models the tested reactions are all essential and the algorithm goes to the
	# next set of reactions (less essential reactions)
	
	for flux_model_ind in range(len(flux_model_list)):
		model_names = list(flux_model_list[flux_model_ind].model_dict["rxn_dict"].keys())
		for rxn in model_names:
			if rxn not in essential_flux_main:
				del flux_model_list[flux_model_ind].model_dict["rxn_dict"][rxn]
	return flux_model_list

def make_uniform_pos(flux_model_list,fva_first = True):
	if fva_first:
		for flux_model in flux_model_list:
			new_bounds = flux_model.fva()
			for rxn_name in flux_model.model_dict["rxn_dict"].keys():
				flux_model.update_reaction_bounds(rxn_name,new_bounds[rxn_name]["lb"],new_bounds[rxn_name]["ub"])
	for flux_model in flux_model_list:
		flux_model.convert_model_to_positive_flux()
	set_of_all_rxns = []
	for flux_model in flux_model_list:
		print(flux_model.model_dict["model_name"], flux_model.test_feasibility())
		set_of_all_rxns += list(flux_model.model_dict["rxn_dict"].keys())
		print(flux_model.model_dict["model_name"], flux_model.test_feasibility())
	set_of_all_rxns = set(set_of_all_rxns)
	
	
	for flux_model in flux_model_list:
		# finds the reactions missing from the current model
		reactions_to_add = set_of_all_rxns - set(flux_model.model_dict["rxn_dict"].keys())
		for rxn_name in reactions_to_add:
			#print(rxn_name in flux_model.rxn_dict.keys())
			# grabs the missing reaction's dictionary value
			rxn_dict_value = {}
			for flux_model_b in flux_model_list:
				if rxn_name in flux_model_b.model_dict["rxn_dict"].keys():
					rxn_dict_value = flux_model_b.model_dict["rxn_dict"][rxn_name]
			# adds the missing reaction's dictionary value
			flux_model.model_dict["rxn_dict"][rxn_name] = copy.deepcopy(rxn_dict_value)
			flux_model.update_reaction_bounds(rxn_name,0,0)
			print(flux_model.test_feasibility())
			if "inverted" in rxn_name and False:
				if rxn_name.split("_")[0] in flux_model.model_dict["rxn_dict"].keys():
					print(rxn_name,flux_model.model_dict["rxn_dict"][rxn_name]["rxn_metabolites"])
					print(rxn_name.split("_")[0] ,flux_model.model_dict["rxn_dict"][rxn_name.split("_")[0]]["rxn_metabolites"])
			else:
				if rxn_name+"_inverted" in flux_model.model_dict["rxn_dict"].keys() and False:
					print(rxn_name,flux_model.model_dict["rxn_dict"][rxn_name]["rxn_metabolites"])
					print(rxn_name+"_inverted" ,flux_model.model_dict["rxn_dict"][rxn_name+"_inverted"]["rxn_metabolites"])
	return flux_model_list

def comp_rxn_dict(rxn_name,model_1,model_2):
	dict_1 = model_1.rxn_dict[rxn_name]
	dict_2 = model_2.rxn_dict[rxn_name]
	dicts_equal = True
	for i in dict_1.keys():
		if dict_1[i]!=dict_2[i]:
			dicts_equal = False
	return dicts_equal