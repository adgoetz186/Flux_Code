import pickle
import sys
from functools import partial
import sympy as sy
from pympler import asizeof
import gurobipy as gp
import copy
import gzip
import re
import copy as cp
import os
from pathlib import Path
import pandas as pd
from gurobipy import GRB
import scipy as sc
import scipy.linalg
import time
import json
import scipy.sparse as sp
import numpy as np
import igraph as ig

print(ig.__version__)
import scipy.io as spio
import scipy.linalg as sla
import random
from matplotlib import pyplot as plt


def loadmat(filename):
	'''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
	data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
	return _check_keys(data)


def _check_keys(dict):
	'''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
	for key in dict:
		if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
			dict[key] = _todict(dict[key])
	return dict


def _todict(matobj):
	'''
    A recursive function which constructs from matobjects nested dictionaries
    '''
	dict = {}
	for strg in matobj._fieldnames:
		elem = matobj.__dict__[strg]
		if isinstance(elem, spio.matlab.mio5_params.mat_struct):
			dict[strg] = _todict(elem)
		else:
			dict[strg] = elem
	return dict


def size_order_key(e):
	return len(e)


class Flux_Balance_Model:
	# Im not adding an add species option, species should all be defined during initialization
	
	model = "none"
	
	def __init__(self, model_dict=None):
		# model_dict should have keys:
		# model_name
		# rxn_dict
		# met_dict
		# gene_dict
		# gurobi_token
		self.model_dict = model_dict
	
	#def add_gp_key_env_to_model(self, TokenServer_name):
		# for hpg, TokenServer_name should be 'grb-ts.ufhpc'
	#	env = gp.Env(empty=True)
	#	env.setParam('TokenServer', TokenServer_name)
	#	env.start()
	#	self.model_dict['gurobi_token'] = env
	
	def add_gp_key_env_to_model(self, TokenServer_name):
		# for hpg, TokenServer_name should be 'grb-ts.ufhpc'
		env = gp.Env({'TokenServer': TokenServer_name,"logfilename":str(time.time()).replace(".","")})
		#env.setParam('TokenServer', TokenServer_name)
		env.start()
		self.model_dict['gurobi_token'] = env
	
	def save_model_as_pkl(self, filepath):
		dict_of_model = {}
		dict_of_model["model_name"] = self.model_name
		dict_of_model["total_flux_limits"] = self.total_flux_limits
		dict_of_model["rxn_dict"] = self.rxn_dict
		dict_of_model["met_dict"] = self.met_dict
		with open(filepath.parent / (filepath.stem + ".pkl"), "wb") as outfile:
			pickle.dump(dict_of_model, outfile)
	
	def save_model_as_mat(self, filename):
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		gr_list = self.get_grRule_list()
		mat_dict = {}
		mat_dict["lb"] = lb
		mat_dict["ub"] = ub
		mat_dict["S"] = S
		mat_dict["b"] = b
		mat_dict["rxn_list"] = rxn_list
		mat_dict["met_list"] = met_list
		spio.savemat(filename, mat_dict)
	
	def create_model_from_components(self, model_name, S, lb, ub, b, rxn_names, met_names,
	                                 additional_rxn_info_dict=None, additional_met_info_dict=None):
		self.model_dict = {}
		self.model_dict["gurobi_token"] = None
		self.model_dict["model_name"] = model_name
		rxn_dict = {}
		for i in range(len(rxn_names)):
			rxn_sub_dict = {}
			rxn_met = [met_names[j] for j in np.nonzero(S[:, i])[0]]
			stochio_values_met = [S[:, i][j] for j in np.nonzero(S[:, i])[0]]
			rxn_sub_dict["rxn_metabolites"] = dict(zip(rxn_met, stochio_values_met))
			rxn_sub_dict["lb"] = lb[i]
			rxn_sub_dict["ub"] = ub[i]
			
			if additional_rxn_info_dict is not None:
				for key in additional_rxn_info_dict.keys():
					if len(additional_rxn_info_dict[key][i]) != 0:
						rxn_sub_dict[key] = additional_rxn_info_dict[key][i]
					else:
						rxn_sub_dict[key] = ""
			
			rxn_dict[rxn_names[i]] = rxn_sub_dict
		
		self.model_dict["rxn_dict"] = rxn_dict
		met_dict = {}
		for i in range(len(met_names)):
			met_sub_dict = {}
			met_sub_dict["b"] = b[i]
			
			if additional_met_info_dict is not None:
				for key in additional_met_info_dict.keys():
					if isinstance(additional_met_info_dict[key][i], np.int32):
						met_sub_dict[key] = str(additional_met_info_dict[key][i])
					elif len(additional_met_info_dict[key][i]) != 0:
						met_sub_dict[key] = additional_met_info_dict[key][i]
					else:
						met_sub_dict[key] = ""
			met_dict[met_names[i]] = met_sub_dict
		self.model_dict["met_dict"] = met_dict
		# this might need some tweaking
		try:
			gene_list = self.get_used_gene_list()
			gene_dict = {}
			for gene in gene_list:
				degen = 0
				for gene_2 in gene_list:
					if gene.split(".")[0] == gene_2.split(".")[0]:
						degen += 1
				gene_dict[gene] = {"entrez": gene.split(".")[0], "degen": degen}
			self.model_dict["gene_dict"] = gene_dict
		except KeyError:
			print(f"No gene data supplied")
	
	def load_mat_model(self, filename, new_model_name=None, read_model_name="", model_comp=None, gz=False):
		# if model_name == "" the program guesses what the model key is by looking at the size of all dict entries
		# the largest entry is likely to be the one that holds the model
		# note sys.getsizeof() is insufficient as the model is often a dictionary which references large objects
		# and sys.getsizeof() will omit those referenced values
		# model_comp: dict which contains the names for each model component with the key being the "standard name"
		#   standard names: Stochiometric matrix = "S", lower bounds = "lb", upper bounds = "ub",
		#   rhs of Sv = b constraint = "b", dictionary of pinched reactions = "pinched_reactions",
		#   list of reaction names = "reaction_names", gene reaction rule list = "grRules",
		#   list of gene names = "genes", list of metabolite names = "metabolite_names"
		#   list of metabolite compositions = "metabolite_comp"
		#   In the case where an optional model component is not expected to be found in the model being loaded
		#   the key should remain unchanged, but the corresponding value should be None
		#   This will be assumed to be the case for pinched_reactions
		if new_model_name is None:
			new_model_name = "Model"
		if gz:
			with gzip.open(filename, 'rb') as f:
				mat = loadmat(f)
		else:
			mat = loadmat(filename)
		if read_model_name == "":
			size_list = [asizeof.asizeof(mat[key]) for key in mat.keys()]
			read_model_name = list(mat.keys())[size_list.index(max([asizeof.asizeof(mat[key]) for key in mat.keys()]))]
		if model_comp == None:
			model_comp = {}
			model_comp["S"] = "S"
			model_comp["lb"] = "lb"
			model_comp["ub"] = "ub"
			model_comp["b"] = "b"
			model_comp["pinched_reactions"] = None
			model_comp["metabolite_names"] = "metabolite_names"
			model_comp["reaction_names"] = "reaction_names"
			model_comp["grRules"] = "grRules"
			model_comp["genes"] = "genes"
			model_comp["metabolite_comp"] = "metabolite_comp"
		try:
			test_key = mat[read_model_name]
		except KeyError:
			print(f"It is likely model name was not valid, valid model name will likely be one of the following:")
			for key in mat.keys():
				print(key)
		req_keys = {"S", "lb", "ub", "b", "pinched_reactions", "reaction_names", "grRules", "genes", "metabolite_names",
		            "metabolite_comp"}
		curr_keys = set(model_comp.keys())
		if len(req_keys - curr_keys) > 0:
			print(f"There are {len(req_keys - curr_keys)} required key(s) missing from model_comp")
			for i in (req_keys - curr_keys):
				print(i)
		print(mat[read_model_name].keys())
		input()
		try:
			model_dict = {}
			tried_key = model_comp["S"]
			S = sp.csr_matrix(mat[read_model_name][tried_key]).astype(np.float)
			S = S.toarray()
			model_dict["S"] = S
			tried_key = model_comp["lb"]
			lb = mat[read_model_name][tried_key].astype(np.float)
			model_dict["lb"] = lb
			tried_key = model_comp["ub"]
			ub = mat[read_model_name][tried_key].astype(np.float)
			model_dict["ub"] = ub
			tried_key = model_comp["b"]
			b = mat[read_model_name][tried_key].astype(np.float)
			model_dict["b"] = b
			
			tried_key = model_comp["reaction_names"]
			if tried_key == None:
				if S != None:
					reaction_names = [f"Reaction_{i}" for i in range(np.shape(self.S)[1])]
				else:
					reaction_names = None
			else:
				reaction_names = list(mat[read_model_name][tried_key])
			model_dict["reaction_names"] = reaction_names
			
			tried_key = model_comp["grRules"]
			if tried_key == None:
				grRules = None
			else:
				grRules = list(mat[read_model_name][tried_key])
			model_dict["grRules"] = grRules
			
			tried_key = model_comp["metabolite_names"]
			if tried_key == None:
				if S != None:
					metabolite_names = [f"metabolite_{i}" for i in range(np.shape(self.S)[0])]
				else:
					metabolite_names = None
			else:
				metabolite_names = list(mat[read_model_name][tried_key])
			model_dict["metabolite_names"] = metabolite_names
			
			tried_key = model_comp["metabolite_comp"]
			if tried_key == None:
				metabolite_comp = None
			else:
				metabolite_comp = list(mat[read_model_name][tried_key])
			model_dict["metabolite_comp"] = metabolite_comp
		
		except KeyError:
			print(f"The model uses a different key for {tried_key}")
			print(f"possible keys: {mat[read_model_name].keys()}")
		self.create_model_from_components(new_model_name, model_dict["S"], model_dict["lb"], model_dict["ub"],
		                                  model_dict["b"], model_dict["reaction_names"], model_dict["metabolite_names"],
		                                  met_comp=model_dict["metabolite_comp"], grRules=model_dict["grRules"])
		return model_dict
	
	def load_recon1_mat_model(self, filename, new_model_name=None):
		# this specifically loads the zeromodel from []
		if new_model_name is None:
			new_model_name = "Model"
		mat = loadmat(filename)
		S = sp.csr_matrix(mat["model"]["S"]).astype(np.float)
		S = S.toarray()
		lb = mat["model"]["lb"].astype(np.float)
		ub = mat["model"]["ub"].astype(np.float)
		b = mat["model"]["b"].astype(np.float)
		reaction_names = list(mat["model"]["rxns"])
		metabolite_names = list(mat["model"]["mets"])
		gene_names = list(mat["model"]["genes"])
		for i in mat["model"].keys():
			if np.size(mat["model"][i]) == len(reaction_names):
				print(i, "rxn")
			if np.size(mat["model"][i]) == len(metabolite_names):
				print(i, "met")
			if np.size(mat["model"][i]) == len(gene_names):
				print(i, "genes")
		
		additional_rxn_info_dict = {}
		
		additional_rxn_info_dict["grRules"] = list(mat["model"]["grRules"])
		additional_rxn_info_dict["subSystems"] = list(mat["model"]["subSystems"])
		additional_rxn_info_dict["rxnReferences"] = list(mat["model"]["rxnReferences"])
		additional_rxn_info_dict["confidenceScores"] = list(mat["model"]["confidenceScores"])
		additional_rxn_info_dict["rxnReferences"] = list(mat["model"]["rxnReferences"])
		additional_rxn_info_dict["rxnECNumbers"] = list(mat["model"]["rxnECNumbers"])
		additional_rxn_info_dict["rxnNotes"] = list(mat["model"]["rxnNotes"])
		additional_rxn_info_dict["rxnNames"] = list(mat["model"]["rxnNames"])
		
		additional_met_info_dict = {}
		additional_met_info_dict["metFormulas"] = list(mat["model"]["metFormulas"])
		additional_met_info_dict["metCharge"] = list(mat["model"]["metCharge"])
		additional_met_info_dict["metNames"] = list(mat["model"]["metNames"])
		additional_met_info_dict["metChEBIID"] = list(mat["model"]["metChEBIID"])
		additional_met_info_dict["metKEGGID"] = list(mat["model"]["metKEGGID"])
		additional_met_info_dict["metPubChemID"] = list(mat["model"]["metPubChemID"])
		additional_met_info_dict["metInChIString"] = list(mat["model"]["metInChIString"])
		
		self.create_model_from_components(new_model_name, S, lb, ub,
		                                  b, reaction_names, metabolite_names,
		                                  additional_rxn_info_dict=additional_rxn_info_dict,
		                                  additional_met_info_dict=additional_met_info_dict)
	
	def load_recon2_2_mat_model(self, filename, new_model_name=None):
		# this specifically loads the zeromodel from []
		if new_model_name is None:
			new_model_name = "Model"
		mat = loadmat(filename)
		print(mat.keys())
		print(mat["MODEL1603150001"].keys())
		print(mat["MODEL1603150001"]["grRules"])
		S = sp.csr_matrix(mat["MODEL1603150001"]["S"]).astype(np.float)
		S = S.toarray()
		lb = mat["model"]["lb"].astype(np.float)
		ub = mat["model"]["ub"].astype(np.float)
		b = mat["model"]["b"].astype(np.float)
		reaction_names = list(mat["model"]["rxns"])
		metabolite_names = list(mat["model"]["mets"])
		gene_names = list(mat["model"]["genes"])
		for i in mat["model"].keys():
			if np.size(mat["model"][i]) == len(reaction_names):
				print(i, "rxn")
			if np.size(mat["model"][i]) == len(metabolite_names):
				print(i, "met")
			if np.size(mat["model"][i]) == len(gene_names):
				print(i, "genes")
		
		additional_rxn_info_dict = {}
		
		additional_rxn_info_dict["grRules"] = list(mat["model"]["grRules"])
		additional_rxn_info_dict["subSystems"] = list(mat["model"]["subSystems"])
		additional_rxn_info_dict["rxnReferences"] = list(mat["model"]["rxnReferences"])
		additional_rxn_info_dict["confidenceScores"] = list(mat["model"]["confidenceScores"])
		additional_rxn_info_dict["rxnReferences"] = list(mat["model"]["rxnReferences"])
		additional_rxn_info_dict["rxnECNumbers"] = list(mat["model"]["rxnECNumbers"])
		additional_rxn_info_dict["rxnNotes"] = list(mat["model"]["rxnNotes"])
		additional_rxn_info_dict["rxnNames"] = list(mat["model"]["rxnNames"])
		
		additional_met_info_dict = {}
		additional_met_info_dict["metFormulas"] = list(mat["model"]["metFormulas"])
		additional_met_info_dict["metCharge"] = list(mat["model"]["metCharge"])
		additional_met_info_dict["metNames"] = list(mat["model"]["metNames"])
		additional_met_info_dict["metChEBIID"] = list(mat["model"]["metChEBIID"])
		additional_met_info_dict["metKEGGID"] = list(mat["model"]["metKEGGID"])
		additional_met_info_dict["metPubChemID"] = list(mat["model"]["metPubChemID"])
		additional_met_info_dict["metInChIString"] = list(mat["model"]["metInChIString"])
		
		self.create_model_from_components(new_model_name, S, lb, ub,
		                                  b, reaction_names, metabolite_names,
		                                  additional_rxn_info_dict=additional_rxn_info_dict,
		                                  additional_met_info_dict=additional_met_info_dict)
	
	def load_bigg_json_model(self, filename, new_model_name=None, model_header=None, rxn_comp=None, met_comp=None,
	                         gene_comp=None, keep_extra=True):
		# if model_name == "" the program guesses what the model key is by looking at the size of all dict entries
		# the largest entry is likely to be the one that holds the model
		# note sys.getsizeof() is insufficient as the model is often a dictionary which references large objects
		# and sys.getsizeof() will omit those referenced values
		# model_comp: dict which contains the names for each model component with the key being the "standard name"
		#   standard names: Stochiometric matrix = "S", lower bounds = "lb", upper bounds = "ub",
		#   rhs of Sv = b constraint = "b", dictionary of pinched reactions = "pinched_reactions",
		#   list of reaction names = "reaction_names", gene reaction rule list = "grRules",
		#   list of gene names = "genes", list of metabolite names = "metabolite_names"
		#   list of metabolite compositions = "metabolite_comp"
		#   In the case where an optional model component is not expected to be found in the model being loaded
		#   the key should remain unchanged, but the corresponding value should be None
		#   This will be assumed to be the case for pinched_reactions
		self.model_dict = {}
		if new_model_name is None:
			self.model_dict["model_name"] = "model"
		else:
			self.model_dict["model_name"] = new_model_name
		
		self.model_dict["gurobi_token"] = None
		
		if model_header is None:
			model_header = {"reactions": "reactions", "metabolites": "metabolites", "genes": "genes"}
		if rxn_comp is None:
			rxn_comp = {"rxn_name": "id", "rxn_metabolites": "metabolites", "lb": "lower_bound", "ub": "upper_bound",
			            "pinched_reactions": None, "subsystem": "subsystem", "grRules": "gene_reaction_rule"}
		if met_comp is None:
			met_comp = {"met_name": "id", "b": "b", "met_composition": "formula", "met_compartment": "compartment"}
		if gene_comp is None:
			gene_comp = {"gene_name": "id", "symbol": "name"}
		with open(filename, ) as json_in:
			model = json.load(json_in)
		self.model_dict["rxn_dict"] = {}
		for i in model[model_header["reactions"]]:
			self.model_dict["rxn_dict"][i[rxn_comp["rxn_name"]]] = {}
			for key in rxn_comp.keys():
				if rxn_comp[key] in i and key != "rxn_name":
					self.model_dict["rxn_dict"][i[rxn_comp["rxn_name"]]][key] = i[rxn_comp[key]]
			extra = {}
			for key_orig in i.keys():
				if key_orig not in [rxn_comp[key_new] for key_new in
				                    self.model_dict["rxn_dict"][i[rxn_comp["rxn_name"]]].keys()] and key_orig != \
						rxn_comp["rxn_name"]:
					extra[key_orig] = i[key_orig]
			if keep_extra:
				self.model_dict["rxn_dict"][i[rxn_comp["rxn_name"]]]["extra"] = extra
		self.model_dict["met_dict"] = {}
		for i in model[model_header["metabolites"]]:
			self.model_dict["met_dict"][i[met_comp["met_name"]]] = {}
			for key in met_comp.keys():
				if met_comp[key] in i and key != "met_name":
					self.model_dict["met_dict"][i[met_comp["met_name"]]][key] = i[met_comp[key]]
			extra = {}
			for key_orig in i.keys():
				if key_orig not in [met_comp[key_new] for key_new in
				                    self.model_dict["met_dict"][i[met_comp["met_name"]]].keys()] and key_orig != \
						met_comp[
							"met_name"]:
					extra[key_orig] = i[key_orig]
			if keep_extra:
				self.model_dict["met_dict"][i[met_comp["met_name"]]]["extra"] = extra
		
		self.model_dict["gene_dict"] = {}
		for i in model[model_header["genes"]]:
			self.model_dict["gene_dict"][i[gene_comp["gene_name"]]] = {}
			for key in gene_comp.keys():
				if gene_comp[key] in i and key != "gene_name":
					self.model_dict["gene_dict"][i[gene_comp["gene_name"]]][key] = i[gene_comp[key]]
			extra = {}
			for key_orig in i.keys():
				if key_orig not in [gene_comp[key_new] for key_new in
				                    self.model_dict["gene_dict"][i[gene_comp["gene_name"]]].keys()] and key_orig != \
						gene_comp[
							"gene_name"]:
					extra[key_orig] = i[key_orig]
			if keep_extra:
				self.model_dict["gene_dict"][i[gene_comp["gene_name"]]]["extra"] = extra
	
	def load_fast_key_json_model(self, filepath):
		with open(filepath, "r") as infile:
			self.model_dict = json.load(infile)
	
	def save_model_as_fast_key_json(self, filepath):
		with open(filepath / f"{self.model_dict['model_name']}.json", "w") as outfile:
			json.dump(self.model_dict, outfile)
	
	def update_reaction_bounds(self, name, new_lb, new_ub):
		if not new_lb == "keep":
			self.model_dict["rxn_dict"][name]["lb"] = new_lb
		if not new_ub == "keep":
			self.model_dict["rxn_dict"][name]["ub"] = new_ub
		if self.model_dict["rxn_dict"][name]["ub"] < self.model_dict["rxn_dict"][name]["lb"]:
			print(new_ub, new_lb)
			print(f"WARNING! Upper bound smaller than lower bound at {name}")
			raise ValueError
	
	def dicts_to_mats(self, alphabetize=True):
		# if omit zero flux is not positive all fluxes are used, if it is positive then if a lb and ub are within
		# that distance the reaction will not be included
		# lists hold teh order of the final metabolites and reactions
		# sometimes this will matter, other times it will not
		# it is better to make it unimportant when possible
		rxn_list = []
		for i in self.model_dict["rxn_dict"].keys():
			rxn_list.append(i)
		met_list = list(self.model_dict["met_dict"].keys())
		if alphabetize:
			rxn_list.sort()
			met_list.sort()
		lb = np.empty((len(rxn_list)))
		ub = np.empty((len(rxn_list)))
		S = np.zeros((len(met_list), len(rxn_list)))
		for rxn_ind in range(len(rxn_list)):
			lb[rxn_ind] = self.model_dict["rxn_dict"][rxn_list[rxn_ind]]['lb']
			ub[rxn_ind] = self.model_dict["rxn_dict"][rxn_list[rxn_ind]]['ub']
			for met in self.model_dict["rxn_dict"][rxn_list[rxn_ind]]['rxn_metabolites'].keys():
				S[met_list.index(met), rxn_ind] = self.model_dict["rxn_dict"][rxn_list[rxn_ind]]['rxn_metabolites'][met]
		b = np.empty((len(met_list)))
		for met_ind in range(len(met_list)):
			if "b" in self.model_dict["met_dict"][met_list[met_ind]]:
				b[met_ind] = self.model_dict["met_dict"][met_list[met_ind]]['b']
			else:
				b[met_ind] = 0
		return lb, ub, S, b, rxn_list, met_list
	
	def get_rxn_comp_flattened(self, comp_name, alphabetize_by_rxn_name=True):
		rxn_list = list(self.rxn_dict.keys())
		if alphabetize_by_rxn_name:
			rxn_list.sort()
		com_list = []
		for rxn_name in rxn_list:
			com_list.append(self.rxn_dict[rxn_name][comp_name])
		return com_list
	
	def limit_flux_bounds(self, absolute_max_value):
		for rxn in self.rxn_dict.keys():
			self.update_reaction_bounds(rxn, -1 * absolute_max_value, absolute_max_value)
	
	def add_reaction(self, rxn_name, S_dict, lb, ub, grRule=None):
		# S_dict must be a dict with the metabolites appearing in the reaction as keys and their coefficients as values
		# Assumes all metabolites of reaction are in model
		if grRule == None:
			print("Gene Reaction Rule should be included, currently assuming no gene association with this reaction")
			self.rxn_dict[rxn_name] = {"S": S_dict, "lb": lb, "ub": ub, "grRule": "[]"}
		else:
			self.rxn_dict[rxn_name] = {"S": S_dict, "lb": lb, "ub": ub, "grRule": grRule}
	
	def convert_scRNAseq_to_RAS_matrix(self, scRNAseq_matrix, gene_name_list):
		# gene_name_list should correspond to the row/column in the scRNAseq_matrix
		RAS_matrix = np.empty((len(self.reaction_names), np.shape(scRNAseq_matrix)[1]))
		for ci in range(np.shape(scRNAseq_matrix)[1]):
			gene_value_dict = dict(zip(gene_name_list, scRNAseq_matrix[:, ci]))
			for ri in range(len(self.reaction_names)):
				# Fixes notation error in recon2.2 model where HGNC:987 and HGNC:2898 appear as HGNC:HGNC:987 and
				# HGNC:HGNC:2898 (respectively)
				print(self.grRules[ri])
				print(self.reaction_names[ri])
				grRule_fixed = self.grRules[ri].replace("HGNC:HGNC:", "HGNC:")
				if self.reaction_names[ri] == "ACCOAC":
					print(grRule_fixed)
					print(self.get_reaction_activity_scores(grRule_fixed, gene_value_dict))
					input()
				RAS_matrix[ri, ci] = self.get_reaction_activity_scores(grRule_fixed, gene_value_dict)
				print(ri, ci)
		return RAS_matrix
	
	def get_reaction_activity_scores(self, string_to_handle, gene_value_dictionary):
		# additional arguements will be needed for different notation systems
		# gene value dictionary should be in the form "gene name" and "abundance"
		if string_to_handle == "[]":
			return 0
		string_to_evaluate = string_to_handle.replace("(", "{").replace(")", "}")
		vart = True
		while vart:
			m = re.search('{([^{}]+)}', string_to_evaluate)
			if m is not None:
				if "or" in m.group(0):
					new_str = m.group(0).replace("or", "+").replace("{", "(").replace("}", ")")
				if "and" in m.group(0):
					new_str = m.group(0).replace("and", ",").replace("{", "min(").replace("}", ")")
				string_to_evaluate = string_to_evaluate.replace(m.group(0), new_str)
			else:
				if "or" in string_to_evaluate:
					string_to_evaluate = "(" + string_to_evaluate.replace("or", "+") + ")"
				if "and" in string_to_evaluate:
					string_to_evaluate = "min(" + string_to_evaluate.replace("and", ",") + ")"
				vart = False
		gene_list = list(gene_value_dictionary.keys())
		
		gene_list.sort(reverse=True, key=size_order_key)
		for i in gene_list:
			string_to_evaluate = string_to_evaluate.replace(i, str(gene_value_dictionary[i]))
		return eval(string_to_evaluate)
	
	def flux_balance_test(self, point, cons_tol):
		# remove false
		if np.sum(np.matmul(self.S.toarray(), np.transpose(point)) - self.b) > -cons_tol and np.sum(
				np.matmul(self.S.toarray(), np.transpose(point)) - self.b) < cons_tol and False:
			print(f"Error of point: {np.abs(np.sum(np.matmul(self.S.toarray(), np.transpose(point))) - self.b)}")
	
	def normalize(self, vec):
		norm = np.linalg.norm(vec)
		if norm == 0:
			return vec
		else:
			return vec / norm
	
	def get_grRule_list(self, alphabetize=True):
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats(alphabetize=alphabetize)
		gr_Rule_list = []
		for i in rxn_list:
			if len(self.model_dict["rxn_dict"][i]["grRules"]) == 0:
				gr_Rule_list.append("")
			else:
				gr_Rule_list.append(self.model_dict["rxn_dict"][i]["grRules"])
		return gr_Rule_list
	
	def update_gene_dict(self):
		gene_list = self.get_used_gene_list()
		for gene in gene_list:
			if gene not in self.model_dict["gene_dict"].keys():
				print(gene)
				degen = 0
				for gene_2 in gene_list:
					if gene.split(".")[0] == gene_2.split(".")[0]:
						degen += 1
				self.model_dict["gene_dict"][gene] = {"entrez": gene.split(".")[0], "degen": degen}
	
	def get_used_gene_list(self):
		gene_list = []
		gr_Rule_list = self.get_grRule_list()
		for i in gr_Rule_list:
			for id in re.findall("\d+\.\d+", i):
				if id not in gene_list:
					gene_list.append(id)
		return gene_list
	
	def grRule_to_RAS(self, grRule, RNA_seq_dict):
		pass
	
	def find_objective_value(self, c):
		if self.model_dict['gurobi_token'] is None:
			env = gp.Env(empty=True)
		else:
			env = self.model_dict['gurobi_token']
		# env.start()
		model = gp.Model("opt_model_" + self.model_name, env=env)
		model.Params.LogToConsole = 0
		react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=self.lb,
		                           ub=self.ub)
		obj = np.reshape(c, (1, -1))
		# print(obj.shape)
		# print(obj @ react_flux)
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		rhs = np.transpose(self.b)
		model.addConstr(self.S @ react_flux == rhs, name="c")
		model.optimize()
		if model.Status == 2:
			return model.objVal
		else:
			print(model.Status)
			return np.nan
	
	def fva(self, lb=None, ub=None, S=None, b=None, rxn_list=None, return_points=False):
		# Returns list of [min,max] elements
		# If reactions are edited this fva no longer applies
		# This is valid for removed reactions, but could lead to issues if 2 reactions are simply switched around
		if lb is None:
			lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		print(np.sum(lb),np.sum(ub),np.sum(S))
		print(rxn_list)
		#input()
		if self.model_dict["gurobi_token"] is not None:
			model = gp.Model("fva_model_" + self.model_dict["model_name"], env=self.model_dict["gurobi_token"])
		else:
			model = gp.Model("fva_model_" + self.model_dict["model_name"])
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		model.addConstr(S @ react_flux == b, name="c")
		model.Params.LogToConsole = 0
		# tag
		# making this 1e-8 causes it
		# this indicates the cause might just be due to little issues in the limits
		# set to 1e-9 fixes it
		model.Params.FeasibilityTol = 1e-9
		model.Params.NumericFocus = 3
		obj = np.zeros((1, S.shape[1]))
		variability_list = [{"lb": 0, "ub": 0} for i in range(len(lb))]
		if return_points:
			point_array = np.zeros((len(lb) * 2, len(lb)))
		for i in range(len(lb)):
			print(len(lb) - i)
			min_max = []
			obj = np.zeros(obj.shape)
			obj[0, i] = 1
			model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
			#if i == 639:
			#	model.write("test_init_1.mps")
			#	model.write("test_init_1.attr")
			#	model.write("test_init_1.prm")
			#	print(model.fingerprint)
			#	input()
			model.optimize()
			#if i == 639:
			#	model.write("test_fin_1.json")
			if return_points:
				point_array[2 * i] = react_flux.X
			min_max.append(model.objVal)
			model.reset()
			obj[0, i] = -1
			model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
			model.optimize()
			min_max.append(model.objVal * -1)
			if return_points:
				point_array[2 * i + 1] = react_flux.X
			# adding a little slack to the range seems to fix it
			variability_list[i]["ub"] = max(min_max) +1e-8
			variability_list[i]["lb"] = min(min_max) -1e-8
		if return_points:
			return dict(zip(rxn_list, variability_list)), point_array
		else:
			return dict(zip(rxn_list, variability_list))
	
	def test_feasibility(self, lb=None, ub=None, S=None, b=None, **kwargs):
		if lb is None:
			lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats(alphabetize=True)
		if self.model_dict["gurobi_token"] != None:
			model = gp.Model("fva_model_" + self.model_dict["model_name"], env=self.model_dict["gurobi_token"])
		else:
			model = gp.Model("fva_model_" + self.model_dict["model_name"])
		if "unique_bounds" in kwargs.keys():
			tflb = kwargs["unique_bounds"][0]
			tfub = kwargs["unique_bounds"][1]
		else:
			tflb = lb
			tfub = ub
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=tflb, ub=tfub)
		model.addConstr(S @ react_flux == b, name="c")
		obj = np.zeros((1, S.shape[1]))
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.LogToConsole = 0
		# model.Params.NumericFocus = 3
		# model.Params.ScaleFlag = 0
		# model.Params.FeasibilityTol = 1e-9
		model.optimize()
		return model.Status == 2
	
	def test_feasibility_pass_model(self, lb, ub, model, react_flux):
		react_flux.lb = lb
		react_flux.ub = ub
		obj = np.zeros((1, self.S.shape[1]))
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.LogToConsole = 0
		model.optimize()
		return model.Status == 2
	
	def initalize_model(self, **kwargs):
		# Performs FVA and updates lb and ub
		varibility_list = self.fva(**kwargs)
		
		for i in range(len(varibility_list)):
			self.update_reaction_bounds(i, varibility_list[i][0], varibility_list[i][1])
		
		if "shrink_S" in kwargs.keys():
			self.pinch_reactions(kwargs["shrink_S"])
	
	def pinch_reactions(self, distance_between_bounds, zero_out=False):
		# This is unstable and should only be used as one of the final steps of model refinement
		# Improvements could be made such as using a dictionary for the metabolite constant storage
		constant_var_index = np.where((self.ub - self.lb) <= distance_between_bounds)[0]
		variable_var_index = np.where((self.ub - self.lb) > distance_between_bounds)[0]
		
		contant_var_val = (self.ub[constant_var_index] + self.lb[constant_var_index]) / 2
		
		cS = self.S.toarray()[:, constant_var_index]
		self.b -= np.matmul(cS, contant_var_val)
		constant_var_index = constant_var_index.tolist()
		constant_var_index.sort(reverse=True)
		for i in constant_var_index:
			metabolite_stochio_storage = {}
			row_values = self.S[:, i].nonzero()[0]
			for j in row_values:
				metabolite_stochio_storage[self.metabolite_names[j]] = self.S[j, i]
			self.pinched_reactions[self.reaction_names[i]] = [(self.lb[i] + self.ub[i]) / 2, metabolite_stochio_storage]
			self.del_reaction(i)
		print(f"Restricted {len(constant_var_index)} of {len(constant_var_index) + len(variable_var_index)} fluxes")
	
	def get_exchanged_reaction_info(self, true_exchange_signifier=""):
		exch_rxn = []
		exch_met = []
		for i in self.model_dict["rxn_dict"].keys():
			if len(self.model_dict["rxn_dict"][i]["rxn_metabolites"]) == 1 and true_exchange_signifier in i:
				exch_rxn.append(i)
				exch_met.append(list(self.model_dict["rxn_dict"][i]["rxn_metabolites"].keys())[0])
		return dict(zip(exch_rxn, exch_met))
	
	def metabolite_info(self, name, flux_dtf=None):
		print(f"Information on {name}:")
		print(f"Molar mass = {np.round(self.get_met_mol_weight(name), 2)} g/mol")
		print(self.model_dict["met_dict"][name])
		print(f"The reactions in which {name} is produced or consumed are:")
		list_of_con = []
		for i in self.model_dict['rxn_dict'].keys():
			if name in self.model_dict['rxn_dict'][i]["rxn_metabolites"].keys():
				print(i)
				if flux_dtf is not None and i in flux_dtf.columns:
					print(
						f"Rxn name: {i}, average met produced: {self.model_dict['rxn_dict'][i]['rxn_metabolites'][name] * np.average(flux_dtf.loc[:, i].to_numpy())}")
					list_of_con.append(self.model_dict['rxn_dict'][i]['rxn_metabolites'][name] * np.average(
						flux_dtf.loc[:, i].to_numpy()))
	
	def fast_make_positive_HR(self, point):
		pass
	
	def generate_mat_wo_exchange(self, prune_specific=[]):
		# convert this back to dict for check to make sure you didnt mess with any reactions
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		rxn_array = np.array(rxn_list)
		iti = [not (np.size(np.nonzero(S[:, rxn_ind])) == 1 or rxn_list[rxn_ind] in prune_specific) for rxn_ind in
		       range(len(rxn_list))]
		internal_rxn = rxn_array[iti]
		internal_lb = lb[iti]
		internal_ub = ub[iti]
		internal_S = S[:, iti]
		return internal_lb, internal_ub, internal_S, b, internal_rxn, met_list
	
	def pos_S(self, inverse_tag="_inverted"):
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		print(rxn_list)
		for i in range(len(rxn_list)):
			if "_inverted" in rxn_list[i]:
				positive_ent = copy.deepcopy(rxn_list[i]).replace("_inverted", "")
				if positive_ent in rxn_list and rxn_list[i - 1] != positive_ent:
					print(positive_ent)
					print(rxn_list[i])
					print(rxn_list[i - 1])
					print(rxn_list[rxn_list.index(positive_ent)])
					input()
		print(np.array(rxn_list))
		input()
	
	def generate_null_space(self, internal_S):
		for i in range(9):
			rounded = np.sum(np.round(internal_S * 10 ** i))
			if rounded == np.sum(internal_S * 10 ** i):
				internal_Sy = sy.matrices.Matrix.zeros(np.shape(internal_S)[0], np.shape(internal_S)[1])
				for row_ind in range(np.shape(internal_S)[0]):
					for col_ind in range(np.shape(internal_S)[1]):
						if internal_S[row_ind, col_ind] != 0:
							internal_Sy[row_ind, col_ind] = sy.S(
								int(np.round(internal_S[row_ind, col_ind] * 10 ** i))) / (sy.S(10) ** sy.S(i))
				break
		Ny = np.hstack([np.array(i) for i in internal_Sy.nullspace()])
		print(Ny)
		return Ny
	
	def generate_null_space_fast(self, internal_S):
		Ny = sla.null_space(internal_S)
		# Ny = np.where(np.abs(Ny) > 1e-10,Ny,np.zeros_like(Ny))
		# print(NzNy)
		flat_array = np.log10(np.abs(Ny.flatten()) + 1e-16)
		return Ny
	
	def generate_therm_model(self, trans_N, point):
		if self.model_dict['gurobi_token'] != None:
			model = gp.Model("therm_model_" + self.model_dict['model_name'], env=self.model_dict['gurobi_token'])
		else:
			model = gp.Model("therm_model_" + self.model_dict['model_name'])
		
		lb_mu = np.where(point > 0, -np.inf, u_bound_tol)
		lb_mu = np.where(point == 0, -np.inf, lb_mu)
		ub_mu = np.where(point > 0, -1 * u_bound_tol, np.inf)
		ub_mu = np.where(point == 0, np.inf, ub_mu)
		# lb_negative_mu = np.where(point >= 0, -1000, tol)
		# ub_positive_mu = np.where(point <= 0, 1000, -1*tol)
		ts1 = trans_N.shape[1]
		mu = model.addMVar(shape=ts1, vtype=GRB.CONTINUOUS, name="mu", lb=lb_mu, ub=ub_mu)
		model.Params.LogToConsole = 0
		model.addConstr(trans_N @ mu == 0, name="c")
		# model.Params.NumericFocus = 3
		model.Params.FeasibilityTol = 1e-9
		obj = np.zeros((1, ts1))
		model.setObjective(obj @ mu, GRB.MAXIMIZE)
		model.optimize()
		return model.Status
	
	def generate_therm_model_new(self, trans_N, point, u_bound_tol=1e-8, slack=1e-6*0):
		if self.model_dict['gurobi_token'] != None:
			model = gp.Model("therm_model_" + self.model_dict['model_name'], env=self.model_dict['gurobi_token'])
		else:
			model = gp.Model("therm_model_" + self.model_dict['model_name'])
		point[np.abs(point) < slack] = 0
		lb_mu = np.where(point > 0, -np.inf, u_bound_tol)
		lb_mu = np.where(point == 0, -np.inf, lb_mu)
		ub_mu = np.where(point > 0, -u_bound_tol, np.inf)
		ub_mu = np.where(point == 0, np.inf, ub_mu)
		# lb_negative_mu = np.where(point >= 0, -1000, tol)
		# ub_positive_mu = np.where(point <= 0, 1000, -1*tol)
		ts1 = trans_N.shape[1]
		mu = model.addMVar(shape=ts1, vtype=GRB.CONTINUOUS, name="mu", lb=lb_mu, ub=ub_mu)
		model.Params.LogToConsole = 0
		model.addConstr(trans_N @ mu == 0, name="c")
		model.Params.NumericFocus = 3
		model.Params.FeasibilityTol = 1e-8
		obj = np.zeros((1, ts1))
		model.setObjective(obj @ mu, GRB.MAXIMIZE)
		model.optimize()
		return model.Status
	
	def test_thermo_feasibilty(self, N, therm_model):
		pass
	
	def element_edge_vert_generator(self, rxn, rxn_dtf=None, element="compound", break_element_list=[]):
		edge_list = []
		weight_list = []
		print(self.model_dict["rxn_dict"][rxn].keys())
		for met in self.model_dict["rxn_dict"][rxn]["rxn_metabolites"].keys():
			
			if element == "compound":
				element_per_stochio = 1
			else:
				element_per_stochio = self.get_met_comp_dict(met)[element]
			if rxn_dtf is not None:
				avg_flux = np.average(rxn_dtf.loc[:, rxn].to_numpy())
			else:
				avg_flux = 1
			print(rxn, avg_flux)
			# input()
			# obtains the amount of element going to a specific product in the reaction
			
			if element_per_stochio * self.model_dict["rxn_dict"][rxn]["rxn_metabolites"][met] * avg_flux < 0:
				if avg_flux < 0:
					if met in break_element_list:
						edge_list.append((met + "_for_" + rxn, "rev-" + rxn + "_RXN"))
					else:
						edge_list.append((met, "rev-" + rxn + "_RXN"))
				else:
					if met in break_element_list:
						edge_list.append((met + "_for_" + rxn, rxn + "_RXN"))
					else:
						edge_list.append((met, rxn + "_RXN"))
				weight_list.append(
					-1 * element_per_stochio * self.model_dict["rxn_dict"][rxn]["rxn_metabolites"][met] * avg_flux)
			elif element_per_stochio * self.model_dict["rxn_dict"][rxn]["rxn_metabolites"][met] * avg_flux > 0:
				if avg_flux < 0:
					if met in break_element_list:
						edge_list.append(("rev-" + rxn + "_RXN", met + "_for_" + rxn))
					else:
						edge_list.append(("rev-" + rxn + "_RXN", met))
				else:
					if met in break_element_list:
						edge_list.append((rxn + "_RXN", met + "_for_" + rxn))
					else:
						edge_list.append((rxn + "_RXN", met))
				weight_list.append(
					element_per_stochio * self.model_dict["rxn_dict"][rxn]["rxn_metabolites"][met] * avg_flux)
		return edge_list, weight_list
	
	# if self.rxn_dict[rxn]["S"][met] < 0:
	#	reactants[met] = np.average(rxn_dtf.loc[:,rxn].to_numpy())*-1*self.met_dict[met]['comp']
	
	def generate_rxn_met_graph(self, edge_list, weight_list, color_edge_of_vert=None, default_edge_color="grey",
	                           vert_leak_tol=None, max_vec_width=5, weight_display="raw"):
		g = ig.Graph(directed=True)
		edge_color_list = []
		vert_list = []
		for i in edge_list:
			if i[0] not in vert_list:
				vert_list.append(i[0])
			if i[1] not in vert_list:
				vert_list.append(i[1])
			if color_edge_of_vert is not None:
				if i[0].replace("_RXN", "").replace("rev-", "") in color_edge_of_vert.keys():
					edge_color_list.append(color_edge_of_vert[i[0].replace("_RXN", "").replace("rev-", "")])
				elif i[1].replace("_RXN", "").replace("rev-", "") in color_edge_of_vert.keys():
					edge_color_list.append(color_edge_of_vert[i[1].replace("_RXN", "").replace("rev-", "")])
				else:
					edge_color_list.append(default_edge_color)
			else:
				edge_color_list.append(default_edge_color)
		
		vert_leak_list_color = []
		for i in vert_list:
			vert_leak = 0
			if vert_leak_tol is not None:
				for j in range(len(edge_list)):
					if i == edge_list[j][0]:
						vert_leak -= weight_list[j]
					elif i == edge_list[j][1]:
						vert_leak += weight_list[j]
				if abs(vert_leak) > vert_leak_tol:
					if vert_leak < 0:
						vert_leak_list_color.append("red")
					else:
						vert_leak_list_color.append("blue")
				else:
					vert_leak_list_color.append("grey")
			else:
				vert_leak_list_color.append("grey")
		
		g.add_vertices(vert_list)
		g.vs["is_RXN"] = ["_RXN" not in i for i in vert_list]
		g.add_edges(edge_list)
		g.es["weight"] = weight_list
		g.es["color"] = edge_color_list
		g.vs["color"] = vert_leak_list_color
		
		print(edge_color_list)
		
		visual_style = {}
		visual_style["vertex_size"] = [float(i) / 5 + .01 for i in g.vs["is_RXN"]]
		visual_style["vertex_color"] = g.vs["color"]
		visual_style["edge_color"] = g.es["color"]
		visual_style["vertex_label"] = [i.replace("_RXN", "") for i in g.vs["name"]]
		visual_style["edge_width"] = [i / max(list(g.es["weight"])) * max_vec_width for i in g.es["weight"]]
		if weight_display == "raw":
			visual_style["edge_label"] = [f"{np.round(i, 3)}" for i in g.es["weight"]]
		elif weight_display == "rel":
			visual_style["edge_label"] = [f"{np.round(i / max(list(g.es['weight'])), 3)}" for i in g.es["weight"]]
		# visual_style["edge_color"] = [color_list[int(i)] for i in g.es["met_created"]]
		visual_style["layout"] = g.layout("kk")
		visual_style["bbox"] = (100000, 100000)
		visual_style["margin"] = [0, 0, 0, 0]
		fig, ax = plt.subplots()
		ig.plot(g, target=ax, **visual_style)
		plt.show()
	
	def reaction_info(self, name):
		# Reaction name can be the name of the reaction or the index location of the reaction
		print(f"Information on {name}:")
		print(
			f"Lower Bound: {self.model_dict['rxn_dict'][name]['lb']}, Upper Bound {self.model_dict['rxn_dict'][name]['ub']}")
		print(f"The relevant chemical species participating in the reaction and their coefficients are:")
		print(self.model_dict['rxn_dict'][name]["rxn_metabolites"])
		print(f"The reaction gene rule is:\n{self.model_dict['rxn_dict'][name]['grRules']}")
	
	def comparable_reaction_report(self, name):
		# Reaction name can be the name of the reaction or the index location of the reaction
		# Returns 2d list, the ith entry being ["metabolite_name_i",reaction_coefficient_for metabolite_i]
		if isinstance(name, str):
			name = self.reaction_names.index(name)
		metabolite_names = []
		row_values = self.S[:, name].nonzero()[0]
		for i in row_values:
			metabolite_names.append([self.metabolite_names[i], self.S[i, name]])
		return metabolite_names
	
	def objective_function_generator(self, list_of_reaction_names, list_of_weights):
		c = np.zeros(np.shape(self.S)[1])
		for i in range(len(list_of_reaction_names)):
			if isinstance(list_of_reaction_names[i], str):
				c[self.reaction_names.index(list_of_reaction_names[i])] = list_of_weights[i]
			else:
				c[list_of_reaction_names[i]] = list_of_weights[i]
		return c
	
	def lin_dep_np_ctp(self, old_on_basis, origin, accepted_points, new_point, old_rank,
	                   unique_vector_threshold_exponent=9):
		# Description:
		#   Takes in a set of basis vectors and a new proposed basis vector and either accepts or rejects
		#       the proposed vector based on whether the addition of the new vector increases the dimensionality of the
		#       basis set
		# Input:
		#   old_basis: 2d numpy array containing points which define a basis wrt the origin point
		#   origin: vectors are defined as the difference between points of interest and the origin
		#       Value should be within feasible space, either the first generated vertex or the center of the point
		#       matrix are valid approaches, we will use the first generated vertex
		#   new_point: proposed point which creates a vector wrt the origin
		# Output:
		#   either the old basis and old rank are returned or, if the proposed vector added to the basis creates
		#       a new larger basis, that basis and rank are returned
		
		# find the dimensionality of the original basis set
		# true_old_rank = np.linalg.matrix_rank(old_basis - origin)
		
		# tests to see if the proposed point is within tolerance of being the origin point and rejects it if so.
		test_unique = np.unique(np.round(np.vstack((new_point, origin)), unique_vector_threshold_exponent), axis=0)
		if np.shape(test_unique)[0] == 1:
			return accepted_points, old_on_basis, old_rank
		
		# make new candidate basis vector
		diff = (new_point - origin)
		# potential_basic_vec = diff
		potential_basic_vec = diff / np.sqrt(diff.dot(diff))
		
		# make new basis set
		if old_rank == 0:
			new_on_basis = np.zeros((1, np.size(origin)))
			new_on_basis[0] = potential_basic_vec
		else:
			new_on_basis = np.vstack((old_on_basis, potential_basic_vec))
		
		# find the dimensionality of the new basis set
		new_rank = np.linalg.matrix_rank((new_on_basis[:, :-1]).T)
		
		# if len(np.shape(old_on_basis)) > 1:
		#	print(np.linalg.svd(old_on_basis)[1][-5:])
		#	print(np.linalg.svd(new_on_basis)[1][-5:])
		
		# if the new basis has a higher dimensionality, the added vector is indeed spanning and is kept
		if new_rank <= old_rank:
			return accepted_points, old_on_basis, old_rank
		else:
			if np.size(accepted_points) == 0:
				accepted_points = np.reshape(new_point, (1, -1))
			else:
				accepted_points = np.vstack([accepted_points, new_point])
			return accepted_points, new_on_basis, new_rank
	
	def min_flux(self):
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		model = gp.Model("warmup_model_" + self.model_name, env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		obj = -1 * np.ones((1, np.size(ub)))
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.optimize()
		min_flux = model.objVal * -1
		
		return min_flux
	
	def generate_warmup_gb_sf(self, max_search_number, sampling_levels, points_save_path="", model_save_path="",
	                          NS=None,
	                          required_rank=np.inf, rxn_penalty_vec=None, unique_vector_threshold_exponent=9, **kwargs):
		# rel_distance_to_min_flux
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		
		# unique_vector_threshold_exponent gives the number of decimal places out the algorithm searches when
		# removing identical vectors
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		
		# If rxn_penalty_vec is None a standard 1 norm is used
		
		# Currently this runs into issues if total flux can be zero as this defies how the 1 norm is handled
		
		# A pass to remove vectors which are identical down to a machine precision might be useful
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		
		if rxn_penalty_vec is None:
			rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
		else:
			rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
		S = np.vstack((S, rxn_penalty_vec))
		S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
		lb = np.hstack((lb, 0))
		ub = np.hstack((ub, np.inf))
		b = np.hstack((b, 0))
		S[-1, -1] = -1
		rxn_list.append("total_flux")
		
		model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		obj = np.zeros((1, np.size(ub)))
		obj[0, -1] = -1
		print(obj)
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.NumericFocus = 3
		model.optimize()
		print(S, lb, ub)
		print("done")
		min_flux = model.objVal * -1
		print(min_flux)
		# input()
		if NS is None:
			NS = sc.linalg.null_space(S)
		
		total_points = np.array([])
		origin = np.array([])
		
		# only 1 origin point should be made
		origin_created = False
		
		for sampling_level in sampling_levels:
			# Defines origin and points as empty
			
			accepted_points = np.array([])
			ortho_norm_basis = np.array([])
			
			ub[-1] = sampling_level * min_flux
			# if you want purely slices use
			# lb[-1] = sampling_level * min_flux
			model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
			
			# Suppress printing to console
			model.Params.LogToConsole = 0
			
			# Initializes number of subsequent vectors which have been rejected
			search_number = 0
			
			#
			rank = 0
			
			print(f"SAMPLING LEVEL {sampling_level}")
			
			while (search_number < max_search_number and rank < required_rank):
				# print(search_number,rank)
				# Finds random vertex
				while model.Status != 2:
					model.reset()
					obj = (2 * np.random.random_sample(S.shape[1]) - 1)
					# for variable objective total flux
					# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
					# model.Params.NumericFocus = 3
					# model.Params.FeasibilityTol = 1e-7
					model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
					model.optimize()
				
				# The below ensures the points fall inside the bounds
				# Does seem to potentially break constraints
				# entries of points are just being moved, the points should be scaled
				# if this becomes a problem, define scaling function
				# point_to_use = react_flux.X
				
				# if np.max(np.abs(np.matmul(S, point_to_use))) > Sv_tol:
				
				point_to_use = react_flux.X
				point_to_use = np.minimum(point_to_use, ub)
				point_to_use = np.maximum(point_to_use, lb)
				point_to_use = np.matmul(NS, np.matmul(np.transpose(NS), point_to_use))
				
				# print(point_to_use)
				# print(np.max(point_to_use[:-1]*rxn_penalty_vec[0]))
				# print(rxn_list[np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				# print(rxn_penalty_vec[0,np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				# input()
				# print(point_to_use[-1])
				# Saves rank to later check for rank change
				old_rank = rank
				if (search_number % 50 == 0 and search_number != 0) or rank == 0:
					print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
				# The first point becomes an origin point
				# The rank is 0 with the origin point since no vectors have been defined
				if origin_created:
					accepted_points, ortho_norm_basis, rank = self.lin_dep_np_ctp(ortho_norm_basis, origin,
					                                                              accepted_points, point_to_use, rank,
					                                                              unique_vector_threshold_exponent=unique_vector_threshold_exponent)
				else:
					origin = point_to_use
					origin_created = True
				# Reset the model to discard the current solution
				model.reset()
				# print(rank)
				# Reset the search_number if a vector is accepted as basis
				# I removed "or np.size(accepted_points) == 0" here which might lead to issues
				if rank != old_rank:
					search_number = 0
				else:
					search_number += 1
			print(rank)
			if np.size(total_points) == 0:
				if np.size(accepted_points) != 0:
					total_points = np.vstack((origin, accepted_points))
				else:
					total_points = origin
			else:
				if np.size(accepted_points) != 0:
					total_points = np.vstack((total_points, accepted_points))
		print(np.average(total_points, axis=0))
		print(1)
		# if multiple 1 norms were used the resulting vectors will span the space, but not be a basis
		# this corrects for that
		new_ortho_norm_basis_sf = np.array([])
		new_accepted_points_sf = np.array([])
		rank_sf = 0
		for i in total_points[1:]:
			new_accepted_points_sf, new_ortho_norm_basis_sf, rank_sf = self.lin_dep_np_ctp(new_ortho_norm_basis_sf,
			                                                                               total_points[0],
			                                                                               new_accepted_points_sf, i,
			                                                                               rank_sf,
			                                                                               unique_vector_threshold_exponent=unique_vector_threshold_exponent)
		
		origin = total_points[0]
		verts = total_points[1:]
		np.random.shuffle(verts)
		
		new_ortho_norm_basis_r = np.array([])
		new_accepted_points_r = np.array([])
		rank_r = 0
		for i in verts:
			new_accepted_points_r, new_ortho_norm_basis_r, rank_r = self.lin_dep_np_ctp(new_ortho_norm_basis_r,
			                                                                            total_points[0],
			                                                                            new_accepted_points_r, i,
			                                                                            rank_r,
			                                                                            unique_vector_threshold_exponent=unique_vector_threshold_exponent)
		print(np.average(new_accepted_points_sf, axis=0))
		print(np.average(new_accepted_points_r, axis=0))
		print(rank_sf)
		print(rank_r)
		print("done")
		plt.hist(new_accepted_points_r[:, -1], alpha=0.5)
		plt.hist(new_accepted_points_sf[:, -1], alpha=0.5)
		plt.show()
		input()
		input()
		
		# This section removes any vectors which are approximately equivalent
		initial_size = np.shape(new_accepted_points)[0]
		# rounds copy of points and uses return_index in np.unique to grab unique original vectors
		# This lets two points be considered too similar without potential loss of machine precision of accepted point
		# This may not actally matter, but is not costly to implement
		new_accepted_points_rnd = np.round(new_accepted_points, unique_vector_threshold_exponent)
		unique_accepted_points_rnd, unique_index = np.unique(new_accepted_points_rnd, axis=0, return_index=True)
		new_accepted_points = new_accepted_points[unique_index]
		new_size = np.shape(new_accepted_points)[0]
		if initial_size != new_size:
			print(
				f"WARNING: There were {initial_size - new_size} repeated vectors in the basis. They have been removed")
			time.sleep(3)
		
		# returns set of points with top point being basis
		total_points = np.vstack((origin, new_accepted_points))
		print(np.average(total_points, axis=0))
		print(3)
		self.model_dict["total_flux_limits"] = [min_flux, max(sampling_levels) * min_flux]
		if points_save_path != "":
			np.save(points_save_path / ("warmup_" + self.model_dict["model_name"]), total_points)
		if model_save_path != "":
			self.save_model_as_fast_key_json(model_save_path)
		return total_points
	
	def generate_warmup_gb_r(self, max_search_number, sampling_levels, points_save_path="", model_save_path="",
	                         required_rank=np.inf, rxn_penalty_vec=None, unique_vector_threshold_exponent=9, **kwargs):
		# rel_distance_to_min_flux
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		
		# unique_vector_threshold_exponent gives the number of decimal places out the algorithm searches when
		# removing identical vectors
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		
		# If rxn_penalty_vec is None a standard 1 norm is used
		
		# Currently this runs into issues if total flux can be zero as this defies how the 1 norm is handled
		
		# A pass to remove vectors which are identical down to a machine precision might be useful
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		
		if rxn_penalty_vec is None:
			rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
		else:
			rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
		S = np.vstack((S, rxn_penalty_vec))
		S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
		lb = np.hstack((lb, 0))
		ub = np.hstack((ub, np.inf))
		b = np.hstack((b, 0))
		S[-1, -1] = -1
		rxn_list.append("total_flux")
		
		model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		obj = np.zeros((1, np.size(ub)))
		obj[0, -1] = -1
		print(obj)
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.NumericFocus = 3
		model.optimize()
		print(S, lb, ub)
		print("done")
		min_flux = model.objVal * -1
		print(min_flux)
		# input()
		N = sc.linalg.null_space(S)
		
		total_points = np.array([])
		origin = np.array([])
		
		# only 1 origin point should be made
		origin_created = False
		
		for sampling_level in sampling_levels:
			# Defines origin and points as empty
			
			accepted_points = np.array([])
			ortho_norm_basis = np.array([])
			
			ub[-1] = sampling_level * min_flux
			# if you want purely slices use
			# lb[-1] = sampling_level * min_flux
			model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
			
			# Suppress printing to console
			model.Params.LogToConsole = 0
			
			# Initializes number of subsequent vectors which have been rejected
			search_number = 0
			
			#
			rank = 0
			
			print(f"SAMPLING LEVEL {sampling_level}")
			
			while (search_number < max_search_number and rank < required_rank):
				# print(search_number,rank)
				# Finds random vertex
				while model.Status != 2:
					model.reset()
					obj = (2 * np.random.random_sample(S.shape[1]) - 1)
					# for variable objective total flux
					# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
					model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
					model.optimize()
				
				# The below ensures the points fall inside the bounds
				# Does seem to potentially break constraints
				# entries of points are just being moved, the points should be scaled
				# if this becomes a problem, define scaling function
				# point_to_use = react_flux.X
				point_to_use = np.matmul(N, np.matmul(np.transpose(N), react_flux.X))
				point_to_use = np.minimum(point_to_use, ub)
				point_to_use = np.maximum(point_to_use, lb)
				
				# print(point_to_use)
				# print(np.max(point_to_use[:-1]*rxn_penalty_vec[0]))
				# print(rxn_list[np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				# print(rxn_penalty_vec[0,np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				# input()
				# print(point_to_use[-1])
				# Saves rank to later check for rank change
				old_rank = rank
				if (search_number % 50 == 0 and search_number != 0) or rank == 0:
					print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
				# The first point becomes an origin point
				# The rank is 0 with the origin point since no vectors have been defined
				if origin_created:
					accepted_points, ortho_norm_basis, rank = self.lin_dep_np_ctp(ortho_norm_basis, origin,
					                                                              accepted_points, point_to_use, rank,
					                                                              unique_vector_threshold_exponent=unique_vector_threshold_exponent)
				else:
					origin = point_to_use
					origin_created = True
				# Reset the model to discard the current solution
				model.reset()
				# print(rank)
				# Reset the search_number if a vector is accepted as basis
				# I removed "or np.size(accepted_points) == 0" here which might lead to issues
				if rank != old_rank:
					search_number = 0
				else:
					search_number += 1
			print(rank)
			if np.size(total_points) == 0:
				if np.size(accepted_points) != 0:
					total_points = np.vstack((origin, accepted_points))
				else:
					total_points = origin
			else:
				if np.size(accepted_points) != 0:
					total_points = np.vstack((total_points, accepted_points))
		print(np.average(total_points, axis=0))
		print(1)
		# if multiple 1 norms were used the resulting vectors will span the space, but not be a basis
		# this corrects for that
		origin = total_points[0]
		verts = total_points[1:]
		np.random.shuffle(verts)
		new_ortho_norm_basis = np.array([])
		new_accepted_points = np.array([])
		rank = 0
		for i in verts:
			new_accepted_points, new_ortho_norm_basis, rank = self.lin_dep_np_ctp(new_ortho_norm_basis, origin,
			                                                                      new_accepted_points, i, rank,
			                                                                      unique_vector_threshold_exponent=unique_vector_threshold_exponent)
		print(np.average(new_accepted_points, axis=0))
		print(2)
		
		# This section removes any vectors which are approximately equivalent
		initial_size = np.shape(new_accepted_points)[0]
		# rounds copy of points and uses return_index in np.unique to grab unique original vectors
		# This lets two points be considered too similar without potential loss of machine precision of accepted point
		# This may not actally matter, but is not costly to implement
		new_accepted_points_rnd = np.round(new_accepted_points, unique_vector_threshold_exponent)
		unique_accepted_points_rnd, unique_index = np.unique(new_accepted_points_rnd, axis=0, return_index=True)
		new_accepted_points = new_accepted_points[unique_index]
		new_size = np.shape(new_accepted_points)[0]
		if initial_size != new_size:
			print(
				f"WARNING: There were {initial_size - new_size} repeated vectors in the basis. They have been removed")
			time.sleep(3)
		
		# returns set of points with top point being basis
		total_points = np.vstack((origin, new_accepted_points))
		print(np.average(total_points, axis=0))
		print(3)
		self.model_dict["total_flux_limits"] = [min_flux, max(sampling_levels) * min_flux]
		if points_save_path != "":
			np.save(points_save_path / ("warmup_" + self.model_dict["model_name"]), total_points)
		if model_save_path != "":
			self.save_model_as_fast_key_json(model_save_path)
		return total_points
	
	def generate_warmup_gb(self, max_search_number, sampling_levels, points_save_path="", model_save_path="",
	                       required_rank=np.inf, rxn_penalty_vec=None, unique_vector_threshold_exponent=9, **kwargs):
		# rel_distance_to_min_flux
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		
		# unique_vector_threshold_exponent gives the number of decimal places out the algorithm searches when
		# removing identical vectors
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		
		# If rxn_penalty_vec is None a standard 1 norm is used
		
		# Currently this runs into issues if total flux can be zero as this defies how the 1 norm is handled
		
		# A pass to remove vectors which are identical down to a machine precision might be useful
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		
		Sv_tol = 1e-9
		
		if rxn_penalty_vec is None:
			rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
		else:
			rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
		S = np.vstack((S, rxn_penalty_vec))
		S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
		lb = np.hstack((lb, 0))
		ub = np.hstack((ub, np.inf))
		b = np.hstack((b, 0))
		S[-1, -1] = -1
		rxn_list.append("total_flux")
		
		model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		obj = np.zeros((1, np.size(ub)))
		obj[0, -1] = -1
		print(obj)
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.NumericFocus = 3
		model.optimize()
		print(S, lb, ub)
		print("done")
		min_flux = model.objVal * -1
		print(min_flux)
		# input()
		N = sc.linalg.null_space(S)
		
		total_points = np.array([])
		origin = np.array([])
		
		# only 1 origin point should be made
		origin_created = False
		
		for sampling_level in sampling_levels:
			# Defines origin and points as empty
			
			accepted_points = np.array([])
			ortho_norm_basis = np.array([])
			
			ub[-1] = sampling_level * min_flux
			# if you want purely slices use
			# lb[-1] = sampling_level * min_flux
			model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
			
			# Suppress printing to console
			model.Params.LogToConsole = 0
			
			# Initializes number of subsequent vectors which have been rejected
			search_number = 0
			
			#
			rank = 0
			
			print(f"SAMPLING LEVEL {sampling_level}")
			
			while (search_number < max_search_number and rank < required_rank):
				# print(search_number,rank)
				# Finds random vertex
				while model.Status != 2:
					model.reset()
					obj = (2 * np.random.random_sample(S.shape[1]) - 1)
					# for variable objective total flux
					# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
					model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
					model.optimize()
				
				# The below ensures the points fall inside the bounds
				# Does seem to potentially break constraints
				# entries of points are just being moved, the points should be scaled
				# if this becomes a problem, define scaling function
				# point_to_use = react_flux.X
				point_to_use = react_flux.X
				point_to_use = np.minimum(point_to_use, ub)
				point_to_use = np.maximum(point_to_use, lb)
				point_to_use = np.matmul(N, np.matmul(np.transpose(N), point_to_use))
				
				if np.max(np.abs(np.matmul(S, point_to_use))) > Sv_tol:
					print(f"Curr WARNING: fidErr {np.max(np.abs(np.matmul(S, point_to_use)))}")
					print(f"Curr Error ub: {max(point_to_use - ub)}, lb: {max(lb - point_to_use)}")
				
				# print(point_to_use)
				# print(np.max(point_to_use[:-1]*rxn_penalty_vec[0]))
				# print(rxn_list[np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				# print(rxn_penalty_vec[0,np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				# input()
				# print(point_to_use[-1])
				# Saves rank to later check for rank change
				old_rank = rank
				if (search_number % 50 == 0 and search_number != 0) or rank == 0:
					print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
				# The first point becomes an origin point
				# The rank is 0 with the origin point since no vectors have been defined
				if origin_created:
					accepted_points, ortho_norm_basis, rank = self.lin_dep_np_ctp(ortho_norm_basis, origin,
					                                                              accepted_points, point_to_use, rank,
					                                                              unique_vector_threshold_exponent=unique_vector_threshold_exponent)
				else:
					origin = point_to_use
					origin_created = True
				# Reset the model to discard the current solution
				model.reset()
				# print(rank)
				# Reset the search_number if a vector is accepted as basis
				# I removed "or np.size(accepted_points) == 0" here which might lead to issues
				if rank != old_rank:
					search_number = 0
				else:
					search_number += 1
			print(rank)
			if np.size(total_points) == 0:
				if np.size(accepted_points) != 0:
					total_points = np.vstack((origin, accepted_points))
				else:
					total_points = origin
			else:
				if np.size(accepted_points) != 0:
					total_points = np.vstack((total_points, accepted_points))
		print(np.average(total_points, axis=0))
		print(1)
		# if multiple 1 norms were used the resulting vectors will span the space, but not be a basis
		# this corrects for that
		
		# This section removes any vectors which are approximately equivalent
		initial_size = np.shape(total_points)[0]
		# rounds copy of points and uses return_index in np.unique to grab unique original vectors
		# This lets two points be considered too similar without potential loss of machine precision of accepted point
		# This may not actally matter, but is not costly to implement
		new_total_points_rng = np.round(total_points, unique_vector_threshold_exponent)
		unique_total_points, unique_index = np.unique(new_total_points_rng, axis=0, return_index=True)
		total_points = total_points[unique_index]
		new_size = np.shape(total_points)[0]
		if initial_size != new_size:
			print(
				f"WARNING: There were {initial_size - new_size} repeated vectors in the basis. They have been removed")
			time.sleep(3)
		
		# returns set of points with top point being basis
		print(np.average(total_points, axis=0))
		print(3)
		input()
		input()
		self.model_dict["total_flux_limits"] = [min_flux, max(sampling_levels) * min_flux]
		if points_save_path != "":
			np.save(points_save_path / ("warmup_" + self.model_dict["model_name"]), total_points)
		if model_save_path != "":
			self.save_model_as_fast_key_json(model_save_path)
		return total_points
	
	def pinch_relative_change_in_bounds(self, max_sampling_level, tol=1e-6, rxn_penalty_vec=None):
		# Generates an array var_ind which specifies index values of variable entries, a nonzero b, a vector v_c which
		# is in the center of all FVA bounds and which can be used to convert from a reduced drawn vector to
		# a full one by v_c[var_ind] = v
		# a recalculation of bounds should not be needed since the true system bounds will shrink from pinching and
		# a point which surpasses the diminished bounds should run afoul of Sv=b constraint, plus lb < 0 is checked
		# and that is the most important constraint
		# fva_dict, v_c = self.find_middle_point(max_sampling_level)
		
		# v_c = v_c[:-1]
		
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		lbc = cp.deepcopy(lb)
		ubc = cp.deepcopy(ub)
		fva_dict = self.fva(lb=lb, ub=ub, S=S, b=b, rxn_list=rxn_list)
		for rxn_ind in range(len(rxn_list)):
			lb[rxn_ind] = fva_dict[rxn_list[rxn_ind]]['lb']
			ub[rxn_ind] = fva_dict[rxn_list[rxn_ind]]['ub']
		
		var_ind = np.arange(np.size(ub))
		var_ind = var_ind[ub - lb > tol]
		
		con_ind = np.arange(np.size(ub))
		con_ind = con_ind[ub - lb <= tol]
		print(np.shape(con_ind), np.shape(var_ind))
		
		min_flux, vec_min_flux = self.get_min_flux(lb, ub, b, S, return_vec=True)
		
		print(min_flux)
		S_con = S[:, con_ind]
		v_con = vec_min_flux[con_ind]
		
		pinch_b = -1 * np.matmul(S_con, v_con)
		min_rxn_names = list(np.array(rxn_list)[var_ind])
		min_flux, vec_min_flux_red = self.get_min_flux(lb[var_ind], ub[var_ind], pinch_b, S[:, var_ind],
		                                               return_vec=True)
		vec_min_flux[var_ind] = vec_min_flux_red[:-1]
		print(min_flux)
		fva_dict_small = self.fva(lb[var_ind], ub[var_ind], S[:, var_ind], pinch_b, min_rxn_names)
		lb_var = np.zeros(len(min_rxn_names))
		ub_var = np.zeros(len(min_rxn_names))
		for rxn_ind in range(len(min_rxn_names)):
			lb_var[rxn_ind] = fva_dict_small[min_rxn_names[rxn_ind]]['lb']
			ub_var[rxn_ind] = fva_dict_small[min_rxn_names[rxn_ind]]['ub']
		
		lb[var_ind] = lb_var
		lb[con_ind] = v_con
		lb[lb < 0] = 0
		
		ub[var_ind] = ub_var
		ub[con_ind] = v_con
		ub[ub < 0] = 0
		
		print(np.max(np.abs(ub - ubc)), np.max(np.abs((ub - ubc) / (ub + ubc + 1e-20))))
		print(np.max(np.abs(lb - lbc)), np.max(np.abs((lb - lbc) / (lb + lbc + 1e-20))))
		print(np.min(lb), np.min(lbc))
		print(np.min(vec_min_flux), np.sum(vec_min_flux))
		print(min_flux)
		input()
		
		# Do not remove this, this section will show the relative change in bounds from pinching the system
		# for rxn in min_rxn_names:
		#	lb_rat = (fva_dict_small[rxn]['lb'] - fva_dict[rxn]['lb'])/(fva_dict_small[rxn]['lb'] + fva_dict[rxn]['lb']+1e-20)
		#	ub_rat = (fva_dict_small[rxn]['ub'] - fva_dict[rxn]['ub'])/(fva_dict_small[rxn]['ub'] + fva_dict[rxn]['ub']+1e-20)
		#	if np.abs(lb_rat) > 1e-5 or np.abs(ub_rat) > 1e-5:
		#		print("rats", lb_rat, ub_rat)
		#		print("lb", fva_dict_small[rxn]['lb'] , fva_dict[rxn]['lb'])
		#		print("ub", fva_dict_small[rxn]['ub'], fva_dict[rxn]['ub'])
		# input()
		
		print(min_flux, np.sum(v_con) + min_flux)
		print(max_flux, np.sum(v_con) + max_flux)
		input()
		fva_dict = self.fva(lb=lb, ub=ub, S=S, b=b, rxn_list=rxn_list)
		for rxn_ind in range(len(rxn_list)):
			lb[rxn_ind] = fva_dict[rxn_list[rxn_ind]]['lb']
			ub[rxn_ind] = fva_dict[rxn_list[rxn_ind]]['ub']
		
		return var_ind, pinch_b, lb, vec_min_flux, ub
	
	def quick_pinch(self, max_sampling_level, tol=1e-6, rxn_penalty_vec=None):
		# Description:
		#   Fully constricts dimensions which are tightly bound, converting S_full*v_full = 0 into S_var * v_var = - S_con * v_con
		#   Where S_full is the stochiometic matrix, v_full is the vector of fluxes,  S_var and v_var contain entries which
		#   correspond to unconstrained fluxes while S_con and v_con correspond to constrained fluxes and as a result
		#   S_con*v_con is constant
		# Inputs:
		#   max_sampling_level - float: the largest multiplicative factor placed on the minimum reaction penalty
		#   tol - float: Any bounds closer together than this value will be pinched
		#   rxn_penalty_vec - numpy array (size shoudl correspond to rxn count): gives weights for 1 norm constraint.
		#       default will weight all fluxes similarly
		# Outputs:
		#   var_ind - numpy array (size will be the number of unpinched reactions): provides indexes for the variable reactions
		#       this assumes alphabetical reaction ordering (which is the defualt order given when dicts_to_mats is called)
		#       This value allows one to convert from full vectors to vectors which correspond to just variable fluxes
		#   pinch_b - numpy array (size will be the number of metabolites): - S_con * v_con which is used in the modified
		#       flux balance constraint, Sv = b
		#   lb - lower bound (size will be number of reactions + 1): gives lower bound values after pinching, IMPORTANTLY
		#       this includes the minimum flux value of the unpinched fluxes (not total flux at min point, but total of variable fluxes)
		#   vec_min_flux - pinch_point (size will be number of reactions + 1): The exact point the bounds are pinched to
		#       contains the full vector with variable fluxes, constant fluxes, and the total flux of the point
		#   ub - upper bound (size will be number of reactions + 1): gives upper bound values after pinching, IMPORTANTLY
		#       this includes the maximum flux value of the unpinched fluxes (not total flux at min point, but total of variable fluxes) as determined by the minimum and the max_sampling_level
		
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		print("keep")
		print(np.sum(lb),np.sum(ub),np.sum(S))

		# uncomment this to see the effect of pinching on the total space bounds
		# lbc = cp.deepcopy(lb)
		# ubc = cp.deepcopy(ub)
		
		fva_dict = self.fva(lb=lb, ub=ub, S=S, b=b, rxn_list=rxn_list)
		print(fva_dict)
		ind = 0
		# 639
		for key in fva_dict:
			print(ind)
			print(fva_dict[key]['lb'],fva_dict[key]['ub'])
			ind += 1
		#input()
		for rxn_ind in range(len(rxn_list)):
			lb[rxn_ind] = fva_dict[rxn_list[rxn_ind]]['lb']
			ub[rxn_ind] = fva_dict[rxn_list[rxn_ind]]['ub']
		
		var_ind = np.arange(np.size(ub))
		var_ind = var_ind[ub - lb > tol]
		
		con_ind = np.arange(np.size(ub))
		con_ind = con_ind[ub - lb <= tol]
		print(np.sum(lb),np.sum(ub),np.sum(b),np.sum(S))
		#input()
		min_flux, vec_min_flux = self.get_min_flux(lb, ub, b, S, return_vec=True)
		print(np.sum(min_flux), np.sum(vec_min_flux))
		#input()
		vec_min_flux = vec_min_flux[:-1]
		S_con = S[:, con_ind]
		v_con = vec_min_flux[con_ind]
		pinch_b = -1 * np.matmul(S_con, v_con)
		min_rxn_names = list(np.array(rxn_list)[var_ind])
		
		S_var = S[:, var_ind]
		lb_var = lb[var_ind]
		ub_var = ub[var_ind]
		# there does seem to be small difference in fva outcomes
		min_flux_pinch, vec_min_flux_red = self.get_min_flux(lb_var, ub_var, pinch_b, S_var,
		                                                     return_vec=True)
		print(np.sum(min_flux_pinch), np.sum(vec_min_flux_red))
		#input()
		min_rxn_names.append("Total_Flux")
		
		if rxn_penalty_vec is None:
			rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
		else:
			rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
		
		rxn_penalty_vec_var = rxn_penalty_vec[0, var_ind]
		S_var_tf = np.vstack((S_var, rxn_penalty_vec_var))
		S_var_tf = np.hstack((S_var_tf, np.zeros((np.shape(S_var_tf)[0], 1))))
		lb_var_tf = np.hstack((lb_var, min_flux_pinch))
		ub_var_tf = np.hstack((ub_var, min_flux_pinch * max_sampling_level))
		pinch_b = np.hstack((pinch_b, 0))
		S_var_tf[-1, -1] = -1
		vec_min_flux[var_ind] = vec_min_flux_red[:-1]
		vec_min_flux = np.hstack((vec_min_flux, vec_min_flux_red[-1]))
		# uncomment this to see the effect of pinching on the total space bounds
		# ub_var_tf[-1] = np.inf
		fva_dict_small = self.fva(lb_var_tf, ub_var_tf, S_var_tf, pinch_b, min_rxn_names)
		lb_var = np.zeros(len(min_rxn_names))
		ub_var = np.zeros(len(min_rxn_names))
		for rxn_ind in range(len(min_rxn_names)):
			lb_var[rxn_ind] = fva_dict_small[min_rxn_names[rxn_ind]]['lb']
			ub_var[rxn_ind] = fva_dict_small[min_rxn_names[rxn_ind]]['ub']
		
		lb[var_ind] = lb_var[:-1]
		lb[con_ind] = v_con
		lb[lb < 0] = 0
		
		ub[var_ind] = ub_var[:-1]
		ub[con_ind] = v_con
		ub[ub < 0] = 0
		# uncomment this to see the effect of pinching on the total space bounds
		# It should be nearly 0
		# print(np.max(np.abs(ub - ubc)), np.max(np.abs((ub - ubc) / (ub + ubc + 1e-20))))
		# print(np.max(np.abs(lb - lbc)), np.max(np.abs((lb - lbc) / (lb + lbc + 1e-20))))
		# input()
		lb = np.hstack((lb, lb_var[-1]))
		ub = np.hstack((ub, ub_var[-1]))
		
		return var_ind, pinch_b[:-1], lb, vec_min_flux, ub
	
	def get_min_flux(self, lb, ub, b, S, env=None, rxn_penalty_vec=None, return_vec=False, invert=False):
		# should be used on lb,ub,b,S without total flux
		if rxn_penalty_vec is None:
			rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
		else:
			rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
		
		# rxn_penalty_vec = rxn_penalty_vec[var_ind]
		S = np.vstack((S, rxn_penalty_vec))
		S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
		lb = np.hstack((lb, 0))
		ub = np.hstack((ub, np.inf))
		b = np.hstack((b, 0))
		S[-1, -1] = -1
		# rxn_list.append("total_flux")
		if self.model_dict["gurobi_token"] is not None:
			model = gp.Model("min_flux_model_" + self.model_dict["model_name"], env=self.model_dict["gurobi_token"])
		else:
			model = gp.Model("min_flux_model_" + self.model_dict["model_name"])
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		# Suppress printing to console
		model.Params.LogToConsole = 0
		
		obj = np.zeros(S.shape[1])
		obj[-1] = -1
		# for variable objective total flux
		# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
		if invert:
			model.setObjective(obj @ react_flux, GRB.MINIMIZE)
		else:
			model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.NumericFocus = 3
		model.Params.FeasibilityTol = 1e-9
		
		model.optimize()
		if not return_vec:
			return model.objVal * -1
		if return_vec:
			return model.objVal * -1, react_flux.X
	
	def correct_point(self, lb, ub, b, S, env, point_to_use):
		test_model = gp.Model("test_model_" + self.model_dict["model_name"], env=env)
		test_react_flux = test_model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="test_react_flux", lb=lb,
		                                     ub=ub)
		# print(react_flux)
		test_model.addConstr(S @ test_react_flux == b, name="test_c")
		# (mf - cpf)^2 = mf^2 - 2*mf*cpf + cpf^2
		Q = np.eye(np.size(point_to_use))
		test_model.setObjective(
			test_react_flux @ Q @ test_react_flux - (2 * point_to_use) @ test_react_flux + np.sum(point_to_use ** 2),
			GRB.MINIMIZE)
		test_model.Params.LogToConsole = 0
		test_model.Params.NumericFocus = 3
		test_model.Params.FeasibilityTol = 1e-9
		
		test_model.optimize()
		return test_react_flux.X
	
	def generate_warmup_gb_pinch(self, max_search_number, sampling_levels, pinch_tol=1e-6,pinch_1_norm = None, repeat_basis=1,
	                             points_save_path="", model_save_path="",
	                             required_rank=np.inf, rxn_penalty_vec=None, unique_vector_threshold_exponent=9, trans_NS = None,
	                             **kwargs):
		# rel_distance_to_min_flux
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		
		# unique_vector_threshold_exponent gives the number of decimal places out the algorithm searches when
		# removing identical vectors
		lb_full, ub_full, S_full, b, rxn_list, met_list = self.dicts_to_mats()
		
		if pinch_1_norm is None:
			pinch_1_norm = max(sampling_levels)
		
		var_ind, pinch_b, lb_full, v_c, ub_full = self.quick_pinch(pinch_1_norm, tol=pinch_tol)
		print(lb_full[-1], ub_full[-1], v_c[-1])
		
		lb_full = lb_full[:-1]
		v_c = v_c[:-1]
		ub_full = ub_full[:-1]
		
		initial_point = cp.deepcopy(v_c)[var_ind]
		
		initial_point = np.hstack((initial_point, np.sum(initial_point)))
		
		v_c[var_ind] = np.ones_like(var_ind) * np.nan
		
		lb = lb_full[var_ind]
		ub = ub_full[var_ind]
		S = S_full[:, var_ind]
		b = pinch_b
		
		# If rxn_penalty_vec is None a standard 1 norm is used
		
		# Currently this runs into issues if total flux can be zero as this defies how the 1 norm is handled
		
		# A pass to remove vectors which are identical down to a machine precision might be useful
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		
		Sv_tol = 1e-9
		lb_ub_tol = 1e-10
		
		if rxn_penalty_vec is None:
			rxn_penalty_vec = np.ones((1, np.shape(S_full)[1]))
		else:
			rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
		
		S_full = np.vstack((S_full, rxn_penalty_vec))
		S_full = np.hstack((S_full, np.zeros((np.shape(S_full)[0], 1))))
		S_full[-1, -1] = -1
		
		rxn_penalty_vec = rxn_penalty_vec[0, var_ind]
		
		S = np.vstack((S, rxn_penalty_vec))
		S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
		lb = np.hstack((lb, 0))
		ub = np.hstack((ub, np.inf))
		b = np.hstack((b, 0))
		S[-1, -1] = -1
		
		model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		obj = np.zeros((1, np.size(ub)))
		obj[0, -1] = -1
		
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.NumericFocus = 3
		model.Params.FeasibilityTol = 1e-9
		model.optimize()
		print(model.objVal)

		
		min_flux = model.objVal * -1
		
		total_points = np.array([])
		origin = np.array([])
		
		# only 1 origin point should be made
		origin_created = False
		
		# add iteration for repeats here
		for repeat in range(repeat_basis):
			
			for sampling_level in sampling_levels:
				# Defines origin and points as empty
				
				accepted_points = np.array([])
				ortho_norm_basis = np.array([])
				
				ub[-1] = sampling_level * min_flux
				# if you want purely slices use
				# lb[-1] = sampling_level * min_flux
				model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
				react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
				# print(react_flux)
				model.addConstr(S @ react_flux == b, name="c")
				
				# Suppress printing to console
				model.Params.LogToConsole = 0
				
				# Initializes number of subsequent vectors which have been rejected
				search_number = 0
				
				#
				rank = 0
				
				print(f"SAMPLING LEVEL {sampling_level}")
				
				while (search_number < max_search_number and rank < required_rank):
					# print(search_number,rank)
					# Finds random vertex
					while model.Status != 2:
						
						model.reset()
						obj = (2 * np.random.random_sample(S.shape[1]) - 1)
						# for variable objective total flux
						# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
						model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
						model.Params.LogToConsole = 0
						model.Params.NumericFocus = 3
						model.Params.FeasibilityTol = 1e-8
						
						model.optimize()
						#print(model.Status)
					point_to_use = react_flux.X
					# point_to_use[0] += 0.3
					
					# The below ensures the points fall inside the bounds
					# Does seem to potentially break constraints
					# entries of points are just being moved, the points should be scaled
					# if this becomes a problem, define scaling function
					# point_to_use = react_flux.X
					
					# point_to_use = np.minimum(point_to_use, ub)
					# point_to_use = np.maximum(point_to_use, lb)
					# point_to_use = np.matmul(N, np.matmul(np.transpose(N), point_to_use))
					
					# print(f"Curr WARNING: fidErr {np.max(np.abs(np.matmul(S, point_to_use)-b))}")
					# print(f"Curr Error ub: {max(point_to_use - ub)}, lb: {max(lb - point_to_use)}")
					
					if (np.max(np.abs(np.matmul(S, point_to_use) - b)) > Sv_tol) or max(
							point_to_use - ub) > lb_ub_tol or max(lb - point_to_use) > lb_ub_tol or min(
							point_to_use) < 0:
						point_to_use = self.correct_point(lb, ub, b, S, env, point_to_use)
						if (np.max(np.abs(np.matmul(S, point_to_use) - b)) > Sv_tol) or max(
								point_to_use - ub) > lb_ub_tol or max(
							lb - point_to_use) > lb_ub_tol or np.min(point_to_use) < 0:
							print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, point_to_use) - b))}")
							print(
								f"WARNING: Error ub: {max(point_to_use - ub)}, lb: {max(lb - point_to_use)}, lb of point = {np.min(point_to_use)}")
					
					# Saves rank to later check for rank change
					old_rank = rank
					#print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
					if (search_number % 50 == 0 and search_number != 0) or rank == 0:
						print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
					# The first point becomes an origin point
					# The rank is 0 with the origin point since no vectors have been defined
					if origin_created:
						accepted_points, ortho_norm_basis, rank = self.lin_dep_np_ctp(ortho_norm_basis, origin,
						                                                              accepted_points, point_to_use,
						                                                              rank,
						                                                              unique_vector_threshold_exponent=unique_vector_threshold_exponent)
					else:
						# origin[-1] = np.sum(origin[:-1])
						point_to_use = initial_point
						if (np.max(np.abs(np.matmul(S, point_to_use) - b)) > Sv_tol) or max(
								point_to_use - ub) > lb_ub_tol or max(lb - point_to_use) > lb_ub_tol or min(
							point_to_use) < 0:
							point_to_use = self.correct_point(lb, ub, b, S, env, point_to_use)
							if (np.max(np.abs(np.matmul(S, point_to_use) - b)) > Sv_tol) or max(
									point_to_use - ub) > lb_ub_tol or max(
								lb - point_to_use) > lb_ub_tol or np.min(point_to_use) < 0:
								print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, point_to_use) - b))}")
								print(
									f"WARNING: Error ub: {max(point_to_use - ub)}, lb: {max(lb - point_to_use)}, lb of point = {np.min(point_to_use)}")
						# input()
						origin = point_to_use
						origin_created = True
					# Reset the model to discard the current solution
					model.reset()
					# print(rank)
					# Reset the search_number if a vector is accepted as basis
					# I removed "or np.size(accepted_points) == 0" here which might lead to issues
					if rank != old_rank:
						search_number = 0
					else:
						search_number += 1
				print(rank)
				if np.size(total_points) == 0:
					if np.size(accepted_points) != 0:
						total_points = np.vstack((origin, accepted_points))
					else:
						total_points = origin
				else:
					if np.size(accepted_points) != 0:
						total_points = np.vstack((total_points, accepted_points))
		print(np.average(total_points, axis=0))
		print(np.shape(total_points))
		print(1)
		# if multiple 1 norms were used the resulting vectors will span the space, but not be a basis
		# this corrects for that
		
		# This section removes any vectors which are approximately equivalent
		initial_size = np.shape(total_points)[0]
		# rounds copy of points and uses return_index in np.unique to grab unique original vectors
		# This lets two points be considered too similar without potential loss of machine precision of accepted point
		# This may not actally matter, but is not costly to implement
		new_total_points_rng = np.round(total_points, unique_vector_threshold_exponent)
		unique_total_points, unique_index = np.unique(new_total_points_rng, axis=0, return_index=True)
		total_points = total_points[unique_index]
		new_size = np.shape(total_points)[0]
		if initial_size != new_size:
			print(
				f"WARNING: There were {initial_size - new_size} repeated vectors in the basis. They have been removed")
			time.sleep(3)
		
		# returns set of points with top point being basis
		print(np.average(total_points, axis=0))
		constant_flux_total = np.nansum(v_c)
		ub_full = np.hstack((ub_full, max(sampling_levels) * min_flux + constant_flux_total))
		lb_full = np.hstack((lb_full, min_flux + constant_flux_total))
		full_total_points = np.zeros((np.shape(total_points)[0], np.size(lb_full)))
		for i in range(np.shape(full_total_points)[0]):
			v_c[var_ind] = total_points[i, :-1]
			full_total_points[i, :-1] = v_c
			full_total_points[i, -1] = np.sum(full_total_points[i, :-1])
			if (np.max(np.abs(np.matmul(S_full, full_total_points[i]))) > Sv_tol) or max(
					full_total_points[i] - ub_full) > lb_ub_tol or max(lb_full - full_total_points[i]) > lb_ub_tol:
				print(f"Curr WARNING: fidErr {np.max(np.abs(np.matmul(S_full, full_total_points[i])))}")
				print(
					f"Curr Error ub: {max(full_total_points[i] - ub_full)}, lb: {max(lb_full - full_total_points[i])}")
				input()
		print(3)
		for i in full_total_points:
			#print(i)
			lb_full, ub_full, S_full, b, rxn_list, met_list = self.dicts_to_mats()
			simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = self.precomp_positive_flux_point_to_bi_dir(
				rxn_list, prune_specific=["biomass_reaction"])

			ent_flux_point = self.positive_flux_point_to_bi_dir(i[:-1], simp_neg_ind,
			                                                    comp_neg_ind, comp_perm,
			                                                    cutter,
			                                                    exchange_cutter=exchange_cutter)
			if self.generate_therm_model_new(trans_NS, ent_flux_point) == 2:
				print("good")
		
		self.model_dict["total_flux_limits"] = [lb_full[-1], ub_full[-1]]
		if points_save_path != "":
			np.save(points_save_path / ("warmup_" + self.model_dict["model_name"]), full_total_points)
		if model_save_path != "":
			self.save_model_as_fast_key_json(model_save_path)
		return total_points
	
	def generate_warmup_gb_qr(self, repeats, samples_per_rep, sampling_levels, points_save_path="", model_save_path="",
	                          required_rank=np.inf, rxn_penalty_vec=None, unique_vector_threshold_exponent=9, **kwargs):
		# rel_distance_to_min_flux
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		
		# unique_vector_threshold_exponent gives the number of decimal places out the algorithm searches when
		# removing identical vectors
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		
		# If rxn_penalty_vec is None a standard 1 norm is used
		
		# Currently this runs into issues if total flux can be zero as this defies how the 1 norm is handled
		
		# A pass to remove vectors which are identical down to a machine precision might be useful
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		
		Sv_tol = 1e-9
		
		if rxn_penalty_vec is None:
			rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
		else:
			rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
		S = np.vstack((S, rxn_penalty_vec))
		S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
		lb = np.hstack((lb, 0))
		ub = np.hstack((ub, np.inf))
		b = np.hstack((b, 0))
		S[-1, -1] = -1
		rxn_list.append("total_flux")
		
		model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		obj = np.zeros((1, np.size(ub)))
		obj[0, -1] = -1
		print(obj)
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.NumericFocus = 3
		model.optimize()
		print(S, lb, ub)
		print("done")
		min_flux = model.objVal * -1
		print(min_flux)
		# input()
		N = sc.linalg.null_space(S)
		
		total_points = np.array([])
		origin = np.array([])
		
		# only 1 origin point should be made
		origin_created = False
		
		for sampling_level in sampling_levels:
			# Defines origin and points as empty
			
			accepted_points = np.array([])
			ortho_norm_basis = np.array([])
			
			ub[-1] = sampling_level * min_flux
			# if you want purely slices use
			# lb[-1] = sampling_level * min_flux
			model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
			
			# Suppress printing to console
			model.Params.LogToConsole = 0
			
			# Initializes number of subsequent vectors which have been rejected
			search_number = 0
			
			#
			rank = 0
			search_matrix = np.zeros((samples_per_rep, np.size(ub)))
			accept_matrix = np.array([])
			print(f"SAMPLING LEVEL {sampling_level}")
			
			while (search_number < repeats and rank < required_rank):
				# print(search_number,rank)
				# Finds random vertex
				for search_vec_ind in range(np.shape(search_matrix)[0]):
					
					while model.Status != 2:
						model.reset()
						obj = (2 * np.random.random_sample(S.shape[1]) - 1)
						# for variable objective total flux
						# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
						model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
						model.optimize()
					
					# The below ensures the points fall inside the bounds
					# Does seem to potentially break constraints
					# entries of points are just being moved, the points should be scaled
					# if this becomes a problem, define scaling function
					# point_to_use = react_flux.X
					point_to_use = react_flux.X
					point_to_use = np.minimum(point_to_use, ub)
					point_to_use = np.maximum(point_to_use, lb)
					point_to_use = np.matmul(N, np.matmul(np.transpose(N), point_to_use))
					
					search_matrix[search_vec_ind] = point_to_use
				Q, R = np.linalg.qr(search_matrix)
				print(Q)
				print(np.shape(Q))
				print(R)
				print(np.shape(R))
				input()
				
				if np.max(np.abs(np.matmul(S, point_to_use))) > Sv_tol:
					print(f"Curr WARNING: fidErr {np.max(np.abs(np.matmul(S, point_to_use)))}")
					print(f"Curr Error ub: {max(point_to_use - ub)}, lb: {max(lb - point_to_use)}")
				
				# print(point_to_use)
				# print(np.max(point_to_use[:-1]*rxn_penalty_vec[0]))
				# print(rxn_list[np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				# print(rxn_penalty_vec[0,np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				# input()
				# print(point_to_use[-1])
				# Saves rank to later check for rank change
				old_rank = rank
				if (search_number % 50 == 0 and search_number != 0) or rank == 0:
					print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
				# The first point becomes an origin point
				# The rank is 0 with the origin point since no vectors have been defined
				if origin_created:
					accepted_points, ortho_norm_basis, rank = self.lin_dep_np_ctp(ortho_norm_basis, origin,
					                                                              accepted_points, point_to_use, rank,
					                                                              unique_vector_threshold_exponent=unique_vector_threshold_exponent)
				else:
					origin = point_to_use
					origin_created = True
				# Reset the model to discard the current solution
				model.reset()
				# print(rank)
				# Reset the search_number if a vector is accepted as basis
				# I removed "or np.size(accepted_points) == 0" here which might lead to issues
				if rank != old_rank:
					search_number = 0
				else:
					search_number += 1
			print(rank)
			if np.size(total_points) == 0:
				if np.size(accepted_points) != 0:
					total_points = np.vstack((origin, accepted_points))
				else:
					total_points = origin
			else:
				if np.size(accepted_points) != 0:
					total_points = np.vstack((total_points, accepted_points))
		print(np.average(total_points, axis=0))
		print(1)
		# if multiple 1 norms were used the resulting vectors will span the space, but not be a basis
		# this corrects for that
		
		# This section removes any vectors which are approximately equivalent
		initial_size = np.shape(total_points)[0]
		# rounds copy of points and uses return_index in np.unique to grab unique original vectors
		# This lets two points be considered too similar without potential loss of machine precision of accepted point
		# This may not actally matter, but is not costly to implement
		new_total_points_rng = np.round(total_points, unique_vector_threshold_exponent)
		unique_total_points, unique_index = np.unique(new_total_points_rng, axis=0, return_index=True)
		total_points = total_points[unique_index]
		new_size = np.shape(total_points)[0]
		if initial_size != new_size:
			print(
				f"WARNING: There were {initial_size - new_size} repeated vectors in the basis. They have been removed")
			time.sleep(3)
		
		# returns set of points with top point being basis
		print(np.average(total_points, axis=0))
		print(3)
		
		self.model_dict["total_flux_limits"] = [min_flux, max(sampling_levels) * min_flux]
		if points_save_path != "":
			np.save(points_save_path / ("warmup_" + self.model_dict["model_name"]), total_points)
		if model_save_path != "":
			self.save_model_as_fast_key_json(model_save_path)
		return total_points
	
	def generate_warmup_gb_cp(self, max_search_number, sampling_levels, points_save_path="", model_save_path="",
	                          required_rank=np.inf, rxn_penalty_vec=None, unique_vector_threshold_exponent=9, **kwargs):
		# rel_distance_to_min_flux
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		
		# unique_vector_threshold_exponent gives the number of decimal places out the algorithm searches when
		# removing identical vectors
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		
		# If rxn_penalty_vec is None a standard 1 norm is used
		
		# Currently this runs into issues if total flux can be zero as this defies how the 1 norm is handled
		
		# A pass to remove vectors which are identical down to a machine precision might be useful
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		
		if rxn_penalty_vec is None:
			rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
		else:
			rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
		S = np.vstack((S, rxn_penalty_vec))
		S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
		lb = np.hstack((lb, 0))
		ub = np.hstack((ub, np.inf))
		b = np.hstack((b, 0))
		S[-1, -1] = -1
		rxn_list.append("total_flux")
		
		model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		obj = np.zeros((1, np.size(ub)))
		obj[0, -1] = -1
		print(obj)
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.NumericFocus = 3
		model.optimize()
		print(S, lb, ub)
		print("done")
		min_flux = model.objVal * -1
		print(min_flux)
		# input()
		N = sc.linalg.null_space(S)
		
		total_points = np.array([])
		origin = np.array([])
		
		# only 1 origin point should be made
		origin_created = False
		
		for sampling_level in sampling_levels:
			# Defines origin and points as empty
			
			accepted_points = np.array([])
			ortho_norm_basis = np.array([])
			
			ub[-1] = sampling_level * min_flux
			# if you want purely slices use
			# lb[-1] = sampling_level * min_flux
			model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
			
			# Suppress printing to console
			model.Params.LogToConsole = 0
			
			# Initializes number of subsequent vectors which have been rejected
			search_number = 0
			
			#
			rank = 0
			
			print(f"SAMPLING LEVEL {sampling_level}")
			
			while (search_number < max_search_number and rank < required_rank):
				# print(search_number,rank)
				# Finds random vertex
				while model.Status != 2:
					model.reset()
					obj = (2 * np.random.random_sample(S.shape[1]) - 1)
					# for variable objective total flux
					# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
					model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
					model.optimize()
				
				# The below ensures the points fall inside the bounds
				# Does seem to potentially break constraints
				# entries of points are just being moved, the points should be scaled
				# if this becomes a problem, define scaling function
				# point_to_use = react_flux.X
				point_to_use = np.matmul(N, np.matmul(np.transpose(N), react_flux.X))
				point_to_use = np.minimum(point_to_use, ub)
				point_to_use = np.maximum(point_to_use, lb)
				
				# print(point_to_use)
				# print(np.max(point_to_use[:-1]*rxn_penalty_vec[0]))
				# print(rxn_list[np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				# print(rxn_penalty_vec[0,np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				# input()
				# print(point_to_use[-1])
				# Saves rank to later check for rank change
				old_rank = rank
				if (search_number % 50 == 0 and search_number != 0) or rank == 0:
					print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
				# The first point becomes an origin point
				# The rank is 0 with the origin point since no vectors have been defined
				if origin_created:
					accepted_points, ortho_norm_basis, rank = self.lin_dep_np_ctp(ortho_norm_basis, origin,
					                                                              accepted_points, point_to_use, rank,
					                                                              unique_vector_threshold_exponent=unique_vector_threshold_exponent)
				else:
					origin = point_to_use
					origin_created = True
				# Reset the model to discard the current solution
				model.reset()
				# print(rank)
				# Reset the search_number if a vector is accepted as basis
				# I removed "or np.size(accepted_points) == 0" here which might lead to issues
				if rank != old_rank:
					search_number = 0
				else:
					search_number += 1
			print(rank)
			if np.size(total_points) == 0:
				if np.size(accepted_points) != 0:
					total_points = np.vstack((origin, accepted_points))
				else:
					total_points = origin
			else:
				if np.size(accepted_points) != 0:
					total_points = np.vstack((total_points, accepted_points))
		print(np.average(total_points, axis=0))
		print(1)
		# if multiple 1 norms were used the resulting vectors will span the space, but not be a basis
		# this corrects for that
		new_ortho_norm_basis = np.array([])
		new_accepted_points = np.array([])
		rank = 0
		for i in total_points[1:]:
			new_accepted_points, new_ortho_norm_basis, rank = self.lin_dep_np_ctp(new_ortho_norm_basis, total_points[0],
			                                                                      new_accepted_points, i, rank,
			                                                                      unique_vector_threshold_exponent=unique_vector_threshold_exponent)
		print(np.average(new_accepted_points, axis=0))
		print(2)
		
		# This section removes any vectors which are approximately equivalent
		initial_size = np.shape(new_accepted_points)[0]
		# rounds copy of points and uses return_index in np.unique to grab unique original vectors
		# This lets two points be considered too similar without potential loss of machine precision of accepted point
		# This may not actally matter, but is not costly to implement
		new_accepted_points_rnd = np.round(new_accepted_points, unique_vector_threshold_exponent)
		unique_accepted_points_rnd, unique_index = np.unique(new_accepted_points_rnd, axis=0, return_index=True)
		new_accepted_points = new_accepted_points[unique_index]
		new_size = np.shape(new_accepted_points)[0]
		if initial_size != new_size:
			print(
				f"WARNING: There were {initial_size - new_size} repeated vectors in the basis. They have been removed")
			time.sleep(3)
		
		# returns set of points with top point being basis
		total_points = np.vstack((origin, new_accepted_points))
		print(np.average(total_points, axis=0))
		print(3)
		self.model_dict["total_flux_limits"] = [min_flux, max(sampling_levels) * min_flux]
		if points_save_path != "":
			np.save(points_save_path / ("warmup_" + self.model_dict["model_name"]), total_points)
		if model_save_path != "":
			self.save_model_as_fast_key_json(model_save_path)
		return total_points
	
	def test_warmup_gb(self, sample_size, save_path, required_rank=np.inf, rel_distance_to_min_flux=-1,
	                   **kwargs):
		# rel_distance_to_min_flux
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		model = gp.Model("warmup_model_" + self.model_name, env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		if rel_distance_to_min_flux > 0:
			obj = -1 * np.ones((1, np.size(ub)))
			model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
			model.optimize()
			min_flux = model.objVal * -1
			max_allowed_flux = rel_distance_to_min_flux * min_flux
			print(max_allowed_flux)
			
			S = np.vstack((S, np.ones((1, np.shape(S)[1]))))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, 0))
			ub = np.hstack((ub, max_allowed_flux))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
			model = gp.Model("warmup_model_" + self.model_name, env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
		
		# Suppress printing to console
		model.Params.LogToConsole = 0
		rxn_count = np.size(ub)
		rxn_ind = np.arange(rxn_count)
		# Initializes rank
		rank = -1
		average_of_average_flux = np.zeros(rxn_count)
		for i in range(rxn_count):
			print(i)
			# Finds random vertex
			
			af = 0
			for j in range(sample_size):
				
				while model.Status != 2:
					np.random.shuffle(rxn_ind)
					obj = (np.random.random_sample(S.shape[1]) - 1)
					obj[rxn_ind[:i]] *= -1
					model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
					model.optimize()
				point_to_use = np.minimum(react_flux.X, ub)
				point_to_use = np.maximum(point_to_use, lb)
				af += point_to_use[-1]
				
				# Reset the model to discard the current solution
				model.reset()
			average_of_average_flux[i] = af / sample_size
			print(af / sample_size)
		
		print(rank)
		# returns set of points with top point being basis
		np.save(save_path / (self.model_name), average_of_average_flux)
	
	def generate_vertex_samples(self, n, save_path, rel_distance_to_min_flux=-1, **kwargs):
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		model = gp.Model("warmup_model_" + self.model_name, env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb,
		                           ub=ub)
		model.addConstr(S @ react_flux == b, name="c")
		if max_combined_flux > 0:
			obj = -1 * np.ones((1, np.size(ub)))
			model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
			model.optimize()
			min_flux = model.objVal * -1
			max_allowed_flux = rel_distance_to_min_flux * min_flux
			S = np.vstack((S, np.ones((1, np.shape(S)[1]))))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, 0))
			ub = np.hstack((ub, max_allowed_flux))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
			model = gp.Model("warmup_model_" + self.model_name, env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
		
		# Suppress printing to console
		model.Params.LogToConsole = 0
		
		# Initializes number of subsequent vectors which have been rejected
		search_number = 0
		
		# Initializes rank
		point_matrix = np.zeros((n, len(lb)))
		for i in range(n):
			if i % 100 == 0:
				print(i)
			# Finds random vertex
			while model.Status == 1:
				obj = (2 * np.random.random_sample(S.shape[1]) - 1)
				model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
				model.optimize()
			point_matrix[i] = react_flux.X
			model.reset()
		np.save(save_path / ("warmup_" + self.model_name), point_matrix)
		np.save(save_path / ("warmup_header_" + self.model_name), rxn_list)
		return point_matrix
	
	def HRSampler(self, origin_and_warmup_points, number_of_points, stepsPerPoint, rxn_penalty_vec=None, **kwargs):
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		Sv_cont_tol = 1e-9
		prev = 0
		
		# find center
		centerPoint = np.mean(origin_and_warmup_points, axis=0)
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = centerPoint
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			prev_point = origin
		
		if "total_flux_limits" in self.model_dict.keys():
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.model_dict["total_flux_limits"])))
			ub = np.hstack((ub, max(self.model_dict["total_flux_limits"])))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		N = sc.linalg.null_space(S)
		
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			randVector = np.random.random_sample(stepsPerPoint)
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				random_row_index = np.random.choice(np.shape(warmup_points)[0])
				rand_point = warmup_points[random_row_index]
				
				# Get a direction from the basis point to the warmup point
				u = self.normalize(rand_point - origin)
				
				# Figure out the distances to upper and lower bounds
				distUb = (ub - prev_point)
				distLb = (prev_point - lb)
				
				# Figure out if we are too close to a boundary
				# validDir = ((distUb > dTol) & (distLb > dTol))
				
				# Finds list of entries of unit vector which are either positive or negative
				# and which fall within tolerance values
				posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
				negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
				
				# Find the largest allowed step size in both the positive and negative directions
				# How many lengths of u can you move in any valid direction before you hit a wall
				pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
				neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
				
				# additional upperbound constraint to prevent exceeding imposed flux limit
				# Cases:
				# sum u is negative
				# sum u is positive
				# sum u is small
				# if max_combined_flux != -1:
				#	current_flux = np.sum(prev_point)
				#	max_add_from_u = max_combined_flux - current_flux
				
				# This might look weird, but we are just seeing how many unit vectors we can move before
				# we hit the flux limit, if that number is positive and if it is smaller than
				# the first positive distance to an upper bound it replaces it as it is more restrictive.
				# later on the code will check to see if it is more restrictive than the other ub distances
				
				# If we have to move a negative amount in the direction of the unit vector the code
				# also just checks to see if that is more restrictive than the first negative distance to an
				# upper bound
				#	rel_u_to_add = max_add_from_u/np.sum(u)
				#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
				#		pos_u_rel_dis_ub[0] = rel_u_to_add
				#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
				#		neg_u_rel_dis_ub[0] = rel_u_to_add
				
				# Find the smallest allowed step size in both the positive and negative directions
				neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u, axis=0)
				pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u, axis=0)
				
				# find the current value of the point
				
				# Find the smallest and largest allowed step sizes
				maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
				minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
				
				if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
					print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
					continue
				
				# grab random value from pre-generated list for the step distance
				stepDist = randVector[step_count - 1] * (maxStep - minStep) + minStep
				curPoint = prev_point + stepDist * u
				
				if total_steps % 100 == 0:
					curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				
				if total_steps % 200000 == 0 and False:
					print(f"Error {max(curPoint - ub)}, {max(lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				if np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0):
					indexer = np.arange(0, np.size(curPoint))
					overInd = indexer[curPoint > ub]
					underInd = indexer[curPoint < lb]
					# print(overInd)
					# print(underInd)
					# scaling might be useful here
					for j in overInd:
						curPoint[j] = ub[j]
					for j in underInd:
						curPoint[j] = lb[j]
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(curPoint[-1])
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	def HRSampler_2(self, origin_and_warmup_points, number_of_points, stepsPerPoint, rxn_penalty_vec=None, **kwargs):
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		Sv_cont_tol = 1e-9
		prev = 0
		
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = origin_and_warmup_points[0]
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			# find center
			prev_point = np.mean(origin_and_warmup_points, axis=0)
		
		if "total_flux_limits" in self.model_dict.keys():
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.model_dict["total_flux_limits"])))
			ub = np.hstack((ub, max(self.model_dict["total_flux_limits"])))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		N = sc.linalg.null_space(S)
		
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			randVector = np.random.random_sample(stepsPerPoint)
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				random_row_index = np.random.choice(np.shape(warmup_points)[0])
				rand_point = warmup_points[random_row_index]
				
				# Get a direction from the basis point to the warmup point
				u = self.normalize(rand_point - origin)
				
				# Figure out the distances to upper and lower bounds
				distUb = (ub - prev_point)
				distLb = (prev_point - lb)
				
				# Figure out if we are too close to a boundary
				# validDir = ((distUb > dTol) & (distLb > dTol))
				
				# Finds list of entries of unit vector which are either positive or negative
				# and which fall within tolerance values
				posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
				negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
				
				# Find the largest allowed step size in both the positive and negative directions
				# How many lengths of u can you move in any valid direction before you hit a wall
				pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
				neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
				
				# additional upperbound constraint to prevent exceeding imposed flux limit
				# Cases:
				# sum u is negative
				# sum u is positive
				# sum u is small
				# if max_combined_flux != -1:
				#	current_flux = np.sum(prev_point)
				#	max_add_from_u = max_combined_flux - current_flux
				
				# This might look weird, but we are just seeing how many unit vectors we can move before
				# we hit the flux limit, if that number is positive and if it is smaller than
				# the first positive distance to an upper bound it replaces it as it is more restrictive.
				# later on the code will check to see if it is more restrictive than the other ub distances
				
				# If we have to move a negative amount in the direction of the unit vector the code
				# also just checks to see if that is more restrictive than the first negative distance to an
				# upper bound
				#	rel_u_to_add = max_add_from_u/np.sum(u)
				#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
				#		pos_u_rel_dis_ub[0] = rel_u_to_add
				#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
				#		neg_u_rel_dis_ub[0] = rel_u_to_add
				
				# Find the smallest allowed step size in both the positive and negative directions
				neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u, axis=0)
				pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u, axis=0)
				
				# find the current value of the point
				
				# Find the smallest and largest allowed step sizes
				maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
				minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
				
				if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
					print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
					continue
				
				# grab random value from pre-generated list for the step distance
				stepDist = randVector[step_count - 1] * (maxStep - minStep) + minStep
				curPoint = prev_point + stepDist * u
				
				if total_steps % 100 == 0:
					curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				
				if total_steps % 200 == 0:
					if max(max(curPoint - ub), max(lb - curPoint)) > 1e-7:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				if np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0):
					indexer = np.arange(0, np.size(curPoint))
					overInd = indexer[curPoint > ub]
					underInd = indexer[curPoint < lb]
					# print(overInd)
					# print(underInd)
					# scaling might be useful here
					for j in overInd:
						curPoint[j] = ub[j]
					for j in underInd:
						curPoint[j] = lb[j]
					curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(curPoint[-1])
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	def HRSampler_gene_bias_single_dfs(self, origin_and_warmup_points, number_of_points, stepsPerPoint, RPS,
	                                   gene_penalty_mod=1, rxn_penalty_vec=None, **kwargs):
		# Note the gene_penalty_mod applies to the value in the exponent
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		Sv_cont_tol = 1e-9
		prev = 0
		
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = origin_and_warmup_points[0]
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			# find center
			prev_point = np.mean(origin_and_warmup_points, axis=0)
		
		if "total_flux_limits" in self.model_dict.keys():
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.model_dict["total_flux_limits"])))
			ub = np.hstack((ub, max(self.model_dict["total_flux_limits"])))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		N = sc.linalg.null_space(S)
		indexer = np.arange(0, np.size(lb))
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		proposed_points = 0
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				passing = True
				while passing:
					randVector = np.random.random_sample(2)
					random_row_index = np.random.choice(np.shape(warmup_points)[0])
					rand_point = warmup_points[random_row_index]
					
					# Get a direction from the basis point to the warmup point
					u = self.normalize(rand_point - origin)
					
					# Figure out the distances to upper and lower bounds
					distUb = (ub - prev_point)
					distLb = (prev_point - lb)
					
					# Figure out if we are too close to a boundary
					# validDir = ((distUb > dTol) & (distLb > dTol))
					
					# Finds list of entries of unit vector which are either positive or negative
					# and which fall within tolerance values
					posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
					negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
					
					# Find the largest allowed step size in both the positive and negative directions
					# How many lengths of u can you move in any valid direction before you hit a wall
					pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
					neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
					
					# additional upperbound constraint to prevent exceeding imposed flux limit
					# Cases:
					# sum u is negative
					# sum u is positive
					# sum u is small
					# if max_combined_flux != -1:
					#	current_flux = np.sum(prev_point)
					#	max_add_from_u = max_combined_flux - current_flux
					
					# This might look weird, but we are just seeing how many unit vectors we can move before
					# we hit the flux limit, if that number is positive and if it is smaller than
					# the first positive distance to an upper bound it replaces it as it is more restrictive.
					# later on the code will check to see if it is more restrictive than the other ub distances
					
					# If we have to move a negative amount in the direction of the unit vector the code
					# also just checks to see if that is more restrictive than the first negative distance to an
					# upper bound
					#	rel_u_to_add = max_add_from_u/np.sum(u)
					#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
					#		pos_u_rel_dis_ub[0] = rel_u_to_add
					#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
					#		neg_u_rel_dis_ub[0] = rel_u_to_add
					
					# Find the smallest allowed step size in both the positive and negative directions
					neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u,
					                                                                               axis=0)
					pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u,
					                                                                               axis=0)
					
					# find the current value of the point
					
					# Find the smallest and largest allowed step sizes
					maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
					minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
					
					if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
						print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
						continue
					
					# grab random value from pre-generated list for the step distance
					stepDist = randVector[0] * (maxStep - minStep) + minStep
					
					curPoint = prev_point + stepDist * u
					
					if np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0):
						
						overInd = indexer[curPoint > ub]
						underInd = indexer[curPoint < lb]
						# print(overInd)
						# print(underInd)
						# scaling might be useful here
						for j in overInd:
							curPoint[j] = ub[j]
						for j in underInd:
							curPoint[j] = lb[j]
						curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
					
					cur_support = self.gene_point_penalty(curPoint[:-1], RPS, gene_penalty_mod)
					prev_support = self.gene_point_penalty(prev_point[:-1], RPS, gene_penalty_mod)
					
					proposed_points += 1
					# print(cur_support,prev_support,np.log(randVector[1]))
					if proposed_points % 500 == 0:
						print(proposed_points, total_steps, (total_steps + 1) / proposed_points)
						print(prev_point[-1], curPoint[-1])
					
					if np.log(randVector[1]) < (cur_support - prev_support):
						passing = False
						continue
					else:
						continue
				
				if total_steps % 200 == 0:
					if max(max(curPoint - ub), max(lb - curPoint)) > 1e-7:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				# curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(curPoint[-1])
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	def HRSampler_gene_bias_single(self, origin_and_warmup_points, number_of_points, stepsPerPoint, RPS,
	                               gene_penalty_mod=1, rxn_penalty_vec=None, **kwargs):
		# Note the gene_penalty_mod applies to the value in the exponent
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		Sv_cont_tol = 1e-7
		step_tol = 1e-9
		prev = 0
		
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = origin_and_warmup_points[0]
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			# find center
			prev_point = np.mean(origin_and_warmup_points, axis=0)
		
		if "total_flux_limits" in self.model_dict.keys():
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.model_dict["total_flux_limits"])))
			ub = np.hstack((ub, max(self.model_dict["total_flux_limits"])))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		N = sc.linalg.null_space(S)
		indexer = np.arange(0, np.size(lb))
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		proposed_points = 0
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				passing = True
				while passing:
					randVector = np.random.random_sample(2)
					random_row_index = np.random.choice(np.shape(warmup_points)[0])
					rand_point = warmup_points[random_row_index]
					
					# Get a direction from the basis point to the warmup point
					u = self.normalize(rand_point - origin)
					
					# Figure out the distances to upper and lower bounds
					distUb = (ub - prev_point)
					distLb = (prev_point - lb)
					
					# Figure out if we are too close to a boundary
					# validDir = ((distUb > dTol) & (distLb > dTol))
					
					# Finds list of entries of unit vector which are either positive or negative
					# and which fall within tolerance values
					posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
					negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
					
					# Find the largest allowed step size in both the positive and negative directions
					# How many lengths of u can you move in any valid direction before you hit a wall
					pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
					neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
					
					# additional upperbound constraint to prevent exceeding imposed flux limit
					# Cases:
					# sum u is negative
					# sum u is positive
					# sum u is small
					# if max_combined_flux != -1:
					#	current_flux = np.sum(prev_point)
					#	max_add_from_u = max_combined_flux - current_flux
					
					# This might look weird, but we are just seeing how many unit vectors we can move before
					# we hit the flux limit, if that number is positive and if it is smaller than
					# the first positive distance to an upper bound it replaces it as it is more restrictive.
					# later on the code will check to see if it is more restrictive than the other ub distances
					
					# If we have to move a negative amount in the direction of the unit vector the code
					# also just checks to see if that is more restrictive than the first negative distance to an
					# upper bound
					#	rel_u_to_add = max_add_from_u/np.sum(u)
					#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
					#		pos_u_rel_dis_ub[0] = rel_u_to_add
					#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
					#		neg_u_rel_dis_ub[0] = rel_u_to_add
					
					# Find the smallest allowed step size in both the positive and negative directions
					neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u,
					                                                                               axis=0)
					pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u,
					                                                                               axis=0)
					
					# find the current value of the point
					
					# Find the smallest and largest allowed step sizes
					maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb))) - step_tol
					maxStep = max(maxStep, 0)
					minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub))) + step_tol
					minStep = min(minStep, 0)
					
					if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
						print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
						continue
					
					# grab random value from pre-generated list for the step distance
					stepDist = randVector[0] * (maxStep - minStep) + minStep
					
					curPoint = prev_point + stepDist * u
					stalls = 0
					while (np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0)) or np.max(
							np.abs(np.matmul(S, curPoint))) > Sv_cont_tol:
						curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
						overInd = indexer[curPoint > ub]
						underInd = indexer[curPoint < lb]
						# print(overInd)
						# print(underInd)
						# scaling might be useful here
						for j in overInd:
							curPoint[j] = ub[j]
						for j in underInd:
							curPoint[j] = lb[j]
						stalls += 1
						if stalls % 100 == 0:
							print(stalls)
							print(f"Curr WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
							print(f"Curr Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
					
					cur_support = self.gene_point_penalty(curPoint[:-1], RPS, gene_penalty_mod)
					prev_support = self.gene_point_penalty(prev_point[:-1], RPS, gene_penalty_mod)
					
					proposed_points += 1
					# print(cur_support,prev_support,np.log(randVector[1]))
					if proposed_points % 500 == 0:
						print(proposed_points, total_steps, (total_steps + 1) / proposed_points)
						print(prev_point[-1], curPoint[-1])
						print(prev_support, cur_support)
						print(f"Curr WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"Curr Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
						print(f"Prev WARNING: fidErr {np.max(np.abs(np.matmul(S, prev_point)))}")
						print(f"Prev Error ub: {max(prev_point - ub)}, lb: {max(lb - prev_point)}")
					
					if np.log(randVector[1]) < (cur_support - prev_support):
						passing = False
						continue
					else:
						continue
				
				if total_steps % 200 == 0:
					if max(max(curPoint - ub), max(lb - curPoint)) > 1e-7:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				# curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(curPoint[-1])
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	def HRSampler_gene_bias_single_ftl(self, origin_and_warmup_points, number_of_points, stepsPerPoint, RPS,
	                                   gene_penalty_mod=1, rxn_penalty_vec=None, **kwargs):
		# Note the gene_penalty_mod applies to the value in the exponent
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		Sv_cont_tol = 1e-9
		prev = 0
		
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = origin_and_warmup_points[0]
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			# find center
			prev_point = np.mean(origin_and_warmup_points, axis=0)
		
		if "total_flux_limits" in self.model_dict.keys():
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.model_dict["total_flux_limits"])))
			ub = np.hstack((ub, max(self.model_dict["total_flux_limits"])))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		N = sc.linalg.null_space(S)
		indexer = np.arange(0, np.size(lb))
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		proposed_points = 0
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				
				random_row_index = np.random.choice(np.shape(warmup_points)[0])
				rand_point = warmup_points[random_row_index]
				
				# Get a direction from the basis point to the warmup point
				u = self.normalize(rand_point - origin)
				
				# Figure out the distances to upper and lower bounds
				distUb = (ub - prev_point)
				distLb = (prev_point - lb)
				
				# Figure out if we are too close to a boundary
				# validDir = ((distUb > dTol) & (distLb > dTol))
				
				# Finds list of entries of unit vector which are either positive or negative
				# and which fall within tolerance values
				posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
				negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
				
				# Find the largest allowed step size in both the positive and negative directions
				# How many lengths of u can you move in any valid direction before you hit a wall
				pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
				neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
				
				# additional upperbound constraint to prevent exceeding imposed flux limit
				# Cases:
				# sum u is negative
				# sum u is positive
				# sum u is small
				# if max_combined_flux != -1:
				#	current_flux = np.sum(prev_point)
				#	max_add_from_u = max_combined_flux - current_flux
				
				# This might look weird, but we are just seeing how many unit vectors we can move before
				# we hit the flux limit, if that number is positive and if it is smaller than
				# the first positive distance to an upper bound it replaces it as it is more restrictive.
				# later on the code will check to see if it is more restrictive than the other ub distances
				
				# If we have to move a negative amount in the direction of the unit vector the code
				# also just checks to see if that is more restrictive than the first negative distance to an
				# upper bound
				#	rel_u_to_add = max_add_from_u/np.sum(u)
				#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
				#		pos_u_rel_dis_ub[0] = rel_u_to_add
				#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
				#		neg_u_rel_dis_ub[0] = rel_u_to_add
				
				# Find the smallest allowed step size in both the positive and negative directions
				neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u,
				                                                                               axis=0)
				pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u,
				                                                                               axis=0)
				
				# find the current value of the point
				
				# Find the smallest and largest allowed step sizes
				maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
				minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
				
				if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
					print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
					continue
				passing = True
				while passing:
					randVector = np.random.random_sample(2)
					# grab random value from pre-generated list for the step distance
					stepDist = randVector[0] * (maxStep - minStep) + minStep
					
					curPoint = prev_point + stepDist * u
					
					if np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0):
						
						overInd = indexer[curPoint > ub]
						underInd = indexer[curPoint < lb]
						# print(overInd)
						# print(underInd)
						# scaling might be useful here
						for j in overInd:
							curPoint[j] = ub[j]
						for j in underInd:
							curPoint[j] = lb[j]
						curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
					
					cur_support = self.gene_point_penalty(curPoint[:-1], RPS, gene_penalty_mod)
					prev_support = self.gene_point_penalty(prev_point[:-1], RPS, gene_penalty_mod)
					
					proposed_points += 1
					# print(cur_support,prev_support,np.log(randVector[1]))
					if proposed_points % 5000 == 0:
						print(proposed_points, total_steps, (total_steps + 1) / proposed_points)
						print(prev_point[-1], curPoint[-1])
					
					if np.log(randVector[1]) < (cur_support - prev_support):
						passing = False
						continue
					else:
						continue
				
				if total_steps % 200 == 0:
					if max(max(curPoint - ub), max(lb - curPoint)) > 1e-7:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				# curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(curPoint[-1])
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	def gene_point_penalty(self, flux_vector, RPS, constant):
		# returns error term in natural log scale
		error_term = np.dot(flux_vector, RPS)
		return -error_term * constant
	
	def RAS_to_RPS_mat(self, RAS, inf_val=1e24, nan_val=None):
		if nan_val is None:
			nan_val = np.nanmedian(RAS)
		RAS = np.nan_to_num(RAS, copy=True, nan=nan_val)
		RAS[RAS == 0] = 1 / inf_val
		return 1 / RAS
	
	def HRSampler_2_ss(self, origin_and_warmup_points, number_of_points, stepsPerPoint, step_size=1,
	                   rxn_penalty_vec=None, **kwargs):
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		Sv_cont_tol = 1e-7
		prev = 0
		
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = origin_and_warmup_points[0]
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			# find center
			prev_point = np.mean(origin_and_warmup_points, axis=0)
		
		if "total_flux_limits" in self.model_dict.keys():
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.model_dict["total_flux_limits"])))
			ub = np.hstack((ub, max(self.model_dict["total_flux_limits"])))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		N = sc.linalg.null_space(S)
		
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			randVector = np.random.random_sample(stepsPerPoint)
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				random_row_index = np.random.choice(np.shape(warmup_points)[0])
				rand_point = warmup_points[random_row_index]
				
				# Get a direction from the basis point to the warmup point
				u = self.normalize(rand_point - origin)
				
				# Figure out the distances to upper and lower bounds
				distUb = (ub - prev_point)
				distLb = (prev_point - lb)
				
				# Figure out if we are too close to a boundary
				# validDir = ((distUb > dTol) & (distLb > dTol))
				
				# Finds list of entries of unit vector which are either positive or negative
				# and which fall within tolerance values
				posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
				negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
				
				# Find the largest allowed step size in both the positive and negative directions
				# How many lengths of u can you move in any valid direction before you hit a wall
				pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
				neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
				
				# additional upperbound constraint to prevent exceeding imposed flux limit
				# Cases:
				# sum u is negative
				# sum u is positive
				# sum u is small
				# if max_combined_flux != -1:
				#	current_flux = np.sum(prev_point)
				#	max_add_from_u = max_combined_flux - current_flux
				
				# This might look weird, but we are just seeing how many unit vectors we can move before
				# we hit the flux limit, if that number is positive and if it is smaller than
				# the first positive distance to an upper bound it replaces it as it is more restrictive.
				# later on the code will check to see if it is more restrictive than the other ub distances
				
				# If we have to move a negative amount in the direction of the unit vector the code
				# also just checks to see if that is more restrictive than the first negative distance to an
				# upper bound
				#	rel_u_to_add = max_add_from_u/np.sum(u)
				#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
				#		pos_u_rel_dis_ub[0] = rel_u_to_add
				#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
				#		neg_u_rel_dis_ub[0] = rel_u_to_add
				
				# Find the smallest allowed step size in both the positive and negative directions
				neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u, axis=0)
				pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u, axis=0)
				
				# find the current value of the point
				
				# Find the smallest and largest allowed step sizes
				maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
				minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
				
				if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
					print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
					continue
				
				# grab random value from pre-generated list for the step distance
				
				stepDist = randVector[step_count - 1] * (maxStep - minStep) + minStep
				curPoint = prev_point + stepDist * u
				
				if np.max(np.abs(np.matmul(S, curPoint))) > Sv_cont_tol:
					curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				
				if total_steps % 200 == 0:
					if max(max(curPoint - ub), max(lb - curPoint)) > 1e-7:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				if np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0):
					indexer = np.arange(0, np.size(curPoint))
					overInd = indexer[curPoint > ub]
					underInd = indexer[curPoint < lb]
					# print(overInd)
					# print(underInd)
					# scaling might be useful here
					for j in overInd:
						curPoint[j] = ub[j]
					for j in underInd:
						curPoint[j] = lb[j]
					curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(curPoint[-1])
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	def HRSampler_gene_bias_lincomb(self, origin_and_warmup_points, number_of_points, stepsPerPoint, RPS,
	                                gene_penalty_mod=1, rxn_penalty_vec=None, **kwargs):
		# Note the gene_penalty_mod applies to the value in the exponent
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		
		# potentially add check for how much the realignment moves a point
		
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 0
		step_tol = 1e-9
		prev = 0
		
		Sv_tol = 1e-7
		lb_ub_tol = 1e-9
		
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = origin_and_warmup_points[0]
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "total_flux_limits" in self.model_dict.keys():
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.model_dict["total_flux_limits"])))
			ub = np.hstack((ub, max(self.model_dict["total_flux_limits"])))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		N = sc.linalg.null_space(S)
		indexer = np.arange(0, np.size(lb))
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			# find center
			prev_point = np.mean(origin_and_warmup_points, axis=0)
			# prev_point = origin_and_warmup_points[np.random.choice(np.shape(warmup_points)[0])]
			print(f"Curr WARNING: fidErr {np.max(np.abs(np.matmul(S, prev_point)))}")
			print(f"Curr Error ub: {max(prev_point - ub)}, lb: {max(lb - prev_point)}")
			
			if (np.max(np.abs(np.matmul(S, prev_point) - b)) > Sv_tol) or max(prev_point - ub) > lb_ub_tol or max(
					lb - prev_point) > lb_ub_tol:
				model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
				react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
				# print(react_flux)
				model.addConstr(S @ react_flux == b, name="c")
				# (mf - cpf)^2 = mf^2 - 2*mf*cpf + cpf^2
				Q = np.eye(np.size(prev_point))
				model.setObjective(
					react_flux @ Q @ react_flux - (2 * prev_point) @ react_flux + np.sum(prev_point ** 2),
					GRB.MINIMIZE)
				model.Params.NumericFocus = 3
				model.Params.FeasibilityTol = 1e-9
				model.optimize()
				prev_point = react_flux.X
				print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, prev_point) - b))}")
				print(f"WARNING: Error ub: {max(prev_point - ub)}, lb: {max(lb - prev_point)}")
			input()
		
		warmup_vecs = np.zeros_like(warmup_points)
		for i in range(np.shape(warmup_points)[0]):
			warmup_vecs[i] = warmup_points[i] - origin
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		proposed_points = 0
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				passing = True
				while passing:
					
					randVector = np.random.random_sample(2)
					
					u = (warmup_vecs.T * np.random.random((np.shape(warmup_points)[0], 1))[:, 0]).T
					u = self.normalize(np.sum(u, axis=0))
					
					# Figure out the distances to upper and lower bounds
					distUb = (ub - prev_point)
					distLb = (prev_point - lb)
					
					# Figure out if we are too close to a boundary
					# validDir = ((distUb > dTol) & (distLb > dTol))
					
					# Finds list of entries of unit vector which are either positive or negative
					# and which fall within tolerance values
					# posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
					# negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
					
					posValidDir = (u > uTol)
					negValidDir = (u < (-1 * uTol))
					
					posValid_UB_dir = (distUb > dTol) & (u > uTol)
					negValid_UB_dir = (distUb > dTol) & (u < (-1 * uTol))
					
					posValid_LB_dir = (distLb > dTol) & (u > uTol)
					negValid_LB_dir = (distLb > dTol) & (u < (-1 * uTol))
					
					# Find the largest allowed step size in both the positive and negative directions
					# How many lengths of u can you move in any valid direction before you hit a wall
					pos_u_rel_dis_ub = np.compress(posValid_UB_dir, distUb, axis=0) / np.compress(posValid_UB_dir, u,
					                                                                              axis=0)
					neg_u_rel_dis_ub = np.compress(negValid_UB_dir, distUb, axis=0) / np.compress(negValid_UB_dir, u,
					                                                                              axis=0)
					
					# additional upperbound constraint to prevent exceeding imposed flux limit
					# Cases:
					# sum u is negative
					# sum u is positive
					# sum u is small
					# if max_combined_flux != -1:
					#	current_flux = np.sum(prev_point)
					#	max_add_from_u = max_combined_flux - current_flux
					
					# This might look weird, but we are just seeing how many unit vectors we can move before
					# we hit the flux limit, if that number is positive and if it is smaller than
					# the first positive distance to an upper bound it replaces it as it is more restrictive.
					# later on the code will check to see if it is more restrictive than the other ub distances
					
					# If we have to move a negative amount in the direction of the unit vector the code
					# also just checks to see if that is more restrictive than the first negative distance to an
					# upper bound
					#	rel_u_to_add = max_add_from_u/np.sum(u)
					#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
					#		pos_u_rel_dis_ub[0] = rel_u_to_add
					#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
					#		neg_u_rel_dis_ub[0] = rel_u_to_add
					
					# Find the smallest allowed step size in both the positive and negative directions
					neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u,
					                                                                               axis=0)
					pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u,
					                                                                               axis=0)
					
					# find the current value of the point
					
					# Find the smallest and largest allowed step sizes
					maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb))) - step_tol
					print(maxStep)
					input()
					maxStep = max(maxStep, 0)
					minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub))) + step_tol
					minStep = min(minStep, 0)
					
					# if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
					#	print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
					#	continue
					
					# grab random value from pre-generated list for the step distance
					stepDist = randVector[0] * (maxStep - minStep) + minStep
					
					curPoint = prev_point + stepDist * u
					
					stalls = 0
					while (np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0)) or np.max(
							np.abs(np.matmul(S, curPoint))) > Sv_cont_tol:
						curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
						overInd = indexer[curPoint > ub]
						underInd = indexer[curPoint < lb]
						# print(overInd)
						# print(underInd)
						# scaling might be useful here
						for j in overInd:
							curPoint[j] = ub[j]
						for j in underInd:
							curPoint[j] = lb[j]
						stalls += 1
						if stalls % 100 == 0:
							print(stalls)
							print(f"Curr WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
							print(f"Curr Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
					
					cur_support = self.gene_point_penalty(curPoint[:-1], RPS, gene_penalty_mod)
					prev_support = self.gene_point_penalty(prev_point[:-1], RPS, gene_penalty_mod)
					
					proposed_points += 1
					# print(cur_support,prev_support,np.log(randVector[1]))
					if proposed_points % 500 == 0:
						print(proposed_points, total_steps, (total_steps + 1) / proposed_points)
						print(prev_point[-1], curPoint[-1])
						print(cur_support, prev_support)
						print(f"CURR WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"CURR Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
						print(f"PREV WARNING: fidErr {np.max(np.abs(np.matmul(S, prev_point)))}")
						print(f"PREV Error ub: {max(prev_point - ub)}, lb: {max(lb - prev_point)}")
					
					if np.log(randVector[1]) < (cur_support - prev_support):
						passing = False
						continue
					else:
						continue
				
				if total_steps % 200 == 0:
					if max(max(curPoint - ub), max(lb - curPoint)) > 1e-7:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
				
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				# curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(curPoint[-1])
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	def HRSampler_gene_bias_lincomb_pinch(self, origin_and_warmup_points, number_of_points, stepsPerPoint, RPS,
	                                      max_sampling_level, pinch_tol=1e-6,
	                                      gene_penalty_mod=1, thermo_const=None, rxn_penalty_vec=None,term_val = 1e7, **kwargs):
		# Note the gene_penalty_mod applies to the value in the exponent
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		
		# potentially add check for how much the realignment moves a point
		
		lb_full, ub_full, S_full, b, rxn_list, met_list = self.dicts_to_mats()
		print(np.sum(lb_full),np.sum(ub_full),np.sum(b),np.sum(S_full))
		var_ind, pinch_b, lb_full, v_c, ub_full = self.quick_pinch(max_sampling_level, tol=pinch_tol)
		print(list(var_ind))
		input()
		print(np.shape(lb_full))
		flux_min = lb_full[-1]
		flux_max = ub_full[-1]
		
		lb_ntf = lb_full[:-1]
		ub_ntf = ub_full[:-1]
		
		print(np.shape(v_c))
		print(np.shape(lb_full))
		v_c = v_c[:-1]
		v_c[var_ind] = np.ones_like(var_ind) * np.nan
		
		lb_full[-1] += np.nansum(v_c)
		ub_full[-1] += np.nansum(v_c)
		
		v_c_fillin = cp.deepcopy(v_c)
		lb = lb_ntf[var_ind]
		ub = ub_ntf[var_ind]
		S = S_full[:, var_ind]
		RPS = RPS[var_ind]
		b = pinch_b
		
		#print(self.model_dict["total_flux_limits"])
		print(-np.nansum(v_c))
		print(v_c)
		print(b)
		#input()
		
		if thermo_const is not None:
			simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = self.precomp_positive_flux_point_to_bi_dir(
				rxn_list, prune_specific=thermo_const["prune_specific"])
			# test_bi_S = self.positive_S_to_bi_dir(S_full, simp_neg_ind, cutter, exchange_cutter=exchange_cutter)
			
			NS = thermo_const["NS_internal"]
			trans_NS = np.transpose(NS)
		
		origin_and_warmup_points_1 = origin_and_warmup_points[:, var_ind]
		origin_and_warmup_points_2 = np.reshape(np.sum(origin_and_warmup_points_1, axis=1), (-1, 1))
		
		origin_and_warmup_points = np.hstack((origin_and_warmup_points_1, origin_and_warmup_points_2))
		
		maxMinTol = 1e-9
		uTol = 1e-7
		dTol = 0
		step_tol = 1e-9
		prev = 0
		
		Sv_tol = 1e-9
		# For warmup slightly going above bounds (~1e-9) does not seem overly pressing
		# For HR a point which sits slightly above an upper bound and slightly below a lower bound
		# has no meaningful way to move. So we must ensure these points sit slightly within bounds.
		# Dimensions with close boundaries have been pinched.
		
		# The below value comes into play when a point sits above a bound, the system attempts to find a point
		# which is close to the current invalid point which sits below upper_bound-lb_ub_tol_barrier and above
		# lower_bound + lb_ub_tol_barrier (the barrier is essential due to numerical instability)
		
		# Numerical issues could still persist, this is not a sure fix, but should hopefully ensure an already
		# rare edge case never practically shows up
		lb_ub_tol_barrier = 1e-10
		
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = origin_and_warmup_points[0]
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "total_flux_limits" in self.model_dict.keys():
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S_full)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			
			S_full = np.vstack((S_full, rxn_penalty_vec))
			S_full = np.hstack((S_full, np.zeros((np.shape(S_full)[0], 1))))
			S_full[-1, -1] = -1
			
			rxn_penalty_vec = rxn_penalty_vec[0, var_ind]
			
			print(np.shape(rxn_penalty_vec))
			# rxn_penalty_vec = rxn_penalty_vec[var_ind]
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, flux_min))
			ub = np.hstack((ub, flux_max + 1e-3))
			print(ub[-1])
			print(lb[-1])
			print("bleh")
			# print(lb_full)
			b = np.hstack((b, 0))
			S[-1, -1] = -1
		
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
			fs = prev_point[-1]
			prev_point = prev_point[var_ind]
			prev_point = np.hstack((prev_point, fs))
		else:
			# find center
			valid_point = True
			move_forward = False
			count = 0
			mid_point = np.mean(origin_and_warmup_points, axis=0)
			# prev_point = origin_and_warmup_points[np.random.choice(np.shape(warmup_points)[0])]
			print(np.sum(S))
			print(np.sum(b))
			print(f"spCurr WARNING: fidErr {np.max(np.abs(np.matmul(S, mid_point) - b))}")
			print(f"spCurr Error ub: {max(mid_point - ub)}, lb: {max(lb - mid_point)}")
			while not move_forward:
				count += 1
				if count % 100 == 0:
					print(count)
				if self.model_dict["gurobi_token"] is not None:
					model = gp.Model("init_model" + self.model_dict["model_name"],
					                 env=self.model_dict["gurobi_token"])
				else:
					model = gp.Model("init_model" + self.model_dict["model_name"])
				react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb,
				                           ub=ub)
				model.addConstr(S @ react_flux == b, name="c")
				model.Params.LogToConsole = 0
				
				n_point_try = 1
				prev_point = np.zeros((n_point_try, S.shape[1]))
				for n_point in range(n_point_try):
					
					accept_point = False
					while not accept_point:
						print(f"spCurr Error ub: {max(mid_point - ub)}, lb: {max(lb - mid_point)}")
						#input()
						obj = (2 * np.random.random_sample(S.shape[1]) - 1)
						# for variable objective total flux
						# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
						model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
						
						model.Params.NumericFocus = 3
						model.Params.FeasibilityTol = 1e-9
						model.optimize()
						if model.Status == 2:
							if np.min(react_flux.X) >= 0:
								accept_point = True
							else:
								model.reset()
						else:
							model.reset()
					# prev_point[n_point] = react_flux.X*(1-1e-15) + mid_point*1e-15
					prev_point[n_point] = react_flux.X
				prev_point = np.average(prev_point, axis=0)
				prev_point = warmup_points[np.random.choice(np.shape(warmup_points)[0])]
				
				# prev_point = mid_point
				# print(np.sum(prev_point**2))
				
				if thermo_const is not None:
					v_c_fillin[var_ind] = prev_point[:-1]
					# print(np.shape(v_c_fillin))
					# print(np.shape(simp_neg_ind))
					# print(np.shape(bi_point))
					ent_flux_point = self.positive_flux_point_to_bi_dir(v_c_fillin, simp_neg_ind,
					                                                    comp_neg_ind, comp_perm,
					                                                    cutter,
					                                                    exchange_cutter=exchange_cutter)
					valid_point = self.generate_therm_model_new(trans_NS, ent_flux_point) == 2
					if valid_point:
						print(1, valid_point, np.min(prev_point) >= 0)
				
				if valid_point and ((np.max(np.abs(np.matmul(S, prev_point) - b)) > Sv_tol) or (
						max(prev_point - ub) > lb_ub_tol_barrier or max(lb - prev_point) > lb_ub_tol_barrier)):
					print(f"mpWARNING: fidErr {np.max(np.abs(np.matmul(S, prev_point) - b))}")
					print(f"mpWARNING: Error ub: {max(prev_point - ub)}, lb: {max(lb - prev_point)}")
					print(np.shape(prev_point))
					print(np.argmin(prev_point))
					if self.model_dict["gurobi_token"] is not None:
						model = gp.Model("correct_model" + self.model_dict["model_name"],
						                 env=self.model_dict["gurobi_token"])
					else:
						model = gp.Model("correct_model" + self.model_dict["model_name"])
					react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
					# print(react_flux)
					model.addConstr(S @ react_flux == b, name="c")
					# (mf - cpf)^2 = mf^2 - 2*mf*cpf + cpf^2
					Q = np.eye(np.size(prev_point))
					model.setObjective(
						react_flux @ Q @ react_flux - (2 * prev_point) @ react_flux + np.sum(prev_point ** 2),
						GRB.MINIMIZE)
					model.Params.LogToConsole = 0
					model.Params.NumericFocus = 3
					model.Params.FeasibilityTol = 1e-9
					
					model.optimize()
					prev_point = react_flux.X
				# print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, prev_point) - b))}")
				# print(f"WARNING: Error ub: {max(prev_point - ub)}, lb: {max(lb - prev_point)}")
				if valid_point and thermo_const is not None:
					v_c_fillin[var_ind] = prev_point[:-1]
					# print(np.shape(v_c_fillin))
					# print(np.shape(simp_neg_ind))
					# print(np.shape(bi_point))
					ent_flux_point = self.positive_flux_point_to_bi_dir(v_c_fillin, simp_neg_ind,
					                                                    comp_neg_ind, comp_perm,
					                                                    cutter,
					                                                    exchange_cutter=exchange_cutter)
					
					if np.min(prev_point) >= 0:
						move_forward = self.generate_therm_model_new(trans_NS, ent_flux_point) == 2
						print(2, move_forward)
				elif valid_point:
					move_forward = True
		
		print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, prev_point) - b))}")
		print(f"WARNING: Error ub: {max(prev_point - ub)}, lb: {max(lb - prev_point)}")
		print(np.min(prev_point))
		print(prev_point)
		print("ope")
		warmup_vecs = np.zeros_like(warmup_points)
		norm_warmup_vecs = np.zeros_like(warmup_points)
		for i in range(np.shape(warmup_points)[0]):
			warmup_vecs[i] = warmup_points[i] - origin
			norm_warmup_vecs[i] = self.normalize(warmup_points[i] - origin)
		# print(list(norm_warmup_vecs[i]), np.shape(warmup_points)[0])
		# input()
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		proposed_points = 0
		points = np.zeros((number_of_points + 1, len(origin)))
		points[0] = prev_point
		point_count = 1
		therm_reject = 0
		gene_reject = 0
		correct_count = 0
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				passing = True
				ts = 0
				vtu_list = []
				while passing:
					prev_point = np.where(prev_point > ub, ub, prev_point)
					prev_point = np.where(prev_point < lb, lb, prev_point)
					
					distUb = (ub - prev_point)
					distLb = (prev_point - lb)
					
					
					
					# test = np.random.multinomial(1+np.random.choice(np.shape(norm_warmup_vecs)[0],1)[0],np.ones(np.shape(norm_warmup_vecs)[0])/np.shape(norm_warmup_vecs)[0])
					
					# u = (warmup_vecs.T * test).T
					# u = self.normalize(np.sum(u, axis=0))
					
					# print(np.shape(u))
					count = 0
					biggest_step = -1
					u_use = norm_warmup_vecs[0]
					max_list = []
					min_list = []
					val_to_use = np.random.choice(np.shape(norm_warmup_vecs)[0])
					u = norm_warmup_vecs[val_to_use]
					posValidDir = (u > uTol)
					negValidDir = (u < (-1 * uTol))
					
					# Find the largest allowed step size in both the positive and negative directions
					# How many lengths of u can you move in any valid direction before you hit a wall
					pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u,
					                                                                          axis=0)
					neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u,
					                                                                          axis=0)
					
					# additional upperbound constraint to prevent exceeding imposed flux limit
					# Cases:
					# sum u is negative
					# sum u is positive
					# sum u is small
					# if max_combined_flux != -1:
					#	current_flux = np.sum(prev_point)
					#	max_add_from_u = max_combined_flux - current_flux
					
					# This might look weird, but we are just seeing how many unit vectors we can move before
					# we hit the flux limit, if that number is positive and if it is smaller than
					# the first positive distance to an upper bound it replaces it as it is more restrictive.
					# later on the code will check to see if it is more restrictive than the other ub distances
					
					# If we have to move a negative amount in the direction of the unit vector the code
					# also just checks to see if that is more restrictive than the first negative distance to an
					# upper bound
					#	rel_u_to_add = max_add_from_u/np.sum(u)
					#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
					#		pos_u_rel_dis_ub[0] = rel_u_to_add
					#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
					#		neg_u_rel_dis_ub[0] = rel_u_to_add
					
					# Find the smallest allowed step size in both the positive and negative directions
					neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u,
					                                                                               axis=0)
					pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u,
					                                                                               axis=0)
					
					# find the current value of the point
					
					# Find the smallest and largest allowed step sizes
					# print(u,test)
					# print(pos_u_rel_dis_ub,neg_u_rel_dis_lb)
					maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
					
					minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
					if maxStep == 0 and minStep == 0:
						if ts%5000 ==0:
							print("0 bounds, count: ", ts)
						ts += 1
						if ts > term_val:
							print("TERMINATING DUE TO NO VALID POINTS")
							return np.array([])
						continue
					if maxStep - minStep <= 0:
						print("WARNING")
						print(maxStep, minStep)
						print(maxStep - minStep)
						min_ind = np.argmin(pos_u_rel_dis_ub)
						print(np.compress(posValidDir, distUb, axis=0)[min_ind])
						print(np.compress(posValidDir, u, axis=0)[min_ind])
						print(np.compress(posValidDir, distUb / u, axis=0)[min_ind])
						# print(np.compress(posValidDir, u, axis=0))
						# print(np.compress(posValidDir, distUb, axis=0)/np.compress(posValidDir, u, axis=0))
						
						min_ind = np.argmin(neg_u_rel_dis_lb)
						print(np.compress(negValidDir, distUb, axis=0)[min_ind])
						print(np.compress(negValidDir, u, axis=0)[min_ind])
						
						max_ind = np.argmax(pos_u_rel_dis_lb)
						print(np.compress(posValidDir, distLb, axis=0)[max_ind])
						print(np.compress(posValidDir, u, axis=0)[max_ind])
						
						max_ind = np.argmax(neg_u_rel_dis_ub)
						print(np.compress(negValidDir, distLb, axis=0)[max_ind])
						print(np.compress(negValidDir, u, axis=0)[max_ind])
					# input()
					randVector = np.random.random_sample(2)
					# if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
					#	print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
					#	continue
					# grab random value from pre-generated list for the step distance
					stepDist = randVector[0] * (maxStep - minStep) + minStep
					
					curPoint = prev_point + stepDist * u
					
					if (np.max(np.abs(np.matmul(S, curPoint) - b)) > Sv_tol) or max(
							curPoint - ub) > lb_ub_tol_barrier or max(
							lb - curPoint) > lb_ub_tol_barrier:
						correct_count += 1
						if correct_count % 50 == 0:
							print(f"ggg WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint) - b))}")
							print(f"ggg WARNING: Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
						if self.model_dict["gurobi_token"] is not None:
							model = gp.Model("correct_model" + self.model_dict["model_name"],
							                 env=self.model_dict["gurobi_token"])
						else:
							model = gp.Model("correct_model" + self.model_dict["model_name"])
						react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux",
						                           lb=lb + lb_ub_tol_barrier, ub=ub - lb_ub_tol_barrier)
						# print(react_flux)
						model.addConstr(S @ react_flux == b, name="c")
						# (mf - cpf)^2 = mf^2 - 2*mf*cpf + cpf^2
						Q = np.eye(np.size(curPoint))
						model.setObjective(
							react_flux @ Q @ react_flux - (2 * curPoint) @ react_flux + np.sum(curPoint ** 2),
							GRB.MINIMIZE)
						model.Params.LogToConsole = 0
						model.Params.NumericFocus = 3
						model.Params.FeasibilityTol = 1e-9
						model.optimize()
						curPoint = react_flux.X
						if correct_count % 50 == 0:
							print(f"rrr WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint) - b))}")
							print(f"rrr WARNING: Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
					
					cur_support = self.gene_point_penalty(curPoint[:-1], RPS, gene_penalty_mod)
					prev_support = self.gene_point_penalty(prev_point[:-1], RPS, gene_penalty_mod)
					
					proposed_points += 1
					# print(cur_support,prev_support,np.log(randVector[1]))
					if proposed_points % 5000 == 0:
						print(proposed_points, total_steps, (total_steps + 1) / proposed_points,
						      therm_reject / proposed_points, gene_reject / proposed_points)
						print(prev_point[-1], curPoint[-1], lb[-1], ub[-1])
						print(cur_support, prev_support)
						print(f"CURR WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint) - b))}")
						print(f"CURR Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
						print(f"PREV WARNING: fidErr {np.max(np.abs(np.matmul(S, prev_point) - b))}")
						print(f"PREV Error ub: {max(prev_point - ub)}, lb: {max(lb - prev_point)}")
					
					if np.min(curPoint) >= 0:
						if np.log(randVector[1]) < (cur_support - prev_support):
							if thermo_const is not None:
								v_c_fillin[var_ind] = curPoint[:-1]
								# print(np.shape(v_c_fillin))
								# print(np.shape(simp_neg_ind))
								# print(np.shape(bi_point))
								ent_flux_point = self.positive_flux_point_to_bi_dir(v_c_fillin, simp_neg_ind,
								                                                    comp_neg_ind, comp_perm,
								                                                    cutter,
								                                                    exchange_cutter=exchange_cutter)
								val = self.generate_therm_model_new(trans_NS, ent_flux_point)
								# print(val)
								if val == 2:
									passing = False
								else:
									therm_reject += 1
									continue
							else:
								passing = False
								continue
						else:
							gene_reject += 1
							continue
					else:
						continue
				
				if total_steps % 200 == 0:
					if max(max(curPoint - ub), max(lb - curPoint)) > 1e-7:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint) - b))}")
						print(f"Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
				
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				# curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint) - b))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(curPoint[-1])
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count] = curPoint
			
			point_count += 1
		print(np.shape(v_c))
		print(v_c)
		full_total_points = np.zeros((np.shape(points)[0], np.size(v_c) + 1))
		for i in range(np.shape(full_total_points)[0]):
			v_c[var_ind] = points[i, :-1]
			full_total_points[i, :-1] = v_c
			full_total_points[i, -1] = np.sum(full_total_points[i, :-1])
			print(np.shape(v_c))
			print(np.shape(S_full), np.shape(full_total_points[i]))
			if (np.max(np.abs(np.matmul(S_full, full_total_points[i]))) > Sv_tol) or max(
					full_total_points[i] - ub_full) > lb_ub_tol_barrier or max(
				lb_full - full_total_points[i]) > lb_ub_tol_barrier:
				print(full_total_points[i, -1], ub_full[-1])
				print(f"ffsCurr WARNING: fidErr {np.max(np.abs(np.matmul(S_full, full_total_points[i])))}")
				# diff = full_total_points[i] - ub_full
				# print(diff)
				# diff = np.sort(diff)
				# print(diff)
				print(
					f"Curr Error ub: {max(full_total_points[i] - ub_full)}, lb: {max(lb_full - full_total_points[i])}")
		print(3)
		return full_total_points
	
	def HRSampler_gene_bias_lincomb_therm(self, origin_and_warmup_points, number_of_points, stepsPerPoint, RPS,
	                                      gene_penalty_mod=1, rxn_penalty_vec=None, **kwargs):
		# Note the gene_penalty_mod applies to the value in the exponent
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-15
		dTol = 1e-9
		Sv_cont_tol = 1e-7
		prev = 0
		
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = origin_and_warmup_points[0]
		
		print(rxn_list)
		print(self.precomp_positive_flux_point_to_bi_dir(rxn_list))
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			# find center
			prev_point = np.mean(origin_and_warmup_points, axis=0)
		
		if "total_flux_limits" in self.model_dict.keys():
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.model_dict["total_flux_limits"])))
			ub = np.hstack((ub, max(self.model_dict["total_flux_limits"])))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		N = sc.linalg.null_space(S)
		indexer = np.arange(0, np.size(lb))
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		warmup_vecs = np.zeros_like(warmup_points)
		for i in range(np.shape(warmup_points)[0]):
			warmup_vecs[i] = warmup_points[i] - origin
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		proposed_points = 0
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				passing = True
				while passing:
					
					randVector = np.random.random_sample(2)
					
					u = (warmup_vecs.T * np.random.random((np.shape(warmup_points)[0], 1))[:, 0]).T
					u = self.normalize(np.sum(u, axis=0))
					
					# Figure out the distances to upper and lower bounds
					distUb = (ub - prev_point)
					distLb = (prev_point - lb)
					
					# Figure out if we are too close to a boundary
					# validDir = ((distUb > dTol) & (distLb > dTol))
					
					# Finds list of entries of unit vector which are either positive or negative
					# and which fall within tolerance values
					posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
					negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
					
					# Find the largest allowed step size in both the positive and negative directions
					# How many lengths of u can you move in any valid direction before you hit a wall
					pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
					neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
					
					# additional upperbound constraint to prevent exceeding imposed flux limit
					# Cases:
					# sum u is negative
					# sum u is positive
					# sum u is small
					# if max_combined_flux != -1:
					#	current_flux = np.sum(prev_point)
					#	max_add_from_u = max_combined_flux - current_flux
					
					# This might look weird, but we are just seeing how many unit vectors we can move before
					# we hit the flux limit, if that number is positive and if it is smaller than
					# the first positive distance to an upper bound it replaces it as it is more restrictive.
					# later on the code will check to see if it is more restrictive than the other ub distances
					
					# If we have to move a negative amount in the direction of the unit vector the code
					# also just checks to see if that is more restrictive than the first negative distance to an
					# upper bound
					#	rel_u_to_add = max_add_from_u/np.sum(u)
					#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
					#		pos_u_rel_dis_ub[0] = rel_u_to_add
					#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
					#		neg_u_rel_dis_ub[0] = rel_u_to_add
					
					# Find the smallest allowed step size in both the positive and negative directions
					neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u,
					                                                                               axis=0)
					pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u,
					                                                                               axis=0)
					
					# find the current value of the point
					
					# Find the smallest and largest allowed step sizes
					maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
					minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
					
					if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
						print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
						continue
					
					# grab random value from pre-generated list for the step distance
					stepDist = randVector[0] * (maxStep - minStep) + minStep
					
					curPoint = prev_point + stepDist * u
					
					stalls = 0
					while (np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0)) or np.max(
							np.abs(np.matmul(S, curPoint))) > Sv_cont_tol:
						curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
						overInd = indexer[curPoint > ub]
						underInd = indexer[curPoint < lb]
						# print(overInd)
						# print(underInd)
						# scaling might be useful here
						for j in overInd:
							curPoint[j] = ub[j]
						for j in underInd:
							curPoint[j] = lb[j]
						stalls += 1
						if stalls % 100 == 0:
							print(stalls)
							print(f"Curr WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
							print(f"Curr Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
					
					cur_support = self.gene_point_penalty(curPoint[:-1], RPS, gene_penalty_mod)
					prev_support = self.gene_point_penalty(prev_point[:-1], RPS, gene_penalty_mod)
					
					proposed_points += 1
					# print(cur_support,prev_support,np.log(randVector[1]))
					if proposed_points % 500 == 0:
						print(proposed_points, total_steps, (total_steps + 1) / proposed_points)
						print(prev_point[-1], curPoint[-1])
						print(cur_support, prev_support)
						print(f"CURR WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"CURR Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
						print(f"PREV WARNING: fidErr {np.max(np.abs(np.matmul(S, prev_point)))}")
						print(f"PREV Error ub: {max(prev_point - ub)}, lb: {max(lb - prev_point)}")
					
					int_point = self.positive_flux_point_to_bi_dir(curPoint[:-1], simp_neg_ind, comp_neg_ind, comp_perm,
					                                               cutter,
					                                               exchange_cutter=exchange_cutter)
					if self.generate_therm_model(NS, int_point) == 2:
						# print("accepted point")
						points[point_count - 1] = curPoint
						point_count += 1
						if np.log(randVector[1]) < (cur_support - prev_support):
							passing = False
							continue
						else:
							continue
					else:
						# print("rejected point")
						continue
				
				if total_steps % 200 == 0:
					if max(max(curPoint - ub), max(lb - curPoint)) > 1e-7:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
				
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				# curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(curPoint[-1])
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	def HRSampler_lincom(self, origin_and_warmup_points, number_of_points, stepsPerPoint, rxn_penalty_vec=None,
	                     **kwargs):
		# Description:
		#   Performs the coordinate direction hit and run algorithm using a linear combination of coordinate vectors
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		
		# Numerical stability is a problem, potentially do a post project purge on bounds and Sv values
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		Sv_error_warning = 1e-6
		prev = 0
		
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		origin = origin_and_warmup_points[0]
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			# Making the origin the center of the warmup points helps to avoid starting it on a boundary
			centerPoint = np.mean(origin_and_warmup_points, axis=0)
			prev_point = centerPoint
		
		if "total_flux_limits" in self.model_dict.keys():
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.model_dict["total_flux_limits"])))
			ub = np.hstack((ub, max(self.model_dict["total_flux_limits"])))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		warmup_vecs = np.zeros_like(warmup_points)
		for i in range(np.shape(warmup_points)[0]):
			warmup_vecs[i] = warmup_points[i] - origin
		
		# Replace this with load call for matlab Null Space
		N = sc.linalg.null_space(S)
		
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			randVector = np.random.random_sample(stepsPerPoint)
			
			step_count = 1
			rand_vec = np.random.random((np.shape(warmup_points)[0], stepsPerPoint))
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				
				# Get a direction from the basis point to the warmup point
				u = (warmup_vecs.T * rand_vec[:, step_count - 1]).T
				u = self.normalize(np.sum(u, axis=0))
				
				# Figure out the distances to upper and lower bounds
				distUb = (ub - prev_point)
				distLb = (prev_point - lb)
				
				# Figure out if we are too close to a boundary
				# validDir = ((distUb > dTol) & (distLb > dTol))
				
				# Finds list of entries of unit vector which are either positive or negative
				# and which fall within tolerance values
				posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
				negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
				
				# Find the largest allowed step size in both the positive and negative directions
				# How many lengths of u can you move in any valid direction before you hit a wall
				pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
				neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
				
				# additional upperbound constraint to prevent exceeding imposed flux limit
				# Cases:
				# sum u is negative
				# sum u is positive
				# sum u is small
				# if max_combined_flux != -1:
				#	current_flux = np.sum(prev_point)
				#	max_add_from_u = max_combined_flux - current_flux
				
				# This might look weird, but we are just seeing how many unit vectors we can move before
				# we hit the flux limit, if that number is positive and if it is smaller than
				# the first positive distance to an upper bound it replaces it as it is more restrictive.
				# later on the code will check to see if it is more restrictive than the other ub distances
				
				# If we have to move a negative amount in the direction of the unit vector the code
				# also just checks to see if that is more restrictive than the first negative distance to an
				# upper bound
				#	rel_u_to_add = max_add_from_u/np.sum(u)
				#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
				#		pos_u_rel_dis_ub[0] = rel_u_to_add
				#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
				#		neg_u_rel_dis_ub[0] = rel_u_to_add
				
				# Find the smallest allowed step size in both the positive and negative directions
				neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u, axis=0)
				pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u, axis=0)
				
				# posValidDir = (u > uTol)
				# negValidDir = (u < (-1 * uTol))
				
				# pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
				# neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
				
				# neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u, axis=0)
				# pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u, axis=0)
				
				# find the current value of the point
				
				# Find the smallest and largest allowed step sizes
				maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
				minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
				
				# if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
				#	print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
				#	continue
				
				# grab random value from pre-generated list for the step distance
				stepDist = randVector[step_count - 1] * (maxStep - minStep) + minStep
				curPoint = prev_point + stepDist * u
				
				if total_steps % 100 == 0:
					if np.max(np.abs(np.matmul(S, curPoint))) > Sv_error_warning:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				
				if total_steps % 200 == 0:
					if max(max(curPoint - ub), max(lb - curPoint)) > 1e-7:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
						print(f"Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
				
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				if (np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0)):
					indexer = np.arange(0, np.size(curPoint))
					overInd = indexer[curPoint > ub]
					underInd = indexer[curPoint < lb]
					# print(overInd)
					# print(underInd)
					# scaling might be useful here
					for j in overInd:
						curPoint[j] = ub[j]
					for j in underInd:
						curPoint[j] = lb[j]
					curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					time_per_step = time_elapsed / total_steps
					prev = total_steps / total_steps_needed * 100
					print(curPoint[-1])
					print(total_steps / total_steps_needed)
					print(time_per_step)
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - total_steps) * time_per_step / 60}")
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	def HRSampler_therm(self, origin_and_warmup_points, number_of_points, stepsPerPoint, skip_therm, NS, simp_neg_ind,
	                    comp_neg_ind, comp_perm, cutter, exchange_cutter=None, **kwargs):
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		Sv_cont_tol = 1e-9
		prev = 0
		
		# find center
		centerPoint = np.mean(origin_and_warmup_points, axis=0)
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = centerPoint
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			prev_point = origin
		
		if self.total_flux_limits is not None:
			S = np.vstack((S, np.ones((1, np.shape(S)[1]))))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.total_flux_limits)))
			ub = np.hstack((ub, max(self.total_flux_limits)))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		N = sc.linalg.null_space(S)
		
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		stepsPerPoint_cp = copy.deepcopy(stepsPerPoint)
		
		first_pass = True
		
		while point_count <= number_of_points:
			if first_pass:
				stepsPerPoint = stepsPerPoint_cp
			else:
				stepsPerPoint = skip_therm
			# draws several random numbers at a time for improved performance
			randVector = np.random.random_sample(stepsPerPoint)
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				random_row_index = np.random.choice(np.shape(warmup_points)[0])
				rand_point = warmup_points[random_row_index]
				
				# Get a direction from the basis point to the warmup point
				u = self.normalize(rand_point - origin)
				
				# Figure out the distances to upper and lower bounds
				distUb = (ub - prev_point)
				distLb = (prev_point - lb)
				
				# Figure out if we are too close to a boundary
				# validDir = ((distUb > dTol) & (distLb > dTol))
				
				# Finds list of entries of unit vector which are either positive or negative
				# and which fall within tolerance values
				posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
				negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
				
				# Find the largest allowed step size in both the positive and negative directions
				# How many lengths of u can you move in any valid direction before you hit a wall
				pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
				neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
				
				# additional upperbound constraint to prevent exceeding imposed flux limit
				# Cases:
				# sum u is negative
				# sum u is positive
				# sum u is small
				# if max_combined_flux != -1:
				#	current_flux = np.sum(prev_point)
				#	max_add_from_u = max_combined_flux - current_flux
				
				# This might look weird, but we are just seeing how many unit vectors we can move before
				# we hit the flux limit, if that number is positive and if it is smaller than
				# the first positive distance to an upper bound it replaces it as it is more restrictive.
				# later on the code will check to see if it is more restrictive than the other ub distances
				
				# If we have to move a negative amount in the direction of the unit vector the code
				# also just checks to see if that is more restrictive than the first negative distance to an
				# upper bound
				#	rel_u_to_add = max_add_from_u/np.sum(u)
				#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
				#		pos_u_rel_dis_ub[0] = rel_u_to_add
				#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
				#		neg_u_rel_dis_ub[0] = rel_u_to_add
				
				# Find the smallest allowed step size in both the positive and negative directions
				neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u, axis=0)
				pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u, axis=0)
				
				# find the current value of the point
				
				# Find the smallest and largest allowed step sizes
				maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
				minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
				
				if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
					print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
					continue
				
				# grab random value from pre-generated list for the step distance
				stepDist = randVector[step_count - 1] * (maxStep - minStep) + minStep
				curPoint = prev_point + stepDist * u
				
				if total_steps % 100 == 0:
					curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				
				if total_steps % 200000 == 0 and False:
					print(f"Error {max(curPoint - ub)}, {max(lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				if np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0):
					indexer = np.arange(0, np.size(curPoint))
					overInd = indexer[curPoint > ub]
					underInd = indexer[curPoint < lb]
					# print(overInd)
					# print(underInd)
					# scaling might be useful here
					for j in overInd:
						curPoint[j] = ub[j]
					for j in underInd:
						curPoint[j] = lb[j]
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			int_point = self.positive_flux_point_to_bi_dir(curPoint[:-1], simp_neg_ind, comp_neg_ind, comp_perm, cutter,
			                                               exchange_cutter=exchange_cutter)
			if self.generate_therm_model(NS, int_point) == 2:
				# print("accepted point")
				points[point_count - 1] = curPoint
				point_count += 1
				first_pass = True
			else:
				# print("rejected point")
				first_pass = False
		
		return points
	
	def precomp_positive_flux_point_to_bi_dir(self, pos_rxn_list, prune_exchange=True, prune_specific=None):
		# pos_rxn_list is used in stead of self.rxn_dict.keys() so that a desired order of keys may be passed
		if prune_specific is None:
			prune_specific = []
		elif not prune_exchange:
			print("prune_exchange is False, but prune specific contains entries, this is a contradiction")
		
		if np.alltrue(np.array(["_inverted" not in i for i in pos_rxn_list])):
			print("WARNING: NO INVERTED RXNS FOUND, WRONG MODEL MAY BE IN USE")
		zeros = np.zeros(len(pos_rxn_list), dtype=float)
		ones = np.ones(len(pos_rxn_list), dtype=float)
		
		cutter = np.nonzero(np.where(
			np.array([("_inverted" in i and i.replace("_inverted", "") in pos_rxn_list) for i in pos_rxn_list]),
			zeros, ones))[0]
		
		simp_neg_ind = np.where(
			np.array([("_inverted" in i and i.replace("_inverted", "") not in pos_rxn_list) for i in pos_rxn_list]),
			-1 * ones, ones)
		comp_neg_ind = []
		comp_perm = []
		for i in range(len(pos_rxn_list)):
			# slowly read through to make sure this is doing what it is supposed to do
			rxn_name = pos_rxn_list[i]
			# if reaction is forward but its inverse is in list
			# make the comp list 0, but tell the comp perm list what the index of the inverted (-1) reaction is
			# in that way if you do comp_neg_ind[comp_perm] the comp perm will request the -1 from the comp neg ind to come to the
			# location of its positive value. Furthermore if you modify comp_neg_ind to have -1* its flux and run comp_neg_ind[comp_perm]
			# you will have the negative fluxes in the right place to simply add the list and apply the cutter
			if "_inverted" not in rxn_name and rxn_name + "_inverted" in pos_rxn_list:
				comp_neg_ind.append(0.0)
				comp_perm.append(pos_rxn_list.index(rxn_name + "_inverted"))
			elif "_inverted" in rxn_name and rxn_name.replace("_inverted", "") in pos_rxn_list:
				comp_neg_ind.append(-1.0)
				comp_perm.append(-5)
			else:
				comp_neg_ind.append(0.0)
				comp_perm.append(-5)
		for i in range(len(comp_perm)):
			if comp_perm[i] == -5:
				comp_perm[i] = np.argmax(comp_neg_ind)
		comp_perm = np.array(comp_perm)
		comp_neg_ind = np.array(comp_neg_ind)
		pos_rxn_array = np.array(pos_rxn_list)
		cut_rxn_array = pos_rxn_array[cutter]
		if prune_exchange:
			exchange_cutter = np.nonzero(np.where(
				np.array(
					[(i in prune_specific or len(self.model_dict['rxn_dict'][i]["rxn_metabolites"].keys()) == 1) for i
					 in cut_rxn_array]),
				zeros[cutter], ones[cutter]))[0]
		else:
			exchange_cutter = None
		
		pos_rxn_array_post_cut = pos_rxn_array[cutter]
		return simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter
	
	def positive_flux_point_to_bi_dir(self, point, simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter=None):
		bi_point = copy.deepcopy(point)
		comp_neg_ind_cp = copy.deepcopy(comp_neg_ind)
		bi_point *= simp_neg_ind
		comp_neg_ind_cp *= bi_point
		comp_neg_ind_cp = comp_neg_ind_cp[comp_perm]
		bi_point += comp_neg_ind_cp
		nbi_point = bi_point[cutter]
		if exchange_cutter is not None:
			nbi_point = nbi_point[exchange_cutter]
		return nbi_point
	
	def positive_S_to_bi_dir(self, S, simp_neg_ind, cutter, exchange_cutter=None):
		bi_dir_S = copy.deepcopy(S)
		for col_ind in range(np.shape(bi_dir_S)[1]):
			bi_dir_S[:, col_ind] *= simp_neg_ind[col_ind]
		bi_dir_S = bi_dir_S[:, cutter]
		if exchange_cutter is not None:
			bi_dir_S = bi_dir_S[:, exchange_cutter]
		return bi_dir_S
	
	def convert_model_to_positive_flux(self):
		# Only works for getting S
		reaction_names = list(self.model_dict["rxn_dict"].keys())
		for rxn_name in reaction_names:
			if self.model_dict["rxn_dict"][rxn_name]['ub'] <= 0:
				new_rxn = rxn_name + "_inverted"
				nrxn_dict = {}
				for key in self.model_dict["rxn_dict"][rxn_name].keys():
					if key == "rxn_metabolites":
						nrxn_dict["rxn_metabolites"] = {
							i: -1 * self.model_dict["rxn_dict"][rxn_name]['rxn_metabolites'][i] for i in
							self.model_dict["rxn_dict"][rxn_name]['rxn_metabolites'].keys()}
					elif key == "lb":
						nrxn_dict["ub"] = self.model_dict["rxn_dict"][rxn_name]["lb"] * -1
					elif key == "ub":
						nrxn_dict["lb"] = self.model_dict["rxn_dict"][rxn_name]["ub"] * -1
					else:
						nrxn_dict[key] = self.model_dict["rxn_dict"][rxn_name][key]
				del self.model_dict["rxn_dict"][rxn_name]
				self.model_dict["rxn_dict"][new_rxn] = copy.deepcopy(nrxn_dict)
			elif self.model_dict["rxn_dict"][rxn_name]['lb'] < 0:
				new_rxn = rxn_name + "_inverted"
				nrxn_dict = {}
				for key in self.model_dict["rxn_dict"][rxn_name].keys():
					if key == "rxn_metabolites":
						nrxn_dict["rxn_metabolites"] = {
							i: -1 * self.model_dict["rxn_dict"][rxn_name]['rxn_metabolites'][i] for i in
							self.model_dict["rxn_dict"][rxn_name]['rxn_metabolites'].keys()}
					elif key == "lb":
						nrxn_dict["ub"] = self.model_dict["rxn_dict"][rxn_name]["lb"] * -1
					elif key == "ub":
						nrxn_dict["lb"] = 0
					else:
						nrxn_dict[key] = self.model_dict["rxn_dict"][rxn_name][key]
				self.model_dict["rxn_dict"][new_rxn] = copy.deepcopy(nrxn_dict)
				self.update_reaction_bounds(rxn_name, 0, 'keep')
	
	def convert_model_to_bidirectional_flux(self):
		# This currently should only be used for S, its not fully tested for other entries
		rxn_names = list(self.model_dict["rxn_dict"].keys())
		for rxn_name in rxn_names:
			if "_inverted" in rxn_name:
				forward_rxn_name = rxn_name.replace("_inverted", "")
				if forward_rxn_name in rxn_names:
					if self.model_dict["rxn_dict"][forward_rxn_name]["ub"] == \
							self.model_dict["rxn_dict"][forward_rxn_name]["lb"] and \
							self.model_dict["rxn_dict"][forward_rxn_name]["lb"] == 0:
						self.update_reaction_bounds(forward_rxn_name, -1 * self.model_dict["rxn_dict"][rxn_name]["ub"],
						                            -1 * self.model_dict["rxn_dict"][rxn_name]["lb"])
					elif self.model_dict["rxn_dict"][rxn_name]["ub"] == self.model_dict["rxn_dict"][rxn_name]["lb"] and \
							self.model_dict["rxn_dict"][rxn_name]["lb"] == 0:
						pass
					else:
						self.update_reaction_bounds(forward_rxn_name, -1 * self.model_dict["rxn_dict"][rxn_name]["ub"],
						                            "keep")
					del self.model_dict["rxn_dict"][rxn_name]
				else:
					
					self.model_dict["rxn_dict"][forward_rxn_name] = self.model_dict["rxn_dict"].pop(rxn_name)
					lb_ph = copy.deepcopy(self.model_dict["rxn_dict"][forward_rxn_name]["lb"])
					ub_ph = copy.deepcopy(self.model_dict["rxn_dict"][forward_rxn_name]["ub"])
					rxn_S = copy.deepcopy(self.model_dict["rxn_dict"][forward_rxn_name]['rxn_metabolites'])
					self.model_dict["rxn_dict"][forward_rxn_name]["rxn_metabolites"] = {i: -1 * rxn_S[i] for i in
					                                                                    rxn_S.keys()}
					self.update_reaction_bounds(forward_rxn_name, -1 * ub_ph, -1 * lb_ph)
	
	def weight_fluxes(self, eq_bound_tol=0.0, new_bound_range_tol=10):
		max_constrained_flux = -1
		max_constrained_flux_id = -1
		for i in range(len(self.reaction_names)):
			if self.lb[i] >= self.ub[i] - eq_bound_tol:
				if abs(self.lb[i]) > max_constrained_flux:
					max_constrained_flux = abs(self.lb[i])
					max_constrained_flux_id = i
		print(f"New extreme bound is {max_constrained_flux * new_bound_range_tol}")
		print(f"Set based on {self.reaction_names[max_constrained_flux_id]}")
		for i in range(len(self.reaction_names)):
			if self.lb[i] < -1 * new_bound_range_tol * max_constrained_flux:
				self.lb[i] = -1 * new_bound_range_tol * max_constrained_flux
			elif self.ub[i] < 1 * new_bound_range_tol * max_constrained_flux:
				self.ub[i] = 1 * new_bound_range_tol * max_constrained_flux
	
	def find_essential_fluxes_fast(self, list_to_return="index", print_progress=False, **kwargs):
		# A more efficient system to find the essential fluxes
		# The old method is kept for comparison
		
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		essential_reactions = []
		true_lb = copy.deepcopy(lb)
		true_ub = copy.deepcopy(ub)
		
		true_lb_2 = copy.deepcopy(lb)
		true_ub_2 = copy.deepcopy(ub)
		if self.model_dict["gurobi_token"] != None:
			model = gp.Model("fva_model_" + self.model_dict["model_name"], env=self.model_dict["gurobi_token"])
		else:
			model = gp.Model("fva_model_" + self.model_dict["model_name"])
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux",
		                           lb=copy.deepcopy(lb),
		                           ub=copy.deepcopy(ub))
		model.addConstr(S @ react_flux == b, name="c")
		# obj = np.zeros((1, self.S.shape[1]))
		
		# model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.LogToConsole = 0
		model.Params.ScaleFlag = 0
		model.Params.NumericFocus = 3
		model.Params.Method = 1
		# model.Params.Presolve = 0
		# model.Params.Aggregate = 0
		model.Params.FeasibilityTol = 1e-9
		
		index_size = np.size(lb)
		index_select = [i for i in range(index_size)]
		if "seed" in kwargs.keys():
			random.Random(kwargs["seed"]).shuffle(index_select)
		else:
			random.shuffle(index_select)
		for i in range(index_size):
			if print_progress and i % 100 == 0:
				print((index_size - i) / index_size)
				print(len(essential_reactions))
			reaction_index_to_check = index_select[i]
			if lb[reaction_index_to_check] > 0 or ub[reaction_index_to_check] < 0:
				# This additional line prevents an issue where a preset bound is broken but the system is still feasible
				# If a metabolite intake can range between -100 and 100 and we force it to be between 9 and 10
				# setting it to 0 and testing for feasibility is not enough
				is_feasible = False
			else:
				save_lb = true_lb_2[reaction_index_to_check]
				save_ub = true_ub_2[reaction_index_to_check]
				true_lb[reaction_index_to_check] = 0
				true_ub[reaction_index_to_check] = 0
				react_flux.lb = true_lb
				react_flux.ub = true_ub
				model.optimize()
				is_feasible = (model.Status == 2)
				if not is_feasible:
					if model.Status != 3:
						print(model.Status)
						input("Error, unexpected model status")
				model.reset()
			react_flux.lb = true_lb_2
			react_flux.ub = true_ub_2
			if is_feasible:
				continue
			else:
				if not (lb[reaction_index_to_check] > 0 or ub[reaction_index_to_check] < 0):
					true_lb[reaction_index_to_check] = save_lb
					true_ub[reaction_index_to_check] = save_ub
				if list_to_return == "index":
					essential_reactions.append(reaction_index_to_check)
				elif list_to_return == "names":
					essential_reactions.append(reaction_names[reaction_index_to_check])
		if "dict_count" in kwargs.keys():
			initial_flux_counts_copy = copy.deepcopy(kwargs["dict_count"])
			for i in essential_reactions:
				initial_flux_counts_copy[i] += 1
			return initial_flux_counts_copy
		else:
			print(len(essential_reactions))
			return essential_reactions
	
	def fit_to_experimental_data(self, experimental_data_file, alpha_array, assign_opt_alpha_fluxes=True, tol=1e-7,
	                             search_save_path=""):
		# experiment vector = ei
		# model vector = vi
		# (vi-aei)^2/(aei)^2 = vi^2/(aei)^2-2vi/aei+1
		# (vi-aei)^2/(aei+vi)^2 = vi^2/(aei)^2-2vi/aei+1
		# returns:
		# model with best experimental fit for given alpha
		# the relative error (objective term)
		# the index values of the experimental data values
		lb, ub, S, b, ordered_rxn, ordered_met = self.dicts_to_mats()
		experimental_fluxes = np.loadtxt(experimental_data_file, delimiter=",", dtype=str)
		model_index_of_experimental_fluxes = []
		e = np.zeros(len(lb))
		amount_to_add_back = 0
		for i in range(np.shape(experimental_fluxes)[0]):
			if i > 0:
				rxn_name = experimental_fluxes[i, 0]
				flux_values = experimental_fluxes[i, 1].astype('float64')
				index_of_reaction = ordered_rxn.index(rxn_name)
				model_index_of_experimental_fluxes.append(index_of_reaction)
				e[index_of_reaction] = flux_values
				amount_to_add_back += 1
		# For very low values of experimental flux gurobi will run into accuracy errors
		# The current solution is to just take the upper limit of the experimental value of the one
		# poorly behaving flux, but for other cases it might be worth establishing a minimal cutoff point
		# And simply do linear minimaztion on the small values error
		array_of_states = np.zeros((np.size(alpha_array), np.size(lb)))
		# objective with and without valine
		array_of_obj_val_no_val = np.zeros((np.size(alpha_array)))
		array_of_obj_val = np.zeros((np.size(alpha_array)))
		
		min_ind = 0
		min_objective = -1
		for alpha_ind in range(np.size(alpha_array)):
			Q = np.zeros((len(e), len(e)))
			for i in range(len(e)):
				if e[i] != 0:
					Q[i, i] = 1 / (alpha_array[alpha_ind] * e[i]) ** 2
				else:
					Q[i, i] = 0
			c = np.zeros(len(e))
			for i in range(len(e)):
				if e[i] != 0:
					c[i] = 2 / (alpha_array[alpha_ind] * e[i])
				else:
					c[i] = 0
			if self.model_dict["gurobi_token"] != None:
				model = gp.Model("fva_model_" + self.model_dict["model_name"], env=self.model_dict["gurobi_token"])
			else:
				model = gp.Model("fva_model_" + self.model_dict["model_name"])
			model.Params.NumericFocus = 3
			# model.Params.ScaleFlag = 0
			model.Params.BarConvTol = 0
			# model.Params.BarHomogeneous = 0
			model.Params.LogToConsole = 0
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb,
			                           ub=ub)
			model.setObjective(react_flux @ Q @ react_flux - c @ react_flux + amount_to_add_back, GRB.MINIMIZE)
			rhs = np.transpose(b)
			model.addConstr(S @ react_flux == rhs, name="c")
			model.optimize()
			if model.Status == 2 or model.Status == 13:
				array_of_states[alpha_ind] = model.X
				ind = copy.deepcopy(model_index_of_experimental_fluxes)
				labels = [ordered_rxn[i] for i in model_index_of_experimental_fluxes]
				# print(error)
				# print(ind)
				del ind[labels.index('EX_val_L(e)')]
				# print(ind)
				expt = np.array([e[i] for i in ind])
				mdl = np.array([model.X[i] for i in ind])
				error = np.sum((mdl - expt * alpha_array[alpha_ind]) ** 2 / (expt * alpha_array[alpha_ind]) ** 2)
				# print(error)
				# print(np.matmul(np.matmul(model.X , Q) , model.X) - np.matmul(c , model.X) + amount_to_add_back)
				# print(model.objVal)
				# input("done")
				print(min_objective, min_ind, error, model.objVal,
				      np.matmul(np.matmul(model.X, Q), model.X) - np.matmul(c, model.X) + amount_to_add_back)
				# time.sleep(1)
				array_of_obj_val_no_val[alpha_ind] = error
				array_of_obj_val[alpha_ind] = model.objVal
				# currently sets min based on no valine
				if array_of_obj_val_no_val[alpha_ind] < min_objective or min_objective == -1:
					print(mdl)
					min_objective = array_of_obj_val_no_val[alpha_ind]
					min_ind = alpha_ind
			else:
				array_of_states[alpha_ind] = np.empty_like(model.X) * np.nan
				array_of_obj_val_no_val[alpha_ind] = np.nan
				array_of_obj_val[alpha_ind] = np.nan
		
		# print(np.matmul(np.matmul(model.X, Q), model.X) - np.matmul(c, model.X) + amount_to_add_back)
		if assign_opt_alpha_fluxes:
			print([array_of_states[min_ind, i] for i in model_index_of_experimental_fluxes])
			for i in model_index_of_experimental_fluxes:
				self.update_reaction_bounds(ordered_rxn[i], array_of_states[min_ind, i] - tol,
				                            array_of_states[min_ind, i] + tol)
		if search_save_path != "":
			result_dict = {}
			result_dict["alpha_list"] = alpha_array
			result_dict["saved_rxn_name_order"] = ordered_rxn
			result_dict["saved_met_order"] = ordered_met
			result_dict["best_alpha_index"] = min_ind
			result_dict["model_alignment"] = array_of_states
			result_dict["model_alignment_performances"] = array_of_obj_val
			result_dict["model_alignment_performances_no_valine"] = array_of_obj_val_no_val
			result_dict["experimental_flux_index"] = model_index_of_experimental_fluxes
			result_dict["experimental_fluxes"] = e
			with open(search_save_path.parent / (search_save_path.stem + ".pkl"), "wb") as outfile:
				pickle.dump(result_dict, outfile)
	
	def fit_to_experimental_data_l1_norm(self, experimental_data_file, alpha_array, assign_opt_alpha_fluxes=True,
	                                     tol=1e-7,
	                                     search_save_path=""):
		# experiment vector = ei
		# model vector = vi
		# |vi-aei|
		# (vi-aei)^2/(aei)^2 = vi^2/(aei)^2-2vi/aei+1
		# returns:
		# model with best experimental fit for given alpha
		# the relative error (objective term)
		# the index values of the experimental data values
		
		experimental_fluxes = np.loadtxt(experimental_data_file, delimiter=",", dtype=str)
		model_index_of_experimental_fluxes = []
		e = np.zeros(len(self.lb))
		amount_to_add_back = 0
		for i in range(np.shape(experimental_fluxes)[0]):
			if i > 0:
				rxn_name = experimental_fluxes[i, 0]
				flux_values = experimental_fluxes[i, 1].astype('float64')
				index_of_reaction = self.reaction_names.index(rxn_name)
				model_index_of_experimental_fluxes.append(index_of_reaction)
				e[index_of_reaction] = flux_values
				amount_to_add_back += 1
		# For very low values of experimental flux gurobi will run into accuracy errors
		# The current solution is to just take the upper limit of the experimental value of the one
		# poorly behaving flux, but for other cases it might be worth establishing a minimal cutoff point
		# And simply do linear minimaztion on the small values error
		array_of_states = np.zeros((np.size(alpha_array), np.size(self.lb)))
		array_of_obj_val = np.zeros((np.size(alpha_array)))
		
		min_ind = 0
		min_objective = -1
		for alpha_ind in range(np.size(alpha_array)):
			if self.env != None:
				model = gp.Model("fva_model_" + self.model_name, env=self.env)
			else:
				model = gp.Model("fva_model_" + self.model_name)
			c = np.zeros(len(e))
			for i in range(len(e)):
				if e[i] != 0:
					c[i] = (alpha_array[alpha_ind] * e[i])
				else:
					c[i] = 0
			model.optimize()
			model.Params.NumericFocus = 3
			model.Params.BarConvTol = 0
			model.Params.BarHomogeneous = 1
			print(gp.sys.version)
			react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=self.lb,
			                           ub=self.ub)
			
			model.setObjective(b, GRB.MINIMIZE)
			rhs = np.transpose(self.b)
			model.addConstr(self.S @ react_flux == rhs, name="c")
			model.optimize()
			input()
			if model.Status == 2 or model.Status == 13:
				array_of_states[alpha_ind] = model.X
				array_of_obj_val[alpha_ind] = model.objVal
				if model.objVal < min_objective or min_objective == -1:
					min_objective = model.objVal
					min_ind = alpha_ind
			else:
				array_of_states[alpha_ind] = np.empty_like(model.X) * np.nan
				array_of_obj_val[alpha_ind] = np.nan
		
		# print(np.matmul(np.matmul(model.X, Q), model.X) - np.matmul(c, model.X) + amount_to_add_back)
		if assign_opt_alpha_fluxes:
			for i in model_index_of_experimental_fluxes:
				self.update_reaction_bounds(i, array_of_states[min_ind, i] - tol, array_of_states[min_ind, i] + tol)
		if search_save_path != "":
			result_dict = {}
			result_dict["alpha_list"] = alpha_array
			result_dict["best_alpha_index"] = min_ind
			result_dict["model_alignment"] = array_of_states
			result_dict["model_alignment_performances"] = array_of_obj_val
			result_dict["experimental_flux_index"] = model_index_of_experimental_fluxes
			result_dict["experimental_fluxes"] = e
			with open(search_save_path.parent / (search_save_path.stem + ".pkl"), "wb") as outfile:
				pickle.dump(result_dict, outfile)
	
	def fit_to_experimental_data_abs_error(self, experimental_data_file, alpha_array, assign_opt_alpha_fluxes=True,
	                                       tol=1e-7,
	                                       search_save_path=""):
		# experiment vector = ei
		# model vector = vi
		# (vi-aei)^2 = vi^2-2vi*ei*a+a^2ei^2
		# returns:
		# model with best experimental fit for given alpha
		# the relative error (objective term)
		# the index values of the experimental data values
		
		experimental_fluxes = np.loadtxt(experimental_data_file, delimiter=",", dtype=str)
		model_index_of_experimental_fluxes = []
		e = np.zeros(len(self.lb))
		for i in range(np.shape(experimental_fluxes)[0]):
			if i > 0:
				rxn_name = experimental_fluxes[i, 0]
				flux_values = experimental_fluxes[i, 1].astype('float64')
				index_of_reaction = self.reaction_names.index(rxn_name)
				model_index_of_experimental_fluxes.append(index_of_reaction)
				e[index_of_reaction] = flux_values
		# For very low values of experimental flux gurobi will run into accuracy errors
		# The current solution is to just take the upper limit of the experimental value of the one
		# poorly behaving flux, but for other cases it might be worth establishing a minimal cutoff point
		# And simply do linear minimaztion on the small values error
		array_of_states = np.zeros((np.size(alpha_array), np.size(self.lb)))
		array_of_obj_val = np.zeros((np.size(alpha_array)))
		
		min_ind = 0
		
		min_objective = -1
		for alpha_ind in range(np.size(alpha_array)):
			Q = np.zeros((len(e), len(e)))
			for i in range(len(e)):
				if e[i] != 0:
					Q[i, i] = 1
				else:
					Q[i, i] = 0
			c = np.zeros(len(e))
			for i in range(len(e)):
				if e[i] != 0:
					c[i] = 2 * e[i] * alpha_array[alpha_ind]
				else:
					c[i] = 0
			if self.env != None:
				model = gp.Model("fva_model_" + self.model_name, env=self.env)
			else:
				model = gp.Model("fva_model_" + self.model_name)
			model.Params.NumericFocus = 3
			# model.Params.ScaleFlag = 0
			model.Params.BarConvTol = 0
			model.Params.BarHomogeneous = 1
			# model.Params.LogToConsole = 0
			react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=self.lb,
			                           ub=self.ub)
			const_add = np.sum((c / 2) ** 2)
			model.setObjective(react_flux @ Q @ react_flux - c @ react_flux + const_add, GRB.MINIMIZE)
			rhs = np.transpose(self.b)
			model.addConstr(self.S @ react_flux == rhs, name="c")
			model.optimize()
			if model.Status == 2 or model.Status == 13:
				array_of_states[alpha_ind] = model.X
				array_of_obj_val[alpha_ind] = model.objVal
				if model.objVal < min_objective or min_objective == -1:
					min_objective = model.objVal
					min_ind = alpha_ind
			else:
				array_of_states[alpha_ind] = np.empty_like(model.X) * np.nan
				array_of_obj_val[alpha_ind] = np.nan
		# print(np.matmul(np.matmul(model.X, Q), model.X) - np.matmul(c, model.X) + amount_to_add_back)
		if assign_opt_alpha_fluxes:
			for i in model_index_of_experimental_fluxes:
				self.update_reaction_bounds(i, array_of_states[min_ind, i] - tol, array_of_states[min_ind, i] + tol)
		if search_save_path != "":
			result_dict = {}
			result_dict["alpha_list"] = alpha_array
			result_dict["best_alpha_index"] = min_ind
			result_dict["model_alignment"] = array_of_states
			result_dict["model_alignment_performances"] = array_of_obj_val
			result_dict["experimental_flux_index"] = model_index_of_experimental_fluxes
			result_dict["experimental_fluxes"] = e
			with open(search_save_path.parent / (search_save_path.stem + ".pkl"), "wb") as outfile:
				pickle.dump(result_dict, outfile)
	
	def pinch_restricted_exchange_reactions(self, allowed_exchange_reaction_file, restore_essential=True,
	                                        true_exchange_signifier="", restore_previously_pinched=False,
	                                        restore_value=1000):
		# Leave option to keep excrete free
		met_rxn_dict = self.get_exchanged_reaction_info(true_exchange_signifier=true_exchange_signifier)
		print(met_rxn_dict)
		reactions_to_allow = list(np.loadtxt(allowed_exchange_reaction_file, dtype="str", delimiter=",")[1:, 0])
		reactions_to_allow = [i.replace(" ", "") for i in reactions_to_allow]
		exchange_type = list(np.loadtxt(allowed_exchange_reaction_file, dtype="str", delimiter=",")[1:, 1])
		exchange_type = [i.replace(" ", "") for i in exchange_type]
		print(len(reactions_to_allow))
		reactions_allowed = []
		count = 0
		rxn_list = list(met_rxn_dict.keys())
		for rxn in rxn_list:
			# print(rxn, self.test_feasibility())
			if rxn not in reactions_to_allow:
				print(rxn)
				# right now this still allows excretion
				if restore_essential:
					ub_save = self.model_dict["rxn_dict"][rxn]['ub']
					lb_save = self.model_dict["rxn_dict"][rxn]['lb']
				self.update_reaction_bounds(rxn, 0, 0)
				if restore_essential:
					if not self.test_feasibility():
						print(f"removing reaction {rxn} breaks feasibility")
						self.update_reaction_bounds(rxn, lb_save, ub_save)
				else:
					del self.model_dict["rxn_dict"][rxn]
			elif rxn in reactions_to_allow:
				reactions_allowed.append(rxn)
				count += 1
				print(rxn, count)
				rxn_index = reactions_to_allow.index(rxn)
				rxn_exchange_type = exchange_type[rxn_index]
				if rxn_exchange_type == "intake":
					if restore_previously_pinched:
						self.update_reaction_bounds(rxn, -1 * restore_value, 0)
					else:
						self.update_reaction_bounds(rxn, "keep", 0)
				elif rxn_exchange_type == "excrete":
					if restore_previously_pinched:
						self.update_reaction_bounds(rxn, 0, restore_value)
					else:
						self.update_reaction_bounds(rxn, 0, "keep")
				elif rxn_exchange_type == "both" and restore_previously_pinched:
					self.update_reaction_bounds(rxn, -1 * restore_value, restore_value)
				elif rxn_exchange_type != "both":
					print("Exchange Reaction direction not understood")
	
	def create_RAS_values(self, gene_array, feature_array, return_dict=True):
		def create_entrez_to_express_dict(gene_names, gene_array, entrez_dict):
			mouse_gene_list = list(gene_names[:, 0])
			
			list_unmatch = []
			list_m_ens_unmatch = []
			
			for i in entrez_dict.keys():
				if "mouse_ortholog_id" in entrez_dict[i].keys():
					if entrez_dict[i]["mouse_ortholog_id"] in mouse_gene_list:
						match_ind = mouse_gene_list.index(entrez_dict[i]["mouse_ortholog_id"])
						entrez_dict[i]["expression_val"] = gene_array[match_ind] / entrez_dict[i]["degen"]
					# print(i, entrez_dict[i], 1)
					# print(avg_exp[mouse_gene_list.index(entrez_dict[i]["mouse_ensemble"])], entrez_dict[i]["degen"],entrez_dict[i]["expression_val"])
					else:
						# print(i, entrez_dict[i], 2)
						entrez_dict[i]["expression_val"] = 0
						list_m_ens_unmatch.append([i, entrez_dict[i]["mouse_ortholog_id"]])
				else:
					entrez_dict[i]["expression_val"] = 0
					# print(i, entrez_dict[i], 3)
					list_unmatch.append(i)
			# print(list_unmatch)
			# print(list_m_ens_unmatch)
			return entrez_dict
		
		def name_to_express(entrez_dict, matchobj):
			# converts from gene name to RNA expression
			entrez_name = matchobj.group(0)
			return str(entrez_dict[entrez_name]["expression_val"])
		
		def strip_pare(matchobj):
			# strips parenthesis which encompass single value
			return matchobj.group(0)[1:-1]
		
		def convert_and(matchobj):
			# converts all pure and sections to value
			ands = matchobj.group(0)
			ands = ands.replace(" ", "")
			ands = np.average(np.array([float(i) for i in ands.split("and")]))
			return str(ands)
		
		def convert_or_internal(matchobj):
			# converts all pure or statements enclosed in parenthesis to value
			ors = matchobj.group(0)[1:-1]
			ors = ors.replace(" ", "")
			if ors.count("and") > 0:
				print(ors)
				print("ERROR")
				input()
			# print(ors.split("or"))
			ors = np.sum(np.array([float(i) for i in ors.split("or")]))
			return str(ors)
		
		def convert_or_external(matchobj):
			# converts final or statements to value, must be ran at very end
			ors = matchobj.group(0)
			ors = ors.replace(" ", "")
			if ors.count("and") > 0:
				print(ors)
				print("ERROR")
				input()
			# print(ors.split("or"))
			ors = np.sum(np.array([float(i) for i in ors.split("or")]))
			return str(ors)
		
		def convert_gr_Rule_to_RAS(gr_rule, entrez_dict):
			# converts name to expression
			name_to_real_express = partial(name_to_express, entrez_dict)
			gr_rule = re.sub("\d+\.\d+", name_to_real_express, gr_rule)
			int_or_altered = True
			while int_or_altered:
				gr_rule_cp_2 = copy.deepcopy(gr_rule)
				and_altered = True
				while and_altered:
					gr_rule_cp = copy.deepcopy(gr_rule)
					gr_rule = re.sub("\(\d+\.?\d*e?-?\d*\)", strip_pare, gr_rule)
					gr_rule = re.sub("(\d+\.?\d*e?-?\d* ?(and) ?)+\d+\.?\d*e?-?\d*", convert_and, gr_rule)
					if gr_rule_cp == gr_rule:
						and_altered = False
				gr_rule = re.sub("\(\d+\.?\d*e?-?\d*\)", strip_pare, gr_rule)
				gr_rule = re.sub("\((\d+\.?\d*e?-?\d* ?(or) ?)+\d+\.?\d*e?-?\d*\)", convert_or_internal, gr_rule)
				if gr_rule_cp_2 == gr_rule:
					int_or_altered = False
			gr_rule = re.sub("(\d+\.?\d*e?-?\d* ?(or) ?)+\d+\.?\d*e?-?\d*", convert_or_external, gr_rule)
			
			if gr_rule == "":
				# print(np.nan)
				return np.nan
			else:
				# print(float(gr_rule))
				return float(gr_rule)
		
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		gr_rule = self.get_grRule_list()
		entrez_express_dict = create_entrez_to_express_dict(feature_array, gene_array, self.model_dict["gene_dict"])
		RAS_list = []
		for i in range(len(gr_rule)):
			RAS = convert_gr_Rule_to_RAS(gr_rule[i], entrez_express_dict)
			RAS_list.append(RAS)
		if return_dict:
			RAS_dict = dict(zip(rxn_list, RAS_list))
			return RAS_dict
		else:
			return RAS_list
	
	def formula_to_dict(self, formula):
		chemical_dict = {}
		initial_list = [i for i in range(len(formula)) if formula[i].isupper()]
		len_list = [formula[initial_list[i]:initial_list[i + 1]] for i in range(len(initial_list) - 1)]
		len_list.append(formula[initial_list[-1]:])
		for i in len_list:
			first_number_location = 0
			for alnum in range(len(i)):
				if i[alnum].isalpha():
					first_number_location += 1
			element = i[:first_number_location]
			if i[first_number_location:] != "":
				number = int(i[first_number_location:])
			else:
				number = 1
			if element in chemical_dict.keys():
				chemical_dict[element] += number
			else:
				chemical_dict[element] = number
		return chemical_dict
	
	def print_nonpinched_exchange(self):
		exch_dict = self.get_exchanged_reaction_info()
		
		for e_rxn in exch_dict.keys():
			if not (self.rxn_dict[e_rxn]["lb"] == 0 and self.rxn_dict[e_rxn]["ub"] == 0):
				if self.rxn_dict[e_rxn]["lb"] == 0:
					print(e_rxn, "excrete")
				elif self.rxn_dict[e_rxn]["ub"] == 0:
					print(e_rxn, "intake")
				else:
					print(e_rxn, "both")
	
	def test_essential_list(self, reaction_index_list, use_gp_env=True):
		# test against set too
		lb, ub, S, b, ordered_rxn, ordered_met = self.dicts_to_mats()
		if isinstance(reaction_index_list[0], str):
			reaction_index_list = [ordered_rxn.index(name) for name in reaction_index_list]
		if len(set(reaction_index_list)) != len(reaction_index_list):
			return False, "Duplicate Essential Flux"
		save_ub_full = copy.deepcopy(ub)
		save_lb_full = copy.deepcopy(lb)
		essential_flux_list = np.zeros(len(ordered_rxn))
		for rxn_ind in reaction_index_list:
			essential_flux_list[rxn_ind] = 1
		for essential_matrix_index in range(len(essential_flux_list)):
			# if the reaction is essential then lb and ub will be the same, if it is not then lb and ub will be 0
			self.update_reaction_bounds(ordered_rxn[essential_matrix_index],
			                            essential_flux_list[essential_matrix_index] * lb[essential_matrix_index],
			                            essential_flux_list[essential_matrix_index] * ub[essential_matrix_index])
			lb[essential_matrix_index] = essential_flux_list[essential_matrix_index] * lb[essential_matrix_index]
			ub[essential_matrix_index] = essential_flux_list[essential_matrix_index] * ub[essential_matrix_index]
		if not self.test_feasibility(lb, ub, S, b):
			return False, "not feasible enough"
		for essential_matrix_index in range(len(essential_flux_list)):
			if essential_matrix_index % 1000 == 0:
				print(essential_matrix_index / len(essential_flux_list))
			if essential_flux_list[essential_matrix_index] == 1:
				if not (lb[essential_matrix_index] > 0 or ub[essential_matrix_index] < 0):
					save_ub = copy.deepcopy(ub[essential_matrix_index])
					save_lb = copy.deepcopy(lb[essential_matrix_index])
					self.update_reaction_bounds(ordered_rxn[essential_matrix_index], 0, 0)
					lb[essential_matrix_index] = 0
					ub[essential_matrix_index] = 0
					if self.test_feasibility(lb, ub, S, b):
						return False, "Too feasible"
					else:
						self.update_reaction_bounds(ordered_rxn[essential_matrix_index], save_lb, save_ub)
						lb[essential_matrix_index] = save_ub
						ub[essential_matrix_index] = save_lb
		for i in range(len(ub)):
			self.update_reaction_bounds(ordered_rxn[i], save_lb_full[i], save_ub_full[i])
		return True, "Just right"
	
	def get_used_genes(self, nomenclature="HGNC"):
		used_genes = []
		for i in self.grRules:
			if nomenclature == "HGNC":
				if i != "[]":
					genes_found = re.findall("(HGNC:\d+)", i)
					for j in genes_found:
						if j not in used_genes:
							used_genes.append(j)
		return used_genes
	
	def generate_essential_flux_dataframe(self, n, output_dir, file_name_base, seed_list=[], print_progress=False,
	                                      seed_name=False, **kwargs):
		lb, ub, S, b, ordered_rxn, ordered_met = self.dicts_to_mats()
		seed_empty = False
		if seed_list == []:
			seed_list = [random.random() for i in range(n)]
			seed_empty = True
		if seed_name:
			file_name = file_name_base + str(seed_list[0]).replace(".", "")
		value_list = np.zeros((n, len(self.model_dict["rxn_dict"])))
		count = 0
		seed_fail = []
		failed = False
		for seed_val in seed_list:
			print(seed_val)
			if seed_empty:
				essential_flux = self.find_essential_fluxes_fast(print_progress=print_progress)
			else:
				essential_flux = self.find_essential_fluxes_fast(print_progress=print_progress, seed=seed_val)
			is_essential, fail_desc = self.test_essential_list(essential_flux)
			if is_essential:
				for i in essential_flux:
					value_list[count, i] += 1
			else:
				failed = True
				seed_fail.append([seed_val, fail_desc])
			count += 1
		if failed:
			print(f"Failed, description of failures:\n{seed_fail}")
		else:
			np.save(output_dir / file_name, pd.DataFrame(value_list).to_numpy())
			np.save(output_dir / (file_name_base + "header"), ordered_rxn)
	
	def get_met_comp_dict(self, met_name):
		chem_count_dict = {"C": 0, "H": 0, "N": 0, "O": 0, "P": 0, "S": 0, "Ca": 0, "Cl": 0, "K": 0, "Na": 0, "Se": 0,
		                   "Co": 0, "I": 0, "Fe": 0, "R": 0, "Ra": 0, "Rb": 0, "Rc": 0, "X": 0, "Y": 0}
		print(self.model_dict["met_dict"][met_name].keys())
		comp = self.model_dict["met_dict"][met_name]['metFormulas']
		flag = False
		if "R" in comp:
			flag = True
		comp = comp.replace("FULLR3", "Rc")
		comp = comp.replace("FULLR2", "Rb")
		comp = comp.replace("FULLR", "Ra")
		search = re.findall("[A-Z][a-z]?\d*", comp)
		for i in search:
			chem = re.findall("[A-Z][a-z]?", i)[0]
			count = 1
			if len(re.findall("\d+", i)) == 1:
				count = re.findall("\d+", i)[0]
			try:
				chem_count_dict[chem] += int(count)
			except KeyError:
				print("Key not found")
				print(chem)
				return {}
		return chem_count_dict
	
	def get_met_mol_weight(self, met_name):
		mol_weight_dict = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "P": 30.97376, "S": 32.06, "Ca": 40.078,
		                   "Cl": 35.45, "K": 39.0983, "Na": 22.98977, "Se": 78.971, "Co": 58.93319, "I": 126.9045,
		                   "Fe": 55.845, "R": 12.011 * 15 + 1.008 * 31, "Ra": (12.011 * 15 + 1.008 * 31) * 1,
		                   "Rb": (12.011 * 15 + 1.008 * 31) * 1, "Rc": (12.011 * 15 + 1.008 * 31) * 1, "X": 0, "Y": 0}
		chem_count_dict = self.get_met_comp_dict(met_name)
		mw = 0
		for chem in chem_count_dict.keys():
			mw += int(chem_count_dict[chem]) * mol_weight_dict[chem]
		return mw
	
	def get_rxn_mass_change(self, rxn_name):
		st = self.model_dict["rxn_dict"][rxn_name]["rxn_metabolites"]
		mass_change = 0
		for met in st.keys():
			mass_change += st[met] * self.get_met_mol_weight(met)
		return mass_change
	
	def get_rxn_element_change(self, rxn_name):
		chem_count_dict = {"C": 0, "H": 0, "N": 0, "O": 0, "P": 0, "S": 0, "Ca": 0, "Cl": 0, "K": 0, "Na": 0, "Se": 0,
		                   "Co": 0, "I": 0, "Fe": 0, "R": 0, "Ra": 0, "Rb": 0, "Rc": 0, "X": 0, "Y": 0}
		st = self.model_dict["rxn_dict"][rxn_name]["rxn_metabolites"]
		for met in st.keys():
			chem_count = self.get_met_comp_dict(met)
			for chem in chem_count.keys():
				chem_count_dict[chem] += st[met] * chem_count[chem]
		return chem_count_dict
	
	def build_essential_flux_matrix_from_dataframes(self, list_of_dataframes):
		
		for dataframe_ind in range(len(list_of_dataframes)):
			if dataframe_ind == 0:
				array = list_of_dataframes[dataframe_ind].to_numpy()
			else:
				array = np.vstack((array, list_of_dataframes[dataframe_ind].to_numpy()))
		return array
	
	def generate_essential_flux_model_from_essential_flux_matrix(self, essential_flux_matrix, use_gp_env=True):
		essential_flux_scores = np.average(essential_flux_matrix, axis=0)
		least_essential_flux_indices = np.argsort(essential_flux_scores)
		
		most_essential_flux_indices = np.flip(least_essential_flux_indices)
		
		saved_ub = copy.deepcopy(self.ub)
		saved_lb = copy.deepcopy(self.lb)
		
		for i in most_essential_flux_indices:
			if (self.ub[i] < 0 or self.lb[i] > 0):
				continue
			else:
				self.update_reaction_bounds(i, 0, 0)
		essential_flux_main = []
		for i in most_essential_flux_indices:
			if (self.ub[i] < 0 or self.lb[i] > 0):
				essential_flux_main.append(i)
				continue
			else:
				self.update_reaction_bounds(i, saved_lb[i], saved_ub[i])
				if self.test_feasibility():
					essential_flux_main.append(i)
					break
				else:
					essential_flux_main.append(i)
		
		most_essential_flux_names = [self.reaction_names[i] for i in most_essential_flux_indices]
		essential_flux_main_names = [self.reaction_names[i] for i in essential_flux_main]
		not_essential_fluxes = list(set(most_essential_flux_names) - set(essential_flux_main_names))
		
		self.del_reaction(not_essential_fluxes)
		self.purge_metabolites()
		flux_limits = self.fva()
		for i in range(len(flux_limits)):
			self.lb[i] = flux_limits[i][0]
			self.ub[i] = flux_limits[i][1]


