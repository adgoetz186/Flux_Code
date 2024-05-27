import pickle

from functools import partial
import sympy as sy
from pympler import asizeof
import gurobipy as gp
import copy
import gzip
import re
import copy as cp
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
	model = "none"

	def __init__(self, model_dict=None):
		"""Description:
			Initializes Flux_Model_Class object

		Input:
			model_dict (dict, optional): contains the keys "model_name, rxn_dict, met_dict, gene_dict, gurobi_token".
			See "Notes" to see the contents of associated values

		Output:
			Flux_Model_Class object

		Notes:
			Often initialization will be empty, model will be loaded in after initialization using
		one of the load functions

			The model dict should contain the following values
			"model_name, rxn_dict, met_dict, gene_dict, gurobi_token"
			with keys given by:
			model_name: a name for the model, used to generate filenames when saving some model related datasets
			rxn_dict: a dict with each reaction name as a key and a dict containing all reaction information as a value
			often reaction info contains: the metabolites involved in the reaction ("rxn_metabolites"),
			the lower bound flux of the rxn ("lb"), the upper bound flux of the rxn ("ub"),
			the subsystem the reaction is associated with ("subsystem"),
			 the gene rule for the reaction ("grRules")
			met_dict: a dict with each metabolite name as a key and a dict containing all metabolite information as a value
			often metabolite info contains: The value of Sv associated with the metabolite ("b"),
			The chemical formula for the metabolite ("met_composition"),
			the compartment of the metabolite ("met_compartment")
			gene_dict: a dict with each gene name as a key and a dict containing all gene information as a value
			often gene info contains: The entrez id ("gene_name") and the gene symbol ("symbol")"""
		self.model_dict = model_dict

	def add_gp_key_env_to_model(self, TokenServer_name):
		"""Description:
			Adds gurobi server token name to the instance for HPC based linear solving.

		Input:
			TokenServer_name (str): contains the token server name, for use with gurobipy

		Output:
			None

		Notes:
			for gurobi documentation see https://www.gurobi.com/documentation/9.0/refman/tokenserver.html
		If this is being used with University of Florida HiPerGator, use 'grb-ts.ufhpc'"""
		env = gp.Env(empty=True)
		env.setParam('TokenServer', TokenServer_name)
		env.start()
		self.model_dict['gurobi_token'] = env

	def save_model_as_mat(self, filename):
		"""Description:
			Saves model as a file which can be read by MATLAB. Contains the following arrays:

			lb: set of lower bounds

			ub: set of upper bounds

			S: stoichiometric matrix

			b: gives right hand side of Sv = b constraint, often a 0 vector except when fluxes are pinched

			rxn_list: list of named reactions

			met_list: list of named metabolites

			gr_list: list of gene rules

		Input:
			filename (str or Path): Path to save file to

		Output:
			None

		Notes:
			rxn_list and met_list are alphabetized. All arrays are aligned with either of these lists when applicable.
		For example if the ith element of rxn_list is A, then the ith element of lb, ub, gr_list,
		and the ith column of S, will correspond to reaction A."""
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		gr_list = self.get_grRule_list()
		mat_dict = {}
		mat_dict["lb"] = lb
		mat_dict["ub"] = ub
		mat_dict["S"] = S
		mat_dict["b"] = b
		mat_dict["rxn_list"] = rxn_list
		mat_dict["met_list"] = met_list
		mat_dict["gr_list"] = gr_list
		spio.savemat(filename, mat_dict)

	def load_recon1_mat_model(self, filename, new_model_name="Model"):
		"""Description:
			Loads the recon1 based zeromodel from []

		Input:
			filename (str or Path): Path to zero model

			new_model_name (str): name of model, e.g. recon1_zero.

		Output:
			None

		Notes:
			None"""
		mat = loadmat(filename)

		S = sp.csr_matrix(mat["model"]["S"]).astype(float)
		S = S.toarray()
		lb = mat["model"]["lb"].astype(float)
		ub = mat["model"]["ub"].astype(float)
		b = mat["model"]["b"].astype(float)
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

	def load_bigg_json_model(self, filename, new_model_name="model", model_header=None, rxn_comp=None, met_comp=None,
							 gene_comp=None, keep_extra=True):
		"""Description:
				Loads model from json file format from BiGG (http://bigg.ucsd.edu/)

			Input:
				filename (str or Path): Path to .json file

				new_model_name (str): name of model, e.g. recon1

				model_header (dict): dict with keys referring to the names of our model components
				(i.e. "reactions", "metabolites", and "genes")
				and values given by the corresponding names of BiGG model components

				rxn_comp (dict): dict with keys referring to the names of our reaction components
				(i.e. "rxn_name", "rxn_metabolites", "lb", "ub", "pinched_reactions", "subsystem", and "grRules")
				and values given by the corresponding names of BiGG model components

				met_comp (dict): dict with keys referring to the names of our metabolite components
				(i.e. "met_name", "b", "met_composition", and "met_compartment")
				and values given by the corresponding names of BiGG model components

				gene_comp (dict): dict with keys referring to the names of our gene components
				(i.e. "gene_name" and "symbol")
				and values given by the corresponding names of BiGG model components

				keep_extra (boolean): If True, saves all data associated with the json file not explicitly stored

			Output:
				None

			Notes:
				None"""
		self.model_dict = {}

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
		"""Description:
			Loads model from json file format based around fast and efficient access to data after loading. See Notes
			for description of this format

		Input:
			filename (str or Path): Path to .json file

		Output:
			None

		Notes:
			json file contains the following 5 keys: "model_name, rxn_dict, met_dict, gene_dict, gurobi_token"
			model_name: a name for the model, used to generate filenames when saving some model related datasets
			rxn_dict: a dict with each reaction name as a key and a dict containing all reaction information as a value
			often reaction info contains: the metabolites involved in the reaction ("rxn_metabolites"),
			the lower bound flux of the rxn ("lb"), the upper bound flux of the rxn ("ub"),
			the subsystem the reaction is associated with ("subsystem"),
			 the gene rule for the reaction ("grRules")
			met_dict: a dict with each metabolite name as a key and a dict containing all metabolite information as a value
			often metabolite info contains: The value of Sv associated with the metabolite ("b"),
			The chemical formula for the metabolite ("met_composition"),
			the compartment of the metabolite ("met_compartment")
			gene_dict: a dict with each gene name as a key and a dict containing all gene information as a value
			often gene info contains: The entrez id ("gene_name") and the gene symbol ("symbol")
			"""
		with open(filepath, "r") as infile:
			self.model_dict = json.load(infile)

	def save_model_as_fast_key_json(self, filepath):
		"""Description:
			Saves model to json file format based around fast and efficient access to data after loading. See Notes
			for description of this format

		Input:
			filename (str or Path): Path to .json file

		Output:
			None

		Notes:
			json file contains the following 5 keys: "model_name, rxn_dict, met_dict, gene_dict, gurobi_token"
			model_name: a name for the model, used to generate filenames when saving some model related datasets
			rxn_dict: a dict with each reaction name as a key and a dict containing all reaction information as a value
			often reaction info contains: the metabolites involved in the reaction ("rxn_metabolites"),
			the lower bound flux of the rxn ("lb"), the upper bound flux of the rxn ("ub"),
			the subsystem the reaction is associated with ("subsystem"),
			 the gene rule for the reaction ("grRules")
			met_dict: a dict with each metabolite name as a key and a dict containing all metabolite information as a value
			often metabolite info contains: The value of Sv associated with the metabolite ("b"),
			The chemical formula for the metabolite ("met_composition"),
			the compartment of the metabolite ("met_compartment")
			gene_dict: a dict with each gene name as a key and a dict containing all gene information as a value
			often gene info contains: The entrez id ("gene_name") and the gene symbol ("symbol")
			"""
		with open(filepath / f"{self.model_dict['model_name']}.json", "w") as outfile:
			json.dump(self.model_dict, outfile)

	def update_reaction_bounds(self, name, new_lb, new_ub):
		"""Description:
			updates reaction bounds

		Input:
			name (str): Name of reaction to update

			new_lb (float or str): if float, updates lb to specified value. If str, should be "keep", doesnt change lb

			new_ub (float or str): if float, updates ub to specified value. If str, should be "keep", doesnt change ub

		Output:
			None

		Notes:
			None"""
		if not new_lb == "keep":
			self.model_dict["rxn_dict"][name]["lb"] = new_lb
		if not new_ub == "keep":
			self.model_dict["rxn_dict"][name]["ub"] = new_ub
		if self.model_dict["rxn_dict"][name]["ub"] < self.model_dict["rxn_dict"][name]["lb"]:
			print(new_ub, new_lb)
			print(f"WARNING! Upper bound smaller than lower bound at {name}")
			raise ValueError

	def dicts_to_mats(self, alphabetize=True):
		"""Description:
			Central to many functions. Generates numpy arrays for S, lb, ub, and b. Also generates lists for
			named reactions and named metabolites.

		Input:
			alphabetize (boolean, optional):  Often assumed to be true, orders named reactions and metabolites

		Output:
			lb (ndarray): array of lower bound values

			ub (ndarray): array of upper bound values

			S (ndarray): array containing stoichiometric matrix

			b (ndarray): array containing the RHS of Sv = b (often 0 vector)

			rxn_list (list): list of reaction names

			met_list (list): list of metabolite names

		Notes:
			All output arrays are aligned with either rxn_list or met_list.
		For example if the ith element of rxn_list is A, then the ith element of lb, ub, gr_list,
		and the ith column of S, will correspond to reaction A. This removes issues which might arise from ordering
		ensuring 2 models with identical components will yield identical S, ub, lb, b, rxn_list, and met_list."""
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

	def normalize(self, vec):
		"""Description:
			Normalizes vector, vec

		Input:
			vec (ndarray):  ray vector

		Output:
			normalized vector (ndarray)

		Notes:
			None"""
		norm = np.linalg.norm(vec)
		if norm == 0:
			return vec
		else:
			return vec / norm

	def get_grRule_list(self, alphabetize=True):
		"""Description:
			obtains a list of gr rules which are ordered in the same rxn order generated by the
			dicts_to_mats function. "alphabetize" is passed to that function and is recommended to be set to True

		Input:
			alphabetize (boolean, optional):  Often assumed to be true, orders named reactions and metabolites

		Output:
			gr_Rule_list (list): list of gene rules

		Notes:
			None"""
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats(alphabetize=alphabetize)
		gr_Rule_list = []
		for i in rxn_list:
			if len(self.model_dict["rxn_dict"][i]["grRules"]) == 0:
				gr_Rule_list.append("")
			else:
				gr_Rule_list.append(self.model_dict["rxn_dict"][i]["grRules"])
		return gr_Rule_list

	def update_gene_dict(self):
		"""Description:
			updates model gene dict, to be run after adding a new reaction with gene rules referencing new genes

		Input:
			None

		Output:
			None

		Notes:
			None"""
		gene_list = self.get_used_gene_list()
		print(gene_list)
		for gene in gene_list:
			if gene not in self.model_dict["gene_dict"].keys():
				print(gene)
				degen = 0
				for gene_2 in gene_list:
					if gene.split(".")[0] == gene_2.split(".")[0]:
						degen += 1
				self.model_dict["gene_dict"][gene] = {"entrez": gene.split(".")[0], "degen": degen}

	def get_used_gene_list(self):
		"""Description:
			Gets list of genes which are used in model, based on all genes which are referenced by a reaction's gene rules

		Input:
			None

		Output:
			gene_list (list): list of gene names in entrez format.

		Notes:
			None"""
		gene_list = []
		gr_Rule_list = self.get_grRule_list()
		for i in gr_Rule_list:
			for id in re.findall("\d+\.\d+", i):
				if id not in gene_list:
					gene_list.append(id)
		return gene_list

	def find_objective_value(self, c, alphabetize=True):
		"""Description:
			Finds objective value for a given linear objective function

		Input:
			c (ndarray):  vector of objective function coefficients

			alphabetize (boolean, optional):  Often assumed to be true, orders named reactions and metabolites

		Output:
			model.objVal (float): value at optimum

			react_flux.X (ndarray): vector of fluxes which optimizes the objective

		Notes:
			None"""
		if self.model_dict['gurobi_token'] is None:
			env = None
		else:
			env = self.model_dict['gurobi_token']
		# env.start()
		if env != None:
			model = gp.Model("opt_model_" + self.model_dict['model_name'], env=env)
		else:
			model = gp.Model("opt_model_" + self.model_dict['model_name'])
		model.Params.LogToConsole = 0
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats(alphabetize=alphabetize)

		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb,
								   ub=ub)
		obj = np.reshape(c, (1, -1))
		# print(obj.shape)
		# print(obj @ react_flux)
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		rhs = np.transpose(b)
		model.addConstr(S @ react_flux == rhs, name="c")
		model.optimize()
		if model.Status == 2:
			return model.objVal, react_flux.X
		else:
			print(model.Status)
			return np.nan

	def fva(self, ignore_pinched=False, lb=None, ub=None, S=None, b=None, rxn_list=None, return_points=False,
			print_progress=True, ub_tol=0):
		"""Description:
			Performs FVA. Obtains min and max possible fluxes for each reaction.
			Min and max values are obtained for each reaction by minimizing and maximizing
			an objective function with a single nonzero objective coefficient belonging to the reaction of interest

		Input:
			ignore_pinched (boolean, optional): If true any fluxes which are constrained to be 0 are not run through
			optimization, saves time when many reactions are pinched to 0.

			lb, ub, S, b, rxn_list (optional):  Output of dicts_to_mats function, very rarely needed, but
			for certain cases it might be advantageous to run outside of function, modify then pass results to fva.

			return_points (boolean, optional):  If false returns both objective value for each reaction.
			If True also returns an array containing the points at which all reaction objective were obtained.

			print_progress (boolean, optional): If true prints number of reactions remaining

			ub_tol (float, optional): Used in certain cases where fva is used to reset bounds. Slight tolerances
			<1e-7 can help avoid numerical issues.

		Output:
			dict(zip(rxn_list, variability_list)) (dict): dictionary with inidividual reaction names as keys and
			the corresponding reaction limits as a value

			point_array (ndarray, optional): Only returned if return_points is True. Array containing the points
			at which all reaction optima were obtained.

		Notes:
			ignore_pinched is the only argument which should be changed without having a strong motivation.
			Other arguments are for niche applications"""
		if lb is None:
			lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		if self.model_dict["gurobi_token"] is not None:
			model = gp.Model("fva_model_" + self.model_dict["model_name"], env=self.model_dict["gurobi_token"])
		else:
			model = gp.Model("fva_model_" + self.model_dict["model_name"])
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		model.addConstr(S @ react_flux == b, name="c")
		model.Params.LogToConsole = 0
		model.Params.FeasibilityTol = 1e-9
		model.Params.NumericFocus = 3
		obj = np.zeros((1, S.shape[1]))
		variability_list = [{"lb": 0, "ub": 0} for i in range(len(lb))]
		if return_points:
			point_array = np.zeros((len(lb) * 2, len(lb)))
		for i in range(len(lb)):
			if ignore_pinched and lb[i] == ub[i] and lb[i] == 0:
				variability_list[i]["ub"] = 0
				variability_list[i]["lb"] = 0
			else:
				if print_progress:
					print(len(lb) - i)
				min_max = []
				obj = np.zeros(obj.shape)
				obj[0, i] = 1
				model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
				model.optimize()
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
				variability_list[i]["ub"] = max(min_max) + ub_tol
				variability_list[i]["lb"] = min(min_max)
		if return_points:
			return dict(zip(rxn_list, variability_list)), point_array
		else:
			return dict(zip(rxn_list, variability_list))

	def test_feasibility(self, lb=None, ub=None, S=None, b=None, **kwargs):
		"""Description:
			Tests a model to see if its feasible

		Input:
			lb, ub, S, b (optional):  Output of dicts_to_mats function, very rarely needed, but
			for certain cases it might be advantageous to run outside of function, modify then pass results to function.

		Output:
			(boolean): model.Status == 2 is true if model is feasible"""
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
		model.optimize()
		return model.Status == 2

	def get_exchanged_reaction_info(self, true_exchange_signifier=""):
		"""Description:
			Generates dict containing all exchange reactions and their corresponding metabolites

		Input:
			true_exchange_signifier (optional):  If a reaction name contains some type of substring which conveys the
			reaction is truly exchange ("EX" for example)

		Output:
			(dict): dict(zip(exch_rxn, exch_met)) is dict containing all exchange reactions and their corresponding metabolites"""
		exch_rxn = []
		exch_met = []
		for i in self.model_dict["rxn_dict"].keys():
			if len(self.model_dict["rxn_dict"][i]["rxn_metabolites"]) == 1 and true_exchange_signifier in i:
				exch_rxn.append(i)
				exch_met.append(list(self.model_dict["rxn_dict"][i]["rxn_metabolites"].keys())[0])
		return dict(zip(exch_rxn, exch_met))

	def metabolite_info(self, name, flux_dtf=None):
		"""Description:
			Prints info about the named metabolite

		Input:
			name (str): name of the metabolite of interest

			flux_dtf (optional):  flux dataframe can be passed, if it is reactions which produce or consume metabolite
			will be printed along with their fluxes"""
		print(f"Information on {name}:")
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

	def generate_mat_wo_exchange(self, prune_specific=[]):
		"""Description:
			Generates important model arrays with exchange reactions removed

		Input:
			prune_specific (optional):  additional reaction names to prune, often reactions which are like exchange
			reactions in that they break mass balance but are not technically exchange reactions,
			for example, the biomass reaction.

		Output:
			internal_lb (ndarray): array of lower bound values with exchange removed

			internal_ub (ndarray): array of upper bound values with exchange removed

			internal_S (ndarray): array containing stoichiometric matrix with exchange removed

			b (ndarray): array containing the RHS of Sv = b (often 0 vector)

			internal_rxn (list): list of reaction names with exchange removed

			met_list (list): list of metabolite names"""
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

	def generate_therm_model_new(self, trans_N, point, u_bound_tol=1e-4, slack=1e-9):
		"""Description:
			Generates a model assigning free energy values to all reactions u_i, such that v_i*u_i < 0
			where v_i is the ith reactions flux. If the model is feasible this is possible. Importantly returns
			the gurobi code with 2 meaning feasible. This is done for better handling possible model outcomes.

		Input:
			trans_N (ndarray):  Null space of the exchange pruned S matrix

			point (ndarray):  Point containing the flux values
			u_bound_tol (optinoal): Assuming v_i is (positive or negative), ensures the max
			value is (-u_bound_tol or u_bound_tol) rather than 0 to avoid numerical issues or trivial solutions


		Output:
			model.Status (int): gurobi model status code, 2 is feasibles"""
		if self.model_dict['gurobi_token'] != None:
			model = gp.Model("therm_model_" + self.model_dict['model_name'], env=self.model_dict['gurobi_token'])
		else:
			model = gp.Model("therm_model_" + self.model_dict['model_name'])
		point[np.abs(point) < slack] = 0
		lb_mu = np.where(point > 0, -np.inf, u_bound_tol)
		lb_mu = np.where(point == 0, -np.inf, lb_mu)
		ub_mu = np.where(point > 0, -u_bound_tol, np.inf)
		ub_mu = np.where(point == 0, np.inf, ub_mu)
		# print(lb_mu)
		# print(ub_mu)
		# lb_negative_mu = np.where(point >= 0, -1000, tol)
		# ub_positive_mu = np.where(point <= 0, 1000, -1*tol)
		ts1 = trans_N.shape[1]
		mu = model.addMVar(shape=ts1, vtype=GRB.CONTINUOUS, name="mu", lb=lb_mu, ub=ub_mu)
		model.Params.LogToConsole = 0
		model.addConstr(trans_N @ mu == 0, name="c")
		model.Params.FeasibilityTol = 1e-9
		obj = np.zeros((1, ts1))
		model.setObjective(obj @ mu, GRB.MAXIMIZE)
		model.optimize()
		return model.Status

	def element_edge_vert_generator(self, rxn, rxn_dtf=None, element="compound", break_element_list=[]):
		"""Description:
			Currently in development, for generating visualizations of fluxes"""
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

	def generate_rxn_met_graph(self, edge_list, weight_list, color_edge_of_vert=None, default_edge_color="grey",
							   vert_leak_tol=None, max_vec_width=5, weight_display="raw"):
		"""Description:
					Currently in development, for generating visualizations of fluxes"""
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
		"""Description:
			Prints information about a given reaction

		Input:
			name: name of reaction to print information on"""
		print(f"Information on {name}:")
		print(
			f"Lower Bound: {self.model_dict['rxn_dict'][name]['lb']}, Upper Bound {self.model_dict['rxn_dict'][name]['ub']}")
		print(f"The relevant chemical species participating in the reaction and their coefficients are:")
		print(self.model_dict['rxn_dict'][name]["rxn_metabolites"])
		print(f"The reaction gene rule is:\n{self.model_dict['rxn_dict'][name]['grRules']}")

	def lin_dep_np_ctp(self, old_on_basis, origin, accepted_points, new_point, old_rank,
					   unique_vector_threshold_exponent=9):
		"""Description:
			Takes in a set of basis vectors and a new proposed basis vector and either accepts or rejects
			the proposed vector based on whether the addition of the new vector increases the dimensionality of the
			basis set

		Input:
			old_on_basis: 2d numpy array containing points which define a basis wrt the origin point

			origin: vectors are defined as the difference between points of interest and the origin
				Value should be within feasible space, either the first generated vertex or the center of the point
				matrix are valid approaches, we will use the first generated vertex

			accepted_points: points which have been accepted

			new_point: proposed point which creates a vector wrt the origin

			old_rank: previous rank


		Output:
			either the old basis and old rank are returned or, if the proposed vector added to the basis creates a
			new larger basis, that basis and rank are returned"""

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

	def quick_pinch(self, max_sampling_level, tol=1e-6, rxn_penalty_vec=None):
		"""Description:
			Fully constricts dimensions which are tightly bound, converting S_full*v_full = 0 into S_var * v_var = - S_con * v_con
			Where S_full is the stochiometic matrix, v_full is the vector of fluxes,  S_var and v_var contain entries which
			correspond to unconstrained fluxes while S_con and v_con correspond to constrained fluxes and as a result
			S_con*v_con is constant

		Input:
			max_sampling_level (float): the largest multiplicative factor placed on the minimum reaction penalty

			tol (float): Any bounds closer together than this value will be pinched

			rxn_penalty_vec (numpy array,optional): gives weights for 1 norm constraint.
				Size shoudl correspond to rxn count, default will weight all fluxes similarly.

		Output:
			var_ind (ndarray, size will be the number of unpinched reactions): provides indexes for the variable reactions
				this assumes alphabetical reaction ordering (which is the defualt order given when dicts_to_mats is called)
				This value allows one to convert from full vectors to vectors which correspond to just variable fluxes
			pinch_b (ndarray, size will be the number of metabolites): - S_con * v_con which is used in the modified
				flux balance constraint, Sv = b
			lb (ndarray, size will be number of reactions + 1): gives lower bound values after pinching, IMPORTANTLY
				this includes the minimum flux value of the unpinched fluxes (not total flux at min point, but total of variable fluxes)
			vec_min_flux (ndarray, size will be number of reactions + 1): The exact point the bounds are pinched to
				contains the full vector with variable fluxes, constant fluxes, and the total flux of the point
			ub (ndarray, size will be number of reactions + 1): gives upper bound values after pinching, IMPORTANTLY
				this includes the maximum flux value of the unpinched fluxes (not total flux at min point, but total of variable fluxes) as determined by the minimum and the max_sampling_level"""
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()

		# uncomment this to see the effect of pinching on the total space bounds
		# lbc = cp.deepcopy(lb)
		# ubc = cp.deepcopy(ub)

		fva_dict = self.fva(lb=lb, ub=ub, S=S, b=b, rxn_list=rxn_list, print_progress=False)
		ind = 0
		for key in fva_dict:
			ind += 1
		for rxn_ind in range(len(rxn_list)):
			lb[rxn_ind] = fva_dict[rxn_list[rxn_ind]]['lb']
			ub[rxn_ind] = fva_dict[rxn_list[rxn_ind]]['ub']
		var_ind = np.arange(np.size(ub))
		var_ind = var_ind[ub - lb > tol]

		con_ind = np.arange(np.size(ub))
		con_ind = con_ind[ub - lb <= tol]

		min_flux, vec_min_flux = self.get_min_flux(lb, ub, b, S, return_vec=True)

		# input()
		# input()
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

		# input()
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
		# fva(self, ignore_pinched=False, lb=None, ub=None, S=None, b=None, rxn_list=None, return_points=False, ub_tol=0):
		fva_dict_small = self.fva(False, lb_var_tf, ub_var_tf, S_var_tf, pinch_b, min_rxn_names, print_progress=False)
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

	def get_min_flux(self, lb, ub, b, S, rxn_penalty_vec=None, return_vec=False, invert=False):
		"""Description:
			Obtains minimum allowed flux

		Input:
			lb, ub, S, b:  Output of dicts_to_mats function

			rxn_penalty_vec (ndarray, optional):  used to weight reaction objectives before minimizing

			return_vec (boolean, optional): If true returns all flux values rather than the minimized objective

			invert (float, optional): If True finds the max flux

		Output:
			(float) and (ndarray,optional): returns min flux model.objVal * -1 and if return_vec is True, returns
			model.objVal * -1, react_flux.X"""
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
		"""Description:
			Corrects a point, aligning it with bounds and flux balance constraint (Sv = b)

		Input:
			lb, ub, S, b:  Output of dicts_to_mats function

			env (gurobi env):  Gurobi environment

			point_to_use (ndarray): Point to correct

		Output:
			(ndarray): corrected point"""
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

	def generate_warmup_gb_pinch(self, max_search_number, sampling_levels, pinch_tol=1e-6, pinch_1_norm=None,
								 repeat_basis=1,
								 points_save_path="", model_save_path="",
								 required_rank=np.inf, rxn_penalty_vec=None, unique_vector_threshold_exponent=9,
								 trans_NS=None, rxn_scan=False,
								 **kwargs):
		"""Description:
			Generates warmup vectors for coordinate direction hit and run

		Input:
			max_search_number (int): integer value which determines how many proposed vectors are subsequently rejected
			before the algorithm terminates.

			sampling_levels (list): Flux is limited by 1 norm which is defined relative to the min possible flux, this
			is a list giving the max flux relative to the min, e.g., [1,1.1,5] would do sampling at the min,
			1.1x the min and 5x the min.

			pinch_tol (float, optional): if lb and ub are under this far apart pinches the reaction to a constant value

			pinch_1_norm (float, optional): the largest multiplicative factor placed on the minimum reaction penalty
			pinching is done with this 1 norm

			repeat_basis (int, optional): repeats sampling this many times, helps to obtain a more complete set of vectors

			points_save_path (path, optional): output path for the warmup points

			points_save_path (path, optional): output path for the output model

			required_rank (int, optional): if given allows model to terminate when a rank is reached

			rxn_penalty_vec (ndarray, optional): Used to weight reactions differently

			unique_vector_threshold_exponent (int, optional): gives the number of decimal places out the algorithm searches when
			removing identical vectors

			trans_NS (ndarray, optional): the exchange removed null space of S, if provided will test for thermo
			feasibility

			rxn_scan (boolean, optional): If true adds points with min/max of each reaction

		Output:
			(ndarray): set of warmup points

		Notes:
			Also creates warmup points in file as well as model to be used with sampler"""

		lb_full, ub_full, S_full, b, rxn_list, met_list = self.dicts_to_mats()

		if pinch_1_norm is None:
			pinch_1_norm = max(sampling_levels)

		var_ind, pinch_b, lb_full, v_c, ub_full = self.quick_pinch(pinch_1_norm, tol=pinch_tol)
		print(lb_full[-1], ub_full[-1], v_c[-1], var_ind)
		print("test3")
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
		print(lb)
		print(lb_full)
		print(var_ind)
		print(np.shape(lb_full), np.shape(lb), np.shape(var_ind))
		# input()
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

					# print(model.Status)

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
					# print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
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
		if rxn_scan:
			ub[-1] = max(sampling_levels) * min_flux
			# if you want purely slices use
			# lb[-1] = sampling_level * min_flux
			model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")

			# Suppress printing to console
			model.Params.LogToConsole = 0
			for i in range(S.shape[1] * 2):
				model.reset()
				obj = np.zeros(S.shape[1])
				if i < S.shape[1]:
					obj[i] -= 1
				else:
					obj[i - S.shape[1]] += 1
				# for variable objective total flux
				# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
				model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
				model.Params.LogToConsole = 0
				model.Params.NumericFocus = 3
				model.Params.FeasibilityTol = 1e-8

				model.optimize()

				# print(model.Status)
				point_to_use = react_flux.X
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
				total_points = np.vstack((total_points, point_to_use))
				# accepted_points, ortho_norm_basis, rank = self.lin_dep_np_ctp(ortho_norm_basis, origin,
				#															  accepted_points, point_to_use,
				#															  rank,
				#															  unique_vector_threshold_exponent=unique_vector_threshold_exponent)
				print(rank, i)
		print('done')
		plt.hist(total_points[:, -2])
		plt.show()
		input()
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
		if trans_NS is not None:
			for i in full_total_points:
				# print(i)
				lb_full, ub_full, S_full, b, rxn_list, met_list = self.dicts_to_mats()
				simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = self.precomp_positive_flux_point_to_bi_dir(
					rxn_list, prune_specific=["biomass_reaction"])

				ent_flux_point = self.positive_flux_point_to_bi_dir(i[:-1], simp_neg_ind,
																	comp_neg_ind, comp_perm,
																	cutter,
																	exchange_cutter=exchange_cutter)
				# print(np.shape(S_full))
				# print(cutter)
				# print(np.shape(ent_flux_point))
				# print(np.shape(trans_NS))
				if self.generate_therm_model_new(trans_NS, ent_flux_point) == 2:
					print("good")

		self.model_dict["total_flux_limits"] = [lb_full[-1], ub_full[-1]]
		if points_save_path != "":
			np.save(points_save_path / ("warmup_" + self.model_dict["model_name"]), full_total_points)
		if model_save_path != "":
			self.save_model_as_fast_key_json(model_save_path)
		return full_total_points

	def generate_warmup_gb_pinch_scan_br(self, max_search_number, sampling_levels, biomass_values,
										 biomass_reaction_name, pinch_tol=1e-6, pinch_1_norm=None,
										 repeat_basis=1,
										 points_save_path="", model_save_path="",
										 required_rank=np.inf, rxn_penalty_vec=None, unique_vector_threshold_exponent=9,
										 trans_NS=None, rxn_scan=False,
										 **kwargs):
		"""Description:
			Generates warmup vectors for coordinate direction hit and run

		Input:
			max_search_number (int): integer value which determines how many proposed vectors are subsequently rejected
			before the algorithm terminates.

			sampling_levels (list): Flux is limited by 1 norm which is defined relative to the min possible flux, this
			is a list giving the max flux relative to the min, e.g., [1,1.1,5] would do sampling at the min,
			1.1x the min and 5x the min.

			pinch_tol (float, optional): if lb and ub are under this far apart pinches the reaction to a constant value

			pinch_1_norm (float, optional): the largest multiplicative factor placed on the minimum reaction penalty
			pinching is done with this 1 norm

			repeat_basis (int, optional): repeats sampling this many times, helps to obtain a more complete set of vectors

			points_save_path (path, optional): output path for the warmup points

			points_save_path (path, optional): output path for the output model

			required_rank (int, optional): if given allows model to terminate when a rank is reached

			rxn_penalty_vec (ndarray, optional): Used to weight reactions differently

			unique_vector_threshold_exponent (int, optional): gives the number of decimal places out the algorithm searches when
			removing identical vectors

			trans_NS (ndarray, optional): the exchange removed null space of S, if provided will test for thermo
			feasibility

			rxn_scan (boolean, optional): If true adds points with min/max of each reaction

		Output:
			(ndarray): set of warmup points

		Notes:
			Also creates warmup points in file as well as model to be used with sampler"""
		lb_full, ub_full, S_full, b, rxn_list, met_list = self.dicts_to_mats()

		if pinch_1_norm is None:
			pinch_1_norm = max(sampling_levels)

		var_ind, pinch_b, lb_full, v_c, ub_full = self.quick_pinch(pinch_1_norm, tol=pinch_tol)
		print(lb_full[-1], ub_full[-1], v_c[-1], var_ind)
		print(ub_full[rxn_list.index(biomass_reaction_name)])
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
		print(lb)
		print(lb_full)
		print(var_ind)
		print(np.shape(lb_full), np.shape(lb), np.shape(var_ind))
		slice_rxn_list = np.array(rxn_list)[var_ind]
		print(ub[list(slice_rxn_list).index(biomass_reaction_name)])
		# input()
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

					# print(model.Status)

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
					# print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
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
		print(ub[list(slice_rxn_list).index(biomass_reaction_name)])
		for biomass_level in np.linspace(lb[list(slice_rxn_list).index(biomass_reaction_name)],
										 ub[list(slice_rxn_list).index(biomass_reaction_name)], biomass_values):
			# Defines origin and points as empty

			accepted_points = np.array([])
			ortho_norm_basis = np.array([])

			ub[-1] = max(sampling_levels) * min_flux
			# if you want purely slices use
			# lb[-1] = sampling_level * min_flux
			model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
			model.addConstr(react_flux[list(slice_rxn_list).index(biomass_reaction_name)] == biomass_level, name="c")

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

				# print(model.Status)

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
				# print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
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
			print(accepted_points[:, list(slice_rxn_list).index(biomass_reaction_name)])
		if rxn_scan:
			ub[-1] = max(sampling_levels) * min_flux
			# if you want purely slices use
			# lb[-1] = sampling_level * min_flux
			model = gp.Model("warmup_model_" + self.model_dict["model_name"], env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")

			# Suppress printing to console
			model.Params.LogToConsole = 0
			for i in range(S.shape[1] * 2):
				model.reset()
				obj = np.zeros(S.shape[1])
				if i < S.shape[1]:
					obj[i] -= 1
				else:
					obj[i - S.shape[1]] += 1
				# for variable objective total flux
				# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
				model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
				model.Params.LogToConsole = 0
				model.Params.NumericFocus = 3
				model.Params.FeasibilityTol = 1e-8

				model.optimize()

				# print(model.Status)
				point_to_use = react_flux.X
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
				total_points = np.vstack((total_points, point_to_use))
				# accepted_points, ortho_norm_basis, rank = self.lin_dep_np_ctp(ortho_norm_basis, origin,
				#															  accepted_points, point_to_use,
				#															  rank,
				#															  unique_vector_threshold_exponent=unique_vector_threshold_exponent)
				print(rank, i)
		print('done')

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
		plt.hist(total_points[:, -2])
		plt.show()
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
		# input()
		print(3)
		if trans_NS is not None:
			for i in full_total_points:
				# print(i)
				lb_full, ub_full, S_full, b, rxn_list, met_list = self.dicts_to_mats()
				simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = self.precomp_positive_flux_point_to_bi_dir(
					rxn_list, prune_specific=["biomass_reaction"])

				ent_flux_point = self.positive_flux_point_to_bi_dir(i[:-1], simp_neg_ind,
																	comp_neg_ind, comp_perm,
																	cutter,
																	exchange_cutter=exchange_cutter)
				# print(np.shape(S_full))
				# print(cutter)
				# print(np.shape(ent_flux_point))
				# print(np.shape(trans_NS))
				if self.generate_therm_model_new(trans_NS, ent_flux_point) == 2:
					print("good")

		self.model_dict["total_flux_limits"] = [lb_full[-1], ub_full[-1]]
		if points_save_path != "":
			np.save(points_save_path / ("warmup_" + self.model_dict["model_name"]), full_total_points)
		if model_save_path != "":
			self.save_model_as_fast_key_json(model_save_path)
		return full_total_points

	def gene_point_penalty(self, flux_vector, RPS, constant):
		"""Description:
			Generates a penalty term for a gene expression given a specified reaction penalty score

		Input:
			flux_vector (ndarray): point in flux space

			RPS (ndarray): reaction penalties (1/RAS), should be as large as the flux_vector

			constant (float, optional): error term is multiplied by this value, increasing its intensity

		Output:
			(float): penalty of the point in flux space"""
		if RPS is None:
			error_term = 0.0
		else:
			error_term = np.dot(flux_vector, RPS)
		# print(RPS)
		# print(flux_vector)
		# input()
		return -error_term * constant

	def RAS_to_RPS_mat(self, RAS, zero_correction_factor=2, nan_val=None):
		"""Description:
			Converts from RAS to RPS

		Input:
			RAS (ndarray): ndarray of RAS values

			zero_correction_factor (float, optional): makes 0 RAS terms zero_correction_factor times greater than the highest
			true penalty score

			nan_val (float, optional): Any nan value in RAS is set to this value, if non the median RAS is used.

		Output:
			(ndarray): ndarray of RPS values"""
		if nan_val is None:
			nan_val = np.nanmedian(RAS)
		RAS = np.nan_to_num(RAS, copy=True, nan=nan_val)
		RAScp = cp.deepcopy(RAS)
		RAScp[RAScp == 0] = np.inf
		minnonzero = np.min(RAScp)
		print(minnonzero)
		RAS[RAS == 0] = (minnonzero / zero_correction_factor)
		return 1 / RAS

	def HRSampler_gene_bias_lincomb_pinch(self, origin_and_warmup_points, number_of_points, stepsPerPoint, RPS,
										  max_sampling_level, pinch_tol=1e-6,
										  gene_penalty_mod=1, thermo_const=None, rxn_penalty_vec=None, term_val=1e5,
										  print_percent=True, print_error=True, **kwargs):
		"""Description:
			Performs coordinate direction hit and run biased by gene data and constrained to obey thermodynamic limitations

		Input:
			origin_and_warmup_points (ndarray): array containing origin and warmup points

			number_of_points (int): number of points to generate

			stepsPerPoint (int): number of points to skip between saved points

			RPS (ndarray): ndarray of RPS values

			max_sampling_level (float): the largest multiplicative factor placed on the minimum reaction penalty

			pinch_tol (float, optional): if lb and ub are under this far apart pinches the reaction to a constant value

			gene_penalty_mod (float, optional): multiplies the RPS score by this value, used to make gene constraints
			more intense

			thermo_const (dict, optional): If a dict is given, ensures each point it thermodynamically feasible.
			Dict should contain "prune_specific" and "NS_internal" keys.
			E.g., thermo_const={"prune_specific": ["biomass_reaction"],"NS_internal": Internal_nullspace_S}

			rxn_penalty_vec (ndarray, optional): provides weights of reactions for minimization

			term_val (int, optional): Number of times a passible initial point is tried before resetting

			print_percent (boolean, optional): If true prints percent completion

			print_error (boolean, optional): If true prints percent error terms

		Output:
			(ndarray): ndarray of flux points"""
		lb_full, ub_full, S_full, b, rxn_list, met_list = self.dicts_to_mats()
		var_ind, pinch_b, lb_full, v_c, ub_full = self.quick_pinch(max_sampling_level, tol=pinch_tol)

		flux_min = lb_full[-1]
		flux_max = ub_full[-1]

		lb_ntf = lb_full[:-1]
		ub_ntf = ub_full[:-1]

		v_c = v_c[:-1]
		v_c[var_ind] = np.ones_like(var_ind) * np.nan

		lb_full[-1] += np.nansum(v_c)
		ub_full[-1] += np.nansum(v_c)

		v_c_fillin = cp.deepcopy(v_c)
		lb = lb_ntf[var_ind]
		ub = ub_ntf[var_ind]
		S = S_full[:, var_ind]
		if RPS is not None:
			RPS = RPS[var_ind]

		# input(2)
		b = pinch_b

		# print(self.model_dict["total_flux_limits"])

		# input()

		if thermo_const is not None:
			simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = self.precomp_positive_flux_point_to_bi_dir(
				rxn_list, prune_specific=thermo_const["prune_specific"])
			# test_bi_S = self.positive_S_to_bi_dir(S_full, simp_neg_ind, cutter, exchange_cutter=exchange_cutter)

			NS = thermo_const["NS_internal"]
			trans_NS = np.transpose(NS)

		origin_and_warmup_points_1 = origin_and_warmup_points[:, var_ind]
		origin_and_warmup_points_2 = np.reshape(np.sum(origin_and_warmup_points_1, axis=1), (-1, 1))

		origin_and_warmup_points = np.hstack((origin_and_warmup_points_1, origin_and_warmup_points_2))

		uTol = 1e-7
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

			# rxn_penalty_vec = rxn_penalty_vec[var_ind]
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, flux_min))
			ub = np.hstack((ub, flux_max + 1e-3))

			# print(lb_full)
			b = np.hstack((b, 0))
			S[-1, -1] = -1

		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		N = sc.linalg.null_space(S)
		indexer = np.arange(0, np.size(lb))
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
			# print(np.sum(S))
			# print(np.sum(b))
			# print(f"spCurr WARNING: fidErr {np.max(np.abs(np.matmul(S, mid_point) - b))}")
			# print(f"spCurr Error ub: {max(mid_point - ub)}, lb: {max(lb - mid_point)}")
			while not move_forward:
				count += 1
				if count % 100 == 0:
					if print_error:
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
						if print_error:
							print(f"spCurr Error ub: {max(mid_point - ub)}, lb: {max(lb - mid_point)}")
						# input()
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
						if print_error:
							print(1, valid_point, np.min(prev_point) >= 0)

				if valid_point and ((np.max(np.abs(np.matmul(S, prev_point) - b)) > Sv_tol) or (
						max(prev_point - ub) > lb_ub_tol_barrier or max(lb - prev_point) > lb_ub_tol_barrier)):
					if print_error:
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
		if print_error:
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
		accepted_point_traj = []
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
					# print(prev_point)
					prev_point = np.where(prev_point > ub, ub, prev_point)
					prev_point = np.where(prev_point < lb, lb, prev_point)

					distUb = (ub - prev_point)
					distLb = (prev_point - lb)

					# test = np.random.multinomial(1+np.random.choice(np.shape(norm_warmup_vecs)[0],1)[0],np.ones(np.shape(norm_warmup_vecs)[0])/np.shape(norm_warmup_vecs)[0])
					# test = np.random.multinomial(3,np.ones(np.shape(norm_warmup_vecs)[0])/np.shape(norm_warmup_vecs)[0])

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
						if ts % 5000 == 0:
							if print_error:
								print("0 bounds, count: ", ts)
						ts += 1
						if ts > term_val:
							if print_error:
								print("TERMINATING DUE TO NO VALID POINTS")
							return np.array([])
						continue
					if maxStep - minStep <= 0:
						if print_error:
							print("WARNING")
							print(maxStep, minStep)
							print(maxStep - minStep)
						min_ind = np.argmin(pos_u_rel_dis_ub)
						if print_error:
							print(np.compress(posValidDir, distUb, axis=0)[min_ind])
							print(np.compress(posValidDir, u, axis=0)[min_ind])
							print(np.compress(posValidDir, distUb / u, axis=0)[min_ind])
						# print(np.compress(posValidDir, u, axis=0))
						# print(np.compress(posValidDir, distUb, axis=0)/np.compress(posValidDir, u, axis=0))

						min_ind = np.argmin(neg_u_rel_dis_lb)
						if print_error:
							print(np.compress(negValidDir, distUb, axis=0)[min_ind])
							print(np.compress(negValidDir, u, axis=0)[min_ind])

						max_ind = np.argmax(pos_u_rel_dis_lb)
						if print_error:
							print(np.compress(posValidDir, distLb, axis=0)[max_ind])
							print(np.compress(posValidDir, u, axis=0)[max_ind])

						max_ind = np.argmax(neg_u_rel_dis_ub)
						if print_error:
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
					if print_error and step_count % 1000 == 0:
						print("start", np.sum(np.abs(np.matmul(S, curPoint) - b)))
					if (np.max(np.abs(np.matmul(S, curPoint) - b)) > Sv_tol) or max(
							curPoint - ub) > lb_ub_tol_barrier or max(
						lb - curPoint) > lb_ub_tol_barrier:
						correct_count += 1
					if correct_count % 50 == 0:
						if print_error:
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

						if print_error and correct_count % 50 == 0:
							print(f"rrr WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint) - b))}")
							print(f"rrr WARNING: Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
					if print_error and step_count % 1000 == 0:
						print(step_count)
						print("stop", np.sum(np.abs(np.matmul(S, curPoint) - b)))
					cur_support = self.gene_point_penalty(curPoint[:-1], RPS, gene_penalty_mod)
					prev_support = self.gene_point_penalty(prev_point[:-1], RPS, gene_penalty_mod)

					proposed_points += 1
					# print(cur_support,prev_support,np.log(randVector[1]))
					if print_error and proposed_points % 5000 == 0:
						print(proposed_points, total_steps, (total_steps + 1) / proposed_points,
							  therm_reject / proposed_points, gene_reject / proposed_points)
						print(prev_point[-1], curPoint[-1], lb[-1], ub[-1])
						print(prev_support)
						print(f"CURR WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint) - b))}")
						print(f"CURR Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")
						print(f"PREV WARNING: fidErr {np.max(np.abs(np.matmul(S, prev_point) - b))}")
						print(f"PREV Error ub: {max(prev_point - ub)}, lb: {max(lb - prev_point)}")

					if np.min(curPoint) >= 0:
						# print(np.log(randVector[1]) , (cur_support - prev_support),cur_support,prev_support)

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
								# print(val)
								# input()
								if val == 2:
									passing = False
								else:
									therm_reject += 1
									continue
							else:
								passing = False
								continue
						else:
							# print(0)
							# input()
							gene_reject += 1
							continue
					else:
						continue

				if print_error and total_steps % 200 == 0:
					if max(max(curPoint - ub), max(lb - curPoint)) > 1e-7:
						print(f"WARNING: fidErr {np.max(np.abs(np.matmul(S, curPoint) - b))}")
						print(f"Error ub: {max(curPoint - ub)}, lb: {max(lb - curPoint)}")

				# Both values should be negative

				time_elapsed = time.time() - t0

				if print_error and step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")

				# curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				if print_error and total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint) - b))}")
				prev_point = curPoint
				step_count += 1

				total_steps += 1
				if print_percent and np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed, flush=True)
					print(curPoint[-1], flush=True)
				# print(accepted_point_traj)
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)

				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count] = curPoint

			point_count += 1
		# input()
		full_total_points = np.zeros((np.shape(points)[0], np.size(v_c) + 1))
		for i in range(np.shape(full_total_points)[0]):
			v_c[var_ind] = points[i, :-1]
			full_total_points[i, :-1] = v_c
			full_total_points[i, -1] = np.sum(full_total_points[i, :-1])
			# print(np.shape(v_c))
			# print(np.shape(S_full), np.shape(full_total_points[i]))
			if print_error and (np.max(np.abs(np.matmul(S_full, full_total_points[i]))) > Sv_tol) or max(
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

		return full_total_points

	def precomp_positive_flux_point_to_bi_dir(self, pos_rxn_list, prune_exchange=True, prune_specific=None):
		"""Description:
			Creates ndarrays of indicies used to rapidly convert from positive fluxes to bidirectional fluxes

		Input:
			pos_rxn_list (list): list of reaction names for model of positive flux points

			prune_exchange (boolean, optional): if true removes exchange reactions

			prune_specific (int, optional): pass specific reactions to prune

		Output:
			(simp_neg_ind), (comp_neg_ind), (comp_perm), (cutter), (exchange_cutter): Returns arrays which are used to
			convert from positive fluxes to bidirectional fluxes"""
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
		"""Description:
			Quickly converts from positive fluxes to bidirectional fluxes

		Input:
			point (ndarray): flux point

			simp_neg_ind (ndarray): array with same length as points with -1 where reaction is inverted and no
			forward reaction exists and 1 everywhere else.

			comp_neg_ind (ndarray): array with same length as points with -1 where reaction is inverted and
			forward reaction exists and 0 everywhere else.

			comp_perm (ndarray): array with same length as points with index of inverted reactions which have non
			inverted reactions also present, all other indicies are placeholders which correspond to 0 values

			cutter (ndarray): array which gives index values for reactions which should be kept due to not being a
			repeat (removing reactions which are reversed and have forward directions in the set)

			exchange_cutter (ndarray): array which gives index values for reactions which should be kept due to not
			being an exchange (removing exchange reactions)

		Output:
			(ndarray): flux of bidirectional flux point"""
		bi_point = copy.deepcopy(point)
		comp_neg_ind_cp = copy.deepcopy(comp_neg_ind)
		bi_point *= simp_neg_ind
		comp_neg_ind_cp *= bi_point
		comp_neg_ind_cp = comp_neg_ind_cp[comp_perm]
		bi_point += comp_neg_ind_cp
		nbi_point = bi_point[cutter]
		print(cutter)
		print(exchange_cutter)
		input(0)
		if exchange_cutter is not None:
			nbi_point = nbi_point[exchange_cutter]
		return nbi_point

	def convert_model_to_positive_flux(self):
		"""Description:
			Converts model to only having positive reactions

		Notes:
			alters the model this function is run on to only have positive fluxes"""
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
		"""Description:
			Converts model from positive flux model to bidirectional flux model

		Notes:
			alters the model this function is run on to have both flux directions"""
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

	def fit_to_experimental_data(self, experimental_data_file, alpha_array, assign_opt_alpha_fluxes=True, tol=1e-7,
								 search_save_path="", method=-1):
		"""Description:
			Fits model to measured fluxes

		Input:
			experimental_data_file (path): path to file containing experimental data

			alpha_array (ndarray): array containing alpha values to use, error for experimental alignment is given by
			(vi-aei)^2/(aei)^2 = vi^2/(aei)^2-2vi/aei+1

			assign_opt_alpha_fluxes (boolean, optional): If true sets the model to have the aligned fluxes after they
			are found

			tol (float): tolerance given to model set experimental fluxes

			search_save_path (path): Gives path to save search data

			method (int): gurobi model method

		Notes:
			Generates model with updated bounds and/our search information"""
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
			if method != -1:
				model.Params.Method = method
			# model.Params.BarHomogeneous = 0
			model.Params.LogToConsole = 0
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb,
									   ub=ub)
			model.setObjective(react_flux @ Q @ react_flux - c @ react_flux + amount_to_add_back, GRB.MINIMIZE)
			rhs = np.transpose(b)
			model.addConstr(S @ react_flux == rhs, name="c")
			model.optimize()
			print(model.Status)
			if model.Status == 2 or model.Status == 13:
				array_of_states[alpha_ind] = model.X
				ind = copy.deepcopy(model_index_of_experimental_fluxes)
				labels = [ordered_rxn[i] for i in model_index_of_experimental_fluxes]
				# print(error)
				# print(ind)
				if 'EX_val_L(e)' in labels:
					del ind[labels.index('EX_val_L(e)')]
				elif 'EX_val__L_e' in labels:
					del ind[labels.index('EX_val__L_e')]
				# print(ind)
				expt = np.array([e[i] for i in ind])
				mdl = np.array([model.X[i] for i in ind])
				print("interest", (mdl - expt * alpha_array[alpha_ind]) ** 2 / (expt * alpha_array[alpha_ind]) ** 2)
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

	def pinch_restricted_exchange_reactions(self, allowed_exchange_reaction_file, restore_essential=True,
											true_exchange_signifier="", restore_previously_pinched=False,
											restore_value=1000):
		"""Description:
			Removes unwanted exchange reactions

		Input:
			allowed_exchange_reaction_file (path): path to file containing allowed exchange reactions, all other
			reactions are removed

			restore_essential (boolean, optional): If true, if removing an exchange makes model unfeasible, re-adds the
			reaction

			true_exchange_signifier (string, optional): If a reaction name contains some type of substring which conveys the
			reaction is truly exchange ("EX" for example).

			restore_previously_pinched (boolean, optional): If true unpinches all exchange reactions at the start

			restore_value (float, optional): If a flux is restored, it is given this value

		Notes:
			alters the model this function is run on to remove unwanted exchange reactions"""
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
						input()
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

	def create_RAS_values(self, gene_array, feature_array, return_dict=True, identifier="entrez"):
		"""Description:
			Creates RAS values

		Input:
			gene_array (ndarray): Gene UMI array

			feature_array (ndarray): array containing gene names, should have format corresponding to identifier

			return_dict (boolean, optional): If true returns dict with keys given as reactions and RAS value as
			dict value

			identifier (string, optional): gene nomenclature

		Output:
			(ndarray): array of RAS values, indices correspond to reaction names
			or (dict) if return_dict=True: dict with keys given as reactions and RAS value as
			dict value"""
		def create_entrez_to_express_dict(gene_names, gene_array, entrez_dict, identifier):
			print(gene_names)
			mouse_gene_list = list(gene_names[:, 0])

			list_unmatch = []
			list_m_ens_unmatch = []
			print(entrez_dict.keys())
			print(entrez_dict)
			for i in entrez_dict.keys():
				if identifier in entrez_dict[i].keys():

					if entrez_dict[i][identifier] in mouse_gene_list:
						match_ind = mouse_gene_list.index(entrez_dict[i][identifier])
						entrez_dict[i]['expression_val'] = gene_array[match_ind] / entrez_dict[i]["degen"]
					# print(i, entrez_dict[i], 1)
					# print(avg_exp[mouse_gene_list.index(entrez_dict[i]["mouse_ensemble"])], entrez_dict[i]["degen"],entrez_dict[i]["expression_val"])
					else:
						# print(i, entrez_dict[i], 2)
						entrez_dict[i]['expression_val'] = 0
						list_m_ens_unmatch.append([i, entrez_dict[i][identifier]])
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

		def convert_and_internal(matchobj):
			# converts all pure and statements enclosed in parenthesis to value
			# print(matchobj)
			# print(1)
			ands = matchobj.group(0)[1:-1]
			# print(2)
			ands = ands.replace(" ", "")
			# print(3)
			if ands.count("or") > 0:
				print(ands)
				print("ERROR")
				input()
			# print(4)
			ands = np.average(np.array([float(i) for i in ands.split("and")]))
			return str(ands)

		def convert_and_external(matchobj):
			# converts all pure and statements enclosed in parenthesis to value
			ands = matchobj.group(0)
			ands = ands.replace(" ", "")
			if ands.count("or") > 0:
				print(ands)
				print("ERROR")
				input()
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
			if gr_rule == "":
				# print(np.nan)
				return np.nan
			name_to_real_express = partial(name_to_express, entrez_dict)
			# print("start")
			# print(gr_rule, 1)
			gr_rule = re.sub("\d+\.\d+", name_to_real_express, gr_rule)
			# print(gr_rule, 2)
			ext_altered = True
			while ext_altered:
				int_or_altered = True
				while int_or_altered:
					gr_rule_cp_2 = copy.deepcopy(gr_rule)
					and_altered = True
					while and_altered:
						gr_rule_cp = copy.deepcopy(gr_rule)
						gr_rule = re.sub("\(\d+\.?\d*e?-?\d*\)", strip_pare, gr_rule)
						# print(gr_rule, 3)
						gr_rule = re.sub(
							"\((?=(?P<tmp>((\d+\.?\d*e?-?\d*)) *))(?P=tmp)((?=(?P<qrs>( *(and)( *\d+\.?\d*e?-?\d*))))(?P=qrs))+\)",
							convert_and_internal, gr_rule)
						# print(gr_rule, 4)
						if gr_rule_cp == gr_rule:
							and_altered = False
					gr_rule = re.sub("\(\d+\.?\d*e?-?\d*\)", strip_pare, gr_rule)
					# print(gr_rule, 5)
					gr_rule = re.sub(
						"\((?=(?P<tmp>((\d+\.?\d*e?-?\d*)) *))(?P=tmp)((?=(?P<qrs>( *(or)( *\d+\.?\d*e?-?\d*))))(?P=qrs))+\)",
						convert_or_internal, gr_rule)
					# print(gr_rule, 6)
					if gr_rule_cp_2 == gr_rule:
						int_or_altered = False
				gr_rule_ext = copy.deepcopy(gr_rule)
				gr_rule = re.sub(
					"(?=(?P<tmp>((\d+\.?\d*e?-?\d*)) *))(?P=tmp)((?=(?P<qrs>( *(and)( *\d+\.?\d*e?-?\d*))))(?P=qrs))+",
					convert_and_external, gr_rule)
				gr_rule = re.sub(
					"(?=(?P<tmp>((\d+\.?\d*e?-?\d*)) *))(?P=tmp)((?=(?P<qrs>( *(or)( *\d+\.?\d*e?-?\d*))))(?P=qrs))+",
					convert_or_external, gr_rule)
				if gr_rule_ext == gr_rule:
					ext_altered = False
			# print(gr_rule, 7)

			return float(gr_rule)

		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		gr_rule = self.get_grRule_list()
		print(gr_rule)
		# print(self.model_dict["gene_dict"])
		print(feature_array)
		print(self.model_dict["gene_dict"])
		entrez_express_dict = create_entrez_to_express_dict(feature_array, gene_array,
															cp.deepcopy(self.model_dict["gene_dict"]), identifier)
		print(entrez_express_dict)
		print(0)
		RAS_list = []
		RAS_nz = []
		test_ras = 0
		used_list = []
		for i in range(len(gr_rule)):

			RAS = convert_gr_Rule_to_RAS(gr_rule[i], entrez_express_dict)
			RAS_list.append(RAS)
			if rxn_list[i].replace("_inverted", "") not in used_list and RAS > 0:
				test_ras += RAS
				used_list.append(rxn_list[i].replace("_inverted", ""))
			if rxn_list[i] == 'PGLYCP' or rxn_list[i] == 'PGLYCP_inverted':
				# print(entrez_express_dict["159.1"])
				print(gr_rule[i], RAS)
		if return_dict:
			RAS_dict = dict(zip(rxn_list, RAS_list))
			return RAS_dict
		else:
			return RAS_list

	def test_essential_list(self, reaction_index_list):
		"""Description:
			Test a list of essential reactions to make sure they are minimal and essential

		Input:
			reaction_index_list (ndarray): list of reactions which should constitute minimal model

		Output:
			(boolean) (string): boolean is true if model is feasible, string contains informative message
			if boolean is false, explains the issue"""
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

	def get_met_comp_dict(self, met_name):
		"""Description:
			Not currently in use, used to get the chemical composition of a metabolite as a dicitonary"""
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
		"""Description:
			Not currently in use, would be used to get metabolite molecular weight"""
		mol_weight_dict = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "P": 30.97376, "S": 32.06, "Ca": 40.078,
						   "Cl": 35.45, "K": 39.0983, "Na": 22.98977, "Se": 78.971, "Co": 58.93319, "I": 126.9045,
						   "Fe": 55.845, "R": 12.011 * 15 + 1.008 * 31, "Ra": (12.011 * 15 + 1.008 * 31) * 1,
						   "Rb": (12.011 * 15 + 1.008 * 31) * 1, "Rc": (12.011 * 15 + 1.008 * 31) * 1, "X": 0, "Y": 0}
		chem_count_dict = self.get_met_comp_dict(met_name)
		mw = 0
		for chem in chem_count_dict.keys():
			mw += int(chem_count_dict[chem]) * mol_weight_dict[chem]
		return mw

	def remove_unused_met(self):
		"""Description:
			Removed metabolites which do not occur in any reactions

		Notes:
			alters the model this function is run on to remove unneeded metabolites"""
		met_names = list(self.model_dict['met_dict'].keys())
		for met_name in met_names:
			# print(recon_flux_model.model_dict['rxn_dict'])
			used = False
			for rxn_name in self.model_dict['rxn_dict'].keys():
				if met_name in self.model_dict['rxn_dict'][rxn_name]['rxn_metabolites'].keys():
					used = True
			if not used:
				del self.model_dict['met_dict'][met_name]
