#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Developed With Python Version 3.5

"""

# Maybe The Max Eigenvalue Being Low is Caused by The FutureWarnings I'm Getting...but Probably Not.
# import warnings
# warnings.simplefilter(action = "ignore", category = FutureWarning)

# Import Core Modules
import csv 					as CSV
import os 					as OS
import sys 					as Sys
import os.path 				as Path
import errno 				as Error
import json 				as JSON
import re 					as Rgx
import pickle 				as Pickle
import random 				as Random

# Import Third-Party Modules
import numpy 				as Npy
import numpy.linalg 		as LinearAlgebra
import networkx 			as NX
import pandas 				as Pandas

# Import Functions
from collections 			import deque
from networkx.algorithms 	import approximation
from pgmpy.models 			import MarkovModel
from pgmpy.factors 			import Factor
from pgmpy.inference 		import VariableElimination
from pgmpy.inference 		import BeliefPropagation

# 
# File Handling Utility Methods
# 

# General File Writing Function
def writeOutput (data, path="output/output.txt", newline = True):
	"""
	Standard File Writing Function
	@params:
		data 		- Required	: data to write (Str, Dict, List, Set)
		path 		- Optional 	: output path (Str)
	Returns: path (Str)
	"""
	with open(path, "w") as file:
		for line in data:
			if newline == True:
				file.write(str(line) + "\n")
			else:
				file.write(str(line))
		file.close()
	return path

# Read Delimited File Columns
def readColsFromXSV(path, colStart = 0, colStop = 0, colStep = 1, keepHeader = True, delimiter = ",", openAs = "r", progress = False):
	"""
	Read Delimited File to
	@params:
		path 		- Required 	: input path (Str)
		delimiter 	- Optional 	: file delimiter (Str)
		openAs 		- Optional 	: file opening flag (Str)
	Returns: Row values (List)
	"""
	rowValues = []
	loadSuffix = "Loaded (" + path + ")"
	if progress == True:
		if keepHeader == False:
			total = open(path).read().count("\n") - 1
		else: 
			total = open(path).read().count("\n")
		with open(path, openAs) as xsv:
			reader = CSV.reader(xsv, delimiter = delimiter)
			if keepHeader == False:
				next(reader)
			i = 0
			for row in reader:
				cols = row[colStart:colStop:colStep]
				if len(cols) > 0:
					rowValues.append(cols)
				i += 1
				printProgress(i, total, prefix = "\t\tLoading File:" , suffix = loadSuffix, decimals = 0, barLength = 50)
				if i == total:
					print("\n")
	else:
		with open(path, openAs) as xsv:
			reader = CSV.reader(xsv, delimiter = delimiter)
			if keepHeader == False:
				next(reader)
			for row in reader:
				cols = row[colStart:colStop:colStep]
				if len(cols) > 0:
					rowValues.append(cols)
	return rowValues

# Prompt user input from command line
def getUserInput (valid, prompt, hint = "", failed = "Error: Invalid input"):
	"""
	Prompts user for and validates input using regular expression
	@params:
		prompt 		- Required 	: verbose user prompt (Str)
		valid 		- Required 	: regex to validate against (Rgx)
	Returns: dicts (List)
	"""
	print(prompt)
	if len(hint) > 0:
		response = input(hint)
	else:
		response = input()
	if Rgx.match(valid, response):
		return response
	else:
		print(failed)
		return getUserInput(valid, prompt, failed)

# Print iterations progress
def printProgress (iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 100):
	"""
	Call in a loop to create terminal progress bar
	@params:
		iterations 	- Required 	: current iteration (Int)
		total 		- Required 	: total iterations (Int)
		prefix 		- Optional	: prefix string (Str)
		suffix 		- Optional 	: suffix string (Str)
	"""
	filledLength 	= int(round(barLength * iteration / float(total)))
	percents 		= round(100.00 * (iteration / float(total)), decimals)
	bar 			= '#' * filledLength + '-' * (barLength - filledLength)
	Sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
	Sys.stdout.flush()

# 
# Graph Utility Methods
# 

def generateRandomComponent (populationGraph, nodes = 1000):
	""" 
	Generate Random Graphs
	@params:
		popualationGraph 	- Required 	: (Dictionary)
		nodes 				- Optional 	: (Int)
	"""
	sampleNodes = populationGraph.nodes()
	sample 	= Random.sample(sampleNodes, nodes)
	graph 	= {node: [n for n in populationGraph[node] if n in sample] for node in sample}
	return NX.Graph(graph)

def generateIntersectionSubgraph (edgeList, nodes):
	""" 
	Generate subgraph from edges with both nodes in node list
	@params:
		edgeList 			- Required 	: (List, Set)
		nodes				- Required 	: (List, Set)
	"""
	componentEdges = set()
	componentNodes = set()
	for a, b in edgeList:
		if a == b:
			continue
		aB, bB = False, False
		for node in nodes:
			if a == node:
				aB = True
			if b == node:
				bB = True
			if aB == True and bB == True:
				componentNodes.add(a), componentNodes.add(b)
				componentEdges.add((a, b))
	return componentEdges, componentNodes

def generateUnionSubgraph (edgeList, nodes):
	""" 
	Generate subgraph from edges with at least one nodes in node list
	@params:
		edgeList 			- Required 	: (List, Set)
		nodes				- Required 	: (List, Set)
	"""
	componentEdges = set()
	componentNodes = set()
	for a, b in edgeList:
		if a == b:
			continue
		for node in nodes:
			if a == node or b == node:
				componentNodes.add(a), componentNodes.add(b)
				componentEdges.add((a, b))
	return componentEdges, componentNodes

# 
# Main Routine
# 

if __name__ == "__main__":
	print("\nNotification: MRF pipeline initialized from commandline\n")
	exitCode 						= 0
	try:
	# Set Meta Data
		print("\tStatus: Setting Global Variables")
		version 						= 0.1
		interactionDataPath 			= "data/homo_sapiens.mitab.interactions.txt"
		oncoNetDataPath 				= "data/Census_allTue Dec  8 23_01_10 2015.csv"
		tcgaNetDataPath 				= "data/tcga-genes.txt"
		outputDirectory  				= "output/"
		objectDirectory 				= "objects/"
		usePickles 						= False
		useRandomComponent 				= False
		usePicklesResponse 				= getUserInput(valid=r"([yYnN]{1}|exit)", prompt="\tPrompt: Use pre-computed data structures when available?", hint = "\t\tHint: enter 'y', 'n', or 'exit' to exit\n\t\tSelection: ", failed = "\t\t\tError: Invalid input")
		if usePicklesResponse.lower() == "y":
			usePickles 					= True
		elif usePicklesResponse.lower() == "n":
			useRandomComponentResponse 		= getUserInput(valid=r"([yYnN]{1}|exit)", prompt="\tPrompt: Generate random component (instead of using big component of network)?", hint = "\t\tHint: enter 'y', 'n', or 'exit' to exit\n\t\tSelection: ", failed = "\t\t\tError: Invalid input")
			if useRandomComponentResponse.lower() == "y":
				useRandomComponent 			= True
				randomComponentSize 			= int(getUserInput(valid=r"[1-9][0-9]{0,6}", prompt="\tPrompt: How many nodes for generated random component?", hint = "\t\tHint: enter integer (n) between [1 - 9999999]\n\t\tSelection: ", failed = "\t\t\tError: Invalid input"))
			elif useRandomComponentResponse == "exit":
				exitCode = 1
				Sys.exit(exitCode)
		elif usePicklesResponse == "exit":
			exitCode = 1
			Sys.exit(exitCode)
	except OSError as error:
		print(error)
		print("Error: unable to configure pipeline.py.\n")
		Sys.exit(exitCode)
	try:

	# 
	# Import Data
	# 

		print("\tStatus: Importing Data")
		allInteractionPairs 			= []
		# Import Interactions Node Data
		for row in readColsFromXSV(interactionDataPath, colStart = 2, colStop = 4, keepHeader = False, delimiter = "\t", progress = True):
			# Convert Row to String (join), Extract Genes Using Rgx Substitution, Split Resulting String Into a List, Coerce into a Tuple
			interaction 				= tuple(Rgx.sub(r"uniprotkb\:[a-zA-Z0-9]{1,8}\_", "", "".join(row)).split("(shortlabel)")[0:2])
			# Remove self-edges, singletons, interactions with empty None values or a length less than or equal to one
			if len(interaction) == 2:
				if interaction[0] != interaction[1]:
					if (interaction[0] != None and len(interaction[0]) > 1) and (interaction[1] != None and len(interaction[1]) > 1):
						allInteractionPairs.append(interaction)
		oncoDataFull 					= set()
		# Import all onco genes as list
		for gene in readColsFromXSV(oncoNetDataPath, colStart = 0, colStop = 1, keepHeader = False, delimiter = ",", openAs = "U", progress = True):
			oncoDataFull.add(gene[0])
		# Import TCGA onco genes as list
		tcgaDataFull					= set()
		for gene in readColsFromXSV(tcgaNetDataPath, colStart = 0, colStop = 1, keepHeader = False, delimiter = ",", openAs = "U", progress = True):
			tcgaDataFull.add(gene[0])
		print("\tStatus: Generating TCGA subgraph")
		# Generate Subgraph
		interactionSubgraph 			= generateUnionSubgraph(allInteractionPairs, tcgaDataFull)[0]



	# 
	# Structure Data
	# 

		print("\tStatus: Retrieving largest TCGA subgraph component")
		fullGraph 						= NX.Graph()
		fullGraph.add_edges_from(interactionSubgraph)
		bigComponent 					= max(NX.connected_component_subgraphs(fullGraph), key = len)
		allInteractionPairs 			= []
		if useRandomComponent == True:
			if bigComponent.number_of_nodes() <= randomComponentSize:
				useRandomComponent 				= False
				print("\tExcept: Random component size larger than big component; proceeding with analysis using big component")
		if useRandomComponent == True:
			targetComponent 				= generateRandomComponent(bigComponent, randomComponentSize)
			print("\t\tRandom Component Metrics - >")
		else:
			targetComponent 				= bigComponent
			print("\t\tBig Component Metrics - >")
		allInteractionPairs 			= targetComponent.edges()
		tcgaGenes 						= set()
		for gene in tcgaDataFull:
			for edge in allInteractionPairs:
				if gene in edge:
					tcgaGenes.add(gene)
		cosmicGenes 					= set()
		for gene in oncoDataFull:
			for edge in allInteractionPairs:
				if gene in edge:
					cosmicGenes.add(gene)
		print("\t\t\tNode Count:\t", targetComponent.number_of_nodes())
		print("\t\t\tEdge Count:\t", targetComponent.number_of_edges())
		print("\t\t\tTCGA Genes:\t", len(tcgaGenes))
		print("\t\t\tCOSMIC Genes:\t", len(cosmicGenes))
		oncoDataFull 					= set()
		# Convert TCGA Genes into String for O(n) Evaluation (in Python, still doing a string search under the hood)
		tcgaGeneStr = "_"	
		for gene in tcgaDataFull:
			tcgaGeneStr 				= tcgaGeneStr + gene + "_"

	# 
	# Generate Factors
	# 

		print("\tStatus: Generating factors")
		interactionFactors  			= []
		for interaction in allInteractionPairs:
			geneA, geneB = "_" + interaction[0] + "_", "_" + interaction[1] + "_"
			if geneA in tcgaGeneStr and geneB in tcgaGeneStr:
				geneA, geneB 			= geneA.replace("_", ""), geneB.replace("_", "")
				phi 	 				= Factor([geneA, geneB], [2, 2], [1, 1, 1, 100])
			elif geneA in tcgaGeneStr:
				geneA, geneB 			= geneA.replace("_", ""), geneB.replace("_", "")
				phi 	 				= Factor([geneA, geneB], [2, 2], [1, 1, 100, 1])
			elif geneB in tcgaGeneStr:
				geneA, geneB 			= geneA.replace("_", ""), geneB.replace("_", "")
				phi 	 				= Factor([geneA, geneB], [2, 2], [1, 100, 1, 1])
			else:
				geneA, geneB 			= geneA.replace("_", ""), geneB.replace("_", "")
				phi 					= Factor([geneA, geneB], [2, 2], [100, 1, 1, 1])
			interactionFactors.append(phi)

	# 
	# Build MRF Model
	# 

		# Generate MRF, Add Edges, Incorporate Factors
		mrf 							= MarkovModel()
		if usePickles == True:
			try:
				with open(objectDirectory + "mrf_model.pickle", "rb") as jar:
					mrf = Pickle.load(jar)
				print("\tStatus: Markov network data structure loaded")
			except:
				print("\tExcept: Unable to locate Markov network data structure...")
				print("\tStatus: Building Markov network")
				mrf.add_edges_from(allInteractionPairs)
				for phi in interactionFactors:
					mrf.add_factors(phi)
				print("\tStatus: Saving Markov network data structure")
				with open(objectDirectory + "mrf_model.pickle", "wb") as jar:
					Pickle.dump(mrf, jar)
		else:
			print("\tStatus: Building Markov network")
			mrf.add_edges_from(allInteractionPairs)
			for phi in interactionFactors:
				mrf.add_factors(phi)
			print("\tStatus: Saving Markov network data structure")
			with open(objectDirectory + "mrf_model.pickle", "wb") as jar:
					Pickle.dump(mrf, jar)
		# Validate that Interaction Genes contain TCGA and all onco-genes
		print("\tStatus: Running Markov network quality check")
		if mrf.check_model():
			print("\tStatus: Markov network passed initial QC")
		else:
			print("\tStatus: Markov network failed initial QC")

	# 
	# MRF Model Inference
	# 
		oncogeneInferenceModel 				= None
		if usePickles == True:
			try:
				with open(objectDirectory + "lbp_inference_model.pickle", "rb") as jar:
					oncogeneInferenceModel 	= Pickle.load(jar)
				print("\tStatus: Belief propagation inference model loaded")
			except:
				print("\tExcept: Unable to locate belief propagation inference model...")
				print("\tStatus: Building belief propagation inference model")
				oncogeneInferenceModel 		= BeliefPropagation(mrf)
				print("\tStatus: Saving belief propagation inference model")
				with open(objectDirectory + "lbp_inference_model.pickle", "wb") as jar:
					Pickle.dump(oncogeneInferenceModel, jar)
				print("\tStatus: Calibrating Belief Propagation Model")
				oncogeneInferenceModel.calibrate()
				print("\tStatus: Saving calibrated belief propagation inference model")
				with open(objectDirectory + "calibrated_lbp_inference_model.pickle", "wb") as jar:
					Pickle.dump(oncogeneInferenceModel, jar)
		else:
			print("\tStatus: Building belief propagation inference model")
			oncogeneInferenceModel 			= BeliefPropagation(mrf)
			print("\tStatus: Saving belief propagation inference model")
			with open(objectDirectory + "lbp_inference_model.pickle", "wb") as jar:
				Pickle.dump(oncogeneInferenceModel, jar)
			print("\tStatus: Calibrating Belief Propagation Model")
			oncogeneInferenceModel.calibrate()
			print("\tStatus: Saving calibrated belief propagation inference model")
			with open(objectDirectory + "calibrated_lbp_inference_model.pickle", "wb") as jar:
				Pickle.dump(oncogeneInferenceModel, jar)

		print("\tStatus: Computing clique beliefs")
		cliqueBeliefs 						= oncogeneInferenceModel.get_clique_beliefs()
		print("\tStatus: Saving clique beliefs")
		with open(objectDirectory + "clique_beliefs.pickle", "wb") as jar:
			Pickle.dump(cliqueBeliefs, jar)
		print("\tStatus: Retrieving Cliques")
		cliques 							= oncogeneInferenceModel.get_cliques()
		print("\tStatus: Saving cliques")
		with open(objectDirectory + "cliques.pickle", "wb") as jar:
			Pickle.dump(cliques, jar)
		print("\tStatus: Retrieving Supset Beliefs")
		supsetBeliefs 						= oncogeneInferenceModel.get_sepset_beliefs()
		print("\tStatus: Saving Supset Beliefs")
		with open(objectDirectory + "supset_beliefs.pickle", "wb") as jar:
			Pickle.dump(supsetBeliefs, jar)

		print("\nStatus: done\n")
	except OSError as error:
		print(error)
		print("Error: unable to run pipeline.py\n")
		Sys.exit(exitCode)
	exitCode = 1
	Sys.exit(exitCode)
else:
	pass