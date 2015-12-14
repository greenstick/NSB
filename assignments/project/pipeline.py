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
from statistics 			import median
from statistics 			import mean

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
def getUserInput (valid, prompt, hint = "", failed =Error: Invalid input"):
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

def bfs (graph, start):
	"""
	Discover edges in a network using a Breadth-First Search
	@params:
		graph 		- Required 	: (Dictionary)
		start 		- Optional 	: (Str, Int)
	"""
	queue, enqueued = deque([(None, start)]), set([start])
	while queue:
		parent, node = queue.popleft()
		yield parent, node
		if node in graph:
			v = set(graph[node]) - enqueued
			enqueued |= v
			queue.extend([(node, child) for child in v])

def generateRandomComponent (graph, nodeCount = 100):
	""" 
	Generate Random Graphs
	@params:
		graph 				- Required 	: (Dictionary, List)
		nodeCount 			- Optional 	: (Int)
	"""
	if type(graph) is list:
		graphEdges 	= graph
		graph 		= convertEdgelistToGraph(graph)
	else:
		graphEdges 	= convertGraphToEdgelist(graph)
	nodes 		= set()
	start 		= Random.sample(graph.keys(), 1)[0]
	edges 		= bfs(graph, start)
	for edge in edges:
		if len(nodes) <= nodeCount:
			nodes.add(edge[0]), nodes.add(edge[1])
		else:
			break
	component = list(generateIntersectionSubgraph(graphEdges, nodes)[0])
	return component

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
# Network Data Structure Conversion
# 

def convertEdgelistToGraph (edgeList):
	"""
	Generate Graph as Dictionary From List of Edges
	@params:
		edgeList 	- Required 	: (Set, List)
	"""
	graph = {}
	for start, end in edgeList: 
		if start in graph:
			graph[start].append(end)
		else:
			graph[start] = [end]
		if end in graph:
			graph[end].append(start)
		else:
			graph[end] = [start]
	return graph

def convertGraphToEdgelist (graph):
	"""
	Generate List of Edges From Graph as Dictionary
	@params:
		edgeList 	- Required 	: (Set, List)
	"""
	return list(set((node, link) for link in list(node for node in graph.keys())))

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
		fiDataPath 					 	= "data/FIsInGene_121514_with_annotations.txt"
		oncoNetDataPath 				= "data/Census_allTue Dec  8 23_01_10 2015.csv"
		tcgaNetDataPath 				= "data/tcga-genes.txt"
		outputDirectory  				= "output/"
		rcsbDirectory 					= "objects/rc-analysis/sepsetbeliefs/"
		rccbDirectory 					= "objects/rc-analysis/cliquebeliefs/"
		rccDirectory 					= "objects/rc-analysis/cliques/"
		bcsbDirectory 					= "objects/bc-analysis/sepsetbeliefs/"
		bccbDirectory 					= "objects/bc-analysis/cliquebeliefs/"
		bccDirectory 					= "objects/bc-analysis/cliques/"
		bcmDirectory 					= "objects/bc-analysis/models/"
		usePickles 						= False
		useRandomComponents 			= False
		randomComponentCount 			= 100
		randomComponentSize 			= 100
		interactionDataResponse 		= "0"
		usePicklesResponse 				= getUserInput(valid=r"([yYnN]{1}|exit)", prompt="\tPrompt: Reproduce last run using pickled data?", hint = "\t\tHint: enter 'y', 'n', or 'exit' to exit\n\t\tSelection: ", failed = "\t\tError: Invalid input")
		if usePicklesResponse.lower() == "y":
			usePickles 					= True
			useRandomComponentsResponse 	= getUserInput(valid=r"([yYnN]{1}|exit)", prompt="\tPrompt: Generate bootstrap analysis (No indicates a big component analysis will be generated)?", hint = "\t\tHint: enter 'y', 'n', or 'exit' to exit\n\t\tSelection: ", failed = "\t\tError: Invalid input")
			if useRandomComponentsResponse.lower() == "y":
				useRandomComponents 		= True
			elif useRandomComponentsResponse == "exit":
				exitCode 				= 1
				Sys.exit(exitCode)
		elif usePicklesResponse.lower() == "n":
			interactionDataResponse 	= getUserInput(valid=r"([0-1]{1}|exit)", prompt="\tPrompt: Use Reactome PSI-MItab [0], FI network [1]?", hint = "\t\tHint: enter 0, 1, or 'exit' to exit\n\t\tSelection: ", failed = "\t\tError: Invalid input")
			if interactionDataResponse == "exit":
				exitCode 				= 1
				Sys.exit(exitCode)
			useRandomComponentsResponse 	= getUserInput(valid=r"([yYnN]{1}|exit)", prompt="\tPrompt: Perform bootstrap sample analysis (no indicates that the network big component will be used; this will take a long time)?", hint = "\t\tHint: enter 'y', 'n', or 'exit' to exit\n\t\tSelection: ", failed = "\t\tError: Invalid input")
			if useRandomComponentsResponse.lower() == "y":
				useRandomComponents 		= True
				print("\tNotification: Lots of components (say, 10,000) with a small sample size (say, 20) tends to offer a good ratio of power to performance.")
				randomComponentCount 	= int(getUserInput(valid=r"[1-9][0-9]{0,5}", prompt="\tPrompt: How many random components would you like to generate?", hint = "\t\tHint: enter integer (n) between [1 - 99999]\n\t\tSelection: ", failed = "\t\tError: Invalid input"))
				randomComponentSize 	= int(getUserInput(valid=r"[1-9][0-9]{0,4}", prompt="\tPrompt: How many nodes per random component should be generated?", hint = "\t\tHint: enter integer (n) between [1 - 9999]\n\t\tSelection: ", failed = "\t\tError: Invalid input"))
			elif useRandomComponentsResponse == "exit":
				exitCode 				= 1
				Sys.exit(exitCode)
		else:
			exitCode 					= 1
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
		if interactionDataResponse == "0":
			# Import Interactions Node Data
			for row in readColsFromXSV(interactionDataPath, colStart = 2, colStop = 4, keepHeader = False, delimiter = "\t", progress = True):
				# Convert Row to String (join), Extract Genes Using Rgx Substitution, Split Resulting String Into a List, Coerce into a Tuple
				interaction 				= tuple(Rgx.sub(r"uniprotkb\:[a-zA-Z0-9]{1,8}\_", "", "".join(row)).split("(shortlabel)")[0:2])
				# Remove self-edges, singletons, interactions with empty None values or a length less than or equal to one
				if len(interaction) == 2:
					if interaction[0] != interaction[1]:
						if (interaction[0] != None and len(interaction[0]) > 1) and (interaction[1] != None and len(interaction[1]) > 1):
							allInteractionPairs.append(interaction)
		else:
			for row in readColsFromXSV(fiDataPath, colStart = 0, colStop = 3, keepHeader = False, delimiter = "\t", progress = True):
				interaction 					= (row[0], row[1])
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
	# Structure Data & Perform Basic EDA
	# 

		print("\tStatus: Retrieving largest TCGA subgraph component")
		fullGraph 						= NX.Graph()
		fullGraph.add_edges_from(interactionSubgraph)
		bigComponent 					= max(NX.connected_component_subgraphs(fullGraph), key = len)
		allInteractionPairs 			= []
		if useRandomComponents == True:
			if bigComponent.number_of_nodes() <= randomComponentSize:
				useRandomComponents 				= False
				print("\tExcept: Random component size larger than big component; proceeding with analysis using big component")

	# 
	# Random Components Analysis
	# 

		if useRandomComponents == True:
			cliqueBeliefFiles 					= []
			cliqueFiles 						= []
			sepsetBeliefFiles 					= []
			if usePickles == True:
				cliqueBeliefFiles 				= OS.listdir(OS.getcwd() + "/" + rccbDirectory)
				cliqueFiles 					= OS.listdir(OS.getcwd() + "/" + rccDirectory)
				sepsetBeliefFiles 				= OS.listdir(OS.getcwd() + "/" + rcsbDirectory)
				if len(cliqueBeliefFiles) < 1 or len(cliqueFiles) < 1 or len(sepsetBeliefFiles) < 1:
					usePickles = False
					print("\tExcept: Unable to locate files from previouse random component generation. Proceeding with random component generation. ")
			if usePickles == False:
				randomComponents 				= []
				nodeCounts 						= []
				edgeCounts 						= []
				tcgaGeneCounts 					= []
				cosmicGeneCounts 				= []
				i 								= 0
				for c in range(randomComponentCount):
					component 						= NX.Graph(generateRandomComponent(bigComponent.edges(), randomComponentSize))
					nodeCounts.append(component.number_of_nodes())
					edgeCounts.append(component.number_of_edges())
					tcgaGenes 						= set()
					for gene in tcgaDataFull:
						for edge in component.edges():
							if gene in edge:
								tcgaGenes.add(gene)
					tcgaGeneCounts.append(len(tcgaGenes))
					cosmicGenes 					= set()
					for gene in oncoDataFull:
						for edge in component.edges():
							if gene in edge:
								cosmicGenes.add(gene)
					cosmicGeneCounts.append(len(cosmicGenes))
					randomComponents.append(component)
					i 							   += 1
					printProgress(i, randomComponentCount, prefix = "\tGenerating Components:" , suffix = "Generated", decimals = 2, barLength = 50)
				print("\n\n\t\tRandom Component Metrics - >")
				print("\n\t\t\tMin Node Count:\t\t\t", min(nodeCounts))
				print("\t\t\tMedian Node Count:\t\t", median(nodeCounts))
				print("\t\t\tMax Node Count:\t\t\t", max(nodeCounts))
				print("\t\t\tMean Node Count:\t\t", mean(nodeCounts))
				print("\n\t\t\tMin Edge Count:\t\t\t", min(edgeCounts))
				print("\t\t\tMedian Edge Count:\t\t", median(edgeCounts))
				print("\t\t\tMax Edge Count:\t\t\t", max(edgeCounts))
				print("\t\t\tMean Edge Count:\t\t", mean(edgeCounts))
				print("\n\t\t\tMin TCGA Gene Count:\t\t", min(tcgaGeneCounts))
				print("\t\t\tMedian TCGA Gene Count:\t\t", median(tcgaGeneCounts))
				print("\t\t\tMax TCGA Gene Count:\t\t", max(tcgaGeneCounts))
				print("\t\t\tMean TCGA Gene Count:\t\t", mean(tcgaGeneCounts))
				print("\n\t\t\tMin COSMIC Gene Count:\t\t", min(cosmicGeneCounts))
				print("\t\t\tMedian COSMIC Gene Count:\t", median(cosmicGeneCounts))
				print("\t\t\tMax COSMIC Gene Count:\t\t", max(cosmicGeneCounts))
				print("\t\t\tMean COSMIC Gene Count:\t\t", mean(cosmicGeneCounts))
				
				print("\n\tStatus: Generating Factors & Building Markov Networks")
				allFactors 			= []
				validModels 		= []
				validModelsCount 	= 0
				completedModels		= 0
				i 					= 0
				tcgaGeneStr = "_"	
				for gene in tcgaDataFull:
					tcgaGeneStr 				= tcgaGeneStr + gene + "_"

				for component in randomComponents:
					interactionFactors  			= []
					for interaction in component.edges():
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
					mrf 							= MarkovModel()
					mrf.add_edges_from(component.edges())
					for phi in interactionFactors:
						mrf.add_factors(phi)
					if mrf.check_model():
						validModels.append(mrf)
						validModelsCount 		   += 1
				print("\tStatus: " + str(round(((float(validModelsCount) / randomComponentCount) * 100), 2)) + "%" + " of Models Pass Quality Check")
				print("\tStatus: Generating & Saving Inferences From Valid Models")
				i = 0
				for model in validModels:
					try:
						inferenceModel 				= BeliefPropagation(model)
						inferenceModel.calibrate()
						cliqueBeliefs 				= inferenceModel.get_clique_beliefs()
						cliques 					= inferenceModel.get_cliques()
						sepsetBeliefs 				= inferenceModel.get_sepset_beliefs()
						completedModels 		   += 1
						i 						   += 1
						with open(rcsbDirectory + str(completedModels) + "_supset_beliefs.pickle", "wb") as jar:
							Pickle.dump(sepsetBeliefs, jar)			
						with open(rccbDirectory + str(completedModels) + "_clique_beliefs.pickle", "wb") as jar:
							Pickle.dump(cliqueBeliefs, jar)
						with open(rccDirectory + str(completedModels) + "_cliques.pickle", "wb") as jar:
							Pickle.dump(cliques, jar)
						printProgress(i, randomComponentCount, prefix = "\tComputing:" , suffix = "Success: (" + str(round(((float(completedModels) / i) * 100), 2)) + "%)" , decimals = 0, barLength = 50)
					except:
						i 						   += 1
						printProgress(i, randomComponentCount, prefix = "\tComputing:" , suffix = "Success: (" + str(round(((float(completedModels) / i) * 100), 2)) + "%)", decimals = 0, barLength = 50)
						continue
				cliqueBeliefFiles 						= OS.listdir(OS.getcwd() + "/" + rccbDirectory)
				cliqueFiles 							= OS.listdir(OS.getcwd() + "/" + rccDirectory)
				sepsetBeliefFiles 						= OS.listdir(OS.getcwd() + "/" + rcsbDirectory)
			for filePath in cliqueBeliefFiles:


	# 
	# Big Component Analysis
	# 

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
					with open(bcmDirectory + "bc_mrf_model.pickle", "rb") as jar:
						mrf = Pickle.load(jar)
					print("\tStatus: Markov network data structure loaded")
				except:
					print("\tExcept: Unable to locate Markov network data structure...")
					print("\tStatus: Building Markov network")
					mrf.add_edges_from(allInteractionPairs)
					for phi in interactionFactors:
						mrf.add_factors(phi)
					print("\tStatus: Saving Markov network data structure")
					with open(bcmDirectory + "bc_mrf_model.pickle", "wb") as jar:
						Pickle.dump(mrf, jar)
			else:
				print("\tStatus: Building Markov network")
				mrf.add_edges_from(allInteractionPairs)
				for phi in interactionFactors:
					mrf.add_factors(phi)
				print("\tStatus: Saving Markov network data structure")
				with open(bcmDirectory + "bc_mrf_model.pickle", "wb") as jar:
						Pickle.dump(mrf, jar)
			# Validate that Interaction Genes contain TCGA and all onco-genes
			print("\tStatus: Running Markov network quality check")
			if mrf.check_model():
				print("\tStatus: Markov network passed initial QC")
			else:
				print("\tStatus: Markov network failed initial QC")

		# 
		# BC MRF Model Inference
		# 
			oncogeneInferenceModel 				= None
			if usePickles == True:
				try:
					with open(bcmDirectory + "bc_lbp_inference_model.pickle", "rb") as jar:
						oncogeneInferenceModel 	= Pickle.load(jar)
					print("\tStatus: Belief propagation inference model loaded")
				except:
					print("\tExcept: Unable to locate belief propagation inference model...")
					print("\tStatus: Building belief propagation inference model")
					oncogeneInferenceModel 		= BeliefPropagation(mrf)
					print("\tStatus: Saving belief propagation inference model")
					with open(bcmDirectory + "bc_lbp_inference_model.pickle", "wb") as jar:
						Pickle.dump(oncogeneInferenceModel, jar)
					print("\tStatus: Calibrating Belief Propagation Model")
					oncogeneInferenceModel.calibrate()
					print("\tStatus: Saving calibrated belief propagation inference model")
					with open(bcmDirectory + "bc_calibrated_lbp_inference_model.pickle", "wb") as jar:
						Pickle.dump(oncogeneInferenceModel, jar)
			else:
				print("\tStatus: Building belief propagation inference model")
				oncogeneInferenceModel 			= BeliefPropagation(mrf)
				print("\tStatus: Saving belief propagation inference model")
				with open(objectDirectory + "bc_lbp_inference_model.pickle", "wb") as jar:
					Pickle.dump(oncogeneInferenceModel, jar)
				print("\tStatus: Calibrating Belief Propagation Model")
				oncogeneInferenceModel.calibrate()
				print("\tStatus: Saving calibrated belief propagation inference model")
				with open(objectDirectory + "calibrated_lbp_inference_model.pickle", "wb") as jar:
					Pickle.dump(oncogeneInferenceModel, jar)

			print("\tStatus: Computing clique beliefs")
			cliqueBeliefs 						= oncogeneInferenceModel.get_clique_beliefs()
			print("\tStatus: Saving clique beliefs")
			with open(bccbDirectory + "bc_clique_beliefs.pickle", "wb") as jar:
				Pickle.dump(cliqueBeliefs, jar)
			print("\tStatus: Retrieving Cliques")
			cliques 							= oncogeneInferenceModel.get_cliques()
			print("\tStatus: Saving cliques")
			with open(bccDirectory + "bc_cliques.pickle", "wb") as jar:
				Pickle.dump(cliques, jar)
			print("\tStatus: Retrieving Supset Beliefs")
			supsetBeliefs 						= oncogeneInferenceModel.get_sepset_beliefs()
			print("\tStatus: Saving Supset Beliefs")
			with open(bcsbDirectory + "bc_supset_beliefs.pickle", "wb") as jar:
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