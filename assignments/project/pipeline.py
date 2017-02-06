#! /usr/bin/env python3

"""
Developed With Python Version 3.5

"""

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
import multiprocessing 		as Multi

# Import Third-Party Modules
import numpy 				as Npy
import numpy.linalg 		as LinearAlgebra
import pandas 				as Pandas
import networkx 			as NX

# Import Functions
from collections 			import deque
from pgmpy.models 			import MarkovModel
from pgmpy.factors 			import Factor
from pgmpy.inference 		import VariableElimination
from pgmpy.inference 		import BeliefPropagation
from statistics 			import median
from statistics 			import mean
from multiprocessing 		import Pool

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

# Delete Files (remove symlinks) in Directory
def unlinkDirectory (path):
	"""
	Delete files in directory (remove symlinks)
	@params
		path  		– Required 	: path
	"""
	import os as OS
	return list(map(OS.unlink, (OS.path.join( path, file) for file in OS.listdir(path))))

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
def getUserInput (valid, prompt, hint = "", failed = "Error: Invalid input", attempts = 3):
	"""
	Prompts user for and validates input using regular expression
	@params:
		valid 		- Required 	: regex to validate against (Rgx)
		prompt 		- Required 	: verbose user prompt (Str)
		hint 		- Optional 	: input hint (Str)
		failed 		- Optional 	: failed input (Str)
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
		if attempts > 0:
			print(failed)
			attempts -= 1
			return getUserInput(valid, prompt, hint = hint, failed = failed, attempts = attempts)
		else:
			print(failed)
			print("Exiting")
			exit()

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
	Sys.stdout.write('\r'),
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
		gmlDirectory  					= "output/gml/"
		rcsbDirectory 					= "objects/rc-analysis/sepsetbeliefs/"
		rccbDirectory 					= "objects/rc-analysis/cliquebeliefs/"
		rccDirectory 					= "objects/rc-analysis/cliques/"
		rcrcDirectory 					= "objects/rc-analysis/components/"
		bcsbDirectory 					= "objects/bc-analysis/sepsetbeliefs/"
		bccbDirectory 					= "objects/bc-analysis/cliquebeliefs/"
		bccDirectory 					= "objects/bc-analysis/cliques/"
		bcmDirectory 					= "objects/bc-analysis/models/"
		usePickles 						= False
		useRandomComponents 			= False
		randomComponentCount 			= 100
		randomComponentSize 			= 100
		interactionDataResponse 		= "0"
		maxCores 						= Multi.cpu_count()

# 
# Command Line Setup 
# 

		nCores 							= getUserInput(valid=r"([1-9][0-9]{0,1}|exit)", prompt="\tPrompt: Set number of cores to use for computations", hint = "\t\tHint: enter a number between 1 and %d\n\t\tSelection: " % (maxCores), failed = "\t\tError: Invalid input, %d cores are available on your system" % (maxCores))
		pValue 							= float(getUserInput(valid=r"([0]\.[0-9]{1,3}|exit)", prompt="\tPrompt: Set rank p-Value...", hint = "\t\tHint: enter a number between 0.001 and 0.999\n\t\tSelection: ", failed = "\t\tError: Invalid input"))
		rankCutoff 					 	= int(getUserInput(valid=r"([1-9][0-9]{1,2}|exit)", prompt="\tPrompt: Number of top ranked novel candidates to return?", hint = "\t\tHint: enter a number between 1 and 999. Top ranks of 25 or 50 are typical.\n\t\tSelection: ", failed = "\t\tError: Invalid input"))
		usePicklesResponse 				= getUserInput(valid=r"([yYnN]{1}|exit)", prompt="\tPrompt: Reproduce last run using pickled data?", hint = "\t\tHint: enter 'y', 'n', or 'exit' to exit\n\t\tSelection: ", failed = "\t\tError: Invalid input")
		if usePicklesResponse.lower() == "y":
			usePickles 					= True
			interactionDataResponse 	= getUserInput(valid=r"([0-1]{1}|exit)", prompt="\tPrompt: Use Reactome PSI-MItab [0], FI network [1]?", hint = "\t\tHint: enter 0, 1, or 'exit' to exit\n\t\tSelection: ", failed = "\t\tError: Invalid input")
			if interactionDataResponse == "exit":
				exitCode 				= 1
				Sys.exit(exitCode)
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

		if usePickles == False:
			print("\n\tStatus: Removing Data From Previous Run...")
			if useRandomComponents == True:
				unlinkDirectory(rcsbDirectory)
				unlinkDirectory(rccbDirectory)
				unlinkDirectory(rccDirectory)
				unlinkDirectory(rcrcDirectory)
			else:
				unlinkDirectory(bcsbDirectory)
				unlinkDirectory(bccbDirectory)
				unlinkDirectory(bccDirectory)
				unlinkDirectory(bcmDirectory)
		print("\tStatus: Loading Data\n")
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
			interactionDataPath = fiDataPath
			for row in readColsFromXSV(interactionDataPath, colStart = 0, colStop = 3, keepHeader = False, delimiter = "\t", progress = True):
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
		allOncogenes = set(tcgaDataFull) | set(oncoDataFull)
		# 619 oncogenes
		# Generate Subgraph
		# Critical error! 
		# interactionSubgraph 			= generateUnionSubgraph(allInteractionPairs, allOncogenes)[0]

	# 
	# Structure Data & Perform Basic EDA
	# 

		print("\tStatus: Generating FI subgraph & Saving to GML")
		fullGraph 						= NX.Graph()
		# Updated to use allInteractionPairs rather than the Subgraph
		fullGraph.add_edges_from(allInteractionPairs)
		NX.write_gml(fullGraph, gmlDirectory + "fi_subgraph.gml")
		print("\tStatus: Retrieving big component & Saving to GML")
		bigComponent 					= max(NX.connected_component_subgraphs(fullGraph), key = len)
		NX.write_gml(bigComponent, gmlDirectory + "fi_subgraph_big_component.gml")
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
			alltcgaGenes 						= set()
			# Get all TCGA and COSMIC genes for downstream analysis
			for gene in tcgaDataFull:
				for edge in bigComponent.edges():
					if gene in edge:
						alltcgaGenes.add(gene)
			allCosmicGenes 						= set()
			for gene in oncoDataFull:
				for edge in bigComponent.edges():
					if gene in edge:
						allCosmicGenes.add(gene)
			# Get list of all oncogenes in big component for downstream analysis
			allBCOncoGenes = alltcgaGenes | allCosmicGenes
			# Get list of all oncogenes (COSMIC & TCGA) that appear in big component
			if usePickles == True:
				cliqueBeliefFiles 				= OS.listdir(OS.getcwd() + "/" + rccbDirectory)
				cliqueFiles 					= OS.listdir(OS.getcwd() + "/" + rccDirectory)
				sepsetBeliefFiles 				= OS.listdir(OS.getcwd() + "/" + rcsbDirectory)
				if len(cliqueBeliefFiles) < 1 or len(cliqueFiles) < 1 or len(sepsetBeliefFiles) < 1:
					usePickles = False
					print("\tExcept: Unable to locate files from previous random component generation. Proceeding with random component generation. ")
			if usePickles == False:
				print("\tStatus: Generating Random Components\n")
				randomComponents 				= []
				nodeCounts 						= []
				edgeCounts 						= []
				tcgaGeneCounts 					= []
				cosmicGeneCounts 				= []
				i 								= 0
				# Parallelize Here
				while i < randomComponentCount:
					component 						= NX.Graph(generateRandomComponent(bigComponent.edges(), randomComponentSize))
					componentEdges 					= component.edges()
					# Get Count of TCGA genes that appear in component (for RC metrics)
					tcgaGenes 						= set()
					for gene in tcgaDataFull:
						for edge in componentEdges:
							if gene in edge:
								tcgaGenes.add(gene)
					# Verify that component has at least one TCGA gene in it
					if len(tcgaGenes) > 0:
						# Get Count of COSMIC genes that appear in component (for RC metrics)
						cosmicGenes 					= set()
						for gene in oncoDataFull:
							for edge in componentEdges:
								if gene in edge:
									cosmicGenes.add(gene)
						# Append component metrics to metric analysis sets
						nodeCounts.append(component.number_of_nodes())
						edgeCounts.append(component.number_of_edges())
						tcgaGeneCounts.append(len(tcgaGenes))
						cosmicGeneCounts.append(len(cosmicGenes))
						# Append component to RC set for training
						randomComponents.append(component)
						i 							   += 1
						printProgress(i, randomComponentCount, prefix = "\t\tRandom Components:" , suffix = "Generated", decimals = 2, barLength = 50)
				with open(rcrcDirectory + str(randomComponentCount) + "_random_components.pickle", "wb") as jar:
					Pickle.dump(randomComponents, jar)
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
				completedModels		= 0
				i 					= 0
				tcgaGeneStr = "~"	
				for gene in tcgaDataFull:
					tcgaGeneStr 				= tcgaGeneStr + gene + "~"
				for component in randomComponents:
					interactionFactors  			= []
					for interaction in component.edges():
						geneA, geneB = "~" + interaction[0] + "~", "~" + interaction[1] + "~"
						if geneA in tcgaGeneStr and geneB in tcgaGeneStr:
							geneA, geneB 			= geneA.replace("~", ""), geneB.replace("~", "")
							phi 	 				= Factor([geneA, geneB], [2, 2], [1, 20, 20, 100])
						elif geneA in tcgaGeneStr:
							geneA, geneB 			= geneA.replace("~", ""), geneB.replace("~", "")
							phi 	 				= Factor([geneA, geneB], [2, 2], [1, 20, 60, 1])
						elif geneB in tcgaGeneStr:
							geneA, geneB 			= geneA.replace("~", ""), geneB.replace("~", "")
							phi 	 				= Factor([geneA, geneB], [2, 2], [1, 60, 20, 1])
						else:
							geneA, geneB 			= geneA.replace("~", ""), geneB.replace("~", "")
							phi 					= Factor([geneA, geneB], [2, 2], [40, 1, 1, 1])
						interactionFactors.append(phi)
					mrf 							= MarkovModel()
					mrf.add_edges_from(component.edges())
					for phi in interactionFactors:
						mrf.add_factors(phi)
					if mrf.check_model():
						validModels.append(mrf)
				# Clear Random Components for Garbage Collection
				randomComponents = None
				print("\tStatus: " + str(round(((float(len(validModels)) / randomComponentCount) * 100), 2)) + "%" + " of Models Pass Quality Check")
				print("\tStatus: Generating & Saving Inferences From Valid Models\n")
				i = 0
				for model in validModels:
					try:
						inferenceModel 				= BeliefPropagation(model)
						inferenceModel.calibrate()
						cliqueBeliefs 				= inferenceModel.get_clique_beliefs()
						cliques 					= inferenceModel.get_cliques()
						sepsetBeliefs 				= inferenceModel.get_sepset_beliefs()
						with open(rcsbDirectory + str(completedModels) + "_sepset_beliefs.pickle", "wb") as jar:
							Pickle.dump(sepsetBeliefs, jar)			
						with open(rccbDirectory + str(completedModels) + "_clique_beliefs.pickle", "wb") as jar:
							Pickle.dump(cliqueBeliefs, jar)
						with open(rccDirectory + str(completedModels) + "_cliques.pickle", "wb") as jar:
							Pickle.dump(cliques, jar)
						completedModels 		   += 1
						i 						   += 1
						printProgress(i, randomComponentCount, prefix = "\t\tComputing:" , suffix = "Success: (" + str(round(((float(completedModels) / i) * 100), 2)) + "%)" , decimals = 2, barLength = 50)
					except:
						i 						   += 1
						printProgress(i, randomComponentCount, prefix = "\t\tComputing:" , suffix = "Success: (" + str(round(((float(completedModels) / i) * 100), 2)) + "%)", decimals = 2, barLength = 50)
						continue
				cliqueBeliefFiles 						= OS.listdir(OS.getcwd() + "/" + rccbDirectory)
				cliqueFiles 							= OS.listdir(OS.getcwd() + "/" + rccDirectory)
				sepsetBeliefFiles 						= OS.listdir(OS.getcwd() + "/" + rcsbDirectory)
				print("\n")

		# 
		# Computing highest scoring genes
		# 
			
			print("\tStatus: Computing Gene Scores\n")
			cbDict = {}
			sbDict = {}
			i = 0
			for cliqueBeliefFile in cliqueBeliefFiles:
				uniqueCBs = set()
				with open(rccbDirectory + cliqueBeliefFile, "rb") as jar:
					cb = Pickle.load(jar)
					for v in cb:
						uniqueCBs.add(frozenset(v))
					for s in uniqueCBs:
						for g in s:
							if g in cbDict:
								cbDict[g] += 1
							else:
								cbDict[g] = 1
				i += 1
				printProgress(i, len(cliqueBeliefFiles), prefix = "\t\tComputing Clique Belief Gene Scores:" , suffix = "of Analysis Complete", decimals = 2, barLength = 50)
			i = 0
			cbList = sorted(cbDict.items(), key = lambda x: x[1], reverse = True)
			# Normalize Scores
			cbListMin 	= cbList[-1][1]
			cbListDenom = cbList[0][1] - cbListMin
			cbList 		= [(value[0], round(((value[1] - cbListMin) / cbListDenom), 4)) for value in cbList]
			print("\n")
			for sepsetBeliefFile in sepsetBeliefFiles:
				uniqueSBs = set()
				with open(rcsbDirectory + sepsetBeliefFile, "rb") as jar:
					sb = Pickle.load(jar)
					for s in sb:
						for v in s:
							uniqueSBs.add(frozenset(v))
					for s in uniqueSBs:
						for g in s:
							if g in sbDict:
								sbDict[g] += 1
							else:
								sbDict[g] = 1
				i += 1
				printProgress(i, len(sepsetBeliefFiles), prefix = "\t\tComputing Separation Set Belief Gene Scores:" , suffix = "of Analysis Complete", decimals = 2, barLength = 50)
			sbList 		= sorted(sbDict.items(), key = lambda x: x[1], reverse = True)
			# Normalize Scores
			sbListMin 	= sbList[-1][1]
			sbListDenom = sbList[0][1] - sbListMin
			sbList 		= [(value[0], round(((value[1] - sbListMin) / sbListDenom), 4)) for value in sbList]
			networkAnalysis = getUserInput(valid=r"([yYnN]{1})", prompt="\n\n\tPrompt: Generate network analysis (in addition to inferences)?", hint = "\t\tHint: enter 'y' or 'n'\n\t\tSelection: ", failed = "\t\tError: Invalid input")
			print("\n\tStatus: Generating Analysis\n")
			# Cross Reference CB Genes With Known
			topCBcosmic = []
			# Are CB Genes in Cosmic?
			for gene in cbList:
				status = False
				for onco in allCosmicGenes:
					if gene[0] == onco:
						status = True
						break
				gene += (status, )
				topCBcosmic.append(gene)
			topCBAll = []
			rank = 1
			# Are CB Genes in TCGA?
			for gene in topCBcosmic:
				status = False
				for onco in alltcgaGenes:
					if gene[0] == onco:
						status = True
						break
				gene += (status, )
				gene += (rank, )
				rank += 1
				topCBAll.append(gene)
			topCBCutoff = topCBAll[0: int(pValue * len(cbList))]
			# Cross Reference SB Genes With Known
			topSB = sbList[0: int(pValue * len(sbList))]
			topSBcosmic = []
			# Are SB Genes in Cosmic?
			for gene in sbList:
				status = False
				for onco in allCosmicGenes:
					if gene[0] == onco:
						status = True
						break
				gene += (status, )
				topSBcosmic.append(gene)
			topSBAll = []
			rank = 1
			# Are SB Genes in TCGA?
			for gene in topSBcosmic:
				status = False
				for onco in alltcgaGenes:
					if gene[0] == onco:
						status = True
						break
				gene += (status, )
				gene += (rank, )
				rank += 1
				topSBAll.append(gene)
			topSBCutoff 				= topSBAll[0: int(pValue * len(sbList))]
			eigenvectorCentralities 	= NX.eigenvector_centrality(bigComponent)

			print("----------------------------------------------------------------------------------------------------------------------------")
			print("--------------------------------------------------------- ANALYSIS ---------------------------------------------------------")
			print("----------------------------------------------------------------------------------------------------------------------------")
			print("\nTop Ranked Genes (Clique Belief)\n")
			for gene in topCBCutoff:
				print("Rank:", gene[4], "\tGene:", gene[0], "\tScore:", gene[1], "\tCOSMIC:", gene[2], "\tTCGA:", gene[3], "\tDegree:", NX.degree(bigComponent, gene[0]), "\tEigenvector:", round(eigenvectorCentralities[gene[0]], 8))
			print("\nAnalysis Statistics:")
			print("p-Value\t\t", "=", pValue)
			print("n\t\t", "=", len(cbList))
			print("----------------------------------------------------------------------------------------------------------------------------")
			print("\nTop Ranked Genes (Sepset Belief)\n")
			for gene in topSBCutoff:
				print("Rank:", gene[4], "\tGene:", gene[0], "\tScore:", gene[1], "\tCOSMIC:", gene[2], "\tTCGA:", gene[3], "\tDegree:", NX.degree(bigComponent, gene[0]), "\tEigenvector:", round(eigenvectorCentralities[gene[0]], 8))
			print("\nAnalysis Statistics:")
			print("p-Value\t\t", "=", pValue)
			print("n\t\t", "=", len(sbList))
			print("----------------------------------------------------------------------------------------------------------------------------")
			print("\nTop", rankCutoff, "Novel Candidate Genes (Clique Belief)\n")
			topRankedCB = []
			count = 0
			for gene in topCBAll:
				if count < rankCutoff:
					if gene[3] is False:
						topRankedCB.append(gene)
						print("Rank:", gene[4], "\tGene:", gene[0], "\tScore:", gene[1], "\tCOSMIC:", gene[2], "\tTCGA:", gene[3], "\tDegree:", NX.degree(bigComponent, gene[0]), "\tEigenvector:", round(eigenvectorCentralities[gene[0]], 8))
						count += 1
			print("----------------------------------------------------------------------------------------------------------------------------")
			print("\nTop", rankCutoff, "Novel Candidate Genes (Separation Set Belief)\n")
			topRankedSB = []
			count = 0
			for gene in topSBAll:
				if count < rankCutoff:
					if gene[3] is False:
						topRankedSB.append(gene)
						print("Rank:", gene[4], "\tGene:", gene[0], "\tScore:", gene[1], "\tCOSMIC:", gene[2], "\tTCGA:", gene[3], "\tDegree:", NX.degree(bigComponent, gene[0]), "\tEigenvector:", round(eigenvectorCentralities[gene[0]], 8))
						count += 1
			print("----------------------------------------------------------------------------------------------------------------------------")
			print("\nNovel Classification Statistics")
			cbNovelValid 	= (float(len([True for gene in topRankedCB if gene[2] is True])) / rankCutoff) * 100
			cbNovel 		= (float(len([True for gene in topRankedCB if gene[2] is False])) / rankCutoff) * 100
			sbNovelValid 	= (float(len([True for gene in topRankedSB if gene[2] is True])) / rankCutoff) * 100
			sbNovel 		= (float(len([True for gene in topRankedSB if gene[2] is False])) / rankCutoff) * 100
			print("\n\tValidated Accurate Novel Classifications")
			print("\t\tClique Belief Classifier\t\t:", str(round(cbNovelValid, 4)) + "%")
			print("\t\tSeparation Set Belief Classifier\t:", str(round(sbNovelValid, 4)) + "%")
			print("\n\tNovel Classifications*")
			print("\t\tClique Belief Classifier\t\t:", str(round(cbNovel, 4)) + "%")
			print("\t\tSeparation Set Belief Classifier\t:", str(round(sbNovel, 4)) + "%")
			print("\n\t*Novel classifications require gene lookup for validation")
			print("----------------------------------------------------------------------------------------------------------------------------")
			if networkAnalysis.lower() == "y":
				print("Generating Network Analysis")
				bcMSP 						= NX.average_shortest_path_length(bigComponent)
				bcDiameter 					= NX.diameter(bigComponent)
				print("\nNetwork (Big Component) Statistics")
				print("Nodes\t\t\t:", bigComponent.number_of_nodes())
				print("Edges\t\t\t:", bigComponent.number_of_edges())
				print("Mean Degree\t\t:", round(float(bigComponent.number_of_edges()) / bigComponent.number_of_nodes(), 8))
				print("Mean Shortest Path\t:", round(bcMSP, 8))
				print("Network Diameter\t:", bcDiameter)
				print("TCGA Nodes\t\t:", len(alltcgaGenes))
				print("COSMIC Nodes\t\t:", len(allCosmicGenes))
				print("All Oncogenes\t\t:", len(allBCOncoGenes))
				print("Data Source\t\t:", interactionDataPath)
				print("----------------------------------------------------------------------------------------------------------------------------")

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
			print("\t\t\tAll Oncogenes: \t", len(allBCOncoGenes))
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
			print("\tStatus: Retrieving Sepset Beliefs")
			sepsetBeliefs 						= oncogeneInferenceModel.get_sepset_beliefs()
			print("\tStatus: Saving Sepset Beliefs")
			with open(bcsbDirectory + "bc_sepset_beliefs.pickle", "wb") as jar:
				Pickle.dump(sepsetBeliefs, jar)


		# Clear Generated Data
		if useRandomComponents == True:
			saveData = getUserInput(valid=r"([yYnN]{1})", prompt="\tPrompt: Save generated data (enables reproduction from pickled data)?", hint = "\t\tHint: enter 'y' or 'n'\n\t\tSelection: ", failed = "\t\tError: Invalid input")
			if saveData.lower == "n":
				unlinkDirectory(rcsbDirectory)
				unlinkDirectory(rccbDirectory)
				unlinkDirectory(rccDirectory)
				unlinkDirectory(rcrcDirectory)
		else:
			saveData = getUserInput(valid=r"([yYnN]{1})", prompt="\tPrompt: Save generated data (enables reproduction from pickled data)?", hint = "\t\tHint: enter 'y' or 'n'\n\t\tSelection: ", failed = "\t\tError: Invalid input")
			if saveData.lower == "n":
				unlinkDirectory(bcsbDirectory)
				unlinkDirectory(bccbDirectory)
				unlinkDirectory(bccDirectory)
				unlinkDirectory(bcmDirectory)
		print("\nStatus: done\n")
	except OSError as error:
		print(error)
		print("Error: unable to run pipeline.py\n")
		Sys.exit(exitCode)
	exitCode = 1
	Sys.exit(exitCode)
else:
	pass