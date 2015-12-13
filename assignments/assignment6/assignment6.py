#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Developed With Python Version 2.7.8

"""
# Maybe The Max Eigenvalue Being Low is Caused by The FutureWarnings I'm Getting...but Probably Not.
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

# Import Dependencies
import matplotlib
import matplotlib.pyplot 	as Plot
import random 				as Random
import csv 					as CSV
import os 					as OS
import sys 					as Sys
import time 				as Time
import os.path 				as Path
import json 				as JSON
import re 					as Rgx
import errno 				as Error
import networkx 			as NX
import numpy 				as Npy
import numpy.linalg 		as LinearAlgebra
import pickle 				as Pickle
from collections 			import deque

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

# Read Delimited File into List of Dictionaries
def readDictsFromXSV(path, delimiter = ",", openAs = "r"):
	"""
	Read Delimited File to Dictionary List
	@params:
		path 		- Required 	: input path (Str)
		delimiter 	- Optional 	: file delimiter (Str)
		openAs 		- Optional 	: file opening flag (Str)
	Returns: dictionaries (List)
	"""
	list = []
	with open(path, openAs) as file: 
		reader = CSV.DictReader(file, delimiter = delimiter)
		for row in reader:
			list.append(row)
	return list

# Read Delimited File Columns
def readColsFromXSV(path, colStart = 0, colStop = 0, colStep = 1, keepHeader = True, delimiter = ",", openAs = "r"):
	"""
	Read Delimited File to
	@params:
		path 		- Required 	: input path (Str)
		delimiter 	- Optional 	: file delimiter (Str)
		openAs 		- Optional 	: file opening flag (Str)
	Returns: Row values (List)
	"""
	rowValues = []
	with open(path, openAs) as xsv:
		reader = CSV.reader(xsv, delimiter = delimiter)
		if keepHeader == False:
			reader.next()
		for row in reader:
			cols = row[colStart:colStop:colStep]
			if len(cols) > 0:
				rowValues.append(cols)
	return rowValues

# Write List of Dictionaries into Delimited File
def writeDictsToCSV(dictList, path = "output/output.txt", delimiter = ",", openAs = "wb"):
	"""
	Write List of Dictionaries to Delimited File
	@params:
		dictList 	- Required	: data to write (Str, Dict, List, Set)
		path 		- Optional 	: directory & file name output path (Str)
		delimiter 	- Optional 	: file delimiter (Str)
		openAs 		- Optional 	: file opening flag (Str)
	Returns: path (Str)
	"""
	with open(path, openAs) as file: 
		writer = CSV.DictWriter(file, delimiter = delimiter, fieldnames = dictList[0].keys())
		writer.writeheader()
		for dict in dictList:
			writer.writerow(dict)
	return path

# Write Python List of Dictionaries to Directory as JSON File
def writeJSON (data, path="output/output.json", openAs = "w", indent = 4, separators = (',', ':'), sortKeys = False):
	"""
	Write JSON output to file.
	@params:
		data 		- Required	: data to write (Dictionary or List)
		path 		- Optional 	: output path (Str)
		openAs 		- Optional 	: file opening flag (Str)
		indent 		- Optional 	: JSON indent size (Int)
		separators 	- Optional 	: JSON separators (Tuple)
		sortKeys 	- Optional 	: Sort JSON keys (Bool)
	Returns: path (Str)
	"""
	with open(path, openAs) as file:
		file.write(JSON.dumps(data, indent = indent, separators=separators, sort_keys = sortKeys))
		file.close()
	return path

# Prompt user input from command line
def getUserInput (valid, prompt, failed = "Error: Invalid input"):
	"""
	Prompts user for and validates input using regular expression
	@params:
		prompt 		- Required 	: verbose user prompt (Str)
		valid 		- Required 	: regex to validate against (Rgx)
	Returns: dicts (List)
	"""
	response = raw_input(prompt)
	if Rgx.match(valid, response):
		return response
	else:
		print failed
		return getUserInput(valid, prompt, failed)

def printProgress (iteration, total, prefix = '', suffix = '', barLength = 100):
	"""
	Call in a loop to create terminal progress bar
	@params:
		iterations 	- Required 	: current iteration (Int)
		total 		- Required 	: total iterations (Int)
		prefix 		- Optional	: prefix string (Str)
		suffix 		- Optional 	: suffix string (Str)
	"""
	filledLength 	= int(round(barLength * iteration / float(total)))
	percents 		= round(100.00 * (iteration / float(total)), 2)
	bar 			= '#' * filledLength + '-' * (barLength - filledLength)
	Sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix))
	Sys.stdout.flush()

# 
# Network Structure Algorithms
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

def getNodesFromEdges (edges):
	"""
	Get nodes from edge list
	@params:
		edges 		- Required 	: (List, Set)
	"""
	nodes = set()
	for edge in edges:
		if len(edge) == 2:
			for node in edge:
				nodes.add(node)
	return list(nodes)

def getComponents (graph):
	"""
	Find network components
	@params:
		graph 		- Required 	: (Dictionary)
	"""
	nodes, seen, components = graph.keys(), set(), []
	for node in nodes:
		if node not in seen:
			seen.add(node)
			search = bfs(graph, node)
			component = []
			for start, end in search:
				seen.add(start)
				seen.add(end)
				component.extend([start, end])
			components.append(component)
	return components

def getShortestPath (graph, start, end):
	"""
	Find shortest path between two vertices 
	@params:
		graph 		- Required 	: (Dictionary)
		start 		- Optional 	: (Str, Int)
		end 		- Optional 	: (Str, Int)
	"""
	parents = {}
	for parent, child in bfs(graph, start):
	    parents[child] = parent
	    if child == end:
	        reversedPath = [end]
	        while True:
	            parent = parents[child]
	            reversedPath.append(parent)
	            if parent == start:
                	break
	            child = parent
	        return list(reversed(reversedPath))
	return []

def getDiameter (graph):
	""" 
	Calculate the Diameter of a Graph
	@params:
		graph 		- Required 	: (Dictionary)
	"""
	nodeList 	= graph.keys()
	l 			= len(nodeList) - 1
	pairs 		= [(nodeList[i], nodeList[j]) for i in range(0, l) for j in range(0, l)]
	shortestPaths = []
	for start, end in pairs:
		if start is not end:
			shortestPaths.append(set(getShortestPath(graph, start, end)))
	shortestPaths.sort(key = len)
	diameter = len(shortestPaths[-1])
	return diameter

# 
# Network Data Quality Control
# 

def cleanSimpleEdges (edges):
	"""
	Remove edges from list that have weights or a length greater or less than two
	@params:
		edges 		- Required 	: (List, Set)
	"""
	clean = set()
	for edge in edges:
		if len(edge) == 2:
			for node in edge:
				if node == None or node == "-" or node == "":
					continue
				else:
					clean.add(edge)
					break
	return list(clean)

def cleanNodes (nodes):
	"""
	Remove nodes from node list that are duplicates, empty, or None values
	@params:
		nodes 		- Required 	: (List, Set)
	"""
	cleaned = set()
	for node in nodes:
		if node == None:
			continue
		elif node == "":
			continue
		else:
			cleaned.add(node)
	return list(cleaned)

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
# Network Generators
# 

def createGeodesicPathMatrix (graph, nodes = False):
	"""
	Generate Matrix of Geodesic Paths in Network
	@params:
		graph 		- Required 	: (Dictionary)
		nodes 		- Optional 	: (List, Set)
	"""
	matrix = []
	if nodes is False:
		nodes 	= graph.keys()
	nodeCount 	= len(nodes)
	for i in range(0, nodeCount):
		row = []
		for j in range(0, nodeCount):
			if nodes[i] == nodes[j]:
				row.append(0)
			else:
				path = getShortestPath(graph, nodes[i], nodes[j])
				row.append(len(path) - 1)
		matrix.append(row)
	return matrix

def generateRandomComponents (populationGraph, count = 10, nodes = 100):
	""" 
	Generate Random Graphs
	@params:
		popualationGraph 	- Required 	: (Dictionary)
		count 				- Optional 	: (Int)
		nodes 				- Optional 	: (Int)
	"""
	sampleNodes = populationGraph.keys()
	for i in range(0, count):
		sample 	= Random.sample(sampleNodes, nodes)
		graph 	= {node: [n for n in populationGraph[node] if n in sample] for node in sample}
		yield graph

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
# Network Statistical Algorithms
# 

def mean (values):
	""" 
	Compute mean of list of integers
	@params:
		values 		- Required 	: (List, Set)
	"""
	return float(sum(values))/len(values) if len(values) > 0 else float('NaN')

def getMeanShortestPath (graph, nodes = False):
	""" 
	Get Mean Shortest Path of Network Graph
	@params:
		graph 		- Required 	: (Dictionary)
	"""
	if nodes is False:
		nodes 	= graph.keys()
	nodeCount 	= len(nodes)
	return sum([sum(l) for l in createGeodesicPathMatrix(graph, nodes)]) / (float(nodeCount * nodeCount) - nodeCount)

# 
# Main Routine
# 

if __name__ == "__main__":
	print "\nStatus: assignment6.py initialized from commandline"
	exitCode 						= 0
	print "\nNotification: Developed with Python Version 2.7.8"
	print "Notification: Dependencies..."
	print "\tmatplotlib:\t", matplotlib.__version__
	print "\tnumpy:\t\t", Npy.__version__
	print "\tnetworkx\t", NX.__version__
	try:
	# Set Meta Data
		if len(Sys.argv) > 1:
			config 							= JSON.loads(Sys.argv[1])
			version 						= config["version"]
			fullNetDataPath 				= config["input"][0]["directory"] + config["input"][0]["file"]
			subNetDataPath 					= config["input"][1]["directory"] + config["input"][1]["file"]
			outputDirectory 				= config["output"]["directory"]
			settings 						= config["settings"]
		else:
			version 						= 0.2
			fullNetDataPath 				= "data/homo_sapiens.mitab.interactions.txt"
			subNetDataPath 					= "data/genes.txt"
			outputDirectory  				= "output/"
			settings 						= []
	except OSError as error:
		print error
		print "Error: unable to configure assignment5.py; try validating the config.json file online at JSONlint\n"
		Sys.exit(exitCode)
	try:

# 
# Assignment 4 Routine
# 

	# Core Logic 
		print "\nStatus: executing Assignment 4 routine\n"
		print "\tStatus: importing data"
		# Import Data
		interactionDataFull	 			= readColsFromXSV(fullNetDataPath, colStart = 2, colStop = 4, keepHeader = False, delimiter = "\t")
		geneDataFull					= set()
		for gene in readColsFromXSV(subNetDataPath, colStart = 0, colStop = 1, keepHeader = False, delimiter = ",", openAs = "U"):
			geneDataFull.add(gene[0])
		allInteractionPairs 			= set()
		matchedInteractions 			= set()
		matchedGenes 					= set()
		edgeData 						= set()
		verticeData 					= set()
		connectedEdgeData 				= set()
		connectedVerticeData 			= set()

	# 
	# Computing Sub Network 
	# 

		# Get Sub Network With At Least One Interaction Gene in TCGA Gene List
		print "\tStatus: mapping subnetwork"
		for row in interactionDataFull:
			# Convert Row to String (join), Extract Genes Using Rgx Substitution & String Replace, Split Resulting String Into a List, Coerce into a Tuple
			rowTuple 						= tuple(Rgx.sub(r"uniprotkb\:[a-zA-Z0-9]{1,8}\_", "", "".join(row)).replace("(shortlabel)", "_").split("_")[0:2])
			allInteractionPairs.add(rowTuple)
			# allEdges.add(rowTuple)
			# for node in rowTuple
			for gene in geneDataFull:
				# Remove Interactions Not Involved With At Least 1 Gene
				if gene in rowTuple:
					matchedInteractions.add(rowTuple)
					matchedGenes.add(gene)
		print "\n\t\tNotification: Simple Graph With at Least One Interaction Gene in TCGA Gene List Contains"
		print "\t\t\tEdges:", len(matchedInteractions)
		print "\t\t\tVertices:", len(matchedGenes)
		# matchedInteractions Now has 7462 Interactions
		# matchedGenes Now has 87 Genes

		print "\t\tNotification: Simple Graph With No Synonymous Interactions or Self Loops Contains"
		edgeData, verticeData = generateIntersectionSubgraph(matchedInteractions, matchedGenes)
		print "\t\t\tEdges:", len(edgeData)
		print "\t\t\tVertices:", len(verticeData)
		# edgeData Now Has 186 Interactions
		# verticeData Now Has 69 Genes

	# 
	# Generating NetworkX Graph & Computing Basic Metrics
	# 

		print "\n\tStatus: generating network graph of largest component & computing metrics\n"
		edgeData = list(edgeData)
		verticeData = list(verticeData)
		# Generate NetworkX Graph Object
		graph 		= NX.Graph()
		# Adding Vertices Is Probably Redundant...
		graph.add_edges_from(edgeData)

		# Saving old Graph for downstream component analysis
		oldGraph 	= graph
		oldNodes 	= oldGraph.nodes()
		oldEdges 	= oldGraph.edges()
		graph 		= max(NX.connected_component_subgraphs(graph), key=len)
		nodeCount 	= graph.number_of_nodes()
		edgeCount 	= graph.number_of_edges()

		print "\t\tNotification: Largest Connected Component"
		print "\t\t\tEdges:", edgeCount
		print "\t\t\tVertices:", nodeCount
		# The Graphs Largest Component Has 163 Interactions
		# The Graphs Largest Component Has 67 Genes
		print "\t\t\tProportion of TCGA Genes =", round(float(nodeCount) / len(geneDataFull), 8)

		# Save Sub Network in Edges File & Vertices File
		# NX.write_edgelist(graph, outputDirectory + "/sub-network-edges.txt", data = False)
		# writeOutput(NX.nodes(graph), outputDirectory + "/sub-network-vertices.txt")

		# Plot Network & Save
		print "\n\tStatus: rendering network graph"
		# NX.draw_spring(graph, with_labels = True, font_size = 4, vertice_size=400, iterations = 1000, alpha = 0.2)
		# Plot.savefig(outputDirectory + "/subnetwork-graph.svg")

	# 
	# Computing Centralities
	# 

		# Compute Degree Centrality & Write Output to JSON
		print "\tStatus: computing degree centralities"
		degreeCentrality = NX.degree_centrality(graph)
		degreeCentralityDicts = []
		for key, value in degreeCentrality.items():
			dict 							= {}
			dict["gene"] 					= key
			dict["degreeCentrality"] 		= value
			degreeCentralityDicts.append(dict)
		# writeDictsToCSV(degreeCentralityDicts, path = outputDirectory + "/sub-network-degree-centrality.txt")
		# writeJSON(degreeCentralityDicts, path = outputDirectory + "/sub-network-degree-centrality.json")
		# Compute Max & Min Degree Centralities
		maxDegreeCentrality 	= max(degreeCentralityDicts, key = lambda d: d["degreeCentrality"])
		minDegreeCentrality 	= min(degreeCentralityDicts, key = lambda d: d["degreeCentrality"])

		# Compute Page Rank & Write Output to JSON
		print "\tStatus: computing pagerank"
		pageRank = NX.pagerank(graph)
		pageRankDicts = []
		for key, value in pageRank.items():
			dict 							= {}
			dict["gene"] 					= key
			dict["pageRank"] 				= value
			pageRankDicts.append(dict)
		# writeDictsToCSV(pageRankDicts, path = outputDirectory + "/sub-network-page-rank.txt")
		# writeJSON(pageRankDicts, path = outputDirectory + "/sub-network-page-rank.json")
		# Compute Max & Min Page Ranks
		maxPagerank 			= max(pageRankDicts, key = lambda d: d["pageRank"])
		minPagerank 			= min(pageRankDicts, key = lambda d: d["pageRank"])

		# Compute Eigenvector Centrality & Write Output to JSON
		print "\tStatus: computing eigenvector centralities"
		eigenvectorCentrality = NX.eigenvector_centrality(graph)
		eigenvectorCentralityDicts = []
		for key, value in eigenvectorCentrality.items():
			dict 							= {}
			dict["gene"] 					= key
			dict["eigenvectorCentrality"] 	= value
			eigenvectorCentralityDicts.append(dict)
		# writeDictsToCSV(eigenvectorCentralityDicts, path = outputDirectory + "/sub-network-eigenvector-centrality.txt")
		# writeJSON(eigenvectorCentralityDicts, path = outputDirectory + "/sub-network-eigenvector-centrality.json")
		# Compute Max & Min Eigenvector Centralities
		maxEigenvectorCentrality= max(eigenvectorCentralityDicts, key = lambda d: d["eigenvectorCentrality"])
		minEigenvectorCentrality= min(eigenvectorCentralityDicts, key = lambda d: d["eigenvectorCentrality"])

		# Get Adjacency Matrix & Compute Eigenvalues for Katz Centrality alpha Parameter
		lMatrix  				= NX.adjacency_matrix(graph)
		eigenVals 				= LinearAlgebra.eigvals(lMatrix.A)
		maxEigenvalueImaginary 	= max(eigenVals)
		maxEigenvalue 			= max(eigenVals).real
		maxEigenvalueInverse 	= 1 / maxEigenvalue
		# Subtracting an additional 0.01 from inverse of max eigenvalue
		katzAlpha 				= maxEigenvalueInverse - 0.0001

		# Compute Katz Centrality & Write Output to JSON
		print "\tStatus: computing katz centralities"
		katzCentrality = NX.katz_centrality_numpy(graph, alpha = katzAlpha)
		katzCentralityDicts = []
		for key, value in katzCentrality.items():
			dict 							= {}
			dict["gene"] 					= key
			dict["katzCentrality"] 			= value
			katzCentralityDicts.append(dict)
		# writeDictsToCSV(katzCentralityDicts, path = outputDirectory + "/sub-network-katz-centrality.txt")
		# writeJSON(katzCentralityDicts, path = outputDirectory + "/sub-network-katz-centrality.json")
		# Compute Max & Min Katz Centralities
		maxKatzCentrality 			= max(katzCentralityDicts, key = lambda d: d["katzCentrality"])
		minKatzCentrality 			= min(katzCentralityDicts, key = lambda d: d["katzCentrality"])
		# Output metrics
		print "\n\t\tNotification: computed metrics..."
		print "\t\t\tMax Degree Centrality \t\t=", 				round(maxDegreeCentrality["degreeCentrality"], 8)			, "\tGene:", maxDegreeCentrality["gene"]
		print "\t\t\tMin Degree Centrality \t\t=", 				round(minDegreeCentrality["degreeCentrality"], 8)			, "\tGene:", minDegreeCentrality["gene"]
		print "\t\t\tMax Pagerank \t\t\t=", 					round(maxPagerank["pageRank"], 8)							, "\tGene:", maxPagerank["gene"]
		print "\t\t\tMin Pagerank \t\t\t=", 					round(minPagerank["pageRank"], 8)							, "\tGene:", minPagerank["gene"]
		print "\t\t\tMax Eigenvector \t\t=", 					round(maxEigenvectorCentrality["eigenvectorCentrality"], 8)	, "\tGene:", maxEigenvectorCentrality["gene"]
		print "\t\t\tMin Eigenvector \t\t=", 					round(minEigenvectorCentrality["eigenvectorCentrality"], 8)	, "\tGene:", minEigenvectorCentrality["gene"]
		print "\t\t\tMax Katz Centrality \t\t=", 				round(maxKatzCentrality["katzCentrality"], 8)				, "\tGene:", maxKatzCentrality["gene"]
		print "\t\t\tMin Katz Centrality \t\t=",				round(minKatzCentrality["katzCentrality"], 8)				, "\tGene:", minKatzCentrality["gene"]
		print "\n\t\tNotification: Katz alpha parameter computations..."
		print "\t\t\tMax Eigenvalue* \t\t=", 					round(maxEigenvalue, 8)
		print "\t\t\tInverse of Max Eigenvalue \t=", 			round(maxEigenvalueInverse, 8)
		print "\t\t\tKatz Centrality alpha Computed Using: (1 / Max Eigenvalue) - 0.0001"
		print "\t\t\tKatz Centrality alpha Parameter Used:", 	round(katzAlpha, 8)
		print "\n\t\t\t*Max Eigenvalue computed is imaginary type", maxEigenvalueImaginary, "only real portion used for computation.\n"

# 
# Assignment 5 Routine
# 

		print "\nStatus: executing Assignment 5 routine..."
		print "\n\tStatus: breadth-first search"
		# Generate Validation Statistics
		print "\n\t\tNotification: NetworkX validation metrics..."
		print "\t\t\tConnected Components:", NX.number_connected_components(oldGraph)
		print "\t\t\tMean Shortest Path:", NX.average_shortest_path_length(graph)
		print "\t\t\tNetwork Diameter:", NX.diameter(graph)

		# Get Components, Paths, & Diameter
		compEdges 					= graph.edges()
		compNodes 					= graph.nodes()
		nodeCount 					= len(compNodes)
		startNode 					= maxDegreeCentrality["gene"]
		graph 						= convertEdgelistToGraph(compEdges)
		oldGraph 					= convertEdgelistToGraph(oldEdges)
		shortestPathMatrix 			= createGeodesicPathMatrix(graph)
		meanShortestPath 			= getMeanShortestPath(graph)
		graphDiameter 				= getDiameter(graph)
		components 					= getComponents(oldGraph)
		compCount 					= len(components)
		# Print Computed Statistics
		print "\n\t\tNotification: Custom function metrics..."
		print "\t\t\tConnected Components:", 	compCount
		print "\t\t\tMean Shortest Path:", 	 	meanShortestPath
		print "\t\t\tNetwork Diameter:", 		graphDiameter

		# Print Shortest Path Matrix
		print "\n\t\tPrompt: Print shortest path matrix?"
		selection 					= getUserInput(valid=r"[yn]{1}", prompt="\t\tHint: enter 'y' or 'n'\n\t\tSelection: ", failed = "\t\t\tError: Invalid input")
		if selection == "y":
			for shortestPaths in shortestPathMatrix:
				print shortestPaths, "\n"
		else:
			pass

# 
# Assignment 6 Routine
# 
		
		print "\nStatus: executing Assignment 6 routine..."
		distribution 				= []
		imgSaveString 				= outputDirectory + 'mean-shortest-path-distribution-tcga-samples.png'
		print "\n\tStatus: finding largest component of full file network"
		allEdges 					= cleanSimpleEdges(allInteractionPairs)
		allNodes 					= cleanNodes(getNodesFromEdges(allEdges))

		# Generate NetworkX Graph Object
		allGraph 					= NX.Graph()
		# Adding Adding Edges; Adding Vertices Too Is Probably Redundant...
		allGraph.add_edges_from(allEdges)

		print "\n\tStatus: generating network graph of largest component & computing metrics"
		components 					= NX.connected_component_subgraphs(allGraph)
		giantComponent 				= max(components, key = len)
		giantEdges 					= giantComponent.edges()
		giantNodes 					= giantComponent.nodes()
		giantEdgeCount 				= len(giantEdges)
		giantNodeCount 				= len(giantNodes)
		giantComponentGraph			= convertEdgelistToGraph(giantEdges)
		tcgaGenes 					= set(giantNodes) & geneDataFull
		tcgaCount 					= len(tcgaGenes)
		print "\n\t\tNotification: giant component has", giantEdgeCount, "edges,", giantNodeCount, "nodes, and includes", tcgaCount, "TCGA genes"
		print "\n\tStatus: computing tcga subgraph mean shortest path (MSP)"
		tcgaMSP 				 	= getMeanShortestPath(giantComponentGraph, list(tcgaGenes))
		print "\n\t\tNotification: TCGA MSP:", tcgaMSP
		print "\n\t\tPrompt: load distribution?"
		loadResponse 				= getUserInput(valid=r"([yn]{1}|exit)", prompt="\t\tHint: enter 'y', 'n', or 'exit'. 'n' will begin distribution generation routine. \n\t\tSelection: ", failed = "\t\t\tError: Invalid input")
		if loadResponse == "exit":
			print "\nStatus: exiting\n"
			exit()
		elif loadResponse == "n":
			print "\n\t\tPrompt: number of samples to be generated for distribution?"
			print "\t\t\tTime Estimations..."
			print "\t\t\t\t10 	~ 10 - 20m"
			print "\t\t\t\t100 	~ 1.5h - 3.5h"
			print "\t\t\t\t1000 	~ 12h - 36h"
			iterationsResponse			= getUserInput(valid=r"([1-9][0-9]{0,4}|exit)", prompt="\t\tHint: enter a number between 0 - 10000 or 'exit' to exit\n\t\tIterations: ", failed = "\t\t\tError: Invalid input")
			if iterationsResponse == "exit":
				print "\nStatus: exiting\n"
				exit()
			print "\n\tStatus: generating distribution"
			iterations 				= int(iterationsResponse)
			# Generate Distribution Samples
			randomComponents 		= generateRandomComponents(giantComponentGraph, iterations, tcgaCount)
			rgIndex 				= 0
			printProgress(rgIndex, iterations, prefix = "\t\tProgress:", suffix = "of Samples Generated", barLength = 50)
			for randomComponent in randomComponents:
				distribution.append(getMeanShortestPath(giantComponentGraph, randomComponent.keys()))
				rgIndex 		   += 1
				printProgress(rgIndex, iterations, prefix = "\t\tProgress:", suffix = "of Samples Generated", barLength = 50)
			with open("distributions/distribution_" + str(iterations) + "-samples.pickle", "wb") as jar:
				Pickle.dump(distribution, jar)
		elif loadResponse == "y":
			print "\n\t\tPrompt: select distribution to load..."
			option = 1
			files 					= OS.listdir(OS.getcwd() + "/distributions")
			for file in files:
				print "\t\t\t[", option, "]", file
				option += 1
			distributionResponse	= int(getUserInput(valid=r"[1-9][0-9]{0,2}", prompt="\t\tHint: enter a number [ n ] to select a distribution file\n\t\tSelected File: ", failed = "\t\t\tError: Invalid input")) - 1
			print "Status: loading distribution\n"
			with open("distributions/" + files[distributionResponse], "rb") as jar:
				distribution = Pickle.load(jar)
			distSamples 		= files[distributionResponse].replace("distribution_", "").replace("-samples.pickle", "")
			iterations 			= int(distSamples)
		else:
			print "\nError: Exit on unexpected input\n"
			exit()
		samples 				= len(distribution)
		distMax 				= max(distribution)
		distMin 				= min(distribution)
		distMean 				= mean(distribution)
		weights 				= Npy.ones_like(distribution)/float(len(distribution))
		pValue 					= 0
		lesserValues 			= [value for value in distribution if value <= tcgaMSP]
		if len(lesserValues) == 0:
			pValue 				= "< 0.01"
		else:
			pValue				= "= " + str(float(len(lesserValues)) / samples)
		print "\n\n\tStatus: distribution generated"
		print "\n\t\tNotification: mean shortest path (MSP) distribution statistics..."
		print "\t\t\tSamples:\t\t", 			samples
		print "\t\t\tMSP Max:\t\t", 			distMax
		print "\t\t\tMSP Min:\t\t", 			distMin
		print "\t\t\tMSP Mean:\t\t", 			distMean
		print "\t\t\t-------------------------------------"
		print "\t\t\tTCGA MSP:\t\t", 			tcgaMSP
		print "\t\t\tNominal p-Value:\t", 		pValue
		Plot.hist(distribution, bins = 30, weights = weights, normed = 0, histtype = 'stepfilled', facecolor = "green", alpha = 0.25, cumulative = False)
		Plot.axvline(tcgaMSP, color = 'red', linestyle = 'dashed', linewidth = 1, label = "TCGA Mean\nShortest Path\n= " + str(tcgaMSP))
		Plot.xlabel('Mean Shortest Path')
		Plot.ylabel('Probability')
		Plot.title('Distribution of Mean Shortest Path from ' + str(iterations) + ' Random Network Samples\nTCGA MSP Nominal p-Value' + pValue)
		Plot.savefig(imgSaveString)
		Plot.show()
		print "\nStatus: done\n"
	except OSError as error:
		print error
		print "Error: unable to run assignment6.py\n"
		Sys.exit(exitCode)
	exitCode = 1
	Sys.exit(exitCode)
else:
	pass