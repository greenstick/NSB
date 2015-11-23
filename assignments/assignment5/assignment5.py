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
import csv 					as CSV
import os 					as OS
import sys 					as Sys
import os.path 				as Path
import json 				as JSON
import re 					as Rgx
import errno 				as Error
import networkx 			as NX
import numpy.linalg 		as LinearAlgebra

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
		data 		- Required	: data to write (Dict or List)
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

# 
# Network Struture Generators
# 

def generateGraphAsDict (edgeList):
	"""
	Generate Graph as Dictionary
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

# 
# Network Analysis Algorithms
# 

def bfs (graph, start):
	"""
	Discover nodes in network using a Breadth-First Search
	@params:
		graph 		- Required 	: (dictionary)
		start 		- Optional 	: (Str, Int)
	"""
	queue, enqueued = deque([(None, start)]), set([start])
	while queue:
		parent, node = queue.popleft()
		yield parent, node
		v = set(graph[node]) - enqueued
		enqueued |= v
		queue.extend([(node, child) for child in v])

def getComponents (graph):
	"""
	Find network components
	@params:
		graph 		- Required 	: (dictionary)
		start 		- Optional 	: (Str, Int)
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

def shortestPath (graph, start, end):
	"""
	Find shortest path between two vertices 
	@params:
		graph 		- Required 	: (dictionary)
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
	return None

def generateGeodesicPathMatrix (graph):
	"""
	Generate Matrix of Geodesic Paths in Network
	@params:
		graph 		- Required 	: (dictionary)
		start 		- Optional 	: (Str, Int)
		end 		- Optional 	: (Str, Int)
	"""
	matrix, nodeList 	= [], graph.keys()
	nodeCount 			= len(nodeList)
	for i in range(0, nodeCount):
		row = []
		for j in range(0, nodeCount):
			if nodeList[i] == nodeList[j]:
				# row.append(0)
				pass
			else:
				path = list(set(shortestPath(graph, nodeList[i], nodeList[j])))
				row.append(len(path) - 1)
		matrix.append(row)
	return matrix

def diameter (graph):
	""" 
	Calculate the Diameter of a Graph
	@params:
		graph 		- Required 	: (dictionary)
	"""
	nodeList 	= graph.keys()
	l 			= len(nodeList) - 1
	pairs 		= [(nodeList[i], nodeList[j]) for i in range(0, l) for j in range(0, l)]
	shortestPaths = []
	for start, end in pairs:
		if start is not end:
			shortestPaths.append(set(shortestPath(graph, start, end)))
	shortestPaths.sort(key = len)
	diameter = len(shortestPaths[-1])
	return diameter

# 
# Main Routine
# 

if __name__ == "__main__":
	print "\nStatus: assignment5.py initialized from commandline\n"
	exitCode 						= 0
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
			outputDirectory  				= "output"
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
		print "Status: executing Assignment 4 routine\n"
		print "\tStatus: importing data"
		# Import Data
		netDataFul	 					= readColsFromXSV(fullNetDataPath, colStart = 2, colStop = 4, keepHeader = False, delimiter = "\t")
		netDataSub						= readColsFromXSV(subNetDataPath, colStart = 0, colStop = 1, keepHeader = False, delimiter = ",", openAs = "U")
		netDataMatch 					= set()
		genDataMatch 					= set()
		synEdgeData 					= set()
		synVerticeData 					= set()
		edgeData 						= set()
		verticeData 					= set()
		connectedEdgeData 				= set()
		connectedVerticeData 			= set()

	# 
	# Computing Sub Network 
	# 

		# Get Sub Network With At Least One Interaction Gene in TCGA Gene List
		print "\tStatus: mapping subnetwork"
		for row in netDataFul:
			# Convert Row to String (join), Extract Genes Using Rgx Substitution & String Replace, Split Resulting String Into a List, Coerce into a Tuple
			rowTuple 						= tuple(Rgx.sub(r"uniprotkb\:[a-zA-Z0-9]{1,8}\_", "", "".join(row)).replace("(shortlabel)", "_").split("_")[0:2])
			for gene in netDataSub:
				# Remove Interactions Not Involved With At Least 1 Gene
				if gene[0] in rowTuple:
					netDataMatch.add(rowTuple)
					genDataMatch.add(gene[0])
		print "\n\t\tNotification: Simple Graph With at Least One Interaction Gene in TCGA Gene List Contains"
		print "\t\t\tEdges:", len(netDataMatch)
		print "\t\t\tVertices:", len(genDataMatch)
		# netDataMatch Now has 7462 Interactions
		# genDataMatch Now has 87 Genes

		# Get Sub Network With Both Interaction Genes in TCGA Gene List
		for geneA, geneB in netDataMatch:
			# Remove Self Edges
			if geneA == geneB:
				continue
			existsGeneA, existsGeneB = False, False
			for gene in genDataMatch:
				if geneA == gene:
					existsGeneA = True
				if geneB == gene:
					existsGeneB = True
				if existsGeneA == True and existsGeneB == True:
					synEdgeData.add((geneA, geneB))
					synVerticeData.add(gene)
		print "\t\tNotification: Simple Graph With Both Interaction Genes in TCGA Gene List Contains"
		print "\t\t\tEdges:", len(synEdgeData)
		print "\t\t\tVertices:", len(synVerticeData)
		# synEdgeData Now Has 186 Interactions
		# synVerticeData Now Has 85 Genes

		print "\t\tNotification: Simple Graph With No Synonymous Interactions or Self Loops Contains"
		# Remove Synonymous Interactions
		for genesA in synEdgeData:
			geneSetA = set(genesA)
			for genesB in synEdgeData:
				geneSetB = set(genesB)
				intersection = set(geneSetA | geneSetB)
				# If Interactions Are Synonymous, the Set Intersection Will Collapse to a Maximum of l = 2
				if len(intersection) > 2:
					edgeData.add(genesA)
		print "\t\t\tEdges:", len(edgeData)
		# edgeData Now Has 186 Interactions

		# Get Unique Vertices From Interactions
		for geneA, geneB in edgeData:
			verticeData.add(geneA)
			verticeData.add(geneB)
		print "\t\t\tVertices:", len(verticeData)
		# verticeData Now Has 69 Genes

		print "\n\tStatus: generating clean network graph\n"
		edgeData = list(edgeData)
		verticeData = list(verticeData)

	# 
	# Generating NetworkX Graph & Computing Basic Metrics
	# 

		# Generate NetworkX Graph Object
		graph 		= NX.Graph()
		# Adding Vertices Is Probably Redundant...
		graph.add_nodes_from(verticeData)
		graph.add_edges_from(edgeData)

		# Saving old Graph for downstream component analysis
		oldGraph 	= graph
		oldNodes 	= oldGraph.nodes()
		oldEdges 	= oldGraph.edges()

		print "\n\tStatus: generating network graph of largest component & computing metrics\n"
		graph 		= max(NX.connected_component_subgraphs(graph), key=len)
		nodeCount 	= graph.number_of_nodes()
		edgeCount 	= graph.number_of_edges()

		print "\t\tNotification: Largest Connected Component"
		print "\t\t\tEdges:", edgeCount
		print "\t\t\tVertices:", nodeCount
		# The Graphs Largest Component Has 163 Interactions
		# The Graphs Largest Component Has 67 Genes
		print "\t\t\tProportion of TCGA Genes =", round(float(nodeCount) / len(netDataSub), 8)

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
		graph 						= generateGraphAsDict(compEdges)
		oldGraph 					= generateGraphAsDict(oldEdges)
		shortestPathMatrix 			= generateGeodesicPathMatrix(graph)
		pathsSum 					= sum([sum(l) for l in shortestPathMatrix])
		meanShortestPath 			= pathsSum / (float(nodeCount * nodeCount) - nodeCount)
		graphDiameter 				= diameter(graph)
		components 					= getComponents(oldGraph)
		compCount 					= len(components)
		# Print Computed Statistics
		print "\n\t\tNotification: Custom function metrics..."
		print "\t\t\tConnected Components:", compCount
		print "\t\t\tMean Shortest Path:", meanShortestPath
		print "\t\t\tNetwork Diameter:", graphDiameter

		# Print Shortest Path Matrix
		print "\n\t\tPrompt: Print shortest path matrix?"
		selection 					= getUserInput(valid=r"[yn]{1}", prompt="\t\t\tHint: enter 'y' or 'n'\n\t\t\tSelection: ", failed = "\t\t\tError: Invalid input")
		if selection is "y":
			for shortestPaths in shortestPathMatrix:
				print shortestPaths, "\n"
		else:
			pass
		print "\nStatus: done\n"
	except OSError as error:
		print error
		print "Error: unable to run assignment5.py\n"
		Sys.exit(exitCode)
	exitCode = 1
	Sys.exit(exitCode)
else:
	pass