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

# 
# Main Routine
# 

if __name__ == "__main__":
	print "\nStatus: assignment4.py initialized from commandline\n"
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
		print "Error: unable to configure assignment4.py; try validating the config.json file online at JSONlint\n"
		Sys.exit(exitCode)
	try:
	# Core Logic 
		print "Status: executing main routine"
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

		print "\n\tStatus: generating clean network graph, computing metrics, and saving files\n"
		edgeData = list(edgeData)
		verticeData = list(verticeData)

	# 
	# Generating NetworkX Graph & Computing Basic Metrics
	# 

		# Generate NetworkX Graph Object
		graph = NX.Graph()
		# Adding Vertices Is Probably Redundant...
		graph.add_nodes_from(verticeData)
		graph.add_edges_from(edgeData)

		# Getting The Largest Component
		graph = max(NX.connected_component_subgraphs(graph), key=len)

		print "\t\tNotification: Largest Connected Graph Component"
		print "\t\t\tEdges:", graph.number_of_edges()
		print "\t\t\tVertices:", graph.number_of_nodes()
		# The Graphs Largest Component Has 163 Interactions
		# The Graphs Largest Component Has 67 Genes
		print "\t\t\tProportion of TCGA Genes =", round(float(graph.number_of_nodes()) / len(netDataSub), 8)

		# Save Sub Network in Edges File & Vertices File
		NX.write_edgelist(graph, outputDirectory + "/sub-network-edges.txt", data = False)
		writeOutput(NX.nodes(graph), outputDirectory + "/sub-network-vertices.txt")

		# Plot Network & Save
		print "\n\tStatus: rendering network graph"
		NX.draw_spring(graph, with_labels = True, font_size = 4, vertice_size=400, iterations = 1000, alpha = 0.2)
		Plot.savefig(outputDirectory + "/subnetwork-graph.svg")

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
		writeDictsToCSV(degreeCentralityDicts, path = outputDirectory + "/sub-network-degree-centrality.txt")
		writeJSON(degreeCentralityDicts, path = outputDirectory + "/sub-network-degree-centrality.json")
		# Compute Max & Min Degree Centralities
		degreeCentralityValues 	= [dict["degreeCentrality"] for dict in degreeCentralityDicts]
		maxDegreeCentrality 	= max(degreeCentralityValues)
		minDegreeCentrality 	= min(degreeCentralityValues)

		# Compute Page Rank & Write Output to JSON
		print "\tStatus: computing pagerank"
		pageRank = NX.pagerank(graph)
		pageRankDicts = []
		for key, value in pageRank.items():
			dict 							= {}
			dict["gene"] 					= key
			dict["pageRank"] 				= value
			pageRankDicts.append(dict)
		writeDictsToCSV(pageRankDicts, path = outputDirectory + "/sub-network-page-rank.txt")
		writeJSON(pageRankDicts, path = outputDirectory + "/sub-network-page-rank.json")
		# Compute Max & Min Page Ranks
		pagerankValues 			= [dict["pageRank"] for dict in pageRankDicts]
		maxPagerank 			= max(pagerankValues)
		minPagerank 			= min(pagerankValues)

		# Compute Eigenvector Centrality & Write Output to JSON
		print "\tStatus: computing eigenvector centralities"
		eigenvectorCentrality = NX.eigenvector_centrality(graph)
		eigenvectorCentralityDicts = []
		for key, value in eigenvectorCentrality.items():
			dict 							= {}
			dict["gene"] 					= key
			dict["eigenvectorCentrality"] 	= value
			eigenvectorCentralityDicts.append(dict)
		writeDictsToCSV(eigenvectorCentralityDicts, path = outputDirectory + "/sub-network-eigenvector-centrality.txt")
		writeJSON(eigenvectorCentralityDicts, path = outputDirectory + "/sub-network-eigenvector-centrality.json")
		# Compute Max & Min Eigenvector Centralities
		eigenvectorValues 		= [dict["eigenvectorCentrality"] for dict in eigenvectorCentralityDicts]
		maxEigenvector 			= max(eigenvectorValues)
		minEigenvector 			= min(eigenvectorValues)

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
		writeDictsToCSV(katzCentralityDicts, path = outputDirectory + "/sub-network-katz-centrality.txt")
		writeJSON(katzCentralityDicts, path = outputDirectory + "/sub-network-katz-centrality.json")
		# Compute Max & Min Katz Centralities
		katzCentralityValues 		= [dict["katzCentrality"] for dict in katzCentralityDicts]
		maxKatzCentrality 			= max(katzCentralityValues)
		minKatzCentrality 			= min(katzCentralityValues)

		print "\n\t\tNotification: computed Metrics:"
		print "\t\t\tMax Degree Centrality \t\t=", round(maxDegreeCentrality, 8)
		print "\t\t\tMin Degree Centrality \t\t=", round(minDegreeCentrality, 8)
		print "\t\t\tMax Pagerank \t\t\t=", round(maxPagerank, 8)
		print "\t\t\tMin Pagerank \t\t\t=", round(minPagerank, 8)
		print "\t\t\tMax Eigenvector \t\t=", round(maxEigenvector, 8)
		print "\t\t\tMin Eigenvector \t\t=", round(minEigenvector, 8)
		print "\t\t\tMax Katz Centrality \t\t=", round(maxKatzCentrality, 8)
		print "\t\t\tMin Katz Centrality \t\t=", round(minKatzCentrality, 8)
		print "\n\t\t\tMax Eigenvalue* \t\t=", round(maxEigenvalue, 8)
		print "\t\t\tInverse of Max Eigenvalue \t=", round(maxEigenvalueInverse, 8)
		print "\t\t\tKatz Centrality alpha Computed Using: (1 / Max Eigenvalue) - 0.0001"
		print "\t\t\tKatz Centrality alpha Parameter Used:", round(katzAlpha, 8)
		print "\n\t\t\t*Max Eigenvalue computed is imaginary", maxEigenvalueImaginary, "only real portion used for computation."

		print "\nStatus: done\n"
	except OSError as error:
		print error
		print "Error: unable to run assignment4.py\n"
		Sys.exit(exitCode)
	exitCode = 1
	Sys.exit(exitCode)
else:
	pass