#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Developed With Python Version 2.7.8

Refer to config.json for setup variables
"""
import matplotlib
matplotlib.use("Svg")
import csv 					as CSV
import os 					as OS
import sys 					as Sys
import os.path 				as Path
import subprocess 			as Subprocess
import json 				as JSON
import re 					as Rgx
import base64 				as B64
import errno 				as Error
import networkx 			as NX
import matplotlib.pyplot 	as Plot

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
	exitCode 					= 0
	try:
	# Parse Config & Set Global Variables
		# Used for assignments frame
		# print "Status: configuring"
		# config 					= JSON.loads(Sys.argv[1])
		# version 				= config["version"]
		# fullNetDataPath 		= config["input"][0]["directory"] + config["input"][0]["file"]
		# subNetDataPath 			= config["input"][1]["directory"] + config["input"][1]["file"]
		# outputDirectory 		= config["output"]["directory"]
		# settings 				= config["settings"]
		# print "Status: done\n"

		# Set Meta Data
		version 				= 0.1
		fullNetDataPath 		= "data/homo_sapiens.mitab.interactions.txt"
		subNetDataPath 			= "data/genes.txt"
		outputDirectory  		= "output"
		settings 				= []
	except OSError as error:
		print error
		print "Error: unable to configure assignment4.py; try validating the config.json file online at JSONlint\n"
		Sys.exit(exitCode)
	try:
	# Core Logic 
		print "Status: executing main routine"
		print "\tStatus: importing data . . ."
		# Import Data
		fullNetData 			= readDictsFromXSV(fullNetDataPath, delimiter = "\t")
		fullNetKeys 			= fullNetData[0].keys()
		subNetData 				= readDictsFromXSV(subNetDataPath, delimiter = ",", openAs = "U")

		# Set List Variables
		subFullData, edgeData, subEdgeTuples = [], [], []

		print "\tStatus: subsetting full network by gene data . . ."
		# Subset Intersection of Full Net Data and Sub Net Data & Write to File (Includes mitab Meta Data) 
		for subDict in subNetData:
			found 	= False
			gene 	= subDict["gene"]
			for fullDict in fullNetData:
				for key in fullNetKeys:
					if fullDict[key].find(gene) > -1:
						subFullData.append(fullDict)
						found = True
						break
				if found == True:
					found = False
					break
		writeDictsToCSV(subFullData, path = outputDirectory + "/full-sub-intersection.txt")
		writeJSON(subFullData, path = outputDirectory + "/full-sub-intersection.json")
		
		print "\tStatus: getting edges . . ."
		# Format into Edges & Write to JSON/CSV
		for subDict in subFullData:
			dict 				= {}
			dict["unique-edge"] = subDict["unique id A"].replace("uniprotkb:", "") + "/" + subDict["unique id B"].replace("uniprotkb:", "")
			# Set Alternative IDs, just in case they'll be needed downstream
			altA 				= subDict["alternative id A"].replace("uniprotkb:", "").replace("(shortlabel)", "").split("_")
			altB 				= subDict["alternative id B"].replace("uniprotkb:", "").replace("(shortlabel)", "").split("_")
			if (len(altA) > 1) and (len(altB) > 1):
				dict["alt-edge"]= altA[1] + "/" + altB[1]
			else:
				dict["alt-edge"]= ""
			edgeData.append(dict)
		writeDictsToCSV(edgeData, path = outputDirectory + "/all-edges.txt")
		writeJSON(edgeData, path = outputDirectory + "/all-edges.json")

		# Write Unique ID Edges to JSON/CSV
		uniqueEdges = []
		for edge in edgeData:
			d = {}
			d["edge"] = edge["unique-edge"]
			uniqueEdges.append(d)
		writeDictsToCSV(uniqueEdges, path = outputDirectory + "/unique-edges.txt")
		writeJSON(uniqueEdges, path = outputDirectory + "/unique-edges.json")

		# Write Alternative ID Edges to JSON/CSV
		altEdges = []
		for edge in edgeData:
			d = {}
			d["edge"] = edge["alt-edge"]
			altEdges.append(d)
		writeDictsToCSV(altEdges, path = outputDirectory + "/alt-edges.txt")
		writeJSON(altEdges, path = outputDirectory + "/alt-edges.json")

		# Convert Dicts into Tuples For Input into NetworkX Module
		for dict in altEdges: 
			genes = dict["edge"].split("/")
			tuple 				= (genes[0], genes[1])
			subEdgeTuples.append(tuple)

		# Generate NetworkX Graph Object
		graph = NX.Graph()
		graph.add_edges_from(subEdgeTuples)

		# Compute Degree Centrality & Write Output to JSON
		print "\tStatus: computing degree centralities . . ."
		degreeCentrality = NX.degree_centrality(graph)
		degreeCentralityDicts = []
		for key, value in degreeCentrality.items():
			dict = {}
			dict["gene"] = key
			dict["degreeCentrality"] = value
			degreeCentralityDicts.append(dict)
		writeDictsToCSV(degreeCentralityDicts, path = outputDirectory + "/sub-network-degree-centrality.txt")
		writeJSON(degreeCentralityDicts, path = outputDirectory + "/sub-network-degree-centrality.json")

		# Compute Katz Centrality & Write Output to JSON
		print "\tStatus: computing katz centralities . . ."
		katzCentrality = NX.katz_centrality(graph)
		katzCentralityDicts = []
		for key, value in katzCentrality.items():
			dict = {}
			dict["gene"] = key
			dict["katzCentrality"] = value
			katzCentralityDicts.append(dict)
		writeDictsToCSV(katzCentralityDicts, path = outputDirectory + "/sub-network-katz-centrality.txt")
		writeJSON(katzCentralityDicts, path = outputDirectory + "/sub-network-katz-centrality.json")

		# Compute Page Rank & Write Output to JSON
		print "\tStatus: computing pagerank . . ."
		pageRank = NX.pagerank(graph)
		pageRankDicts = []
		for key, value in pageRank.items():
			dict = {}
			dict["gene"] = key
			dict["pageRank"] = value
			pageRankDicts.append(dict)
		writeDictsToCSV(pageRankDicts, path = outputDirectory + "/sub-network-page-rank.txt")
		writeJSON(pageRankDicts, path = outputDirectory + "/sub-network-page-rank.json")

		# Compute Eigenvector Centrality & Write Output to JSON
		print "\tStatus: computing eigenvector centralities . . ."
		eigenvectorCentrality = NX.eigenvector_centrality_numpy(graph)
		eigenvectorCentralityDicts = []
		for key, value in eigenvectorCentrality.items():
			dict = {}
			dict["gene"] = key
			dict["eigenvectorCentrality"] = value
			eigenvectorCentralityDicts.append(dict)
		writeDictsToCSV(eigenvectorCentralityDicts, path = outputDirectory + "/sub-network-eigenvector-centrality.txt")
		writeJSON(eigenvectorCentralityDicts, path = outputDirectory + "/sub-network-eigenvector-centrality.json")

		# Plot Network
		# print "\tStatus: rendering network graph . . ."
		# NX.draw_spring(graph, node_color = "lightblue", node_size=100, iterations = 1000, alpha = 0.4)
		# Plot.savefig(outputDirectory + "/full-graph.svg")
		print "Status: done\n"
	except OSError as error:
		print error
		print "Error: unable to run assignment4.py\n"
		Sys.exit(exitCode)
	exitCode = 1
	Sys.exit(exitCode)
else:
	pass