#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Developed With Python Version 2.7.8

Refer to config.json for setup variables
"""

import csv 			as CSV
import os 			as OS
import sys 			as Sys
import os.path 		as Path
import subprocess 	as Subprocess
import json 		as JSON
import re 			as Rgx
import errno 		as Error

# 
# Utility Methods
# 

def rnd (value, precision = 0.1, rounding = "ROUND_HALF_UP"):
    """
    An accurate, predictable, no monkey business rounding function.
    @params:
        value     - Required - Numeric (String or Float)
        precision - Optional - Numeric (String or Float)
        rounding  - Optional - String  (see https://docs.python.org/3/library/decimal.html for acceptable values)
    """
    from decimal import Decimal
    return Decimal(str(value)).quantize(Decimal(str(precision)), rounding = rounding)

# Parse JSON Objects to Python Dictionaries
def parseJSONToDicts (path):
	"""
	Parse from strings to dictionaries
	@params:
		path 		- Required 	: input path (Str)
	Returns: dicts (List)
	"""
	with open(path, "r") as file:
		dicts = JSON.load(file)
	return dicts

# Make Directories From Path
def makeDirectories (path):
	"""
	Create directories from path
	@params:
		path 		- Required 	: directory path (Str)
	Returns: path (Str)
	"""
	try:
		OS.makedirs(path)
	except:
		try: 
			import os as OS
			OS.makedirs(path)
		except OSError as exception:
			if exception.errno != Error.EEXIST:
				raise
	return path

# Prompt user input from command line
def getUserInput (valid, prompt):
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
		print "Error: Invalid input"
		getUserInput(valid, prompt)

# 
# Main Routine
# 

if __name__ == "__main__": 
	print "\nStatus: core.py initialized from commandline\n"
	exitCode 					= 0
	try:
	# Parse Config & Set Global Variables
		print "Status: configuring"
		config 					= parseJSONToDicts("config.json")
		startupNotification 	= config["startupNotification"]
		assignmentConfigs 		= config["assignments"]
		validAssignments 		= 0
		currentAssignment 		= 0
		voidAssignments 		= []
		option 					= 1
		print "Status: done\n"
		if startupNotification != False:
			print "Notification:", startupNotification, "\n"
	except:
		print "Error: unable to configure core.py; try validating the config.json file online at JSONlint\n"
		Sys.exit(exitCode)
	print "Status: verifying assignments"
	try: 
		# Assignment Verification Routine
		for assignment in assignmentConfigs: 
			if Path.isfile(assignment["source"]) == False:
				print "Warning: Unable to find", assignment["name"]
				voidAssignments.append(currentAssignment)
				currentAssignment 	+= 1
			else:
				validAssignments 	+= 1
				currentAssignment 	+= 1
		# Remove Voided Assignments From Assignment Config Dictionary
		for assignment in voidAssignments:
			del assignmentConfigs[assignment]
		print "Status: done\n"
	except:
		print "Error: unable to verify"
		Sys.exit(exitCode)
	print "Input: select assignment . . ."
	try:
		# User Selection of Assignment
		for assignment in assignmentConfigs:
			print "\t[", option, "]", assignment["name"]
			option += 1
		selection 					= int(getUserInput(valid=r"[0-9]{1,2}", prompt="Hint: enter [ n ] to select the appropriate assignment\nSelection: "))
		print "Status: input received\n"
	except:
		print "Error: unable to receive user input"
		Sys.exit(exitCode)
	print "Status: retrieving assignment & scaffolding output directory"
	try: 
	# Scaffold Assignment Output Directory
		assignmentConfig 			= assignmentConfigs[selection - 1]
		assignmentName 				= assignmentConfig["name"]
		assignmentSource 			= assignmentConfig["source"]
		assignmentVersion 			= assignmentConfig["version"]
		assignmentOutputDirectory 	= assignmentConfig["output"]["directory"]
		assignmentInitArgs 			= [assignmentSource, JSON.dumps(assignmentConfig)]
		makeDirectories(assignmentOutputDirectory)
		print "Status: done\n"
	except:
		print "Error: unable to setup assignment\n"
		Sys.exit(exitCode)
	print "Status: running", assignmentName, "- Version:", assignmentVersion
	# Run Selected Assignment
	try: 
		returnCode = Subprocess.call(assignmentInitArgs, close_fds = True)
		if returnCode == 0:
			raise
		print "Status:", assignmentName, "ran successfully\n"
	except: 
		print "Error: unable to execute", assignmentName, "\n"
		Sys.exit(exitCode)
	exitCode = 1
	Sys.exit(exitCode)
else:
	pass