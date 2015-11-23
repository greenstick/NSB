#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Developed With Python Version 2.7.8

Refer to config.json for setup variables
"""

# 
# Markov Random Field (MRF) Class Definition
# 

class MarkovRandomField(object):

	# Import Dependencies
	import csv 			as CSV
	import os 			as OS
	import sys 			as Sys
	import os.path 		as Path
	import errno 		as Error
	import MarkovChain 	as markov


	# Define Default Private Variables
	__iterations = 		100

	def __init__(self, transitionMatrix = [], iterations = __iterations):
		"""
		Initialize Makov Random Field Class
		"""
		self.transitionMatrix 	= transitionMatrix
		self.iterations 		= iterations

	def solve (self, transitionMatrix = None):
		# 
		transitionMatrix = setDefault(transitionMatrix, self.transitionMatrix)
		print transitionMatrix
		for r in transitionMatrix:
			for c in r:
				print c
		return transitionMatrix

matrix 	= [
	[0, 1, 2],
	[2, 3, 4],
	[1, 5, 6]
]
mrf = MarkovRandomField(matrix)
mrf.solve()