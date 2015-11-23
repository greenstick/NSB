#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Developed With Python Version 2.7.8

Refer to config.json for setup variables
"""

# 
# Markov Random Field (MRF) Class Definition
# 

class MarkovChain(object):

	# Import Dependencies
	import random as Random

	def __init__(self, data):
		self.cache 				= {}
		self.data 				= data
		self.nodes 				= self.getNodes()
		self.pathLength 		= len(self.nodes)
		self.__buildCache()

	def __buildCache (self):
		for n1, n2, n3 in self.triplets():
			key 				= (n1, n2)
			if key in self.cache:
				self.cache[key].append(n3)
			else:
				self.cache[key] = [n3]

	def getNodes (self):
		self.data.seek(0)
		__data 					= self.data.read()
		nodes 					= __data.split()
		return nodes

	def triplets(self):
		""" 
		Generates node triplets from network. 

		For example, if we have a pathway represented by [A, B, C, D]
		The generator will yield:
			[A, B, C]
			[B, C, D]
		"""
		if len(self.nodes) < 3:
			return
		for i in range(len(self.nodes) - 2):
			yield (self.nodes[i], self.nodes[i + 1], self.nodes[i + 2])

	def generatePathway (self, size = 25):
		seed 					= random.randint(0, self.pathLength - 3)
		seedNode, nextNode 		= self.nodes[seed], self.nodes[seed + 1]
		n1, n2, genNodes 		= seedNode, nextNode, []
		for i in xrange(size):
			genNodes.append(n1)
			n1, n2 				= n2, random.choice(self.cache[(n1, n2)])
			genNodes.append(n2)
		return genNodes

mc = MarkovChain(["a", "b", "c", "d", "e"])
