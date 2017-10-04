# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 08:56:36 2017

@author: Brittany
"""

##########################################################################
# Inputs for the calculate_events function:
#		N: number of neurons in the network
#		results: a python dictionary containing "PRM time", "spikemon times", "spikemon indices"
#		neuron_bin_size: how many neurons to use per neuron bin (default = 100)
#
# Calls the analyze_results.py script
# 		uses functions: create_subPR, get_thresholds, and get_events
#
# Returns a numpy array/list of event "objects" as tuples:
#		(start_neuron_bin, end_neuron_bin, start_time, end_time)
#
#########################################################################
def calculate_events(N, results, neuron_bin_size=100): 
	import numpy as np
	import analyze_results as ar
 	import math

	num_neuron_bins = math.ceil(N/neuron_bin_size)

	(subPR, time_bin_size) = ar.create_subPR(results, neuron_bin_size, num_neuron_bins)
	thresholds = ar.get_thresholds(subPR, num_neuron_bins)
	events = ar.get_events(subPR, thresholds, num_neuron_bins, time_bin_size) # a numpy array of tuples

	return(events)

def analyze_events(N, events, simulation_time):
	import numpy as np
	
