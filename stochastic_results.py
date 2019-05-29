# -*- coding: utf-8 -*-
"""
Created on Wed May 29 9:23 2019

@author: Brittany
"""
## This script uses a given adjacency matrix and runs a stochastic simulation
## then produces a "rastor" plot of the results.

## other scripts needed to run this code: 
## 	stochastic_model.py


import sys
import time


# data_dir = 'matrices/N1000_Linf_recurr_alpha_div_rand/'
# data_dir = 'matrices/N3000_LL100_LR0_ff_alpha_div_rand/'
# data_dir = 'matrices/N1000_erdos_renyi/'
# data_dir = "matrices/"
# data_dir = 'matrices/N3000_Linf_homogeneous_alphas_all_zero/'
data_dir = 'matrices/N3000_Linf_homogeneous_alpha_div_half/'
# data_dir = 'matrices/test/'
print("data_dir: {0}".format(data_dir))
sys.stdout.flush()


try:
	import cPickle as pickle ## used to store python data types
except:
	import pickle

import math
import numpy as np
import scipy as sp
from scipy import sparse

try:
	del input ## brian simulator overwrites this, we want to reset it to the python default
except:
	pass
input_orig = input

import stochastic_model as sm

import matplotlib.pyplot as plt


N = 3000 ## Number of excitatory neurons
p = 50/N ## average probability of connectivity between neurons
L_left = float("inf")
# L_left = 100

coupling_strength = 0.03

reload=False
# reload=True

summary_filename = "{0}StochasticSummary_W_N{1}_p{2}_coupling{3}.pickle".format(data_dir,N,p,coupling_strength) 


if reload:
	with open(summary_filename, "rb") as rf:
		results_summary = pickle.load(rf)

else:
	## which matrices do you want to use:
	if len(sys.argv) >= 3:
		start_index = int(sys.argv[1])
		end_index = int(sys.argv[2])
	else:
		start_index = int(input_orig("enter starting index: "))
		end_index = int(input_orig("enter end index: "))

	n_indices = end_index-start_index+1

	results_summary = dict([('N', np.zeros(n_indices)), 
		('L_left', np.zeros(n_indices)), 
		('L_right', np.zeros(n_indices)), 
		('p_AVG', np.zeros(n_indices)), 
		('alpha_recip', np.zeros(n_indices)), 
		('alpha_conv', np.zeros(n_indices)), 
		('alpha_div', np.zeros(n_indices)), 
		('alpha_chain', np.zeros(n_indices)), 
		('p_hat', np.zeros(n_indices)), 
		('alpha_recip_hat', np.zeros(n_indices)), 
		('alpha_conv_hat', np.zeros(n_indices)), 
		('alpha_div_hat', np.zeros(n_indices)), 
		('alpha_chain_hat', np.zeros(n_indices)),
		('num_sim', np.zeros(n_indices)),
		('coupling_strength', np.zeros(n_indices)), 
		('event_time excess_kurtosis', np.zeros(n_indices)), 
		('event_time skew', np.zeros(n_indices)),
		('event_rate', np.zeros(n_indices)),
		('total_time', np.zeros(n_indices))])
		# ('plateau', np.zeros(n_indices)),
		# ('threshold', np.zeros(n_indices)),
		# ('event_times', np.zeros(n_indices)),
		# ('runtime', np.zeros(n_indices))])

	actual_index = -1

	## now, let's access the results for the chosen indices
	for w_index in range(start_index, end_index+1):	
		## load simulation results
		result_filename = "{0}StochasticResult_W_N{1}_p{2}_coupling{3}_index{4}.pickle".format(data_dir,N,p,coupling_strength,w_index)
		try:
			with open(result_filename, "rb") as rf:
				results = pickle.load(rf)
		except:
			print("couldn't load {0}".format(result_filename))
			continue
		
		actual_index += 1

		## add simulation data to results
		for key in results_summary:
			results_summary[key][actual_index]=results[key]

	## resize each result vector to fit the actual number of entries for the index
	for key in results_summary:
		results_summary[key].resize(actual_index+1, refcheck=False)

	## Save results summary
	with open(summary_filename, "wb") as rf:
		pickle.dump(results_summary, rf)
	print("result summary pickled")


# print("data_dir: {0}".format(data_dir))
print("\nevent_rates: {0}".format(results_summary['event_rate']))
print("mean event rate: {0}".format(np.mean(results_summary['event_rate'])))
print("\nevent_times skew: {0}".format(results_summary["event_time skew"]))
print("mean skew: {0}".format(np.mean(results_summary['event_time skew'])))