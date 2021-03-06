# -*- coding: utf-8 -*-
"""
Created on Wed May 29 9:49 2019

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
from scipy import stats

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
nswitch = 25000
num_sim = 200

if len(sys.argv) >= 3:
	start_index = int(sys.argv[1])
	end_index = int(sys.argv[2])
else:
	start_index = int(input_orig("enter network starting index: "))
	end_index = int(input_orig("enter network end index: "))
	
w_indices = end_index-start_index+1

## now, let's access the results for the chosen indices
for w_index in range(start_index, end_index+1):
	print("\n\ndata_dir: {0}".format(data_dir))
	print("\nBeginning network {0}".format(w_index))	
	sys.stdout.flush()

	results = dict([
		('num_switches', np.zeros(num_sim)),
		('simulation_index', np.zeros(num_sim)),
		('time_bin_size', np.zeros(num_sim)),
		('num_time_bins', np.zeros(num_sim)),
		('tmax', np.zeros(num_sim)),
		('time80percent', np.zeros(num_sim)),
		('plateau', np.zeros(num_sim)),
		# ('plateau_time', np.zeros(num_sim)),
		('threshold', np.zeros(num_sim)),
		('event_time_bin', np.zeros(num_sim)),
		('event_time', np.zeros(num_sim)),
		('runtime', np.zeros(num_sim))])

	actual_index = -1
	skipped = 0

	for sim in range(num_sim):
		if sim%30 == 0:
			print("index {0} simulation {1}".format(w_index, sim))
			sys.stdout.flush()

		## load simulation data
		stochastic_filename = "{0}Stochastic_Results_N{1}_p{2}_coupling{3}_index{4}_simulation{5}".format(
			data_dir, N, p, coupling_strength, w_index, sim)
		try:
			with open(stochastic_filename+".pickle", "rb") as stochf:
				samp_results = pickle.load(stochf)
		except:
			print("couldn't load {0}.pickle".format(stochastic_filename))
			print("Error: %s" % e)
			sys.stdout.flush()
			continue
		
		actual_index += 1

		## add simulation data to results
		for key in results:
			results[key][actual_index]=samp_results[key]

	## resize each result vector to fit the actual number of entries for the index
	for key in results:
		results[key].resize(actual_index+1, refcheck=False)

	## load network statistics
	stats_filename = "{0}Stats_W_N{1}_p{2}_L{3}_{4}.pickle".format(data_dir, N, p, L_left, w_index)
	with open(stats_filename, 'rb') as statf:
		try:
			Wstats = pickle.load(statf)
		except (EOFError):
			print("unpickling stats error")
			sys.stdout.flush()

	## Add network (index) specific results
	results['N'] = Wstats['N']
	results['L_left'] = Wstats['L_left']  
	results['L_right'] = Wstats['L_right']
	results['p_AVG'] = p
	results['alpha_recip'] = Wstats['alpha_recip']
	results['alpha_conv'] = Wstats['alpha_conv']
	results['alpha_div'] = Wstats['alpha_div']
	results['alpha_chain'] = Wstats['alpha_chain']
	results['p_hat'] = Wstats['p_hat']
	results['alpha_recip_hat'] = Wstats['alpha_recip_hat']
	results['alpha_conv_hat'] = Wstats['alpha_conv_hat']
	results['alpha_div_hat'] = Wstats['alpha_div_hat']
	results['alpha_chain_hat'] = Wstats['alpha_chain_hat']

	results['skipped'] = skipped
	results['num_sim'] = num_sim
	results['coupling_strength'] = coupling_strength


	## calculate skew and excess kurtosis of event times (IEIs)
	excess_kurtosis = stats.kurtosis(results['event_time'], bias=False)
	skew = stats.skew(results['event_time'], bias=False)

	results['event_time skew'] = skew
	results['event_time excess_kurtosis'] = excess_kurtosis

	## Calculate num events per unit 'time'
	# total_time = np.sum(results['time80percent'])
	# plateau_time_bin = float(np.argwhere(samp_results['active_count']>=samp_results['plateau'])[0])
	# plateau_time = samp_results['time_bin_size']*plateau_time_bin
	# total_time = np.sum(plateau_time)
	total_time = np.sum(results['event_time'])
	event_rate = np.true_divide(num_sim,total_time)


	results['total_time'] = total_time
	results['event_rate'] = event_rate

	## Save index simulation results
	result_filename = "{0}StochasticResult_W_N{1}_p{2}_coupling{3}_index{4}.pickle".format(data_dir,N,p,coupling_strength,w_index)
	with open(result_filename, "wb") as rf:
		pickle.dump(results, rf)

	print("results for network {0} pickled".format(w_index))
	sys.stdout.flush()