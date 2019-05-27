# -*- coding: utf-8 -*-
"""
Created on Son Jan 6 14:19 2019

@author: Brittany
"""
## This script uses a given adjacency matrix and runs a stochastic simulation
## then produces a "rastor" plot of the results.

## other scripts needed to run this code: 
## 	stochastic_model.py


import sys
import time
start = time.time()
# print("start time: {0}\n".format(start))
# sys.stdout.flush()


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
nswitch = 25000

print("coupling strength: {0}".format(coupling_strength))
# print("index: {0} \ncoupling strength: {1}".format(w_index, coupling_strength))

## which matrices do you want to use:
if len(sys.argv) >= 3:
   start_index = int(sys.argv[1])
   end_index = int(sys.argv[2])
else:
   start_index = int(input_orig("enter starting index: "))
   end_index = int(input_orig("enter end index: "))

## now, let's access the matrices with chosen indices
for w_index in range(start_index, end_index+1):
	np.random.seed(w_index)

	print("\n\nNow running simulations on network {0}\n".format(w_index))

	## Load the W matrix
	W_filename = "{0}Wsparse_N{1}_p{2}_L{3}_{4}".format(data_dir, N, p, L_left, w_index)
	with open(W_filename+'.pickle', 'rb') as wf:
	    try:
	        Wsparse = pickle.load(wf)
	    except (EOFError):
	        print("unpickling W error")
	        sys.stdout.flush()
	W = np.array(Wsparse.todense())

	## Load the stats
	stats_filename = "{0}Stats_W_N{1}_p{2}_L{3}_{4}.pickle".format(data_dir, N, p, L_left, w_index)
	with open(stats_filename, 'rb') as statf:
		try:
			stats = pickle.load(statf)
		except (EOFError):
			print("unpickling stats error")
			sys.stdout.flush()

	## Print the stats
	print('Matrix {0}'.format(w_index))
	for k,v in sorted(stats.items()):
		if (not isinstance(v,np.ndarray)) and (not isinstance(v, list)):
			print(k+":{0}".format(v))
		if k=="events":
			print(k+"{0}".format(v[0:20]))
	# print("\n")


    ## Initialize simulation 
	time_bin_size = 1/N
	event_times = []
	num_sim = 200

	## run simulations
	for sim in range(num_sim):
		sim_start = time.time()
		print("\n\nSimulation {0} of {1} on network {2} of ({3} to {4})".format(sim+1, num_sim, w_index, start_index, end_index))
		(times, neurons, tmax, time80percent, fired_neurons, on_times, off_times) = sm.stochastic_model(W, N, coupling_strength, nswitch)

		print("data_dir: {0}".format(data_dir))
		print("alpha_div_hat: {0}".format(stats['alpha_div_hat']))
		
		print("time of 80% of switches: {0}".format(time80percent))
		print("tmax: {0}".format(tmax))
		sys.stdout.flush()


		## Plot simulation
		# sm.stochastic_raster_plot(N, fired_neurons, on_times, off_times, tmax)


		## Now make a list of active neurons at beginning of each time bin
		num_time_bins = math.ceil(tmax/time_bin_size)

		active_count = []
		
		## Fill the active_count arry
		for time_bin_index in range(num_time_bins):
			t = time_bin_size*time_bin_index
			## define num_active = the number switched on before time t - number that switched off before time t)
			num_active = len(on_times[np.where(on_times <= t)]) - len(off_times[np.where(off_times <= t)]) 
			active_count.append(num_active)


		## find plateau and threshold values & add to plot
		time_bin_80percent = math.floor(time80percent/time_bin_size)
		plateau = float(np.array(active_count[time_bin_80percent:]).mean())
		threshold = 0.5*plateau
		print("plateau = {0}".format(plateau))
		sys.stdout.flush()


		## plot num active neurons vs time
		# sm.num_active_plot(sim, active_count, plateau, threshold, time_bin_size, tmax)


		## Check that an event occurred
		if plateau < (.1*N):
			print("\nsimulation {0} may not have an event\n".format(sim))


		## define event time = first time num active neurons > threshold
		active_count = np.asarray(active_count)
		event_time_bin = float(np.argwhere(active_count >= threshold)[0])
		event_time = time_bin_size*event_time_bin

		event_times.append(event_time)
		print("event time: {0}".format(event_time))

		## record short details
		stats['coupling_strength'] = coupling_strength
		stats['num_switches'] = nswitch
		stats['simulation_index'] = sim
		stats['time_bin_size'] = time_bin_size
		stats['num_time_bins'] = num_time_bins

		stats['tmax'] = tmax
		stats['time80percent'] = time80percent

		stats['plateau'] = plateau
		stats['threshold'] = threshold
		stats['active_count'] = active_count
		stats['event_time_bin'] = event_time_bin
		stats['event_time'] = event_time

		stats['runtime'] = time.time()-sim_start

		## save stochastic stats
		stochastic_filename = "{0}Stochastic_Results_N{1}_p{2}_coupling{3}_index{4}_simulation{5}_Clean".format(
			data_dir, N, p, coupling_strength, w_index, sim)
		with open(stochastic_filename+".pickle", "wb") as stochf:
				pickle.dump(stats, stochf)


		## Add in more details
		stats['times'] = times
		stats['neurons'] = neurons
		stats['fired_neurons'] = fired_neurons
		stats['on_times'] = on_times
		stats['off_times'] = off_times

		## save stochastic stats
		stochastic_filename = "{0}Stochastic_Results_N{1}_p{2}_coupling{3}_index{4}_simulation{5}".format(
			data_dir, N, p, coupling_strength, w_index, sim)
		with open(stochastic_filename+".pickle", "wb") as stochf:
				pickle.dump(stats, stochf)


	print("\ntotal runtime for {0} simulations: {1} seconds\n".format(num_sim, time.time()-start))
	sys.stdout.flush()

	# plt.show()


	## now update the event_times and num_sim lists
	# try:
	# 	stats['event_times'] += event_times ## this is a list (not numpy array)
	# 	stats['num_sim'] += num_sim ## this is a number (integer)
	# except:
	# 	stats['event_times'] = event_times
	# 	stats['num_sim'] = num_sim



	# ## calculate (and update) event time stats
	# mean_event_time = np.mean(stats['event_times'])
	# median_event_time = np.median(stats['event_times'])
	# std_event_time = np.std(stats['event_times'])
	# stats['mean_event_time'] = mean_event_time
	# stats['median_event_time'] = median_event_time
	# stats['std_event_time'] = std_event_time

	# print("\ntotal number of simulations: {0}".format(stats['num_sim']))
	# print("mean event time: {0}".format(mean_event_time))
	# print("median event time: {0}".format(median_event_time))
	# print("std event time: {0}".format(std_event_time))
	# sys.stdout.flush()

	# plt.show()


	## save stochastic stats
	# with open(stochastic_filename+".pickle", "wb") as stochf:
	# 		pickle.dump(stats, stochf)
