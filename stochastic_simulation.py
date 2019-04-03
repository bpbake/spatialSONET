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


data_dir = 'matrices/N1000_Linf_recurr_alpha_div_rand/'
# data_dir = 'matrices/N3000_LL100_LR0_ff_alpha_conv_div_rand/'
# data_dir = 'matrices/N1000_erdos_renyi/'
# data_dir = "matrices/"
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


N = 1000 ## Number of excitatory neurons
p = 50/N ## average probability of connectivity between neurons
# w_index = 4
L_left = float("inf")

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

	print("\n\n\nNow running simulations on network {0}\n\n".format(w_index))

	W_filename = "{0}Wsparse_N{1}_p{2}_L{3}_{4}".format(data_dir, N, p, L_left, w_index)
	with open(W_filename+'.pickle', 'rb') as wf:
	    try:
	        Wsparse = pickle.load(wf) ## load in W matrix
	    except (EOFError):
	        print("unpickling W error")
	        sys.stdout.flush()
	W = np.array(Wsparse.todense())

	stochastic_filename = "{0}Stochastic_Results_N{1}_p{2}_{3}_{4}".format(data_dir, N, p, w_index, coupling_strength)
	stats_filename = "{0}Stats_W_N{1}_p{2}_L{3}_{4}.pickle".format(data_dir, N, p, L_left, w_index)
	try:
		with open(stochastic_filename, 'rb') as stochf:
			try:
				stats = pickle.load(stochf)
			except (EOFError):
				print("unpickling stochastic stats error")
				sys.stdout.flush()
	except: 
		with open(stats_filename, 'rb') as statf:
			try:
				stats = pickle.load(statf)
			except (EOFError):
				print("unpickling stats error")
				sys.stdout.flush()



	time_bin_size = 1/N
	# event = False
	event_times = []
	num_sim = 500

	for sim in range(num_sim):
		sim_start = time.time()
		print("\n\nSimulation {0} of {1} on network {2} of ({3} to {4})\n".format(sim+1, num_sim, w_index, start_index, end_index))
		(times, neurons, tmax, time80percent, fired_neurons, on_times, off_times) = sm.stochastic_model(W, N, coupling_strength, nswitch)
		# print("number of switches: {0}".format(len(times)))
		print("time of 80%% of switches: {0}".format(time80percent))
		# print("last time: {0}".format(times[-1]))
		print("tmax: {0}".format(tmax))
		sys.stdout.flush()

		# now_time = time.time()
		# print("\nfinished with stochastic simulation.  \nruntime: {0}\n".format(now_time-sim_start))


		## Plot simulation
		# plt.figure()
		# plt.suptitle("rastor plot of simulation {0}".format(sim))
		# sm.stochastic_plot(N, fired_neurons, on_times, off_times, tmax)


		## Now make a list of active neurons at beginning of each time bin
		num_time_bins = math.ceil(tmax/time_bin_size)
		# print("num time bins: {0}".format(num_time_bins))
		# sys.stdout.flush()

		active_count = []
		
		for time_bin_index in range(num_time_bins):
			t = time_bin_size*time_bin_index
			num_active = len(on_times[np.where(on_times <= t)]) - len(off_times[np.where(off_times <= t)]) ## the number switched on before time t - number that switched off before time t)
			active_count.append(num_active)

		# print("\nfinished calculating num active neurons per time bin")#. \nruntime: {0}".format(time.time()-sim_start))
		# print("time spent calculating num active neurons: {0}\n".format(time.time()-now_time))
		# now_time = time.time()
		# sys.stdout.flush()


		## find plateau and threshold values & add to plot
		time_bin_80percent = math.floor(time80percent/time_bin_size)
		plateau = float(np.array(active_count[time_bin_80percent:]).mean())
		threshold = 0.5*plateau
		print("plateau = {0}".format(plateau))
		# print("\nready to plot num active neurons vs. time.  \nruntime: {0}\n".format(time.time()-sim_start))
		sys.stdout.flush()


		## plot num active neurons vs time
		# sm.num_active_plot(sim, active_count, plateau, threshold, time_bin_size, tmax)


		## Check that an event occurred
		# if plateau > (.1*N):
		# 	event=True
		# else:
		# 	break
		if plateau < (.1*N):
			print("\n\nsimulation {0} may not have an event\n\n".format(sim))


		## event time = first time num active neurons > threshold
		# print("now calculating event time")
		# sys.stdout.flush()
		active_count = np.asarray(active_count)
		event_time_bin = float(np.argwhere(active_count >= threshold)[0])
		event_time = time_bin_size*event_time_bin

		event_times.append(event_time)
		print("event time: {0}".format(event_time))


	print("\ntotal runtime for {0} simulations: {1} seconds\n".format(num_sim, time.time()-start))
	sys.stdout.flush()

	# plt.show()


	## now update the event_times and num_sim lists
	try:
		stats['event_times'] += event_times ## this is a list (not numpy array)
		stats['num_sim'] += num_sim ## this is a number (integer)
	except:
		stats['event_times'] = event_times
		stats['num_sim'] = num_sim



	## calculate (and update) event time stats
	mean_event_time = np.mean(stats['event_times'])
	median_event_time = np.median(stats['event_times'])
	std_event_time = np.std(stats['event_times'])
	stats['mean_event_time'] = mean_event_time
	stats['median_event_time'] = median_event_time
	stats['std_event_time'] = std_event_time

	print("\ntotal number of simulations: {0}".format(stats['num_sim']))
	print("mean event time: {0}".format(mean_event_time))
	print("median event time: {0}".format(median_event_time))
	print("std event time: {0}".format(std_event_time))
	sys.stdout.flush()



	# stochastic_filename = "{0}Stochastic_Results_N{1}_p{2}_{3}_{4}".format(data_dir, N, p, w_index, coupling_strength)
	with open(stochastic_filename+".pickle", "wb") as stochf:
			pickle.dump(stats, stochf)
