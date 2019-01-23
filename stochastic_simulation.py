# -*- coding: utf-8 -*-
"""
Created on Son Jan 6 14:19 2019

@author: Brittany
"""
## This script uses a given adjacency matrix and runs a stochastic simulation
## then produces a "rastor" plot of the results.

## other scripts needed to run this code: 
## 	stochastic_model.py


# data_dir = 'matrices/N3000_LL50_LR50_recurr_alphas_all_rand/'
# data_dir = 'matrices/N3000_LL100_LR0_ff_alpha_conv_div_rand/'
data_dir = 'matrices/N1000_erdos_renyi/'
# data_dir = "matrices/"
print("data_dir: {0}".format(data_dir))

import sys

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
p_AVG = 50/N ## average probability of connectivity between neurons
w_index = 1
L_left = float("inf")
tmax = 50
nswitch = 25000


W_filename = "{0}Wsparse_N{1}_p{2}_L{3}_{4}".format(data_dir, N, p_AVG, L_left, w_index)
with open(W_filename+'.pickle', 'rb') as wf:
    try:
        Wsparse = pickle.load(wf) ## load in W matrix
    except (EOFError):
        print("unpickling W error")
W = np.array(Wsparse.todense())

stats_filename = "{0}Stats_W_N{1}_p{2}_L{3}_{4}.pickle".format(data_dir, N, p_AVG, L_left, w_index) 
with open(stats_filename, 'rb') as statf:
	try:
		stats = pickle.load(statf)
	except (EOFError):
		print("unpickling stats error")

# active_count_matrix = []

time_bin_size = 1/N
event = False
event_times = []
for sim in range(1):

	# (times,neurons,active_neurons) = sm.stochastic_model(W, N, tmax)
	(times,neurons,time80percent) = sm.stochastic_model(W, N, nswitch, tmax)
	print("number of switches: {0}".format(len(times)))
	print("time of 80%% of switches: {0}".format(time80percent))
	print("last time: {0}".format(times[-1]))
	tmax = times[-1]


	## First find these corresponding lists for easy plotting
	fired_neurons = []
	on_times = []
	off_times = []

	states = np.zeros((N,2))
 
	for i in range(len(times)):
		n = neurons[i]
		if states[n-1,0] == 0:
			states[n-1,0] = 1
			states[n-1,1] = times[i]
		else:
			fired_neurons.append(n)
			on_times.append(states[n-1,1])
			off_times.append(times[i])

			states[n-1,0] = 0

	for n in range(1,N+1):
		if states[n-1,0]==1:
			fired_neurons.append(n)
			on_times.append(states[n-1,1])
			off_times.append(tmax)


	## Now make a list of active neurons at beginning of each time bin
	num_time_bins = math.ceil(tmax/time_bin_size)

	active_count = []
	
	for time_bin_index in range(num_time_bins):
		t = time_bin_size*time_bin_index
		num_active = sum(i <= t for i in on_times) - sum(i <= t for i in off_times)
		## the number switched on before time t - number that switched off before time t)
		active_count.append(num_active)

	# active_count_matrix.append(active_count)

	## plot num active neurons vs time
	time_bins = time_bin_size*np.arange(num_time_bins)
	plt.plot(active_count,time_bins,linestyle="-")

	## find plateau and threshold values
	time_bin_80percent = math.floor(time80percent/time_bin_size)
	plateau = float(np.array(active_count[time_bin_80percent:]).mean())
	threshold = 0.5*plateau
	print("plateau = "+plateau)
	
	plt.show()

	# if plateau > (.1*N):
	# 	event=True
	# else:
	# 	break

	## event time = first time num active neurons > threshold
	event_time = 0
	for i in range(num_time_bins):
		c = active_count[i]
		while c < threshold:
			event_time = i*time_bin_size

	if event==False:
		event_time = float(nan)

	event_times.append(event_time)





# stats['times'] = times
# stats['fired_neurons'] = neurons
# # stats['active_neurons'] = active_neurons

# stochastic_filename = "{0}stoch_results_N{1}_p{2}_{3}".format(data_dir, N, p_AVG, w_index)
# with open(stochastic_filename+".pickle", "wb") as stochf:
# 		pickle.dump(stats, stochf)

	# sm.stochastic_plot(N, fired_neurons, on_times, off_times)