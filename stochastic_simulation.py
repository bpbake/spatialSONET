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
print("start time: {0}\n".format(start))
sys.stdout.flush()


# data_dir = 'matrices/N3000_LL50_LR50_recurr_alphas_all_rand/'
# data_dir = 'matrices/N3000_LL100_LR0_ff_alpha_conv_div_rand/'
data_dir = 'matrices/N1000_erdos_renyi/'
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
p_AVG = 50/N ## average probability of connectivity between neurons
w_index = 1
L_left = float("inf")

nswitch = 25000


W_filename = "{0}Wsparse_N{1}_p{2}_L{3}_{4}".format(data_dir, N, p_AVG, L_left, w_index)
with open(W_filename+'.pickle', 'rb') as wf:
    try:
        Wsparse = pickle.load(wf) ## load in W matrix
    except (EOFError):
        print("unpickling W error")
        sys.stdout.flush()
W = np.array(Wsparse.todense())

stats_filename = "{0}Stats_W_N{1}_p{2}_L{3}_{4}.pickle".format(data_dir, N, p_AVG, L_left, w_index) 
with open(stats_filename, 'rb') as statf:
	try:
		stats = pickle.load(statf)
	except (EOFError):
		print("unpickling stats error")
		sys.stdout.flush()



time_bin_size = 1/N
event = False
event_times = []

for sim in range(1):
	(times, neurons, tmax, time80percent, fired_neurons, on_times, off_times) = sm.stochastic_model(W, N, nswitch)
	print("number of switches: {0}".format(len(times)))
	print("time of 80%% of switches: {0}".format(time80percent))
	print("last time: {0}".format(times[-1]))
	print("tmax: {0}".format(tmax))
	sys.stdout.flush()
	# tmax = times[-1]

	on_times = np.asarray(on_times)
	off_times = np.asarray(off_times)

	now_time = time.time()
	print("\nfinished with stochastic simulation.  \nruntime: {0}\n".format(now_time-start))


	## Now make a list of active neurons at beginning of each time bin
	num_time_bins = math.ceil(tmax/time_bin_size)
	print("num time bins: {0}".format(num_time_bins))
	sys.stdout.flush()

	active_count = []
	
	for time_bin_index in range(num_time_bins):
		t = time_bin_size*time_bin_index
		num_active = len(on_times[np.where(on_times <= t)]) - len(off_times[np.where(off_times <= t)]) ## the number switched on before time t - number that switched off before time t)
		active_count.append(num_active)

	print("\nfinished calculating num active neurons per time bin. \nruntime: {0}".format(time.time()-start))
	print("time spent calculating num active neurons: {0}\n".format(time.time()-now_time))
	now_time = time.time()
	sys.stdout.flush()


	## plot num active neurons vs time
	time_bins = time_bin_size*np.arange(num_time_bins)
	plt.plot(time_bins, active_count, "b-")


	## find plateau and threshold values & add to plot
	time_bin_80percent = math.floor(time80percent/time_bin_size)
	plateau = float(np.array(active_count[time_bin_80percent:]).mean())
	threshold = 0.5*plateau
	plt.plot(time_bins, plateau*np.ones(num_time_bins), "r-")
	plt.plot(time_bins, threshold*np.ones(num_time_bins), "r-")
	print("plateau = {0}".format(plateau))
	print("\nready to plot num active neurons vs. time.  \nruntime: {0}\n".format(time.time()-start))
	sys.stdout.flush()

	# if plateau > (.1*N):
	# 	event=True
	# else:
	# 	break

	print("now calculating event time")
	sys.stdout.flush()
	## event time = first time num active neurons > threshold
	# event_time = 0
	# for i in range(num_time_bins):
	# 	c = active_count[i]
	# 	if c >= threshold:
	# 		event_time = i*time_bin_size
	# 		break

	active_count = np.asarray(active_count)
	event_time_bin = float(np.argwhere(active_count >= threshold)[0])
	event_time = time_bin_size*event_time_bin

	# if event==False:
	# 	event_time = float(nan)

	event_times.append(event_time)
	print("event time: {0}".format(event_time))


print("total runtime: {0}".format(time.time()-start))

plt.show()
sys.stdout.flush()


# stats['times'] = times
# stats['fired_neurons'] = neurons
# # stats['active_neurons'] = active_neurons

# stochastic_filename = "{0}stoch_results_N{1}_p{2}_{3}".format(data_dir, N, p_AVG, w_index)
# with open(stochastic_filename+".pickle", "wb") as stochf:
# 		pickle.dump(stats, stochf)

	# stochastic_plot(N, fired_neurons, on_times, off_times)