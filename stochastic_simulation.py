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
tmax = 500


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

# (times,neurons,active_neurons) = sm.stochastic_model(W, N, tmax)
(times,neurons) = sm.stochastic_model(W, N, tmax)

stats['times'] = times
stats['fired_neurons'] = neurons
# stats['active_neurons'] = active_neurons

stochastic_filename = "{0}stoch_results_N{1}_p{2}_{3}".format(data_dir, N, p_AVG, w_index)
with open(stochastic_filename+".pickle", "wb") as stochf:
	pickle.dump(stats, stochf)

sm.stochastic_plot(N, tmax, times, neurons)