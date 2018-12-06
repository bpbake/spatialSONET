## -*- coding: utf-8 -*-
"""
Created on Tues Dec 4 10:11 2018

@author: Brittany
"""

## Simulation of a stochastic process

def stochastic_model(W, N, tmax):
	import numpy as np	

	t = 0
	dt = 0.01 #ms
	death = 1 ## constant "death rate" (rate at which active neurons become inactive)
	birth = np.sum(W,axis=1) #

	times = []
	neurons = []

	active_N = np.zeros(N) ## active_N[i] = 1 if node i is active(fired)

	while t<tmax:
		## calculate the current total rate
		rate_vector = death*active_N + np.multiply((1-active_N), birth)
		total_rate =  np.sum(rate_vector)

		## generated an exponential random variable 
		delta_t = np.random.exponential(1/total_rate)
		t += delta_t
		times.append(t) ## update the list of times

		## determine which neuron flipped
		X_flip = np.random.uniform(0,1) 
		N_flip = #use X_flip to define which neuron flipped

		## update neurons list and active neuron vector
		neurons.appennd(N_flip)
		active_N[N_flip] = 1-active_N[N_flip]
