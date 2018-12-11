## -*- coding: utf-8 -*-
"""
Created on Tues Dec 4 10:11 2018

@author: Brittany
"""

## Simulation of a stochastic process

def stochastic_model(W, N, tmax):
	import numpy as np	

	t = 0
	coupling_strength = 0.01
	ext_input = 0.001 ## constant rate of external input
	death = 1 ## constant "death rate" (rate at which active neurons become inactive)
	nodes = np.arange(1,N+1).flatten()

	times = []
	neurons = []

	active_N = np.zeros(N) ## active_N[i] = 1 if node i is active(fired)

	while t<tmax:
		## calculate the current total rate
		birth = ext_input + np.squeeze(np.asarray(np.multiply(coupling_strength, np.matmul(W,active_N)))) #
		# print("birth vector shape: {0}".format(birth.shape))
		rate_vector = death*active_N + np.multiply((1-active_N), birth)
		total_rate =  np.sum(rate_vector)

		## generated an exponential random variable 
		x_t = np.random.exponential(1/total_rate)
		t += x_t
		times.append(t) ## update the list of times

		## determine which neuron flipped_rate
		probabilities = np.true_divide(rate_vector, total_rate).flatten()
		N_flip = int(np.random.choice(nodes, 1, p=probabilities))

		## update neurons list and active neuron vector
		neurons.append(N_flip)
		active_N[N_flip-1] = 1-active_N[N_flip-1]

	return(times,neurons)
