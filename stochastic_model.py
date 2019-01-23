## -*- coding: utf-8 -*-
"""
Created on Tues Dec 4 10:11 2018

@author: Brittany
"""

## Simulation of a stochastic process

def stochastic_model(W, N, nswitch=25000, tmax=100):
	import numpy as np	

	t = 0
	coupling_strength = 0.03
	ext_input = 0.00005 ## constant rate of external input
	death = 1 ## constant "death rate" (rate at which active neurons become inactive)
	nodes = np.arange(1,N+1)

	times = []
	neurons = []

	active_N = np.zeros(N) ## active_N[i] = 1 if node i is active

	time20percent = 0

	fired_neurons = []
	on_times = []
	off_times = []

	states = np.zeros((N,2))

	for i in range(nswitch):
		## calculate the current total rate
		birth = ext_input + np.squeeze(np.asarray(np.multiply(coupling_strength, np.matmul(W,active_N))))
		rate_vector = death*active_N + np.multiply((1-active_N), birth)
		total_rate =  np.sum(rate_vector)

		## generated an exponential random variable 
		x_t = np.random.exponential(1/total_rate)
		t += x_t
		times.append(t) ## update the list of times

		## determine which neuron flipped_rate
		probabilities = np.true_divide(rate_vector, total_rate)
		n = np.random.choice(nodes, p=probabilities)

		## update neurons list and active neuron vector
		neurons.append(n)
		active_N[n-1] = 1-active_N[n-1]


		## update lists
		if states[n-1,0] == 0:
			states[n-1,0] = 1
			states[n-1,1] = times[i]
		else:
			fired_neurons.append(n)
			on_times.append(states[n-1,1])
			off_times.append(times[i])

			states[n-1,0] = 0

		if i <= (.8*nswitch):
			time80percent = t ## time when 80% of switches have occurred

	return(times,neurons,time80percent)


def stochastic_plot(N, fired_neurons, on_times, off_times):
	import matplotlib.pyplot as plt

	plt.xlim(0, tmax+5)
	plt.ylim(0,N)
	plt.hlines(fired_neurons, on_times, off_times)
	plt.show()

# def active_neuron_plot(num_active_neurons, plateau, threshold, time_bin_size, tmax):