## -*- coding: utf-8 -*-
"""
Created on Tues Dec 4 10:11 2018

@author: Brittany
"""

## Simulation of a stochastic process

def stochastic_model(W, N, nswitch=25000):
	import numpy as np	

	t = 0
	tmax = 0 ## the time of the last switch
	coupling_strength = 0.03
	ext_input = 0.00005 ## constant rate of external input
	death = 1 ## constant "death rate" (rate at which active neurons become inactive)
	nodes = np.arange(1,N+1)

	times = np.zeros(nswitch)
	neurons = np.zeros(nswitch)

	active_N = np.zeros(N) ## active_N[i] = 1 if node i is active

	time20percent = 0

	fired_neurons = np.zeros(nswitch)
	on_times = np.zeros(nswitch)
	off_times = np.zeros(nswitch)
	on_count = 0

	states = np.zeros((N,2))

	for i in range(nswitch):
		## calculate the current total rate
		birth = ext_input + np.squeeze(np.asarray(np.multiply(coupling_strength, np.matmul(W,active_N))))
		rate_vector = death*active_N + np.multiply((1-active_N), birth)
		total_rate =  np.sum(rate_vector)

		## generate an exponential random variable 
		x_t = np.random.exponential(1/total_rate)
		t += x_t
		times[i] = t ## update the list of times

		## determine which neuron flipped_rate
		probabilities = np.true_divide(rate_vector, total_rate)
		n = np.random.choice(nodes, p=probabilities)

		## update neurons list and active neuron vector
		neurons[i] = n
		active_N[n-1] = 1-active_N[n-1]


		## update lists
		if states[n-1,0] == 0:
			states[n-1,0] = 1
			states[n-1,1] = t
		else:
			fired_neurons[on_count] = n
			on_times[on_count] = states[n-1,1]
			off_times[on_count]= t
			on_count += 1

			states[n-1,0] = 0

		if i <= (.8*nswitch):
			time80percent = t ## time when 80% of switches have occurred

		tmax = t ## update tmax to be the time of the last switch.

	## turn off any neurons that are still on at the end of the simulation
	for n in range(1,N+1):
		if states[n-1,0]==1:
			fired_neurons[on_count] = n
			on_times[on_count] = states[n-1,1]
			off_times[on_count] = tmax
			on_count += 1

	## keep only the added values to fired_neurons, on_times, and off_times arrays
	fired_neurons = fired_neurons[:on_count]
	on_times = on_times[:on_count]
	off_times = off_times[:on_count]


	return(times, neurons, tmax, time80percent, fired_neurons, on_times, off_times)


def stochastic_plot(N, fired_neurons, on_times, off_times):
	import matplotlib.pyplot as plt

	plt.xlim(0, tmax+5)
	plt.ylim(0,N)
	plt.hlines(fired_neurons, on_times, off_times)
	plt.show()

# def active_neuron_plot(num_active_neurons, plateau, threshold, time_bin_size, tmax):