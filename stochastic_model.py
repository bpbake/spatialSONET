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
	# active_neurons = [] ## neurons are labeled 1 through N

	active_N = np.zeros(N) ## active_N[i] = 1 if node i is active(fired)

	time20percent = 0


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
		N_flip = np.random.choice(nodes, p=probabilities)

		## update neurons list and active neuron vector
		neurons.append(N_flip)
		active_N[N_flip-1] = 1-active_N[N_flip-1]
		# active_neurons.append(np.copy(active_N))

		if i <= (.8*nswitch):
			time80percent = t ## time when 80% of switches have occurred

	# active_neurons = np.asarray(active_neurons)
	# return(times,neurons,active_neurons)
	return(times,neurons,time80percent)


def stochastic_plot(N, tmax, times, neurons):
	import matplotlib.pyplot as plt
	import numpy as np

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

	plt.xlim(0, tmax+5)
	plt.hlines(fired_neurons, on_times, off_times)
	plt.show()