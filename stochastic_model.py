## -*- coding: utf-8 -*-
"""
Created on Tues Dec 4 10:11 2018

@author: Brittany
"""

## Simulation of a stochastic process

#######################################################################################
##
##
#######################################################################################
def stochastic_model(W, N, coupling_strength, nswitch=25000):
	import numpy as np	

	t = 0
	tmax = 0 ## the time of the last switch
	# coupling_strength = 0.03
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



#######################################################################################
##
##
#######################################################################################
def load_sim_results(N, p, w_index, coupling_strength=0.03, data_dir='matrices/'):
	import pickle
	import numpy as np

	result_filename = "{0}StochasticResult_W_N{1}_p{2}_coupling{3}_index{4}.pickle".format(data_dir,N,p,coupling_strength,w_index)
	try:
		with open(result_filename, "rb") as rf:
			results = pickle.load(rf)
	except:
		print("couldn't load {0}".format(result_filename))

	return(results)



#######################################################################################
##
##
#######################################################################################
def load_single_sim_results(N, p, w_index, sim_index, coupling_strength=0.03, data_dir='matrices/'):
	import pickle
	import numpy as np

	stochastic_filename = "{0}Stochastic_Results_N{1}_p{2}_coupling{3}_index{4}_simulation{5}".format(
			data_dir, N, p, coupling_strength, w_index, sim_index)
	try:
		with open(stochastic_filename+".pickle", "rb") as stochf:
			results = pickle.load(stochf)
	except:
		print("couldn't load {0}.pickle".format(stochastic_filename))

	return(results)


#######################################################################################
##
##
#######################################################################################
def stochastic_plot_multiple_results(N, p, coupling_strength=0.03, data_dir1='matrices/', data_dir2='matrices/'):
	import pickle
	import numpy as np
	import matplotlib.pyplot as plt

	summary_filename1 = "{0}StochasticSummary_W_N{1}_p{2}_coupling{3}.pickle".format(data_dir1,N,p,coupling_strength)
	with open(summary_filename1, "rb") as rf:
		results1 = pickle.load(rf)

	summary_filename2 = "{0}StochasticSummary_W_N{1}_p{2}_coupling{3}.pickle".format(data_dir2,N,p,coupling_strength)
	with open(summary_filename2, "rb") as rf:
		results2 = pickle.load(rf)

	# alpha_div_hats = np.concatenate(results1['alpha_div_hat'], results2['alpha_div_hat'])
	# event_rates = np.concatenate(results1['event_rate'], results2['event_rate'])
	# results = load_sim_results(N, p, w_index, coupling_strength=0.03, data_dir='matrices/')

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif', size=80)
	plt.rc('xtick', labelsize=50)
	plt.rc('ytick', labelsize=50)
	plt.figure(figsize=(12,10))
	# plt.plot(alpha_div_hats, event_rates, 'o', color='b', markersize=30)
	plt.plot(results1['alpha_div_hat'], results1['event_rate'], 'o', color='C0', markersize=30)
	plt.plot(results2['alpha_div_hat'], results2['event_rate'], 'o', color='C0', markersize=30)
	plt.xlabel(r'$\hat\alpha_{\mathrm{div}}$')
	plt.ylabel('event rate')
	plt.xticks(np.arange(0.1,1,0.1))
	plt.xlim(0.02, 0.82)
	plt.tight_layout()
	mng = plt.get_current_fig_manager()
	mng.window.showMaximized()



#######################################################################################
##
##
#######################################################################################
def stochastic_raster_plot(N, fired_neurons, on_times, off_times, tmax):
	import matplotlib.pyplot as plt
	import numpy as np

	plt.figure()

	plt.rc('font', family='serif', size=60)
	plt.rc('xtick', labelsize=50)
	plt.rc('ytick', labelsize=50)
	plt.rc('lines', markersize=2, linewidth=2)

	plt.xlabel('Time')
	plt.ylabel('Neuron index')

	plt.xlim(0, tmax)
	plt.yticks(np.arange(0,3100,500))
	plt.ylim(-10,3010)

	plt.hlines(fired_neurons, on_times, off_times)
	# plt.show()



#######################################################################################
##
##
#######################################################################################
def num_active_plot(sim_index, active_count, plateau, threshold, time_bin_size, tmax):
	import math
	import numpy as np
	import matplotlib.pyplot as plt

	num_time_bins = math.ceil(tmax/time_bin_size)
	time_bins = time_bin_size*np.arange(num_time_bins)

	## plot num active neurons vs time
	plt.figure()

	plt.rc('font', family='serif', size=60)
	plt.rc('xtick', labelsize=50)
	plt.rc('ytick', labelsize=50)
	plt.rc('lines', markersize=5, linewidth=5)

	plt.xlabel('Time')
	plt.ylabel('Number of active neurons')
	# plt.suptitle("number active neurons in simulation index {0}".format(sim_index))
	plt.plot(time_bins, active_count, "b-")
	plt.plot(time_bins, plateau*np.ones(num_time_bins), "r--", label = "plateau")
	plt.plot(time_bins, threshold*np.ones(num_time_bins), "--", color = 'gold', label = "threshold")
	plt.legend()
	plt.show()



#######################################################################################
##
##
#######################################################################################
def event_times_hist_plot(N, p, i, coupling_strength, data_dir="matrices/", bins=50):
	import numpy as np
	import matplotlib.pyplot as plt
	import pickle

	stochastic_filename = "{0}Stochastic_Results_N{1}_p{2}_{3}_{4}".format(data_dir, N, p, i, coupling_strength)
	with open(stochastic_filename+".pickle", "rb") as stochf:
		stats = pickle.load(stochf)
	
	print("alpha_div: {0}".format(stats['alpha_div_hat']))
	print("coupling strength: {0}".format(coupling_strength))
	print("number of simulations: {0}".format(stats['num_sim']))
	print("mean event time: {0}".format(stats["mean_event_time"]))
	print("median event time: {0}".format(stats["median_event_time"]))
	print("std event time: {0}".format(stats["std_event_time"]))

	plt.hist(stats["event_times"],bins)
	plt.show()




#######################################################################################
##
##
#######################################################################################
def mean_event_time_plot(N, p, coupling_strength, indices, data_dir="matrices/", reload = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import pickle

	stoch_summary_filename = "{0}Stochastic_Summary_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p,coupling_strength)

	if reload:
		with open(stoch_summary_filename, "rb") as rf:
			results = pickle.load(rf)
	
	else: 
		n_indices = len(indices)

		results = dict([('N', N*np.ones(n_indices)), ('L_left', np.zeros(n_indices)), ('L_right', np.zeros(n_indices)), 
			('p', np.zeros(n_indices)), ('alpha_recip', np.zeros(n_indices)), ('alpha_conv', np.zeros(n_indices)), 
			('alpha_div', np.zeros(n_indices)), ('alpha_chain', np.zeros(n_indices)), ('p_hat', np.zeros(n_indices)), 
			('alpha_recip_hat', np.zeros(n_indices)), ('alpha_conv_hat', np.zeros(n_indices)), 
			('alpha_div_hat', np.zeros(n_indices)), ('alpha_chain_hat', np.zeros(n_indices)),
			('mean_event_time', np.zeros(n_indices)), ('median_event_time', np.zeros(n_indices)), 
			('std_event_time', np.zeros(n_indices)), ('num_sim', np.zeros(n_indices))])
			# ('average firing rate', np.zeros(n_indices)), ('event_rate', np.zeros(n_indices)), 
			# ('IEI excess_kurtosis', np.zeros(n_indices)), ('IEI skew', np.zeros(n_indices))])

		for i in range(n_indices):
			index = indices[i]
			stochastic_filename = "{0}Stochastic_Results_N{1}_p{2}_{3}_{4}".format(data_dir, N, p, index, coupling_strength)
			with open(stochastic_filename+".pickle", "rb") as stochf:
				stoch_stats = pickle.load(stochf)

			for key in results:
				try:
					results[key][i]=stoch_stats[key]
				except:
					print('No {0} in stochastic simulation {1}.  Value is set to zero'.format(key, i))

		with open(stoch_summary_filename, "wb") as rf:
			pickle.dump(results, rf)



	# plt.rc('text', usetex=True)
	# plt.rc('font', family='serif', size=80)
	# plt.rc('xtick', labelsize=50)
	# plt.rc('ytick', labelsize=50)


	plt.figure()
	plt.plot(results['alpha_conv_hat'], np.divide(np.ones(n_indices),results['mean_event_time']), 'o', markersize=30)
	plt.xlabel(r'$\alpha_{conv}$')
	plt.ylabel('mean event frequency')
	plt.xticks(np.arange(0,0.5,0.1))

	plt.figure()
	plt.plot(results['alpha_div_hat'], np.divide(np.ones(n_indices),results['mean_event_time']), 'o', markersize=30)
	plt.xlabel(r'$\alpha_{div}$')
	plt.ylabel('mean event frequency')
	plt.xticks(np.arange(0,0.5,0.1))

	plt.figure()
	plt.plot(results['alpha_chain_hat'], np.divide(np.ones(n_indices),results['mean_event_time']), 'o', markersize=30)
	plt.xlabel(r'$\alpha_{chain}$')
	plt.ylabel('mean event frequency')
	plt.xticks(np.arange(0,0.5,0.1))


	plt.figure()
	plt.scatter(results['alpha_conv_hat'], results['alpha_chain_hat'], c=np.divide(np.ones(n_indices),results['mean_event_time']),s=500)
	plt.xlabel(r'$\alpha_{conv}$')
	plt.xticks(np.arange(0,0.5,0.1))
	plt.ylabel(r'$\alpha_{chain}$')
	cb=plt.colorbar()
	cb.set_label('mean event frequency')

	plt.figure()
	plt.scatter(results['alpha_div_hat'], results['alpha_conv_hat'], c=np.divide(np.ones(n_indices),results['mean_event_time']),s=500)
	plt.xlabel(r'$\alpha_{div}$')
	plt.ylabel(r'$\alpha_{conv}$')
	cb=plt.colorbar()
	cb.set_label('mean event frequency')


	plt.figure()
	plt.scatter(results['alpha_div_hat'], results['alpha_chain_hat'], c=np.divide(np.ones(n_indices),results['mean_event_time']),s=500)
	plt.xlabel(r'$\alpha_{div}$')
	plt.ylabel(r'$\alpha_{chain}$')
	cb=plt.colorbar()
	cb.set_label('mean event frequency')

	plt.show()