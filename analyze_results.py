# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 09:46:52 2017

@author: Brittany
"""
def plot_results(N,p,i, data_dir='matrices/'):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
 	  import pickle
    
    result_filename = "{0}Results_W_N{1}_p{2}_slower{3}.pickle".format(data_dir,N,p,i) 
    with open(result_filename, "rb") as rf:
       	results = pickle.load(rf)

    print('Matrix {0}'.format(i))
    for k,v in sorted(results.items()):
        if not isinstance(v,np.ndarray):
            print(k+":{0}".format(v))
    print("\n")


    matplotlib.rcParams.update({'font.size': 20})
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    # #plot the results of the simulation
    plt.figure(figsize=(20,7))
    # #subplot(122)
    # #plot(simulation_statemon.t/ms, simulation_statemon.v[0])
    # #xlabel('Time (ms)')
    # #ylabel('v')
    
    mintime=500
    maxtime=3500
    plt.subplot(211)

    plt.suptitle('Matrix {0}'.format(i))    

    inds = np.logical_and(results['spikemon times']>mintime, results['spikemon times'] < maxtime)
    plt.plot(results['spikemon times'][inds],results['spikemon indices'][inds], '.k', markersize=1)
    
    #axis([mintime, maxtime, 1, N])
    plt.xlabel('Time (ms)')
    #plt.xticks([])
    plt.ylabel('Neuron index')
    plt.tight_layout()
    
    plt.subplot(212)
    ind1=np.min(np.where(results['PRM time']>mintime))
    ind2=np.max(np.where(results['PRM time']<maxtime))
    plt.plot(results['PRM time'][ind1:ind2],results['PRM rate'][ind1:ind2])
    plt.xlabel('Time (ms)')
    plt.ylabel('Population rate (Hz)')
    #axis([1500, 2000, 0, ceil(max(results['PRM rate'])/100)*100])
    plt.tight_layout()


def create_subPR(N,p,i, data_dir='matrices/', neuron_bin_size):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import pickle
    import math

    
    result_filename = "{0}Results_W_N{1}_p{2}_slower{3}.pickle".format(data_dir,N,p,i) 
    with open(result_filename, "rb") as rf:
       	results = pickle.load(rf)

    print('Matrix {0}'.format(i))
    for k,v in sorted(results.items()):
        if not isinstance(v,np.ndarray):
            print(k+":{0}".format(v))
    print("\n")

   	# subspike_indices = []
   	# subspike_times = []
   	# # First attempt method... it definitely works, but...
   	# # takes a long time to run since we're looping through 'spikemon indices'(len ~ 30,000) N/100 = 1000 times
   	# subi = 0
   	# while subi < N:
   	# 	si = []
   	# 	st = []
   	# 	for n in range(len(results['spikemon indices'])):
   	# 		index = results['spikemon indices'][n]
   	# 		time = results['spikemon times'][n]
   	# 		if index in range(subi, subi+100):
   	# 			si.append(index)
   	# 			st.append(time)
   	# 	subspike_indices.append(si)
   	# 	subspike_times.append(st)
   	# 	subi = subi + 100


    num_neuron_bins = math.ceil(N/neuron_bin_size)
    time_bin_size = results['PRM time'][1] - results['PRM time'][0]
    num_time_bins = len(results['PRM time'])

    subPR = np.zeros((num_neuron_bins, num_time_bins))

    for n in range(len(results['spikemon indices'])):
        neuron_bin_index = math.floor((results['spikemon indices'][n])/neuron_bin_size)
        time_bin_index = results['spikemon times'][n]
        
        subPR[neuron_bin_index, time_bin_index] += 1

    return(np.divide(subPR,(np.multiply(neuron_bin_size,time_bin_size))))