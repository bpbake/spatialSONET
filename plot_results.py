#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 10:25:06 2017

@author: nykamp
"""

def plot_results(N,p,i, data_dir='matrices/'):
    
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import pickle
    
    result_filename = "{0}Results_W_N{1}_p{2}_slower{3}.pickle".format(data_dir,N,p,i) 
    with open(result_filename, "rb") as rf:
       results = pickle.load(rf)

    for k,v in sorted(results.items()):
        if not isinstance(v,np.ndarray):
            print(k+":{0}".format(v))
    print("\n")


    matplotlib.rcParams.update({'font.size': 20})
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    # #plot the results of the simulation
    plt.figure(figsize=(20,10))
    # #subplot(122)
    # #plot(simulation_statemon.t/ms, simulation_statemon.v[0])
    # #xlabel('Time (ms)')
    # #ylabel('v')
    
      
    mintime=500
    maxtime=2500
    plt.subplot(211)
    

    inds = np.logical_and(results['spikemon times']>mintime, results['spikemon times'] < maxtime)
    plt.plot(results['spikemon times'][inds],results['spikemon indices'][inds], '.k', markersize=1)
    
    #axis([mintime, maxtime, 1, N])
    #xlabel('Time (ms)')
    plt.xticks([])
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