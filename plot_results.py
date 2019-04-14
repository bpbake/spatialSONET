#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 10:25:06 2017

@author: nykamp
"""

def plot_results(N,p,i, style, data_dir='matrices/'):
    
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import pickle
    
    result_filename = "{0}Results_W_N{1}_p{2}_{3}{4}.pickle".format(data_dir,N,p,style, i) 
    #result_filename = "{0}Results_W_N{1}_p{2}_slower{3}.pickle".format(data_dir,N,p,i) 
    with open(result_filename, "rb") as rf:
       results = pickle.load(rf)

    print(data_dir)

    print('Matrix {0}, style {1}'.format(i, style))
    for k,v in sorted(results.items()):
        if (not isinstance(v,np.ndarray)) and (not isinstance(v, list)):
            print(k+":{0}".format(v))
        if k=="events":
            print(k+"{0}".format(v[0:20]))
    print("\n")


    matplotlib.rcParams.update({'font.size': 60})
    plt.rc('font', family='serif', size=60)
    plt.rc('xtick', labelsize=50)
    plt.rc('ytick', labelsize=50)

    # #plot the results of the simulation
    # plt.figure(figsize=(20,7))
    # #subplot(122)
    # #plot(simulation_statemon.t/ms, simulation_statemon.v[0])
    # #xlabel('Time (ms)')
    # #ylabel('v')
    plt.figure()

      
    mintime=500
    maxtime=results['simulation_time']+mintime

    ## Raster plot
    # plt.subplot(211)
    # plt.suptitle('Matrix: {0}, style: {1} \n data_dir: {2}'.format(i, style, data_dir))    

    inds = np.logical_and(results['spikemon times']>mintime, results['spikemon times'] < maxtime)
    plt.plot(results['spikemon times'][inds],results['spikemon indices'][inds], '.k', markersize=1)
    
    #axis([mintime, maxtime, 1, N])
    plt.xlabel('Time (ms)')
    plt.xticks(np.arange(600,1501,300))
    plt.xlim(499,1550)
    plt.ylabel('Neuron index')
    plt.yticks(np.arange(0,3001,500))
    plt.ylim(-1,3001)

    # plt.tight_layout()
    # plt.grid(True)
    

    ## PRM Plot: 
    # sigma = 100
    # gaussian_kernel = np.exp(-((np.arange(-4*sigma, 4*sigma+1, 1))**2)/(2*sigma))
    # gaussian_kernel = gaussian_kernel/np.sum(gaussian_kernel)
    # # print(gaussian_kernel)
    # plt.subplot(212)
    # ind1=np.min(np.where(results['PRM time']>mintime))
    # ind2=np.max(np.where(results['PRM time']<maxtime))
    # plt.plot(results['PRM time'][ind1:ind2],np.convolve(results['PRM rate'][ind1:ind2],gaussian_kernel,mode="same"))
    # plt.xlabel('Time (ms)')
    # plt.ylabel('Population rate (Hz)')
    # # plt.xlim(10000,15000)

    plt.show()



def plot_hist(N, p, i, style, data_dir):
    import analyze_results as ar
    import matplotlib.pyplot as plt
    import numpy as np


    print(data_dir)

    print('Matrix {0}, style {1}'.format(i, style))

    ## Load in results
    results = ar.load_results(N, p, i, style, data_dir)
    print("mean IEI:{0}\nevent_rate:{1}".format(np.mean(results['IEIs']), results['event_rate']))

    ### Plot histograms of IEIs:
    plt.figure()
    # matplotlib.rcParams.update({'font.size': 60})
    plt.rc('font', family='serif', size=60)
    plt.rc('xtick', labelsize=50)
    plt.rc('ytick', labelsize=50)
    plt.hist(results['IEIs'],50) #data , number of bins
    plt.xlabel('Inter-event interval (ms)')
    plt.ylabel('Count')
    # plt.tight_layout()

    plt.show()


#### Plot histogram of start_neuron_bin:
# start_neurons = []
# for event in results["events"]:
#   start_neurons.append(event["start_neuron_bin"])
# plt.rc('font', family='serif', size=80)
# plt.rc('xtick', labelsize=70)
# plt.rc('ytick', labelsize=70)
# plt.hist(start_neurons,50) #data , number of bins
# plt.xlabel('starting neuron bin')
# plt.ylabel('count')
# plt.suptitle("long irregular")
