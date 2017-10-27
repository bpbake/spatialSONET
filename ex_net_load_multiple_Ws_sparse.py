# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:06:45 2017

@author: rhino
"""

data_dir = 'matrices/N10000_LL70_LR0_ff_alpha_div_rand/'
# data_dir = 'matrices/N10000_LL70_LR0_ff_alpha_chain_zero/'
# data_dir = 'matrices/N10000_LL70_LR0_ff_alphas_all_rand/'
import sys

try:
   import cPickle as pickle # used to store python data types
except:
   import pickle
#import dill #pickle works fine

#from analyze import analyze_autocor # used to analyze synchrony of networks
from analyze_results import calculate_events, analyze_events

try:
    del input # brian overwrites this, we want to reset it to the python default
except:
    pass
input_orig = input # rename the python default for input (brian will overwrite it when imported)

from brian2 import * #scary, but there are so many things that we need from brian2 for simulation that it would be a pain to import them all individually


import numpy as np
np.set_printoptions(threshold=np.nan)
import scipy as sp
from scipy import sparse

import matplotlib.pyplot as plt
import math

start_scope() # start fresh with magic settings

N = 10000 # Number of excitatory neurons
p_AVG = 50/N # average probability of connectivity between neurons
neuron_bin_size = 100 # number of neurons in each neuron bin (for analysis of network simulation)

# variables and equation for voltage decay back to equilibrium (-60) for firing potential
tau = 10*ms 
Er = -60*mV # rest potential (equilibrium voltage)
tauS = 5*ms # used for inhibitory only (not used in this program)
Ee = 0*mV  #used for inhibitory only (not used in this program)
eqs= '''
dv/dt = (-v+Er)/tau: volt
''' # leaky integrate and fire model

vthreshold = -55*mV # if voltage passes this, it counts as a spike
vreset = -65*mV # reset voltage
refract = 1*ms # "cool down" time between spikes (after a spike, it can't spike again for this amount of time)

transienttime = 500*ms # getting the network into place (the start bit of the simulation)
simulationtime = 12000*ms # the part of the simulation we care about


#Set up the Neuron Groups for simulation
G = NeuronGroup(N, eqs, threshold='v>-55*mV', reset='v=-65*mV', refractory='refract', method='euler') 
# we had trouble using the variables for threshold and reset in the G setup... could probably be improved
G.v='vreset+(vthreshold-vreset)*rand()' # sets voltage dip below reset after spike

# variables that control the PoissonGroup
ext_rate = 118*Hz # rate of external input (how often input happens)
ext_mag = 1.5*mV # how much the voltage gets affected by the external input

P = PoissonGroup(N, ext_rate) # adds noise to the simulation
Sp = Synapses(P,G, on_pre="v+=ext_mag") # synapes P onto G
Sp.connect(j='i') # where to connect P and G


j = 0.18*mV # coupling strength
# Weight of neuron connection (when neuron j fires, and is connected to neuron i, this is how much voltage is passed from j to i)

S = Synapses(G, G,"w:volt",on_pre='v_post +=w') # connects G onto itself.  
S.connect() # no specificications of where connections are made... W will be used for this later


# initialize monitors for simulation... Look up in Brian documentation!!!
transient_statemon = StateMonitor(G, 'v', record=0) # 
transient_spikemon = SpikeMonitor(G) # recording times of spikes
transient_PRM = PopulationRateMonitor(G) # records voltage (of some kind)
 
store() # record state of simulation for future reference

######Load in matrices one at a time and simulate!
if len(sys.argv) >= 3:
   start_index = int(sys.argv[1])
   end_index = int(sys.argv[2])
else:
   start_index = int(input_orig("enter starting index: "))
   end_index = int(input_orig("enter end index: "))

for w_index in range(start_index, end_index+1):
    
    print("\n\nw_index = {0}".format(w_index))
    
    restore() # set the state back to what it was when the store() command was called
    
    Wsparse_filename = "{0}Wsparse_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,w_index)
    with open(Wsparse_filename, 'rb') as wf:
        try:
            Wsparse = pickle.load(wf) # load in W matrix
        except (EOFError):
            print("unpickling error")
    W = np.array(Wsparse.todense())
    
    
    stats_filename = "{0}Stats_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,w_index)
    with open(stats_filename, 'rb') as sf:
        try:
            stats = pickle.load(sf) # load in the stats for the W matrix (L, p_hat, alpha values, alpha_hat values)
        except (EOFError):
            print("unpickling error")
    
    # translate the pickled matrix to be useful with brian, scale to appropriate connection weight (j)
    S.w=W.transpose().flatten()*j

    # run to setup the simulation 
    run(transienttime)
    
    # if the setup yeilded either of these problems (too many spikes or not enough), we go to next iteration of loop (next W)
    # are the threshold values for these if statements appropriate?
    if transient_spikemon.num_spikes > (transienttime*N/refract*0.5): # if the number of spikes it too large, assume it's saturated
        print("\nnetwork saturated, skipping matrix {0}\n".format(w_index))
        stats['saturated'] = True # add to the stats dict
        result_filename = "{0}Results_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,w_index) 
        with open(result_filename, "wb") as rf:
            pickle.dump(stats, rf) #pickle the new stats dict 
        continue # go to next matrix
        
    if transient_spikemon.num_spikes < (2*N): # if the number of spikes is too small, we assume it's not spiking
        print("\nnetwork not spiking, skipping matrix {0}\n".format(w_index))
        stats['not spiking'] = True #add to the stats dict
        result_filename = "{0}Results_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,w_index) 
        with open(result_filename, "wb") as rf:
            pickle.dump(stats, rf) # pickle the new stats file
        continue # go to next matrix

    print("\nnumber of spikes in transient: {0}\n".format(transient_spikemon.num_spikes))
    
    # print the stats for W
    for k,v in sorted(stats.items()):
        print(k+":{0}".format(v))
    print("\n")
        
    # reset monitors before run(simulationtime)
    simulation_statemon = StateMonitor(G, 'v', record=0)
    simulation_spikemon = SpikeMonitor(G)
    simulation_PRM = PopulationRateMonitor(G)   
    run(simulationtime)
    
    print("\nnumber of spikes in simulation: {0}".format(simulation_spikemon.num_spikes))


    #synchrony = analyze_autocor(simulation_PRM.rate) # monotonic measure of synchrony
    #print("\nExcitatory Synchrony = {0}\n".format(synchrony))
    
    # update the stats dict
    stats['j'] = j
    stats['ext_rate'] = ext_rate
    stats['ext_mag'] = ext_mag
    stats['simulation_time'] = simulationtime/ms
    
    #stats['synchrony'] = synchrony
    stats['PRM rate'] = simulation_PRM.rate/hertz
    stats['PRM time'] = simulation_PRM.t/ms
    stats['spikemon times'] = simulation_spikemon.t/ms
    stats['spikemon indices'] = simulation_spikemon.i/1

    events, simulation_time = calculate_events(N, stats, neuron_bin_size) # numpy array of tuples representing events
    stats['events'] = events
    stats['num events'] = len(events)
    print("\nnumber of events: {0}\n".format(len(events)))

    (event_rate, event_mag, IEIs, excess_kurtosis, skew) = analyze_events(N, events, simulationtime/ms, neuron_bin_size)
    stats['event_rate'] = event_rate
    stats['event_mag'] = event_mag
    stats['IEIs'] = IEIs 
    stats['IEI excess_kurtosis'] = excess_kurtosis
    stats['IEI skew'] = skew
    
    try:
        # #plot the results of the simulation
        figure(figsize=(20,10))
        # #subplot(122)
        # #plot(simulation_statemon.t/ms, simulation_statemon.v[0]) # plot of voltage of one node vs. time
        # #xlabel('Time (ms)')
        # #ylabel('v')
        
          
        subplot(211)
        plot(simulation_spikemon.t/ms,simulation_spikemon.i, '.k') # raster plot: y-axis = neuron index, x-axis = time, dot at (t,i) if neuron i fires at time t
        #axis([simulationtime/ms-500, simulationtime/ms, 1, N])
        xlabel('Time (ms)')
        ylabel('Neuron index')
        plt.tight_layout()
        
        subplot(212)
        plot(simulation_PRM.t/ms,simulation_PRM.smooth_rate(window='flat', width=0.5*ms)/Hz) # smooth curve of how many neurons fire at each time
    except Exception as e:
        print("Error: %s" % e)

    #delete monitors so they don't cause restore() to get an error when looped through    
    del simulation_statemon 
    del simulation_spikemon 
    del simulation_PRM 

    # save results (pickle new stats dictionary)
    # result_filename = "{0}Results_W_N{1}_p{2}_slower{3}.pickle".format(data_dir,N,p_AVG,w_index) 
    result_filename = "{0}Results_W_N{1}_p{2}_long{3}.pickle".format(data_dir,N,p_AVG,w_index) 
    with open(result_filename, "wb") as rf:
       pickle.dump(stats, rf)