# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:06:45 2017

@author: rhino
"""

import sys

# data_dir = 'matrices/N10000_LL70_LR0_ff_alphas_all_zero/'
# data_dir = 'matrices/N10000_LL70_LR0_ff_alpha_div_half/'
# data_dir = 'matrices/N10000_LL70_LR0_ff_alpha_div_rand/'
# data_dir = 'matrices/N10000_LL70_LR0_ff_alpha_chain_zero/'
# data_dir = 'matrices/N10000_LL70_LR0_ff_alphas_all_rand/'
data_dir = "matrices/layer_test/"
print("data_dir: "+data_dir)
sys.stdout.flush()

try:
   import cPickle as pickle # used to store python data types
except:
   import pickle
#import dill #pickle works fine

#from analyze import analyze_autocor # used to analyze synchrony of networks
import analyze_results as ar

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

n = 100 #number of neurons in each layer
neuron_bin_size = 100 # number of neurons in each neuron bin (for analysis of network simulation)

sigma_square = 1 #variance of neurons in first layer
rho = 0 #rho*sigma_squre = cov of any pair of neurons in first layer

if len(sys.argv) >= 4:
   start_index = int(sys.argv[1])
   end_index = int(sys.argv[2])
   num_layers = int(sys.argv[3])
else:
   start_index = int(input_orig("enter starting index: "))
   end_index = int(input_orig("enter end index: "))
   num_layers = int(input_orig("how many layers? "))

N = n*num_layers # total number neurons in network
p = 10/n #50/N

# variables and equation for voltage decay back to equilibrium (-60) for firing potential
tau = 10*ms 
Er = -60*mV # rest potential (equilibrium voltage)
tauE = 2*ms # excitatory synaptic time constant
tauS = 5*ms # used for inhibitory only (not used in this program)
Ee = 0*mV  #used for inhibitory only (not used in this program)
eqs= '''
dv/dt = (-v+Er+s)/tau: volt
ds/dt = -s/tauE: volt
''' # leaky integrate and fire model

vthreshold = -55*mV # if voltage passes this, it counts as a spike
vreset = -65*mV # reset voltage
refract = 1*ms # "cool down" time between spikes (after a spike, it can't spike again for this amount of time)

transienttime = 500*ms # getting the network into place (the start bit of the simulation)
simulationtime = 2500*ms # the part of the simulation we care about


#Set up the Neuron Groups and subgroups for simulation
G = NeuronGroup(N, eqs, threshold='v>-55*mV', reset='v=-65*mV', refractory='refract', method='euler') 
# we had trouble using the variables for threshold and reset in the G setup... could probably be improved
G.v='vreset+(vthreshold-vreset)*rand()' # sets voltage dip below reset after spike

layers = []
for index in range(num_layers):
    layers.append(G[index*n:(index+1)*n])

# variables that control the PoissonGroup
ext_rate = 150*Hz # rate of external input (how often input happens)
ext_mag = 10*mV # how much the voltage gets affected by the external input

P = PoissonGroup(N, ext_rate) # adds noise to the simulation
Sp = Synapses(P,G, on_pre="s+=ext_mag") # synapes P onto G
Sp.connect(j='i') # where to connect P and G


j = 3.8*mV # coupling strength
# Weight of neuron connection (when neuron j fires, and is connected to neuron i, this is how much voltage is passed from j to i)

S = Synapses(G, G,"w:volt",on_pre='s_post+=w') # connects G onto itself.  
S.connect() # no specificications of where connections are made... W will be used for this later


# initialize monitors for simulation... Look up in Brian documentation!!!
transient_statemon = StateMonitor(G, 'v', record=0) # records voltage of some kind
transient_spikemon = SpikeMonitor(G) # recording times of spikes
transient_PRM = PopulationRateMonitor(G) # records "average" firing rate
 
store() # record state of simulation for future reference

######Load in matrices one at a time and simulate!
for index in range(start_index, end_index+1):
    
    print("\n\nindex = {0}".format(index))
    sys.stdout.flush()
    
    restore() # set the state back to what it was when the store() command was called
    
    W_filename = "{0}W_N{1}_p{2}_numLay{3}_{4}".format(data_dir, N, p, num_layers, index)
    with open(W_filename+'.pickle', 'rb') as wf:
        try:
            Wsparse = pickle.load(wf) # load in W matrix
        except (EOFError):
            print("unpickling error")
            sys.stdout.flush()
    W = np.array(Wsparse.todense())
    
    
    stat_filename = "{0}Stats_N{1}_p{2}_numLay{3}_rho{4}_{5}.pickle".format(data_dir, N, p, num_layers, rho, index)
    with open(stat_filename, 'rb') as sf:
        try:
            stats = pickle.load(sf) # load in the stats for the W matrix (L, p_hat, alpha values, alpha_hat values)
        except (EOFError):
            print("unpickling error")
            sys.stdout.flush()
    
    # translate the pickled matrix to be useful with brian, scale to appropriate connection weight (j)
    S.w=W.transpose().flatten()*j

    # run to setup the simulation 
    run(transienttime)
    
    # if the setup yeilded either of these problems (too many spikes or not enough), we go to next iteration of loop (next W)
    # are the threshold values for these if statements appropriate?
    if transient_spikemon.num_spikes > (transienttime*N/refract*0.5): # if the number of spikes it too large, assume it's saturated
        print("\nnetwork saturated, skipping matrix {0}\n".format(index))
        sys.stdout.flush()
        stats['saturated'] = True # add to the stats dict
        result_filename = "{0}Results_W_N{1}_p{2}_saturated_{3}.pickle".format(data_dir,N,p,index) 
        with open(result_filename, "wb") as rf:
            pickle.dump(stats, rf) #pickle the new stats dict 
        continue # go to next matrix
        
    if transient_spikemon.num_spikes < (2*N): # if the number of spikes is too small, we assume it's not spiking
        print("\nnetwork not spiking, skipping matrix {0}\n".format(index))
        sys.stdout.flush()
        stats['not spiking'] = True #add to the stats dict
        result_filename = "{0}Results_W_N{1}_p{2}_noSpike_{3}.pickle".format(data_dir,N,p,index) 
        with open(result_filename, "wb") as rf:
            pickle.dump(stats, rf) # pickle the new stats file
        continue # go to next matrix

    print("\nnumber of spikes in transient: {0}\n".format(transient_spikemon.num_spikes))
    sys.stdout.flush()
        
    # reset monitors before run(simulationtime)
    simulation_statemon = StateMonitor(G, 'v', record=0)
    simulation_spikemon = SpikeMonitor(G)
    simulation_PRM = PopulationRateMonitor(G) 
    # layers_simulation_PRM0 = PopulationRateMonitor(layers[0])
    # layers_simulation_PRM1 = PopulationRateMonitor(layers[1])
    # layers_simulation_PRM2 = PopulationRateMonitor(layers[2])
    # layers_simulation_PRM3 = PopulationRateMonitor(layers[3])
    # layers_simulation_PRM4 = PopulationRateMonitor(layers[4])
    # layers_simulation_PRM5 = PopulationRateMonitor(layers[5])
    # layers_simulation_PRM6 = PopulationRateMonitor(layers[6])
    # layers_simulation_PRM7 = PopulationRateMonitor(layers[7])
    # layers_simulation_PRM8 = PopulationRateMonitor(layers[8])
    # layers_simulation_PRM9 = PopulationRateMonitor(layers[9]) 
    layer_simulation_PRM_last = PopulationRateMonitor(layers[num_layers-1])
    
    run(simulationtime)

    print("\nnumber of spikes in simulation: {0}".format(simulation_spikemon.num_spikes))
    sys.stdout.flush()


    #synchrony = analyze_autocor(simulation_PRM.rate) # monotonic measure of synchrony
    
    # update the stats dict
    stats['j'] = j
    stats['ext_rate'] = ext_rate
    stats['ext_mag'] = ext_mag
    stats['simulation_time'] = simulationtime/ms
    
    #stats['synchrony'] = synchrony
    stats['PRM rate'] = simulation_PRM.rate/hertz
    stats['PRM time'] = simulation_PRM.t/ms
    stats['last layer PRM rate'] = layer_simulation_PRM_last.rate/hertz
    stats['last layer PRM time'] = layer_simulation_PRM_last.t/ms

    # layers_PRM_rate = []
    # layers_PRM_time = []
    # layers_PRM_rate.append(layers_simulation_PRM0.rate/hertz)
    # layers_PRM_rate.append(layers_simulation_PRM1.rate/hertz)
    # layers_PRM_rate.append(layers_simulation_PRM2.rate/hertz)
    # layers_PRM_rate.append(layers_simulation_PRM3.rate/hertz)
    # layers_PRM_rate.append(layers_simulation_PRM4.rate/hertz)
    # layers_PRM_rate.append(layers_simulation_PRM5.rate/hertz)
    # layers_PRM_rate.append(layers_simulation_PRM6.rate/hertz)
    # layers_PRM_rate.append(layers_simulation_PRM7.rate/hertz)
    # layers_PRM_rate.append(layers_simulation_PRM8.rate/hertz)
    # layers_PRM_rate.append(layers_simulation_PRM9.rate/hertz)
    # layers_PRM_time.append(layers_simulation_PRM0.t/ms)
    # layers_PRM_time.append(layers_simulation_PRM1.t/ms)
    # layers_PRM_time.append(layers_simulation_PRM2.t/ms)
    # layers_PRM_time.append(layers_simulation_PRM3.t/ms)
    # layers_PRM_time.append(layers_simulation_PRM4.t/ms)
    # layers_PRM_time.append(layers_simulation_PRM5.t/ms)
    # layers_PRM_time.append(layers_simulation_PRM6.t/ms)
    # layers_PRM_time.append(layers_simulation_PRM7.t/ms)
    # layers_PRM_time.append(layers_simulation_PRM8.t/ms)
    # layers_PRM_time.append(layers_simulation_PRM9.t/ms)
    # stats['layers PRM rate'] = layers_PRM_rate
    # stats['layers PRM time'] = layers_PRM_time

    stats['spikemon times'] = simulation_spikemon.t/ms
    stats['spikemon indices'] = simulation_spikemon.i/1
    stats['average firing rate'] = simulation_spikemon.num_spikes/(N*simulationtime/second)

    events, simulation_time = ar.calculate_events(N, stats, neuron_bin_size) # numpy array of tuples representing events
    stats['events'] = events
    stats['num events'] = len(events)
    print("\nnumber of events: {0}\n".format(len(events)))
    sys.stdout.flush()

    (event_rate, event_mag, IEIs, excess_kurtosis, skew) = ar.analyze_events(N, events, simulationtime/ms, neuron_bin_size)
    stats['event_rate'] = event_rate
    stats['event_mag'] = event_mag
    stats['IEIs'] = IEIs 
    stats['IEI excess_kurtosis'] = excess_kurtosis
    stats['IEI skew'] = skew

    # save results (pickle new stats dictionary)
    style = "numLay{0}_".format(num_layers)
    ar.save_results(N, p, index, stats, style, data_dir)
    ar.clean_results(N, p, index, style, data_dir)
    
    # print the stats for W
    for key,value in sorted(stats.items()):
        if not isinstance(value,np.ndarray) and not isinstance(value,list):
            print(key+":{0}".format(value))
            sys.stdout.flush()
    print("\n")
    
    try:
        # plot the results of the simulation
        figure(figsize=(20,10))
        try: 
            numPRMshow = 3
            for index in range(num_layers):
                if index < numPRMshow:
                    subplot(numPRMshow+1,1,index+2)
                    # plot(layers_simulation_PRM[index].t/ms,layers_simulation_PRM[index].smooth_rate(window='flat', width=0.5*ms)/Hz) # smooth curve of how many neurons fire at each time
                    plot(layers_PRM_time[num_layers-index-1],layers_PRM_rate[num_layers-index-1]) # smooth curve of how many neurons fire at each time
        
        except Exception as e:
            numPRMshow = 1
            subplot(numPRMshow+1, 1, 2)
            plot(layer_simulation_PRM_last.t/ms, layer_simulation_PRM_last.rate/hertz)
            xlabel('Time (ms)')
            ylabel('Last layer PRM rate')
            # plot(simulation_PRM.t/ms, simulation_PRM.rate/hertz)
            # xlabel('Time (ms)')
            # ylabel('Network PRM rate')

            # print("Error: %s" %e)
            # sys.stdout.flush()
          
        subplot(numPRMshow+1, 1, 1)
        plot(simulation_spikemon.t/ms,simulation_spikemon.i, '.k') # raster plot: y-axis = neuron index, x-axis = time, dot at (t,i) if neuron i fires at time t
        xlabel('Time (ms)')
        ylabel('Neuron index')
        plt.tight_layout()
        # show()
    except Exception as e:
        # raise e
        print("Error: %s" % e)

    #delete monitors so they don't cause restore() to get an error when looped through    
    del simulation_statemon 
    del simulation_spikemon 
    del simulation_PRM 
    del layer_simulation_PRM_last
    # del layers_simulation_PRM0
    # del layers_simulation_PRM1
    # del layers_simulation_PRM2
    # del layers_simulation_PRM3
    # del layers_simulation_PRM4
    # del layers_simulation_PRM5
    # del layers_simulation_PRM6
    # del layers_simulation_PRM7
    # del layers_simulation_PRM8
    # del layers_simulation_PRM9

    # save results (pickle new stats dictionary)
    style = "numLay{0}_".format(num_layers)
    ar.save_results(N, p, index, stats, style, data_dir)
    ar.clean_results(N, p, index, style, data_dir)

show()