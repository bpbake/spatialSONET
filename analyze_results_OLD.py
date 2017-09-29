# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 09:46:52 2017

@author: Brittany
"""

# functions in this script are:
# plot_results(N,p,i, data_dir='matrices/') 
#     -- plots raster plot and PRM for matrix i
# create_subPR(N,p,i, neuron_bin_size=100, data_dir='matrices/') 
#     -- returns numpy array: subPR and float: time_bin_size
# get_thresholds(N, subPR, neuron_bin_size=100) 
#     -- returns numpy array: thresholds
# get_events(N, subPR, thresholds, neuron_bin_size = 100, time_bin_size = .1, consecutive_count = 1, group_spacing = 1, consecutive_bin = 1) 
#     -- returns numpy array of event "objects"
# plot_subPR_thresh(nbin, subPR, thresholds) 
#     -- plots subPR and threshold value for a given neuron bin
#
# to reimport a module, use importlib.reload(module), you must first import importlib
#
#####


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


###############################################################################################################
def create_subPR(N,p,i, neuron_bin_size=100, data_dir='matrices/'):
  import numpy as np
  import matplotlib
  import matplotlib.pyplot as plt
  import pickle
  import math

  
  result_filename = "{0}Results_W_N{1}_p{2}_slower{3}.pickle".format(data_dir,N,p,i) 
  with open(result_filename, "rb") as rf:
   	results = pickle.load(rf)

  # print('Matrix {0}'.format(i))
  # for k,v in sorted(results.items()):
  #   if not isinstance(v,np.ndarray):
  #     print(k+":{0}".format(v))
  # print("\n")

  num_neuron_bins = math.ceil(N/neuron_bin_size)
  time_bin_size = results['PRM time'][1] - results['PRM time'][0] # return this too, to save for future use
  num_time_bins = len(results['PRM time'])

  subPR = np.zeros((num_neuron_bins, num_time_bins))

  for n in range(len(results['spikemon indices'])):
    neuron_bin_index = int(math.floor((results['spikemon indices'][n])/neuron_bin_size))
    time_bin_index = int(np.round((results['spikemon times'][n] - results['spikemon times'][0])/time_bin_size))
    subPR[neuron_bin_index, time_bin_index] += 1

  scale_factor = np.multiply(neuron_bin_size, time_bin_size)

  return(np.true_divide(subPR, scale_factor), time_bin_size)



###############################################################################################################
def get_thresholds(N, subPR, neuron_bin_size=100):
  import numpy as np
  import math

  num_neuron_bins = math.ceil(N/neuron_bin_size)
  tempThresh = np.zeros(num_neuron_bins)

  for i in range(num_neuron_bins):
    std = np.std(subPR[i])
    median = np.median(subPR[i])
    tempThresh[i] = median + (6*std)
 
  tempSubPR = np.minimum(subPR, tempThresh.reshape((num_neuron_bins,1)))

  thresholds = np.zeros(num_neuron_bins)
  for i in range(num_neuron_bins):
    std = np.std(tempSubPR[i])
    median = np.median(tempSubPR[i])
    thresholds[i] = median + (15*std)

  return(thresholds.reshape((num_neuron_bins,1)))


###############################################################################################################
def get_events(N, subPR, thresholds, neuron_bin_size = 100, time_bin_size = .1, consecutive_time = 1, time_spacing = 1, consecutive_bin = 1):
  import numpy as np
  import math

  above_lower_thresh = np.greater_equal(subPR, thresholds)#_matrix)
  above_higher_thresh = np.greater_equal(subPR, 3*thresholds)#_matrix)

  events_by_bin = [] # a list of event tuples

  num_neuron_bins = math.ceil(N/neuron_bin_size)
  for i in range(num_neuron_bins):
    count = 0
    success = False
    for j in range(len(subPR[i])):
      if above_lower_thresh[i,j]: # add to count (number consecutive)
        count += 1
        if count >= consecutive_time or above_higher_thresh[i,j]:
          success = True
      else:
        if success: # save time and bin, reset count
          events_by_bin.append((i, i, (j-count)*time_bin_size, (j-1)*time_bin_size))
        count = 0
        success = False

  #return(events_by_bin)

  dtype = [('start_neuron_bin', int), ('end_neuron_bin', int), ('start_time', float), ('end_time', float)] # these are labels for the event tuples
  sorted_events_by_bin_array = np.sort(np.array(events_by_bin, dtype=dtype), order=['start_time', 'start_neuron_bin'])
  #return(sorted_events_by_bin_array)

  events_list = [] # each entry will be tuple: (start_neuron_bin, end_neuron_bin, start_time, end_time)
  starting_events_list = []

  for newe in sorted_events_by_bin_array:
    used = False
    for event in events_list:
      if (newe['start_time'] >= event[2]) and (newe['start_time'] <= (event[3]+time_spacing)):
        if (event[0] <= event[1]) and (newe['start_neuron_bin'] == (event[1]+1)): # event is moving forward through neuron bins
          if used:
            print("error: new event fits with multiple existing events")
          else: # update event in events_list
            event[1] = newe['end_neuron_bin'] # new end_neuron_bin
            event[3] = newe['end_time']
            used = True
            break

        elif (event[0] >= event[1]) and (newe['start_neuron_bin'] == (event[1]-1)): # event is moving backwards through neuron bins
          if used:
            print("error: new event fits with multiple existing events")
          else: # update event in events_list
            event[1] = newe['end_neuron_bin'] # new end_neuron_bin
            event[3] = newe['end_time']
            used = True
            break

    if not used: # if it doesn't work with any current event in events_list, add it to the list
      events_list.append(newe)
      events_list.append(newe.copy()) # add twice because it could travel in both directions
      starting_events_list.append(newe) # record info about first event in chain
      used = True

  for event in events_list: # now remove all events that cover fewer than a certain number (consecutive_bin) of neuron bins
    if abs(event[0] - event[1]) <= consecutive_bin:
      events_list.remove(event) 

  events = np.array(events_list, dtype=dtype)
  return(events)


###############################################################################################################
def plot_subPR_thresh(nbin, subPR, thresholds): # hard coded ranges for specific scenerios
    import numpy as np
    import matplotlib.pyplot as plt

    plt.plot(subPR[nbin])
    plt.hold(True)
    plt.plot([0,30000], [thresholds[nbin], thresholds[nbin]])
    plt.plot([0,30000], [3*thresholds[nbin], 3*thresholds[nbin]])
    plt.hold(False)