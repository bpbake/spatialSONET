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

  print('Matrix {0}'.format(i))
  for k,v in sorted(results.items()):
    if not isinstance(v,np.ndarray):
      print(k+":{0}".format(v))
  print("\n")

  num_neuron_bins = math.ceil(N/neuron_bin_size)
  time_bin_size = results['PRM time'][1] - results['PRM time'][0]
  num_time_bins = len(results['PRM time'])

  subPR = np.zeros((num_neuron_bins, num_time_bins))

  for n in range(len(results['spikemon indices'])):
    neuron_bin_index = int((results['spikemon indices'][n])//neuron_bin_size)
    time_bin_index = int(10*(results['spikemon times'][n] - results['spikemon times'][0]))
    subPR[neuron_bin_index, time_bin_index] += 1

  scale_factor = np.multiply(neuron_bin_size, time_bin_size)

  return(np.true_divide(subPR, scale_factor))



###############################################################################################################
def get_thresholds(N, subPR, neuron_bin_size=100):
  import numpy as np
  import math

  num_neuron_bins = math.ceil(N/neuron_bin_size)
  thresholds = np.zeros(num_neuron_bins)

  for i in range(num_neuron_bins):
    #thresholds[i] = np.percentile(subPR[i],90)
    std = np.std(subPR[i])
    mean = np.median(subPR[i])
    thresholds[i] = mean + (2.5*std)

  return(thresholds.reshape((num_neuron_bins,1)))


###############################################################################################################
def get_events(N, subPR, thresholds, consecutive_count = 1, group_spacing = 1, neuron_bin_size = 100):
  import numpy as np
  import math

  #thresh_matrix = np.broadcast_to(thresholds, np.shape(subPR))
  above_thresh = np.greater_equal(subPR, thresholds)#_matrix)

  events_by_bin = [] # a list of event tuples

  num_neuron_bins = math.ceil(N/neuron_bin_size)
  for i in range(num_neuron_bins):
    count = 0
    for j in range(len(subPR[i])):
      if above_thresh[i,j]: # add to count (number consecutive)
        count += 1
      elif count >= consecutive_count: # save time and bin, reset count
        events_by_bin.append((i, i, (j-count)*.1, j*.1))
        count = 0

  dtype = [('start_neuron_bin', int), ('end_neuron_bin', int), ('start_time', float), ('end_time', float)] # these are labels for the event tuples
  sorted_events_by_bin_array = np.sort(np.array(events_by_bin, dtype=dtype), order=['start_time', 'start_neuron_bin'])

  events_list = [(-1, -1, -1, -1)] # each entry is (start_neuron_bin, end_neuron_bin, start_time, end_time)
  starting_events_list = []

  for newe in sorted_events_by_bin_array:
    used = False
    for event in events_list:
      #while used == False:
      # if (newe['start_time'] >= event['start_time']) and (newe['start_time'] <= (event['end_time']+group_spacing)):
      #   if (event['start_neuron_bin'] <= event['end_neuron_bin']) and (newe['start_neuron_bin'] == (event['end_neuron_bin']+1)):
      if (newe['start_time'] >= event[2]) and (newe['start_time'] <= (event[3]+group_spacing)):
        if (event[0] <= event[1]) and (newe['start_neuron_bin'] == (event[1]+1)): # event is moving forward through neuron bins
          # update event in events_list
          if used:
            print("error: new event fits with multiple existing events")
          else:
            event[1] = newe['end_neuron_bin'] # new end_neuron_bin
            event[3] = newe['end_time']
            used = True

        # elif (event['start_neuron_bin'] >= event['end_neuron_bin']) and (newe['start_neuron_bin'] == (event['end_neuron_bin']-1)):
        elif (event[0] >= event[1]) and (newe['start_neuron_bin'] == (event[1]-1)): # event is moving backwards through neuron bins
          # update event in events_list
          if used:
            print("error: new event fits with multiple existing events")
          else:
            event[1] = newe['end_neuron_bin'] # new end_neuron_bin
            event[3] = newe['end_time']
            used = True

    if not used: # if it doesn't work with any current event in events_list, add it to the list
      events_list.append(newe)
      events_list.append(newe) # add twice because it could travel in both directions
      starting_events_list.append(newe) # record info about first event in chain
      used = True

      if events_list[0] == (-1, -1, -1, -1):
        events_list.remove((-1, -1, -1, -1)) 

  events = np.array(events_list, dtype=dtype)

  return(events)