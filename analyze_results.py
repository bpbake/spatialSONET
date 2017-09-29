# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 09:46:52 2017

@author: Brittany
"""

# functions in this script are:
# create_subPR(N, results, neuron_bin_size=100) 
#     -- returns numpy array: subPR and float: time_bin_size
# get_thresholds(subPR, num_neuron_bins) 
#     -- returns numpy array: thresholds
# get_events(subPR, thresholds, num_neuron_bins, time_bin_size = .1, consecutive_count = 1, group_spacing = 1, consecutive_bin = 1) 
#     -- returns numpy array of event "objects"
#
# to reimport a module, use importlib.reload(module), you must first import importlib
#
#####

###############################################################################################################
def create_subPR(results, neuron_bin_size=100, num_neuron_bins=100):
  import numpy as np
  import math

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
def get_thresholds(subPR, num_neuron_bins):
  import numpy as np
  import math

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
def get_events(N, subPR, thresholds, num_neuron_bins, time_bin_size, consecutive_time = 1, time_spacing = 1, consecutive_bin = 1):
  import numpy as np
  import math

  above_lower_thresh = np.greater_equal(subPR, thresholds)#_matrix)
  above_higher_thresh = np.greater_equal(subPR, 3*thresholds)#_matrix)

  events_by_bin = [] # a list of event tuples

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
def calculate_events(N, results, neuron_bin_size=100):
  import numpy as np
  import math

  num_neuron_bins = math.ceil(N/neuron_bin_size)

  (subPR, time_bin_size) = create_subPR(results, neuron_bin_size, num_neuron_bins)
  thresholds = get_thresholds(subPR, num_neuron_bins)
  events = get_events(subPR, thresholds, num_neuron_bins, time_bin_size)

  return(events)


###############################################################################################################
def analyze_events(events, ):