# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 09:46:52 2017

@author: Brittany
"""

# functions in this script are:
# load_results(N, p, i, data_dir='matrices/')
#     -- returns a python dict: results
# save_results(N, p, i, results, data_dir='matrices/')
#     -- saves results as pickle file (ovewrites existing pickle file if one exists), 
#     -- returns nothing
# create_subPR(N, results, neuron_bin_size=100, num_neuron_bins=100) 
#     -- returns numpy array: subPR; float: time_bin_size; float: simulation_time
# get_thresholds(subPR, num_neuron_bins) 
#     -- returns numpy array: thresholds
# get_events(subPR, thresholds, num_neuron_bins, time_bin_size = .1, consecutive_count = 1, 
#             group_spacing = 1, consecutive_bin = 1) 
#     -- returns numpy array of event tuples
# calculate_events(N, results, neuron_bin_size=100)
#     -- calls other funtions in this script and 
#     -- returns numpy array of event tuples; float: simulation_time 
# analyze_events(N, events, simulation_time, neuron_bin_size=100)
#     -- returns event_rate, event_mag, IEIs, excess_kurtosis(of IEIs), and skew (of IEIs)
# update_results(N, p, i, data_dir='matrices/', neuron_bin_size=100, num_neuron_bins=100, 
#                  time_bin_size = .1, consecutive_count = 1, group_spacing = 1, consecutive_bin = 1)
#     -- calls load_results(), calculate_events(), analyze_events(), and save_results()
#     -- updates the results pickle file of a given matrix with events and outputs from analyze_events
#     -- also prints new results
#
# to reimport a module, use importlib.reload(module), you must first import importlib
#
#####

#######################################################################################################
def load_results(N, p, i, style, data_dir='matrices/'):
  import pickle

  result_filename = "{0}Results_W_N{1}_p{2}_{3}{4}.pickle".format(data_dir,N,p,style,i) 
  with open(result_filename, "rb") as rf:
    results = pickle.load(rf)

  return(results)

#######################################################################################################
def save_results(N, p, i, results, style, data_dir='matrices/'):
  import pickle

  result_filename = "{0}Results_W_N{1}_p{2}_{3}{4}.pickle".format(data_dir,N,p,style,i) 
  with open(result_filename, "wb") as rf:
    pickle.dump(results, rf)

#######################################################################################################
def clean_results(N, p, i, style, data_dir='matrices/'):
  import pickle
  import numpy as np

  results = load_results(N, p, i, style, data_dir)
  cleanResults = dict()

  for key,value in results.items():
        if not isinstance(value,np.ndarray):
            cleanResults[key] = results[key]
  try:  
    cleanResults['largest eigenvalue'] = results['largest eigenvalue']
  except:
    pass

  save_results(N, p, i, cleanResults, style+"Clean_", data_dir)


###############################################################################################################
def update_results(N, p, i, style, data_dir='matrices/', neuron_bin_size=100):
  import numpy as np
  import sys
  
  results = load_results(N, p, i, style, data_dir)
  events, simulation_time = calculate_events(N, results, neuron_bin_size)
  try:
    simulation_time = results['simulation_time']
  except:
    pass

  event_rate, event_mag, IEIs, excess_kurtosis, skew = analyze_events(N, events, simulation_time, neuron_bin_size)

  results['event_rate'] = event_rate
  results['event_mag'] = event_mag
  results['IEIs'] = IEIs 
  results['IEI excess_kurtosis'] = excess_kurtosis
  results['IEI skew'] = skew

  save_results(N, p, i, results, data_dir)

  for k,v in sorted(results.items()):
    if not isinstance(v,np.ndarray):
      print(k+":{0}".format(v))
  print("\n")
  sys.stdout.flush()

#######################################################################################################
#######################################################################################################
def create_subPR(results, neuron_bin_size=100, num_neuron_bins=100):
  import numpy as np
  import math

  time_bin_size = 5*(results['PRM time'][1] - results['PRM time'][0]) # return this too, to save for future use
  num_time_bins = int(math.ceil(len(results['PRM time'])/5))+1
  simulation_time = results['PRM time'][-1] - results['PRM time'][0]

  subPR = np.zeros((num_neuron_bins, num_time_bins)) # subPR is a num_neuron_bin x num_time_bin matrix

  for n in range(len(results['spikemon indices'])):
    neuron_bin_index = int(math.floor((results['spikemon indices'][n])/neuron_bin_size))
    time_bin_index = int(np.round((results['spikemon times'][n] - results['spikemon times'][0])
      /time_bin_size))
    subPR[neuron_bin_index, time_bin_index] += 1 
    # subPR[i,j] = the number of neurons in bin i that fired in time bin j

  scale_factor = np.multiply(neuron_bin_size, time_bin_size) # subPR will be devided by this to scale appropriately

  return(np.true_divide(subPR, scale_factor), time_bin_size, simulation_time)



#######################################################################################################
def get_thresholds(subPR, num_neuron_bins):
  import numpy as np
  import math

  tempThresh = np.zeros(num_neuron_bins) # each neuron bin will have a threshold value

  for i in range(num_neuron_bins):
    std = np.std(subPR[i])
    median = np.median(subPR[i])
    tempThresh[i] = median + (6*std) # set the tempThresh to be 6 standard deviations above the median subPR value for each neuron bin
 
  tempSubPR = np.minimum(subPR, tempThresh.reshape((num_neuron_bins,1)))
  # if subPR[i,j] > tempThresh[i], then tempSubPR[i,j] = tempThresh[i], otherwise tempSubPR[i,j] = subPR[i,j]
  # This is creating a new subPR matrix which replaces the extremely high values in subPR with the tempThresh value for that neuron bin

  thresholds = np.zeros(num_neuron_bins)
  for i in range(num_neuron_bins):
    std = np.std(tempSubPR[i])
    median = np.median(tempSubPR[i])
    thresholds[i] = median + (15*std) # the actual threshold values are computed using the new tempSubPR matrix

  return(thresholds.reshape((num_neuron_bins,1)))


#######################################################################################################
def get_events(N, subPR, thresholds, num_neuron_bins, time_bin_size, 
  consecutive_time = 1, time_spacing = 1, consecutive_bin = 2):
  import numpy as np
  import math

  above_lower_thresh = np.greater_equal(subPR, thresholds) # matrix
  above_higher_thresh = np.greater_equal(subPR, 3*thresholds) # matrix

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

  dtype = [('start_neuron_bin', int), ('end_neuron_bin', int), ('start_time', float), 
    ('end_time', float)] # these are labels for the event tuples
  sorted_events_by_bin_array = np.sort(np.array(events_by_bin, dtype=dtype), 
    order=['start_time', 'start_neuron_bin'])
  #return(sorted_events_by_bin_array)

  events_list = [] # each entry will be tuple: 
  #  (start_neuron_bin, end_neuron_bin, start_time, end_time)
  # starting_events_list = []

  for newe in sorted_events_by_bin_array:
    used = False
    for event in events_list:
      if (newe['start_time'] >= event[2]) and (newe['start_time'] <= (event[3]+time_spacing)):
        if (event[0] <= event[1]) and (newe['start_neuron_bin'] == (event[1]+1)): 
          # event is moving forward through neuron bins
          if used:
            print("error: new event fits with multiple existing events")
          else: # update event in events_list
            event[1] = newe['end_neuron_bin'] # new end_neuron_bin
            event[3] = newe['end_time']
            used = True
            break

        elif (event[0] >= event[1]) and (newe['start_neuron_bin'] == (event[1]-1)): 
          # event is moving backwards through neuron bins
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
      # starting_events_list.append(newe) # record info about first event in chain
      used = True

  final_events_list = []
  for event in events_list: 
    # now remove all events that cover fewer than a certain number (consecutive_bin) of neuron bins
    if abs(event[1] - event[0]) >= consecutive_bin:
      final_events_list.append(event)

  events = np.array(final_events_list, dtype=dtype)
  return(events)


####################################################################################################### 
# Inputs for the calculate_events function:
#   N: number of neurons in the network
#   results: a python dictionary containing "PRM time", "spikemon times", "spikemon indices"
#   neuron_bin_size: how many neurons to use per neuron bin (default = 100)
#
# Uses functions: create_subPR, get_thresholds, and get_events
#
# Returns a numpy array/list of event "objects" as tuples:
#   (start_neuron_bin, end_neuron_bin, start_time, end_time)
#########################################################################
def calculate_events(N, results, neuron_bin_size=100):
  import numpy as np
  import math

  num_neuron_bins = math.ceil(N/neuron_bin_size)

  subPR, time_bin_size, simulation_time = create_subPR(results, neuron_bin_size, num_neuron_bins)
  thresholds = get_thresholds(subPR, num_neuron_bins)
  events = get_events(N, subPR, thresholds, num_neuron_bins, time_bin_size) 

  return(events, simulation_time) # a numpy array of tuples


###############################################################################################################
def analyze_events(N, events, simulation_time=3000, neuron_bin_size=100):
  import numpy as np
  from scipy import stats

  event_rate = (len(events)/simulation_time)*1000 # number of events per second (simulation time units: ms)

  events_length = 0 # will be the total number of neurons that are included in all events

  IEIs = [] # will be a list of inter-event intervals (times)
  stime = None # we don't care about the time before the first event

  for event in events:
    events_length += neuron_bin_size*(abs(event['end_neuron_bin'] - event['start_neuron_bin']) + 1) 
    # update num neurons covered by event
    
    if stime is not None:
      # update IEIs list with time beteween previous event and current event
      IEIs.append(event['start_time']-stime) 
    
    stime = event['start_time'] # update with start time of current event

  event_mag = events_length/(N*simulation_time) # normalized magnitude of events 
  # (number of neurons covered by events, proportional to num neurons in network and simulation time)

  try:
    excess_kurtosis = stats.kurtosis(IEIs, bias=False)
    skew = stats.skew(IEIs, bias=False)
  except:
    excess_kurtosis = float('nan')
    skew = float('nan')

  return(event_rate, event_mag, IEIs, excess_kurtosis, skew)