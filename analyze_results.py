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
  simulation_time = results['simulation_time']

  events = calculate_events(N, results, neuron_bin_size)
  results['events'] = events
  results['num events'] = len(events)

  event_rate, event_mag, IEIs, excess_kurtosis, skew = analyze_events(N, events, simulation_time, neuron_bin_size)

  results['event_rate'] = event_rate
  results['event_mag'] = event_mag
  results['IEIs'] = IEIs 
  results['IEI excess_kurtosis'] = excess_kurtosis
  results['IEI skew'] = skew

  save_results(N, p, i, results, style, data_dir)

  # for k,v in sorted(results.items()):
  #   if not isinstance(v,np.ndarray):
  #     print(k+":{0}".format(v))
  # print("\n")
  # sys.stdout.flush()

#######################################################################################################
#######################################################################################################
## called by calculate_events
def create_subPR(results, neuron_bin_size=100, num_neuron_bins=100): #, time_units_per_bin=1):
  import numpy as np
  import math

  time_bin_size = (results['PRM time'][1] - results['PRM time'][0]) # return this too, to save for future use
  num_time_bins = len(results['PRM time'])
  # time_bin_size = time_units_per_bin*(results['PRM time'][1] - results['PRM time'][0]) # return this too, to save for future use
  # num_time_bins = int(math.ceil(len(results['PRM time'])/time_units_per_bin))
  # simulation_time = results['PRM time'][-1] - results['PRM time'][0]
  # print("PRM simulation time: {0}".format(simulation_time))

  subPR = np.zeros((num_neuron_bins, num_time_bins)) # subPR is a num_neuron_bin x num_time_bin matrix

  for n in range(len(results['spikemon indices'])):
    neuron_bin_index = int(math.floor((results['spikemon indices'][n])/neuron_bin_size))
    time_bin_index = int(np.floor((results['spikemon times'][n] - results['spikemon times'][0])
      /time_bin_size))
    subPR[neuron_bin_index, time_bin_index] += 1 
    # subPR[i,j] = the number of neurons in bin i that fired in time bin j

  scale_factor = np.multiply(neuron_bin_size, time_bin_size) # subPR will be devided by this to scale appropriately

  return(np.true_divide(subPR, scale_factor), time_bin_size)#, simulation_time)



#######################################################################################################
## called by calculate_events
def get_thresholds(subPR, num_neuron_bins):
  import numpy as np
  import math

  tempThresh = np.zeros(num_neuron_bins) # each neuron bin will have a threshold value

  for i in range(num_neuron_bins):
    std = np.std(subPR[i])
    median = np.median(subPR[i])
    # tempThresh[i] = median + (5*std) 
    # set the tempThresh to be 6 standard deviations above the median subPR value for each neuron bin
    # tempThresh[i] = np.percentile(subPR[i],90)
    for temp in range(80, 100, 1):
      tempThresh[i] = np.percentile(subPR[i],temp)
      if tempThresh[i] > 0:
        break
 
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
## called by calculate_events
def get_events(N, subPR, thresholds, num_neuron_bins=1, time_bin_size=.1, 
  consecutive_time = 1, time_spacing = 1, consecutive_bin = 2): 
  import numpy as np
  import math

  above_lower_thresh = np.greater_equal(subPR, thresholds) # matrix
  above_higher_thresh = np.greater_equal(subPR, 3*thresholds) # matrix

  events_by_bin = [] # a list of event tuples

  for i in range(num_neuron_bins):
    count = 0
    success = False
    for j in range(len(subPR[i])): # subPR is a num_neuron_bin x num_time_bin matrix of "population rates"
      if above_lower_thresh[i,j]: # add to count (number consecutive time bins above thresh)
        count += 1 # count is the number of time bins that the event takes up (for one neuron bin)
        if count >= consecutive_time or above_higher_thresh[i,j]:
          success = True # this really is an event... 
          # keep increasing the time bin until the end of the event
          # then record this as an event
      else:
        if success: # save time and bin, reset count
          events_by_bin.append((i, i, (j-count)*time_bin_size, (j-1)*time_bin_size, None))
          # append event tuple: (start_neuron_bin, end_neuron_bin, start_time, end_time, direction)
        count = 0
        success = False

  # print(events_by_bin)

  dtype = [('start_neuron_bin', int), ('end_neuron_bin', int), ('start_time', float), 
    ('end_time', float), ('direction', str)] # these are labels for the event tuples
  sorted_events_by_bin_array = np.sort(np.array(events_by_bin, dtype=dtype), 
    order=['start_time', 'start_neuron_bin']) # sort events_by_bin by start time (then by start neuron bin)
  # print(sorted_events_by_bin_array[0:100])
  print("dtype: {0}".format(dtype))

  events_list = [] # each entry will be tuple: 
  #  (start_neuron_bin, end_neuron_bin, start_time, end_time, direction)
  # starting_events_list = []

  num_events = 0

  for newe in sorted_events_by_bin_array:
    used = False
    for event in events_list:
      # if (newe['start_time'] >= event[2]) and (newe['start_time'] <= (event[3]+(time_spacing*time_bin_size))):
      if (newe['start_time'] >= event[2]) and (newe['start_time'] <= (event[3]+time_spacing)): 
        # if newe starts not before start time of existing event
        # and if newe starts before end time of existing event, 
        # then concatinate the events:

        # if (event[0] <= event[1]) and (newe['start_neuron_bin'] == (event[1]+1)): 
        if (event[4] is 'up') and ((newe['start_neuron_bin'] == (event[1]+1)) or ((newe['start_neuron_bin']==0) and (event[1]==num_neuron_bins-1))):
          # if existing event start neuron bin is before end neuron bin
          # and if newe start neuron bin is the bin after the existing event end neuron bin,
          # then event is moving forward through neuron bins (perhaps wrap-around)
          if used:
            print("error: new event fits with multiple existing events")
          else: # update event in events_list
            event[1] = newe['end_neuron_bin'] # new end_neuron_bin
            event[3] = newe['end_time']
            # event[4] = 'up'
            used = True # newe has been added to the event list (concatinated with an existing event)
            break

        # elif (event[0] >= event[1]) and (newe['start_neuron_bin'] == (event[1]-1)): 
        elif (event[4] is 'down') and ((newe['start_neuron_bin'] == (event[1]-1)) or ((newe['start_neuron_bin']==num_neuron_bins-1) and (event[1]==0))): 
          # event is moving backwards through neuron bins (maybe wrap-around)
          if used:
            print("error: new event fits with multiple existing events")
          else: # update event in events_list
            event[1] = newe['end_neuron_bin'] # new end_neuron_bin
            event[3] = newe['end_time']
            # event[4] = 'down'
            used = True # newe has now been added to the event list (concatinated with an existing event)
            break

    if not used: # if it doesn't work with any current event in events_list, add it to the list
      newe['direction'] = 'up'
      events_list.append(newe)
      print("newe: {0}".format(newe))
      copy_newe = newe.copy()
      copy_newe['direction'] = 'down'
      events_list.append(copy_newe) # add twice because it could travel in both directions
      print("newe.copy: {0}".format(copy_newe))
      num_events += 1 # but it's only one event, so add only once
      # starting_events_list.append(newe) # record info about first event in chain
      used = True # newe has now been added to the event list

  print("events list created: {0}".format(events_list))

  sorted_events_list = np.sort(np.array(events_list, dtype=dtype), order=['start_time', 'start_neuron_bin'])
  print("events list sorted, but needs to be cleaned by event length")

  final_events_list = []
  prev_start_bin = None
  prev_start_time = None
  for event in sorted_events_list: 
    # now remove all events that cover fewer than a certain number (consecutive_bin) of neuron bins
    if abs(event['start_neuron_bin'] - event['end_neuron_bin']) >= consecutive_bin:
      final_events_list.append(event) 
    elif (event['start_neuron_bin'] == prev_start_bin) and  (event['start_time'] == prev_start_time):
      num_events -= 1
      # if both directions of an event are not long enough, then we removed the entire event
      # this will work because the events are sorted by start_time then by start_neuron_bin
    prev_start_bin = event['start_neuron_bin']
    prev_start_time = event['start_time']

  events = np.array(final_events_list, dtype=dtype)
  return(events, num_events)


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

  subPR, time_bin_size = create_subPR(results, neuron_bin_size, num_neuron_bins)
  thresholds = get_thresholds(subPR, num_neuron_bins)
  events, num_events = get_events(N, subPR, thresholds, num_neuron_bins, time_bin_size) 
  print("events calculated")

  return(events, num_events) # a numpy array of tuples


###############################################################################################################
def analyze_events(N, events, num_events, simulation_time, neuron_bin_size=100):
  import numpy as np
  from scipy import stats

  event_rate = (num_events/simulation_time)*1000 # number of events per second (simulation time units: ms)

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


###############################################################################################################
###############################################################################################################
def new_measure(N, results, neuron_bin_size=100, num_time_bins=300, simulation_time=3000):
  import numpy as np
  import math
  from scipy import linalg

  num_neuron_bins = math.ceil(N/neuron_bin_size)
  simulation_time = results['PRM time'][-1] - results['PRM time'][0]
  time_bin_size = math.ceil(simulation_time/num_time_bins)

  R_spikes = np.zeros((N, num_time_bins))
  # R_spikes = np.zeros((num_neuron_bins, num_time_bins))

  for n in range(len(results['spikemon indices'])):
    neuron_index = int(results['spikemon indices'][n])
    neuron_bin_index = int(math.floor(neuron_index/neuron_bin_size))
    time_bin_index = int(np.floor((results['spikemon times'][n] - results['spikemon times'][0])/time_bin_size))
    R_spikes[neuron_index, time_bin_index] = 1
    # R_spikes[neuron_bin_index, time_bin_index] = 1

  c = 1 # max "speed" of spike wave from neuron layer j to i (number of time bins?)
  # constant for all i,j pairs?  Maybe need to make c an array (matrix)???

  dist = linalg.toeplitz(np.floor(np.linspace(0,num_neuron_bins,num=N,endpoint=False)).astype(int)) 
  # dist is an N by N matrix with entry (i,j) being the distance between the layers of neurons i and j

  k = np.multiply(c,dist)

  cov_R=np.zeros((N,N))
  # MUST Figure out a better way to do this triple loop!!!
  for i in range(N):
    for j in range(N):
      for t in range(num_time_bins):
        if t+k[i,j] < num_time_bins:
          cov_R[i,j] += R_spikes[i,t]*R_spikes[j,(t+k[i,j])]

  return np.mean(cov_R)