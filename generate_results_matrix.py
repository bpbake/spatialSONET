# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 00:02:45 2017

@author: rhino
"""
try:
   import cPickle as pickle # used to store python data types
except:
   import pickle
   
import numpy as np

try:
    del input # brian overwrites this, we want to reset it to the python default
except:
    pass
input_orig = input

data_dir = 'matrices/N10000_LL70_LR0_ff/'

N = 10000
p_AVG = 50/N

results_header = ['w_index', 'N', 'L_left', 'L_right', #'largest eigenvalue', 
'alpha_chain', 'alpha_chain_hat', 'alpha_conv', 'alpha_conv_hat', 'alpha_div', 'alpha_div_hat', 
'alpha_recip', 'alpha_recip_hat', 'p_AVG', 'p_hat', 
'simulation_time', 'event_rate', 'event_mag', 'IEI excess_kurtosis', 'IEI skew']

results_type = [('w_index', int), ('N', int), ('L_left', float), ('L_right', float), #('largest eigenvalue', complex list), 
('alpha_chain', float), ('alpha_chain_hat', float), ('alpha_conv', float), ('alpha_conv_hat', float), 
('alpha_div', float), ('alpha_div_hat', float), ('alpha_recip', float), ('alpha_recip_hat', float), 
('p_AVG', float), ('p_hat', float), ('simulation_time', float), #('IEIs', list),
('event_rate', float), ('event_mag', float), ('IEI excess_kurtosis', float), ('IEI skew', float)]

start_index = int(input_orig("enter a starting index: "))
end_index = int(input_orig("enter end index: "))

results_matrix = np.zeros((end_index-start_index+1, len(results_header)))
result_list = []
                                                  
for w_index in range(start_index, end_index+1):
   
    results_filename = "{0}Results_W_N{1}_p{2}_slower{3}.pickle".format(data_dir,N,p_AVG,w_index)
    with open(results_filename, 'rb') as rf:
        try:
            results = pickle.load(rf) 
        except (EOFError):
            print("unpickling error")
          
    results_matrix[w_index-start_index, 0] = w_index
    result = [w_index]

    for k in range(1, len(results_header)):    
        results_matrix[w_index-start_index, k] = results[results_header[k]]
        result.append(results[results_header[k]])

    result_list.append(result)


r_filename = "{0}Result_Matrices_{1}_{2}.pickle".format(data_dir, start_index, end_index) 
with open(r_filename, "wb") as rf:
    pickle.dump(np.array(result_list, results_type), rf)

np.savetxt('{0}Result_Matrices_{1}_{2}.csv'.format(data_dir, start_index, end_index), 
    results_matrix, delimiter=',', header=str(results_header)) 