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

start_index = int(input_orig("enter a starting index: "))
end_index = int(input_orig("enter end index: "))

result_matrix = np.zeros((end_index-start_index+1, len(results_header)))
                                                  
for w_index in range(start_index, end_index+1):
   
    results_filename = "{0}Results_W_N{1}_p{2}_slower{3}.pickle".format(data_dir,N,p_AVG,w_index)
    with open(results_filename, 'rb') as rf:
        try:
            results = pickle.load(rf) 
        except (EOFError):
            print("unpickling error")
          
    result_matrix[w_index-start_index, 0] = w_index
    
    for k in range(1, len(results_header)):    
        result_matrix[w_index-start_index, k] = results[results_header[k]]


np.savetxt('{0}Result_Matrices_{1}_{2}.csv'.format(data_dir, start_index, end_index), 
    result_matrix, delimiter=',', header=str(results_header)) 
    
    


