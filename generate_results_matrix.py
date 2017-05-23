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

N = 3000
p_AVG = 50/N

results_header = ['w_index', 'synchrony', 'N', 'L_left', 'L_right', 'largest eigenvalue', 'alpha_chain', 'alpha_chain_hat', 'alpha_conv', 'alpha_conv_hat', 'alpha_div', 'alpha_div_hat', 'alpha_recip', 'alpha_recip_hat', 'p_AVG', 'p_hat']

start_index = int(input_orig("enter a starting index: "))
end_index = int(input_orig("enter end index: "))

result_matrix = np.zeros((end_index-start_index+1, len(results_header)))
                                                  
for w_index in range(start_index, end_index+1):
   
    results_filename = "matrices\Results_W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,w_index)
    with open(results_filename, 'rb') as rf:
        try:
            results = pickle.load(rf) 
        except (EOFError):
            print("unpickling error")
          
    result_matrix[w_index-start_index, 0] = w_index
    
    for k in range(1, len(results_header)):    
        result_matrix[w_index-start_index, k] = results[results_header[k]]
#    result_matrix[w_index, 2] = results['N']
#    result_matrix[w_index, 2] = results['L_left']
#    result_matrix[w_index, 2] = results['L_right']
#    result_matrix[w_index, 2] = results['largest eigenvalue']
#    result_matrix[w_index, 3] = results['alpha_chain_hat']
#    result_matrix[w_index, 4] = results['alpha_conv_hat']
#    result_matrix[w_index, 5] = results['alpha_div_hat']
#    result_matrix[w_index, 6]= results['alpha_recip_hat']
#    result_matrix[w_index, 7] = results['alpha_chain']
#    result_matrix[w_index, 8] = results['alpha_conv']
#    result_matrix[w_index, 9] = results['alpha_div']
#    result_matrix[w_index, 10] = results['alpha_recip']
#    result_matrix[w_index, 11] = results['p_hat']
#    result_matrix[w_index, 12] = results['p_AVG']


np.savetxt('result_matrix.csv', result_matrix, delimiter=',', header=str(results_header)) 
    
    


