# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:06:45 2017

@author: rhino
"""

try:
   import cPickle as pickle # used to store python data types
except:
   import pickle
#import dill #pickle works fine

from analyze import analyze_autocor # used to analyze synchrony of networks

try:
    del input # brian overwrites this, we want to reset it to the python default
except:
    pass
input_orig = input # rename the python default for input (brian will overwrite it when imported)

import numpy as np
np.set_printoptions(threshold=np.nan)

import matplotlib.pyplot as plt
import math

N = 3000 # Number of excitatory neurons
p_AVG =50/N # average probability of connectivity between neurons


######Load in matrices one at a time and simulate!
import sys
if len(sys.argv) >= 3:
   start_index = int(sys.argv[1])
   end_index = int(sys.argv[2])
else:
   start_index = int(input_orig("enter a starting index: "))
   end_index = int(input_orig("enter end index: "))

n_indices = end_index-start_index+1


results = dict([('N', np.zeros(n_indices)), ('L_left', np.zeros(n_indices)), ('L_right', np.zeros(n_indices)), ('p_AVG', np.zeros(n_indices)), ('alpha_recip', np.zeros(n_indices)), ('alpha_conv', np.zeros(n_indices)), ('alpha_div', np.zeros(n_indices)), ('alpha_chain', np.zeros(n_indices)), ('p_hat', np.zeros(n_indices)), ('alpha_recip_hat', np.zeros(n_indices)), ('alpha_conv_hat', np.zeros(n_indices)), ('alpha_div_hat', np.zeros(n_indices)), ('alpha_chain_hat', np.zeros(n_indices)), ('largest eigenvalue', np.zeros(n_indices)), ('synchrony', np.zeros(n_indices)), ('j', np.zeros(n_indices)),('ext_rate', np.zeros(n_indices)),('ext_mag', np.zeros(n_indices))])

        

for w_index in range(start_index, end_index+1):
    
    print("w_index = {0}".format(w_index))
    
    results_filename = "matrices/Results_W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,w_index) 

    with open(results_filename, 'rb') as sf:
        try:
            stats = pickle.load(sf) # load in the stats for the W matrix (L, p_hat, alpha values, alpha_hat values)
        except (EOFError):
            print("unpickling error")
    
    ind = w_index-start_index
    
    for key in results:
        results[key][ind]=stats[key]

# save results (pickle new stats dictionary)
result_filename = "matrices/Summary_W_N{0}_p{1}.pickle".format(N,p_AVG) 
with open(result_filename, "wb") as rf:
    pickle.dump(results, rf)



figure()

plot(results['alpha_chain'], results['synchrony'], 'o')
