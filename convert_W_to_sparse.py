# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 10:01:50 2017

@author: Brittany
"""

data_dir = 'matrices/N3000_LL90_LR0_W/'

try:
   import cPickle as pickle # used to store python data types
except:
   import pickle
   
try:
    del input # brian overwrites this, we want to reset it to the python default
except:
    pass
input_orig = input
   
import numpy as np
np.set_printoptions(threshold=np.nan)

import scipy as sp

import matplotlib.pyplot as plt
import math


N = 3000
p_AVG = 50/N

#if len(sys.argv) >= 3:
#   start_index = int(sys.argv[1])
#   end_index = int(sys.argv[2])
#else:
start_index = int(input_orig("enter a starting index: "))
end_index = int(input_orig("enter end index: "))

for w_index in range(start_index, end_index+1):
    print("w_index = {0}".format(w_index))    
    
    W_filename = "{0}W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,w_index)
    with open(W_filename, 'rb') as wf:
        try:
            W = pickle.load(wf) # load in W matrix
        except (EOFError):
            print("unpickling error")
    
    Wsparse = sp.sparse.csr_matrix(W)
    plt.matshow(Wsparse.toarray())
    
    Wsparse_filename = "{0}Wsparse_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,w_index)
    with open(Wsparse_filename, 'wb') as fp:
        pickle.dump(Wsparse, fp)
    
#    stats_filename = "{0}Stats_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,w_index)
#    with open(stats_filename, 'rb') as sf:
#        try:
#            stats = pickle.load(sf) # load in the stats for the W matrix (L, p_hat, alpha values, alpha_hat values)
#        except (EOFError):
#            print("unpickling error")