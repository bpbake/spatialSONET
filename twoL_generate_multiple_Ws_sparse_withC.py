# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 16:41:12 2017

@author: Brittany
"""

# Other scripts needed to run this code:  
#   useP_to_createW_withC.py
#   create_P.py
#   produceW.pyx converted into produceW.c


data_dir = 'matrices/N10000_LL70_LR0_alphaRandom/'
import os
try:
   os.mkdir(data_dir)
except FileExistsError:
   pass

import math
import numpy as np
import scipy as sp
from scipy import sparse

from useP_to_createW_withC import *
from create_P import *

try:
   import cPickle as pickle
except:
   import pickle

import sys

try:
    del input # brian overwrites this, we want to reset it to the python default
except:
    pass
input_orig = input


N = 10000 # Number of excitatory neurons
p_AVG = 50/N # average probability of connectivity between neurons

if len(sys.argv) >= 3:
   start_index = int(sys.argv[1])
   end_index = int(sys.argv[2])
else:
   start_index = int(input_orig("enter starting index: "))
   end_index = int(input_orig("enter end index: "))

for w_index in range(start_index, end_index+1): #so i=start_index, start_index+1,start_index+2,...,start_index+num_matrices-1
    np.random.seed(w_index)
    
    trying = True
    while trying:  
        try:
            print("\nmaking matrix {0}".format(w_index))
            sys.stdout.flush()
            
            #generate Ls, alphas
#            L_left = math.exp(np.random.uniform(math.log(45), math.log(10000)))# L=[90,22000]ish
#            L_right = math.exp(np.random.uniform(math.log(45), math.log(10000)))# L=[90,22000]ish  L_left #math.exp(np.random.uniform(4.5, 10))
            alpha_recip = np.random.uniform(0, 0.3)
            alpha_conv = np.random.uniform(0, 0.3)
            alpha_div = np.random.uniform(0, 0.3)
            alpha_chain = np.random.uniform(-0.4, 0.3)
            L_left = 70 #float("inf")# math.inf
            L_right = 0 #float("inf")#math.inf
            # alpha_recip = 0.3
            # alpha_conv = 0.3
            # alpha_div = 0.3
            # alpha_chain = 0.3
            

            P = create_P(N, L_left, L_right, p_AVG)
            print("P has been created \n")
            sys.stdout.flush()
            
            #call other program to create W (and return W)
            W = create_W(N, P, alpha_recip, alpha_conv, alpha_div, alpha_chain)
            print("W has been created \n")
            sys.stdout.flush()
            Wsparse = sparse.csr_matrix(W)
            
            #save the W
            W_filename = "{0}Wsparse_N{1}_p{2}_{3}.pickle".format(data_dir, N, p_AVG, w_index)
            with open(W_filename, 'wb') as fp:
                pickle.dump(Wsparse, fp)    
            print("W has been pickled.")
            sys.stdout.flush()
                
            # generate statistics of W
            print("generating stats")
            sys.stdout.flush()
            WL1 = sparse.csr_matrix.sum(Wsparse)
            Wsquare = Wsparse**2
            p_hat = WL1/(N*(N-1))
            alpha_recip_hat = (np.trace(Wsquare.toarray())/(N*(N-1)*math.pow(p_hat,2)))-1
            alpha_conv_hat = ((sparse.csr_matrix.sum((Wsparse.transpose())*Wsparse) -WL1) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
            alpha_div_hat = ((sparse.csr_matrix.sum(Wsparse*(Wsparse.transpose())) - WL1) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
            alpha_chain_hat = ((sparse.csr_matrix.sum(Wsquare) - np.trace(Wsquare.toarray())) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
            
            #evals = np.linalg.eigvals(W)
            #largest_eigenval = np.ndarray.max(evals)
            max_eigenval = sparse.linalg.eigs(Wsparse, k=1, which='LR', return_eigenvectors=False, ncv = 500)
            print("stats have been calculated. Now saving as a dict")
            sys.stdout.flush()
            
            #if largest_eigenval != max_eigenval:
            #    print("different max eigenvalues. \n np largest = {0} \n sp largest = {1}".format(largest_eigenval, max_eigenval))
                
            #create dictionary of stats and save
            stats = dict([('N', N), ('L_left', L_left), ('L_right', L_right), ('p_AVG', p_AVG), ('alpha_recip', alpha_recip), ('alpha_conv', alpha_conv), ('alpha_div', alpha_div), ('alpha_chain', alpha_chain), ('p_hat', p_hat), ('alpha_recip_hat', alpha_recip_hat), ('alpha_conv_hat', alpha_conv_hat), ('alpha_div_hat', alpha_div_hat), ('alpha_chain_hat', alpha_chain_hat), ('largest eigenvalue', max_eigenval)])
            
            stat_filename = "{0}Stats_W_N{1}_p{2}_{3}.pickle".format(data_dir, N,p_AVG, w_index) #pickle the dictionary of stats for each W
            with open(stat_filename, "wb") as f:
                pickle.dump(stats, f)
            print("stats have been pickled")
            sys.stdout.flush()

            # print the stats
            for k,v in sorted(stats.items()):
                print(k+":{0}".format(v))
                sys.stdout.flush()

            trying = False

            # plot the W matrix
            #plt.matshow(W)
            #plt.show()
                
                
        except Exception as error:
            print(error)
            continue
        
        
