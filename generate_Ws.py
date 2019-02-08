# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 16:41:12 2017

@author: Brittany
"""
## This script creates multiple adjacency matrices (Ws) 
## and saves them as sparse matrices in python pickle files

## Other scripts needed to run this code:  
##   useP_to_createW.py
##   create_P.py
##   produceW.pyx converted into produceW.c

## Define a directory for adjacency matrices (saved as python pickle files)
# data_dir = 'matrices/N3000_LL50_LR50_recurr_alphas_all_rand/'
data_dir = 'matrices/N1000_LL100_LR0_ff_alpha_conv_div_rand/'
# data_dir = 'matrices/N1000_erdos_renyi/'
# data_dir = 'matrices/N1000_Linf_recurr_alphas_all_rand/'
# data_dir = "matrices/"
print("data_dir: {0}".format(data_dir))

import os
try:
   os.mkdir(data_dir)
except FileExistsError:
   pass

import math
import numpy as np
import scipy as sp
from scipy import linalg, special, sparse
from scipy.sparse import linalg as splinalg

## all functions and objects from useP_to_createW.py and create_p.py are needed
from useP_to_createW import *
from create_P import *

try:
   import cPickle as pickle
except:
   import pickle

import sys

try:
    del input ## brian simulator overwrites this, we want to reset it to the python default
except:
    pass
input_orig = input


N = 1000 ## Number of excitatory neurons
p_AVG = 50/N ## average probability of connectivity between neurons

## what indices do you want to use for these matrices:
if len(sys.argv) >= 3:
   start_index = int(sys.argv[1])
   end_index = int(sys.argv[2])
else:
   start_index = int(input_orig("enter starting index: "))
   end_index = int(input_orig("enter end index: "))

## now, let's create the matrices with chosen indices
for w_index in range(start_index, end_index+1):
    np.random.seed(w_index) ## randomness can be used in many places
    ## randomness is always used in the matrix creation
    
    trying = True
    while trying: ## sometimes the parameters won't work, so we keep trying
    ## watch out: this could be an infinite loop with some network parameters 
        try:
            print("\nmaking matrix {0}".format(w_index))
            sys.stdout.flush()
            
            ## Define the L_left, L_right, and alpha values here:
            # L_left = math.exp(np.random.uniform(math.log(45), math.log(10000)))# L=[90,22000]ish
            # L_right = math.exp(np.random.uniform(math.log(45), math.log(10000)))# L=[90,22000]ish
            L_left = 100 
            # L_left = float("inf") ## for homogeneous networks (as in Zhao et al.)
            L_right = 0 
            # L_right = L_left ## for symmetric/recurrent networks

            # alpha_recip = np.random.uniform(-0.5, 1)
            alpha_conv = np.random.uniform(0, 0.5)
            alpha_div = np.random.uniform(0, 0.5)
            # alpha_chain = np.random.uniform(-0.5, 0.5)
            alpha_recip = 0
            # alpha_conv = 0
            # alpha_div = 0
            alpha_chain = 0

            # print('alpha_recip={0}'.format(alpha_recip))
            # print('alpha_conv={0}'.format(alpha_conv))
            # print('alpha_div={0}'.format(alpha_div))
            # print('alpha_chain={0}'.format(alpha_chain))
            

            ## call the function in create_P.py to create the matrix of connection probabilities
            P = create_P(N, L_left, L_right, p_AVG)
            print("P has been created \n")
            sys.stdout.flush()
            
            ## call the function in useP_to_createW_withC.py to generate the adjacency matrix
            W = create_W(N, P, alpha_recip, alpha_conv, alpha_div, alpha_chain)
            print("W has been created \n")
            sys.stdout.flush()

            ## the matrix is sparse, so let's save it that way
            Wsparse = sparse.csr_matrix(W)

            ## truncate W to make it a lower triangular matrix (used in the Feed Forward case)
            # W_lowerTri = np.tril(W) 
            # Wsparse = sparse.csr_matrix(W_lowerTri) 
            

            ## save the adjacency matrix W as a python pickle file
            W_filename = "{0}Wsparse_N{1}_p{2}_L{3}_{4}".format(data_dir, N, p_AVG, L_left, w_index)
            with open(W_filename+'.pickle', 'wb') as fp:
                pickle.dump(Wsparse, fp)
            print("W has been pickled.")
            sys.stdout.flush()

            ## or save the full adjacency matrix W as a csv file 
            # np.savetxt(W_filename+'.csv', W, delimiter=',')  
            # print("W has been saved as csv file.")
            # sys.stdout.flush()

                
            ## generate statistics of W
            print("generating stats")
            sys.stdout.flush()

            WL1 = sparse.csr_matrix.sum(Wsparse)
            Wsquare = Wsparse**2

            p_hat = WL1/(N*(N-1))
            alpha_recip_hat = (np.trace(Wsquare.toarray())/(N*(N-1)*math.pow(p_hat,2)))-1
            alpha_conv_hat = ((sparse.csr_matrix.sum((Wsparse.transpose())*Wsparse) -WL1) 
                / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
            alpha_div_hat = ((sparse.csr_matrix.sum(Wsparse*(Wsparse.transpose())) - WL1) 
                / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
            alpha_chain_hat = ((sparse.csr_matrix.sum(Wsquare) - np.trace(Wsquare.toarray())) 
                / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
            
            max_eigenval = splinalg.eigs(Wsparse, k=1, which='LR', return_eigenvectors=False, ncv = 500)

            print("stats have been calculated. Now saving as a dict")
            sys.stdout.flush()
            

            ## create a python dictionary of the stats and save as python pickle file
            stats = dict([('N', N), ('L_left', L_left), ('L_right', L_right), ('p_AVG', p_AVG),  
                ('alpha_recip', alpha_recip), ('alpha_conv', alpha_conv), ('alpha_div', alpha_div), 
                ('alpha_chain', alpha_chain), ('p_hat', p_hat), ('alpha_recip_hat', alpha_recip_hat), 
                ('alpha_conv_hat', alpha_conv_hat), ('alpha_div_hat', alpha_div_hat), ('alpha_chain_hat', alpha_chain_hat), 
                ('largest eigenvalue', max_eigenval)])
            
            stat_filename = "{0}Stats_W_N{1}_p{2}_L{3}_{4}.pickle".format(data_dir, N, p_AVG, L_left, w_index) 
            with open(stat_filename, "wb") as f:
                pickle.dump(stats, f)

            print("stats have been pickled")
            sys.stdout.flush()


            ## print the stats
            for k,v in sorted(stats.items()):
                print(k+":{0}".format(v))
                sys.stdout.flush()

            trying = False

            ## plot the W matrix
            # plt.matshow(W)
            # plt.show()
                
        
        ## if there's an error (probably from create_W), just print the error and keep trying        
        except Exception as error:
            print(error)
            continue