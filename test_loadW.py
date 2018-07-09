# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:24:58 2017

@author: Brittany
"""
try:
   import cPickle as pickle
except:
   import pickle
   
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import sparse
import math

data_dir = 'matrices/N10000_LL70_LR0_ff_alpha_div_rand/'
# data_dir = 'matrices/N10000_LL70_LR0_ff_alpha_chain_zero/'
# data_dir = 'matrices/N10000_LL70_LR0_ff_alphas_all_rand/'
# data_dir = 'matrices/N10000_LL70_LR0_ff/'
# data_dir = 'matrices/CNS18/'

N = 10000 #minimum of 1000
p_AVG = 50/N
print("N={0}".format(N))
L_left=45

w_index = 1

#read in pickled W matrix
# W_filename = "{0}W_N{1}_p{2}_L{3}_achain0.3_{4}".format(data_dir, N, p_AVG, L_left, w_index)
W_filename = "{0}Wsparse_N{1}_p{2}_{3}".format(data_dir, N, p_AVG, w_index)

with open(W_filename+".pickle", 'rb') as wf:
    try:
        Wsparse = pickle.load(wf)
    except (EOFError):
        print("unpickling error")
        
W = Wsparse.todense()
#print("W= {0}".format(W))


#eigenvals = sparse.linalg.eigs(Wsparse, k=10, which='LR', return_eigenvectors=False, ncv = 500)
#plt.scatter(np.real(eigenvals), np.imag(eigenvals))

# reading in pickled stats file:        
# stat_filename = "{0}Stats_W_N{1}_p{2}_L{3}_achain0.3_{4}.pickle".format(data_dir, N, p_AVG, L_left, w_index)
stat_filename = "{0}Stats_W_N{1}_p{2}_{3}.pickle".format(data_dir, N, p_AVG, w_index)
with open(stat_filename, 'rb') as sf:
    try:
        stats = pickle.load(sf)
    except (EOFError):
        print("unpickling error")

for k,v in sorted(stats.items()):
    print(k+":{0}".format(v))

# add a new entry to the stats dictionary
# stats["hi"] = 5 #key= "hi", value=5

# print(stats)

# save updated stats dictionary as a picle file
# stats_with_results_filename = "matrices\StatsWithSync_W_N{0}_p{1}_{2}.pickle".format(N,p,i)
# with open(stats_with_results_filename, 'wb') as f:
#     pickle.dump(stats, f)

# plot the W matrixplt.rc('text', usetex=True)
plt.rc('font', family='serif', size=100)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)

plt.figure()
plt.spy(W, markersize=1)
plt.xticks(np.arange(0,10001,2000))
plt.yticks(np.arange(0,10001,2000))
plt.show()
    
# generate statistics of matrix
#p_hat = np.sum(W)/(N*(N-1))
#alpha_recip_hat = (np.trace(np.matmul(W,W))/(N*(N-1)*math.pow(p_hat,2)))-1
#alpha_conv_hat = ((np.sum(np.matmul(np.transpose(W),W)) - np.sum(W)) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
#alpha_div_hat = ((np.sum(np.matmul(W,np.transpose(W))) - np.sum(W)) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
#alpha_chain_hat = ((np.sum(np.matmul(W,W)) - np.trace(np.matmul(W,W))) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1