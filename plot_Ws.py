# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:24:58 2017

@author: Brittany
"""
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
   
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import sparse
import math
from loadW import loadW

# data_dir = 'matrices/N3000_LL100_LR0_ff_alpha_div_rand/'
# data_dir = 'matrices/N3000_LL50_LR50_recurr_alpha_div_rand/'
data_dir = 'matrices/test/'
print("data_dir: {0}".format(data_dir))

N = 3000 ## minimum of 1000
p = 30/N
print("N={0}".format(N))
L_left=float('inf')
# L_left=50
print("L_left:{0}".format(L_left))


## what matrices do we want to plot:
if len(sys.argv) >= 3:
   start_index = int(sys.argv[1])
   end_index = int(sys.argv[2])
else:
   start_index = int(input_orig("enter starting index: "))
   end_index = int(input_orig("enter end index: "))

## now, let's load those matrices
for w_index in range(start_index, end_index+1):
    print("\nmatrix {0}".format(w_index))
    sys.stdout.flush()

    ## load in pickled W matrix and stats
    Wsparse, stats = load_W(N,p,w_index,L_left,data_dir)
            
    W = Wsparse.todense()

    ## plot the W matrix
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=100)
    plt.rc('xtick', labelsize=50)
    plt.rc('ytick', labelsize=50)
    plt.figure(figsize=(11,10))

    # plt.figure()
    # plt.suptitle('Matrix: {0}'.format(w_index))
    plt.spy(W, markersize=0.5)
    plt.xticks(np.arange(0,N+1,1000))
    plt.yticks(np.arange(0,N+1,1000))
    # mng = plt.get_current_fig_manager()
    # mng.window.state('zoomed')
    # mng.window.showMaximized()

plt.show()