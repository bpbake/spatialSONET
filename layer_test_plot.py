# -*- coding: utf-8 -*-
"""
Created on Mon Mar 5 2018

@author: Brittany
"""

import sys
import numpy as np
import scipy as sp
from scipy import linalg, special, sparse
from scipy.sparse import linalg as splinalg
import math
import matplotlib
import matplotlib.pyplot as plt
try:
   import cPickle as pickle
except:
   import pickle

data_dir = "matrices/layer_test/"

n = 1000 #number of neurons in each layer
num_samp = 10000 #number of samples to consider (necessary to be able to calculate cov matrix)

sigma_square = 1 #variance of neurons in first layer
rho = 0.05 #rho*sigma_squre = cov of any pair of neurons in first layer 
#(rho is the correleation coeff - it's a value between 0 & 1)

alpha_divs = []
alpha_chains = []
alpha_recips = []
alpha_convs = []
y_corr_coeffs = []

if len(sys.argv) >= 4:
   start_index = int(sys.argv[1])
   end_index = int(sys.argv[2])
   num_layers = int(sys.argv[3])
else:
   start_index = int(input("enter starting index: "))
   end_index = int(input("enter end index: "))
   num_layers = int(input("how many layers? "))

N = n*num_layers

for index in range(start_index, end_index+1):

	stat_filename = "{0}Stats_N{1}_numLay{2}_rho{3}_{4}.pickle".format(data_dir, N, num_layers, rho, index)
	with open(stat_filename, "rb") as sf:
		stats = pickle.load(sf)

	alpha_divs.append(stats['alpha_div_hat'])
	alpha_chains.append(stats['alpha_chain_hat'])
	alpha_recips.append(stats['alpha_recip_hat'])
	alpha_convs.append(stats['alpha_conv_hat'])
	y_corr_coeffs.append(stats['y_corr_coeff'])

# configure the plots
matplotlib.rcParams.update({'font.size': 20})
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

plt.subplot(221)
plt.plot(alpha_convs, y_corr_coeffs, 'o')
plt.ylabel('y_corr_coeff')
# plt.ylim(0,0.05)
plt.xlabel('alpha_conv_hat')
plt.grid(True)

plt.subplot(222)
plt.plot(alpha_chains, y_corr_coeffs, 'o')
plt.ylabel('y_corr_coeff')
# plt.ylim(0,0.05)
plt.xlabel('alpha_chain_hat')
plt.grid(True)

plt.subplot(223)
plt.plot(alpha_divs, y_corr_coeffs, 'o')
plt.ylabel('y_corr_coeff')
# plt.ylim(0,0.05)
plt.xlabel('alpha_div_hat')
plt.grid(True)

plt.subplot(224)
plt.plot(alpha_recips, y_corr_coeffs, 'o')
plt.ylabel('y_corr_coeff')
# plt.ylim(0,0.05)
plt.xlabel('alpha_recip_hat')
plt.grid(True)

plt.show()