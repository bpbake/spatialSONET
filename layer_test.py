# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 2018

@author: Brittany
"""
data_dir = "matrices/layer_test/"
import os
try:
   os.mkdir(data_dir)
except FileExistsError:
   pass
from useP_to_createW import *
import numpy as np
import scipy as sp
from scipy import linalg, special, sparse
from scipy.sparse import linalg as splinalg
import math
try:
   import cPickle as pickle
except:
   import pickle
import sys
import produceW
import matplotlib.pyplot as plt


if len(sys.argv) >= 4:
   start_index = int(sys.argv[1])
   end_index = int(sys.argv[2])
   num_layers = int(sys.argv[3])
else:
   start_index = int(input("enter starting index: "))
   end_index = int(input("enter end index: "))
   num_layers = int(input("how many layers? "))

for index in range(start_index, end_index+1):
	np.random.seed(index)
	print("\n\n\nindex = {0}".format(index))
	sys.stdout.flush()

	trying = True
	while trying:
		try:
			n = 100 #number of neurons in each layer
			num_samp = 10000 #number of samples to consider (necessary to be able to calculate cov matrix)

			sigma_square = 1 #variance of neurons in first layer
			rho = 0.15 #rho*sigma_squre = cov of any pair of neurons in first layer 
			#(rho is the correleation coeff - it's a value between 0 & 1)

			cov_x = (sigma_square*np.identity(n)) + (rho*sigma_square*(np.ones((n,n))-np.identity(n)))

			L = np.linalg.cholesky(cov_x) #cov_x = L * (L conjugate transpose)

			Z = np.random.standard_normal((n, num_samp)) # array of iid standard normals (each column reps one sample vector)
			X = np.matmul(L, Z) # sample of r.v. representing neurons in first layer

			N = num_layers*n

			p = 10/n #constant probability of connection between any two neurons
			base_P = np.zeros((N, N))
			for i in range(num_layers):
				base_P[(((i+1)%num_layers)*n):((((i+1)%num_layers)+1)*n), (i*n):((i+1)*(n))] = np.ones((n,n))
			P = p*base_P # probability matrix


			alpha_recip = np.random.uniform(-0.5, 1)
			alpha_conv = np.random.uniform(0, 0.5)
			alpha_div = np.random.uniform(0, 0.5)
			# alpha_chain = np.random.uniform(-0.4, 0.5)
			alpha_chain = 0
			# alpha_conv = 0
			# alpha_div = 0
			# alpha_recip = 0
			# print("\nalphas: \nchain {0} \nconv {1} \ndiv {2} \nrecip {3}\n".format(
			# 	alpha_chain, alpha_conv, alpha_div, alpha_recip))

			W = create_W(N, P, alpha_recip, alpha_conv, alpha_div, alpha_chain)
			W_lowerTri = np.tril(W) # truncates W to make it a lower triangular matrix... Feed Forward case
			Wsparse = sparse.csr_matrix(W_lowerTri)
			# Wsparse = sparse.csr_matrix(W)
			# plt.matshow(W_lowerTri)
			# plt.show()

			W_filename = "{0}W_N{1}_p{2}_numLay{3}_{4}".format(data_dir, N, p, num_layers, index)
			with open(W_filename+'.pickle', 'wb') as fp:
			    pickle.dump(Wsparse, fp)  

			WL1 = sparse.csr_matrix.sum(Wsparse)
			Wsquare = Wsparse**2
			p_hat = WL1/(N*(N-1))
			alpha_recip_hat = (np.trace(Wsquare.toarray())/(N*(N-1)*math.pow(p_hat,2)))-1
			alpha_conv_hat = ((sparse.csr_matrix.sum((Wsparse.transpose())*Wsparse) -WL1) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
			alpha_div_hat = ((sparse.csr_matrix.sum(Wsparse*(Wsparse.transpose())) - WL1) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
			alpha_chain_hat = ((sparse.csr_matrix.sum(Wsquare) - np.trace(Wsquare.toarray())) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1

			Y = X
			for i in range(num_layers-1):
				Y = np.matmul(W[((i+1)*n):((i+2)*n), (i*n):((i+1)*(n))], Y)

			Y_cov = np.cov(Y)

			y_avg_cov = (np.sum(Y_cov)-np.trace(Y_cov))/(n*(n-1)) 
			y_avg_var = np.trace(Y_cov)/n

			y_corr_coeff = y_avg_cov/y_avg_var

			stats = dict([('index', index), ('N', N), ('num_layers', num_layers), ('n', n), ('num_samp', num_samp), 
				('rho', rho), ('sigma_square', sigma_square), ('p', p), ('p_hat', p_hat),
			    ('alpha_recip', alpha_recip), ('alpha_recip_hat', alpha_recip_hat), 
			    ('alpha_conv', alpha_conv), ('alpha_conv_hat', alpha_conv_hat), 
			  	('alpha_div', alpha_div), ('alpha_div_hat', alpha_div_hat), 
			  	('alpha_chain', alpha_chain), ('alpha_chain_hat', alpha_chain_hat), 
			    ('y_avg_cov', y_avg_cov), ('y_avg_var', y_avg_var), ('y_corr_coeff', y_corr_coeff)])


			stat_filename = "{0}Stats_N{1}_p{2}_numLay{3}_rho{4}_{5}.pickle".format(data_dir, N, p, num_layers, rho, index) #pickle the dictionary of stats for each W
			with open(stat_filename, "wb") as f:
			    pickle.dump(stats, f) # write the python pickle file for stats

			# print the stats
			for k,v in sorted(stats.items()):
				print(k+":{0}".format(v))
				sys.stdout.flush()

			trying = False

		except Exception as error:
			print(error)
			sys.stdout.flush()
			continue