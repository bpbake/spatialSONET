# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 2018

@author: Brittany
"""

from useP_to_createW import *
import numpy as np

w_index = int(input_orig("enter W index: "))
np.random.seed(w_index)

num_x = 100 #number of neurons in first layer
num_y = 100 #number of neurons in second layer
num_samp = 10000 #number of samples to consider (necessary to be able to calculate cov matrix)

sigma_square = 1 #variance of neurons in first layer
rho = .7 #rho*sigma_squre = cov of any pair of neurons in first layer 
#(rho is the correleation coeff - it's a value between 0 & 1)

cov_x = (sigma_square*np.identity(num_x)) + (rho*sigma_square*(np.ones((num_x,num_x))-np.identity(num_x)))

L = np.linalg.cholesky(cov_x) #cov_x = L * (L conjugate transpose)

Z = np.random.standard_normal((num_x, num_samp)) # array of iid standard normals (each column reps one sample vector)
X = np.matmul(L, Z) # sample of r.v. representing neurons in first layer
print("samples of X have been calcuated.  Now generating W.")
sys.stdout.flush()

#N = max(num_x, num_y)
N = num_x+num_y

p = .1 #constant probability of connection between any two neurons
P = p*np.ones((N,N)) # probability matrix

alpha_chain = 0
alpha_conv = 0
alpha_div = 0
alpha_recip = 0
print("\nalphas: \nchain {0} \nconv {1} \ndiv {2} \nrecip {3}\n".format(alpha_chain, alpha_conv, alpha_div, alpha_recip))

W = create_W(N, P, alpha_recip, alpha_conv, alpha_div, alpha_chain)
Wsparse = sparse.csr_matrix(W)
W_star = W[num_x:, :num_x] # the connectivity matrix for x's to y's (W_star[i,j] = 1 if x_j to y_i, 0 o.w.)

data_dir = "layer_test"
W_filename = "{0}W_numX{1}_numY{2}_rho{3}_{4}".format(data_dir, num_x, num_y, rho, w_index)
with open(W_filename+'.pickle', 'wb') as fp:
    pickle.dump(W_star, fp)  
print("W has been generated and pickled.")
sys.stdout.flush()

WL1 = sparse.csr_matrix.sum(Wsparse)
Wsquare = Wsparse**2
p_hat = WL1/(N*(N-1))
alpha_recip_hat = (np.trace(Wsquare.toarray())/(N*(N-1)*math.pow(p_hat,2)))-1
alpha_conv_hat = ((sparse.csr_matrix.sum((Wsparse.transpose())*Wsparse) -WL1) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
alpha_div_hat = ((sparse.csr_matrix.sum(Wsparse*(Wsparse.transpose())) - WL1) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
alpha_chain_hat = ((sparse.csr_matrix.sum(Wsquare) - np.trace(Wsquare.toarray())) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
print("hat values for W have been calculated. Now calculating Y.")
sys.stdout.flush()


Y = np.matmul(W_star, X)

Y_cov = np.cov(Y)

y_avg_cov = (np.sum(Y_cov)-np.trace(Y_cov))/(num_y*(num_y-1)) 
y_avg_var = np.trace(Y_cov)/num_y

y_corr_coeff = y_avg_cov/y_avg_var

print("Y and wants have been calculated. Now saving all stats as a dict")
sys.stdout.flush()

stats = dict([('w_index', w_index), ('num_x', num_x), ('num_y', num_y), ('num_samp', num_samp), 
	('rho', rho), ('sigma_square', sigma_square), ('p', p), ('p_hat', p_hat),
    ('alpha_recip', alpha_recip), ('alpha_recip_hat', alpha_recip_hat), 
    ('alpha_conv', alpha_conv), ('alpha_conv_hat', alpha_conv_hat), 
  	('alpha_div', alpha_div), ('alpha_div_hat', alpha_div_hat), 
  	('alpha_chain', alpha_chain), ('alpha_chain_hat', alpha_chain_hat), 
    ('y_avg_cov', y_avg_cov), ('y_avg_var', y_avg_var), ('y_corr_coeff', y_corr_coeff)])


stat_filename = "{0}Stats_numX{1}_numY{2}_rho{3}_{4}.pickle".format(data_dir, num_x, num_y, rho, w_index) #pickle the dictionary of stats for each W
with open(stat_filename, "wb") as f:
    pickle.dump(stats, f) # write the python pickle file for stats
print("stats have been pickled")
sys.stdout.flush()

# print the stats
for k,v in sorted(stats.items()):
    print(k+":{0}".format(v))
    sys.stdout.flush()