# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 12:32:09 2017

@author: Brittany
"""

## to create produceW.c use the command:  python setup.py build_ext --inplace


import numpy as np
import sys


###########################################
## this program uses C (runs faster than python) to create the adjacency matrix W
## parameters:
##   A, B, M_tilde, c, d, e: the components of S, the approximate square root of Sigma 
##      (Sigma is the covariance matrix of the Gaussian r.v.'s Z for dichotomization)
##   M_theta: thresholds for dichotomization (W_ij = 1 if Z_ij > M_theta_ij, W_ij = 0 otherwise)
##   N: the number of neurons in the network
## returns: 
##   adjacency matrix W
##########################################
def produceW(double[:,:] A, double[:,:] B, double[:,:] M_tilde, double[:,:] M_theta, double c, double d, double e, int N):
    cdef int i, j, n
    
    cdef double[:,:] X = np.random.standard_normal((N,N)) ## random matrix of iid standard normal entries
    ## We don't actually need to create/initialize a matrix Z, but here's how it's used
    # cdef double[:,:] Z = np.zeros((N,N)) ## Z = SX, a random vector with entries iid standard normal
    ## Z is Gaussian with covariance matrix approximately Sigma
    cdef double[:,:] W = np.zeros((N,N)) # dichomatized Z

    cdef double[:,:] MX = np.multiply(M_tilde, X) ## entry MX[i,j] is sigma_tilde_ij*X_ij
    cdef double[:] MXrowsum = np.sum(MX, axis=1) ## vector of sums of rows of MX
    cdef double[:] MXcolsum = np.sum(MX, axis=0) ## vector of sums of columns of MX

    ## initialize some variables 
    ##    these will be defined in the double for loop where W is created
    cdef double cMij, dMij, eMij
    cdef double sumCDEij
    cdef double Zij 
    
    ## Now: create W
    for i in range(N):
        for j in range(N):
            cMij = c*M_tilde[i,j]
            dMij = d*M_tilde[i,j]
            eMij = e*M_tilde[i,j]

            sumCDEij = cMij*MXrowsum[i] + dMij*MXcolsum[j] + eMij*(MXrowsum[j] + MXcolsum[i])

            Zij = A[i,j]*X[i,j] + B[i,j]*X[j,i] - (cMij*MX[i,j] + dMij*MX[i,j] + 2*eMij*MX[j,i]) + sumCDEij

            if Zij >= M_theta[i,j]:
                W[i,j]=1
            ## W was initialized as a matrix of zeros, so we don't really need:
            #else:
            #   W[i,j]=0

        ## Let's keep track of our progress through this Matrix creation:
        if (i%500 == 0):
            print("row={0} out of N={1}".format(i,N))
            sys.stdout.flush()
    print("\n")
    return W