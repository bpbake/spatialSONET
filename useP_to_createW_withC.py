# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 12:37:05 2017

@author: Brittany
"""

import numpy as np
import scipy as sp
from scipy import linalg, special
import math
import matplotlib.pyplot as plt
import warnings
import sys
import produceW ## this is the cython script (python using C)
## you may get an error if you don't compile produceW.pyx to create produceW.c on your on machine

class MyException(Exception):
    pass
###################################
## How W is created (overview):
## We will genenrate Gaussian r.v.'s (Z) 
##    then threshold on each Z to create the adjacency matrix W
##    the Z'a will have a certain covariance matrix 
##          so that W has prescribed first and second order statistics
##    note: the covariance matrix of Z (and it's square root)
##          are too large for memory, so only the components are recorded
##
## To do this (algorithm):
## 1. Define threshold values (M_theta)
## 2. Define approximate covariance matrix of Z
## 3. Approximate the squre root of the covariance matrix of Z (call this S)
## 4. Call produceW (cython script) to Generate iid standard normal r.v.s X
## 5. Within produceW (using C), calculate Z's by Z=SX
## 6. Within produceW (using C), threshold Z's to create W 
##      W_ij = 1 if Z_ij > M_theata_ij, W_ij = 0 otherwise
## 7. Return W
##################################
def create_W(N, P, alpha_recip, alpha_conv, alpha_div, alpha_chain):
    print("creating W")
    sys.stdout.flush()
    ## Generate the sonet matrix with spatial structure from P
    c1=0.2658  ## coefficients for best fit equation of rho as a function of p1, p2, alpha
    c2=0.0470  ## rho(p1,p2,alpha)~c3*(p1+c1)*(p2+c1)*(alpha+c2)
    c3=1.9531
    

    ## Generate the matrix of thresholds.  If Z_ij > M_theta(i,j) then W_ij = 1 
    ## this will be used much later (when we call C)       
    M_theta = math.sqrt(2)*special.erfinv(1 - (2*P))
    print("M_theta has been created. Line 51")
    sys.stdout.flush()
    
    ## Make the matrices of sigma and sigma_tilde
    Msigma = np.ones((N,N)) - np.identity(N) ## used in the covariance matrix of gaussian Z's
    M_tilde = np.multiply(Msigma, np.add(P,c1)) ## used in the approximation of the covariance matrix of Z's
    ## we multiply sigma_ij by (p_ij + c1) to get sigma_ij^tilde

    m = np.mean(np.sum(np.square(M_tilde), axis=0)) ## average row/column sum of M_tilde (the sigma tilde squared)
    ## Each row/column sum should be the same, but we're averaging just in case they're not
    

    ## now, we will solve for c, d, and e using (AA)*(AA)^T = E
    ## c, d, and e are parameters in our approximation of the square root of the covariance matrix of Z
    ## E is a 2x2 matrix:
    E = np.array([[(c3*(alpha_conv +c2)), (c3*(alpha_chain +c2))], [(c3*(alpha_chain +c2)), (c3*(alpha_div +c2))]])

    print("Now finding eigenvalues of E to calculate c, d, e.  Line 67")
    sys.stdout.flush()
    D, V = np.linalg.eig(E)
    ## D is a vector of the eigenvalues, each repeated according to its multiplicity. The eigenvalues are not necessarily ordered. 
    ##   The resulting array will be of complex type, unless the imaginary part is zero in which case it will be cast to a real type. 
    ##   When a is real the resulting eigenvalues will be real (0 imaginary part) or occur in conjugate pairs
    ## V is an array of the corresponding normalized (unit length) eigenvectors, 
    ##   such that the column V[:,i] is the eigenvector corresponding to the eigenvalue D[i]

    print("Eigenvalues have been calculated.  Now checking for negatives. Line 77")
    sys.stdout.flush()
    
    ## If an eigenvalue of E is really small and negative, we can assume there's a roundoff error and it should be zero
    for i in range(2):
        if (D[i] < 0) and (D[i] > -1e-12):
            D[i] = 0
    ## Store eigenvalues of E in diagonal matrix BB
    BB = np.zeros((2, 2))
    if (D[0] >= 0) and (D[1] >= 0):
        BB[0,0] = math.sqrt(D[0])
        BB[1,1] = math.sqrt(D[1])
    # If an eigenvalue is smaller (more negative) than -1e-12, then that's a problem (they should be non-negative)
    else: 
        raise MyException("Some eigenvalues in D are negative.") 
    

    ## Calculate AA, reacall (AA)*(AA)^T = E
    AA = np.matmul(np.matmul(V,BB), np.linalg.inv(V)) 
    c = AA[0,0]/np.sqrt(m)
    d = AA[1,1]/np.sqrt(m)
    e = AA[0,1]/np.sqrt(m)
    print("c,d,e have been calculated. Line 99")
    sys.stdout.flush()
    

    ## Now we'll define F, G, and H (this is how we store "components" of S)
    print("Next creating F, G, H. Line 104")
    sys.stdout.flush()
    ## We'll make this computation easier by defining these:
    M_square = np.square(M_tilde) ## M_square[i,j] = sigma_ij^tilde^2
    m_sq = np.subtract(m,M_square) ## m_sq[i,j] = m - sigma_ij^tilde^2
    m_sqT = np.subtract(m,np.transpose(M_square)) ## m_sqT[i,j] = m - sigma_ji^tilde^2
    
    F = np.subtract(np.square(Msigma), np.multiply(M_square,np.add(np.multiply(math.pow(c,2), m_sq),np.add(np.multiply(math.pow(d,2), m_sq), np.multiply(2*math.pow(e,2),m_sqT)))))   
    ## F[i,j] = A[i,j]^2+B[i,j]^2 = Msigma[i,j]^2 - M[i,j]^2*((c^2 + d^2)*(m-M[i,j]^2) + 2*e^2*(m-M[j,i]^2));
    ## note this requires F[i,j] >=0 
    print("matrix F has been created. Line 114")
    sys.stdout.flush()   
    if (F < np.zeros(N)).any():
        raise MyException("F had some negative entries. Break at line 117")

    G = np.transpose(F)## G[i,j] = A[j,i]^2 + B[i,j]^2 = Msigma[j,i]^2 - M[j,i]^2*((c^2 + d^2)*(m-M[j,i]^2) + 2*e^2*(m - M[i,j]^2));
    ## note this requires G[i,j]>=0
    print("matrix G has been created. Line 121")
    sys.stdout.flush()
    
    H = np.multiply(np.multiply(M_square, np.transpose(M_square)),np.subtract(c3*(alpha_recip +c2), np.multiply(e*(c+d),np.add(m_sq, m_sqT))))
    ## H[i,j] = 2*A[i,j]*B[i,j] = M[i,j]*M[j,i]*rho_recip - M[i,j]*M[j,i]*(c*e*(m-M[i,j])^2) + c*e*(m-M[j,i]^2) + d*e*(m-M[i,j]^2) + d*e*(m-M[j,i]^2);
    ## note sign(H[i,j]) = sign(B[i,j]) since A[i,j]>=0
    print("matrix H has been created. Line 127")
    sys.stdout.flush()
    
    ## note that F[i,j], G[i,j], and H[i,j] do not depend on m because of how c, d, and e depend on m
    ## note that F[i,j]=G[j,i] and H[i,j]=H[j,i]
    
    print("Now creating A, B. Line 133")
    sys.stdout.flush()
    A = np.zeros((N, N)) ## A[i,j] = a_ij
    B = np.zeros((N, N)) ## B[i,j] = b_ij
    
    ## We use the quadratic forumula.  Here are the components:        
    coeff_a = np.square(np.subtract(F,G)) + 4*np.square(H)
    coeff_b = -2*(np.multiply(np.square(H),np.subtract(F,G))) - 2*(np.multiply(F,np.square(np.subtract(F,G)))) - 4*(np.multiply(np.square(H),F))
    coeff_c = np.square(np.square(H)) + 2*(np.multiply(np.multiply(np.square(H),F),np.subtract(F,G))) + np.multiply(np.square(np.subtract(F,G)),np.square(F))
    Disc = np.square(coeff_b) - 4*np.multiply(coeff_a, coeff_c)
    print("Discriminant calculated to find A and B. Line 143")
    sys.stdout.flush()

    Asq = np.divide(np.add(-1*coeff_b, np.sqrt(Disc)),2*coeff_a)
    A = np.sqrt(Asq)
    B = np.divide(H, np.add(A,np.transpose(A)))
    print("Initial A and B have been created.  Next making zeros along diagonals. Line 149")
    sys.stdout.flush()
        
    np.fill_diagonal(A,0)
    np.fill_diagonal(B,0)

    ## Check for non-negative Disc 
    ## it's okay if it is negative, but then W won't accurately represent the prescribed alpha values    
    if  (Disc < np.zeros(N)).any():
        warnings.warn("discriminant<=0 for some i,j.  alpha_recip not accounted for.")
        sys.stdout.flush()
        for i in range(N):
            for j in range(i):
                if Disc[i,j] < 0:
                    B[i,j] = 0
                    A[i,j] = math.sqrt(F[i,j])
                if Disc[j,i] < 0:
                    B[j,i] = 0
                    A[j,i] = math.sqrt(F[j,i])
    print("Values in A have been edited to be sqrt(F) if discriminant was negative. Line 168")
    sys.stdout.flush()
    
    ## A^2 + B^2 should be very close to F                
    checkFmat = np.amax(np.absolute(np.subtract(np.add(np.square(A),np.square(B)),F)))
    print("|A^2 + B^2 - F| leq {0}".format(checkFmat))
    sys.stdout.flush()


    ## We've finished computing the components of S
    ## Now, use produceW (cython) to calculate W
    print("Now calling cython function to create W. Line 179")
    sys.stdout.flush()
    W = np.array(produceW.produceW(A, B, M_tilde, M_theta, c, d, e, N))
    print("W creation is complete. Line 182")
    sys.stdout.flush()
    
    return(W)