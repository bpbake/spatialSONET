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
import produceW
import sys

class MyException(Exception):
    pass

def create_P(N, L_left, L_right, p_AVG):
    print("creating P")
    sys.stdout.flush()
    
    ## set a threshold for the exponent to ensure tightness
    threshold = 1e100 ## this is so extreme, that it's not actually doing anything
    
    p = np.zeros(N)
    ## p= first row of the matrix P, when i=1
    ## P(i, j)= probability of connection from j to i.  
    ##  This depends on the value of i-j mod N, 
    ##  since we are using a ring model to define the p's
    ## Then P is circulant, so we can construct P from its first row.
    ## We'll set the max possible value of any P(i, j) to be p_max (see below)
    
    for i in range(1,N):
        if L_right == 0:
            # v = (N-i)/L_left # #exponential kernel
            v = math.pow((N-i), 2)/(2*math.pow(L_left,2)) ## gaussian kernel
        elif L_left == 0:
            # v = i/L_right ## exponential kernel
            v = math.pow(i, 2)/(2*math.pow(L_right,2)) ## gaussian kernel
        else:
            # v = min((N-i)/L_left, i/L_right)
            v = min(math.pow((N-i), 2)/(2*math.pow(L_left,2)), math.pow(i, 2)/(2*math.pow(L_right,2)))
            
        if v <= threshold: ## threshold may be so extreme, that this will always be true whenever v is nonzero (in float)
            p[i] = math.exp(-v)
    
    p_avg_current = (1/(N-1))*np.sum(p) ## current average value for p
    
    ## scale values to make average p value be p_AVG
    for i in range(N):
        p[i] = (p_AVG*p[i])/p_avg_current
    
    ## check that new average p value is about p_AVG
    p_avg_check = (1/(N-1))*np.sum(p)     
    if p_avg_check>1:
        raise ValueError("Average p value is greater than 1")
    
    ## find the max and min of the p's, and make sure they are all in [0,1]
    p_max = np.amax(p)
    p_min = np.amin(p)
    print("p_max = {0}".format(p_max))
    p_valid = True 
    if (p_max>1) or (p_min<0):
        valid_p = False
        raise ValueError("Some probability ouside of [0,1]")
    
    
    ## Finally, make the P matrix from the p vector
    p_short = p[1::]
    p_new = np.concatenate(([p[0]], p_short[::-1]), 0)
    P = linalg.toeplitz(p_new, p)
    
    
    # plt.matshow(P)
    # plt.show()
    
    if p_valid:
        return(P)
    else:
        raise MyException("P was not created.")