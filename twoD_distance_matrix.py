# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 11:13 2017

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

def create_distance(N, A):
	# N = 10000 # perfect square would be nice
	# A = 100 # a divisor of N, the number of x-coordinates in the grid

	#location = np.random.uniform(0, 10, (N,2)) # for random locations

	location = np.zeros((N,2)) # for a grid of locations with A x-coordinates
	for i in range(0,N):
		location[i,0] = i%A
		location[i,1] = i//A

	distance = np.zeros((N,N))
	for i in range(0,N):
		for j in range(0, N):
			distance[i,j] = np.linalg.norm(location[i]-location[j])

	return(location, distance)



def create_P(N, distance, L_left, L_right, p_AVG):
    print("creating P")
    sys.stdout.flush()
    
    threshold = 1e100


    P = np.zeros((N,N))

    for i in range(0,N):
    	for j in range(0,i):
	        if L_right == 0:
	            #v = distance[i,j]/L_left #exponential kernel
	            v = math.pow(distance[i,j], 2)/(2*math.pow(L_left,2)) #gaussian kernel
	        elif L_left == 0:
	            #v = distance[i,j]/L_right #exponential kernel
	            v = math.pow(distance[i,j], 2)/(2*math.pow(L_right,2)) #gaussian kernel
	        else:
	            #v = min(distance[i,j]/L_left, distance[i,j]/L_right)
	            v = min(math.pow(distance[i,j], 2)/(2*math.pow(L_left,2)), math.pow(distance[i,j], 2)/(2*math.pow(L_right,2)))
	            
	        if v <= threshold:
	            P[i,j] = math.exp(-v)
	            P[j,i] = P[i,j]

	#######################################################
	#### Fix below here!
    
    p_avg_current = (1/(N*N))*np.sum(P) #current average value for p
    
    # #scale values to make average p value be p_AVG
    # for i in range(0,N):
    # 	for j in range(0,i):
    #     	P[i,j] = (p_AVG*P[i,j])/p_avg_current
    #     	P[j,i] = (p_AVG*P[j,i])/p_avg_current

    finalP = np.divide(np.multiply(p_AVG,P),p_avg_current)
    	
    p_avg_check = (1/(N*N))*np.sum(finalP) #check that new average p value is about p_AVG
    
    if p_avg_check>1:
        raise ValueError("Average p value is greater than 1")
    
    p_max = np.amax(finalP)
    p_min = np.amin(finalP)
    print("p_max = {0}".format(p_max))
    
    p_valid = True #make sure our entries in P are between 0 and 1
    if (p_max>1) or (p_min<0):
        valid_p = False
        raise ValueError("Some probability ouside of [0,1]")

    if p_valid:
        return(finalP)
    else:
        raise MyException("P was not created.")