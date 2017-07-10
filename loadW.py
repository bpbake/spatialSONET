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
import math


def loadW(N,p,i, data_dir='matrices/') :


   #read in pickled W matrix
   W_filename = "{0}Wsparse_N{1}_p{2}_{3}.pickle".format(data_dir,N,p,i)

   with open(W_filename, 'rb') as wf:
       try:
           W = pickle.load(wf)
       except (EOFError):
           print("unpickling error")

   # reading in pickled stats file:        
   stat_filename = "{0}Stats_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p,i)
   with open(stat_filename, 'rb') as sf:
       try:
           stats = pickle.load(sf)
       except (EOFError):
           print("unpickling error")

   for k,v in sorted(stats.items()):
       print(k+":{0}".format(v))
   
   return W    
