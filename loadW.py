# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:24:58 2017

@author: Brittany
"""
try:
   import cPickle as pickle
except:
   import pickle


def loadW(N,p,i,L, data_dir='matrices/') :


   #read in pickled W matrix
   W_filename = "{0}Wsparse_N{1}_p{2}_L{3}_{4}.pickle".format(data_dir, N, p, L, i)

   with open(W_filename, 'rb') as wf:
       try:
           Wsparse = pickle.load(wf)
       except (EOFError):
           print("unpickling error")

   # reading in pickled stats file:        
   stat_filename = "{0}Stats_W_N{1}_p{2}_L{3}_{4}.pickle".format(data_dir, N, p, L, i)
   with open(stat_filename, 'rb') as sf:
       try:
           stats = pickle.load(sf)
       except (EOFError):
           print("unpickling error")

   for k,v in sorted(stats.items()):
       print(k+":{0}".format(v))
   
   return(Wsparse, stats)  
