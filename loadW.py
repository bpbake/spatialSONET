# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:24:58 2017

@author: Brittany
"""
try:
   import cPickle as pickle
except:
   import pickle


def load_W(N,p,i,L, data_dir='matrices/') :


   ## read in pickled W matrix
   W_filename = "{0}Wsparse_N{1}_p{2}_L{3}_{4}.pickle".format(data_dir, N, p, L, i)

   with open(W_filename, 'rb') as wf:
       try:
           Wsparse = pickle.load(wf)
       except (EOFError):
           print("unpickling error")

   ## read in pickled stats file:        
   stat_filename = "{0}Stats_W_N{1}_p{2}_L{3}_{4}.pickle".format(data_dir, N, p, L, i)
   with open(stat_filename, 'rb') as sf:
       try:
           stats = pickle.load(sf)
       except (EOFError):
           print("unpickling error")

   ## print stats
   for k,v in sorted(stats.items()):
       print(k+":{0}".format(v))
   
   return(Wsparse, stats)  


def plot_degree_hist(N,p,i,L,data_dir='matrices/test/'):
  import matplotlib.pyplot as plt
  import numpy as np
  import scipy as sp
  from scipy import sparse
  import matplotlib.colors as colors

  ## Load W
  Wsparse, stats = load_W(N,p,i,L, data_dir)

  W = np.array(Wsparse.todense())

  indeg = np.sum(W,axis=1) ## sum for each row
  outdeg = np.sum(W,axis=0) ## sum for each column


  print("mean in-degree:{0}\nmean out-degree:{1}".format(np.mean(indeg), np.mean(outdeg)))

  ## Plot histograms:
  # plt.figure()
  # # matplotlib.rcParams.update({'font.size': 60})
  # plt.rc('font', family='serif', size=60)
  # plt.rc('xtick', labelsize=50)
  # plt.rc('ytick', labelsize=50)
  # plt.hist(indeg,40) ## data , number of bins
  # plt.xlabel('in-degree')
  # plt.ylabel('Count')
  # mng = plt.get_current_fig_manager()
  # # mng.window.state('zoomed')
  # mng.window.showMaximized()
  # # plt.tight_layout()


  # plt.figure()
  # # matplotlib.rcParams.update({'font.size': 60})
  # plt.rc('font', family='serif', size=60)
  # plt.rc('xtick', labelsize=50)
  # plt.rc('ytick', labelsize=50)
  # plt.hist(outdeg,40) ## data , number of bins
  # plt.xlabel('out-degree')
  # plt.ylabel('Count')
  # mng = plt.get_current_fig_manager()
  # # mng.window.state('zoomed')
  # mng.window.showMaximized()
  # plt.tight_layout()

  ## 2-d historam
  plt.figure()
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif', size=80)
  plt.rc('xtick', labelsize=50)
  plt.rc('ytick', labelsize=50)
  plt.hist2d(outdeg, indeg, bins=80, norm=colors.SymLogNorm(linthresh=1))
  plt.xlabel('out-degree')
  plt.ylabel('in-degree')
  # plt.clim(0,10)
  cb=plt.colorbar()
  cb.set_label('Count')
  mng = plt.get_current_fig_manager()
  # mng.window.state('zoomed')
  mng.window.showMaximized()
  # plt.tight_layout()




  plt.show()