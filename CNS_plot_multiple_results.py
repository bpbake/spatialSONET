# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:06:45 2017

@author: rhino
"""

try:
   import cPickle as pickle # used to store python data types
except:
   import pickle
#import dill #pickle works fine


try:
    del input # brian overwrites this, we want to reset it to the python default
except:
    pass
input_orig = input # rename the python default for input (brian will overwrite it when imported)

import numpy as np
np.set_printoptions(threshold=np.nan)

import matplotlib.pyplot as plt

import analyze_results as ar

N = 3000 # Number of excitatory neurons
p_AVG =50/N # average probability of connectivity between neurons
# data_dir = 'matrices/CNS18/'
# data_dir = "matrices/N10000_LL70_LR0_ff_alpha_chain_zero/"
# data_dir = "matrices/N3000_LL70_LR70_sym_alpha_div_rand/"
data_dir = "matrices/N3000_LL70_LR0_ff_alphas_all_rand/"
# data_dir = "matrices/N3000_LL70_LR70_sym_alphas_all_rand/"
# Style = "FF_L70"
# Style = "Regular120s_Clean"
# Style = "CnsLClean" #Cns = Irregular, CnsL = Regular
Style = "Regular5s_Clean_"
# Style = "ttLongClean"

# reload=False
reload=True

summary_filename = "{0}Summary_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,Style) 


if reload:
    with open(summary_filename, "rb") as rf:
        results = pickle.load(rf)

else:
    ######Load in matrices one at a time and simulate!
    import sys
    if len(sys.argv) >= 3:
       start_index = int(sys.argv[1])
       end_index = int(sys.argv[2])
    else:
       start_index = int(input_orig("enter a starting index: "))
       end_index = int(input_orig("enter end index: "))
    
    n_indices = end_index-start_index+1
    
    

    results = dict([('N', np.zeros(n_indices)), ('L_left', np.zeros(n_indices)), ('L_right', np.zeros(n_indices)), 
        ('p_AVG', np.zeros(n_indices)), ('alpha_recip', np.zeros(n_indices)), ('alpha_conv', np.zeros(n_indices)), 
        ('alpha_div', np.zeros(n_indices)), ('alpha_chain', np.zeros(n_indices)), ('p_hat', np.zeros(n_indices)), 
        ('alpha_recip_hat', np.zeros(n_indices)), ('alpha_conv_hat', np.zeros(n_indices)), ('alpha_div_hat', np.zeros(n_indices)), 
        ('alpha_chain_hat', np.zeros(n_indices)), #('largest eigenvalue', np.zeros(n_indices)), 
        ('event_rate', np.zeros(n_indices)), ('IEI excess_kurtosis', np.zeros(n_indices)), ('IEI skew', np.zeros(n_indices)),]) 
        #('j', np.zeros(n_indices)),('ext_rate', np.zeros(n_indices)),('ext_mag', np.zeros(n_indices))])
    actual_index = -1

    for w_index in range(start_index, end_index+1):
            
        print("w_index = {0}".format(w_index))

        try:
            samp_results = ar.load_results(N, p_AVG, w_index, Style, data_dir)
        except:
            print("couldn't load {0}Results_W_N{1}_p{2}_{3}{4}.pickle".format(data_dir,N,p_AVG,Style,w_index))
            continue

        if samp_results['average firing rate'] > 100:
            print("skipped network {0} because average firing rate {1} > 100".format(w_index, samp_results['average firing rate']))
            continue
        # print("average firing rate is {0}".format(samp_results['average firing rate']))
        # refractory = 1
        # spike_indices = samp_results['spikemon indices']
        # unique, counts = np.unique(spike_indices, return_counts=True)
        # if np.maximum(counts) > (samp_results['simulation_time']/refractory*0.5):
        #     print("skipped network {0} because a neuron was firing too much".format(w_index))
        #     continue
        
        # results_filename = "{0}Results_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,w_index) 
        # with open(results_filename, 'rb') as sf:
        #     try:
        #         stats = pickle.load(sf) # load in the stats for the W matrix (L, p_hat, alpha values, alpha_hat values)
        #     except (EOFError):
        #         print("unpickling error")
        
        actual_index += 1
        # ind = w_index-start_index
        
        for key in results:
            results[key][actual_index]=samp_results[key]

        if samp_results["alpha_chain_hat"]>=0 and samp_results["event_rate"]<=200 and samp_results["alpha_conv_hat"]>=0.3:
            print("\nindex {0} has \nlow event_rate {1} \nand high alpha_chain {2}\n".format(w_index, samp_results["event_rate"], samp_results["alpha_chain_hat"]))

        # EXCLUSIONS 
        # if w_index not in [23,27,39,42,71,77,95]: #for Regular sym_alphas_all_rand
        # if w_index not in [23,27,39]: #for Irregular sym_alphas_all_rand
        # if w_index not in [8]:
            
        #     print("w_index = {0}".format(w_index))

        #     samp_results = ar.load_results(N, p_AVG, w_index, Style, data_dir)
            
        #     ind = w_index-start_index
            
        #     for key in results:
        #         results[key][ind]=samp_results[key]

        #     if samp_results["alpha_chain_hat"]>=0 and samp_results["event_rate"]<=100:
        #         print("\nindex {0} has \nlow event_rate {1} \nand high alpha_chain {2}\n".format(w_index, samp_results["event_rate"], samp_results["alpha_chain_hat"]))
        # END EXCLUSIONS


    for key in results:
        results[key].resize(actual_index+1, refcheck=False)
    # save results (pickle new stats dictionary)
    with open(summary_filename, "wb") as rf:
        pickle.dump(results, rf)

kurtosis_mean = np.nanmean(results['IEI excess_kurtosis'])
print("mean kurtosis: {0}".format(kurtosis_mean))



# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', size=80)
# plt.rc('xtick', labelsize=50)
# plt.rc('ytick', labelsize=50)
# plt.figure(figsize=(12,10))
# plt.hist(results['IEI excess_kurtosis'], 50) # data, number of bins
# plt.xlabel('IEI excess kurtosis')
# plt.ylabel('event rate')



# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', size=80)
# plt.rc('xtick', labelsize=50)
# plt.rc('ytick', labelsize=50)
# plt.figure(figsize=(12,10))
# plt.hist(results['IEI skew'], 50) # data, number of bins
# plt.xlabel('IEI skew')
# plt.ylabel('event rate')



# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', size=60)
# plt.rc('xtick', labelsize=50)
# plt.rc('ytick', labelsize=50)
# plt.figure()
# # plt.subplot(222)
# plt.semilogx(results['L_left'], results['event_rate'], 'o', markersize=30)
# plt.xlabel(r'$\sigma$')
# plt.ylabel('event rate')
# # plt.tight_layout()


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=80)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure(figsize=(12,10))
plt.plot(results['alpha_chain_hat'], results['event_rate'], 'o', markersize=30)
plt.xlabel(r'$\alpha_{chain}$')
plt.ylabel('event rate')
plt.xticks(np.arange(-0.5,0.5,0.2))
# plt.tight_layout()



plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=80)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure(figsize=(12,10))
plt.plot(results['alpha_conv_hat'], results['event_rate'], 'o', markersize=30)
plt.xlabel(r'$\alpha_{conv}$')
plt.ylabel('event rate')
plt.xticks(np.arange(0,0.5,0.1))
# plt.tight_layout()




plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=80)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure(figsize=(12,10))
plt.plot(results['alpha_div_hat'], results['event_rate'], 'o', markersize=30)
plt.xlabel(r'$\alpha_{div}$')
plt.ylabel('event rate')
plt.xticks(np.arange(0,0.5,0.1))
# plt.tight_layout()



plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=100)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure(figsize=(12,10))
plt.scatter(results['alpha_conv_hat'], results['alpha_chain_hat'], c=results['event_rate'],s=400)
plt.xlabel(r'$\alpha_{conv}$')
plt.ylabel(r'$\alpha_{chain}$')
cb=plt.colorbar()
# plt.clim(0,3)
cb.set_label('event rate')


plt.figure()
plt.scatter(results['alpha_div_hat'], results['alpha_conv_hat'], c=results['event_rate'],s=400)
plt.xlabel(r'$\alpha_{div}$')
plt.ylabel(r'$\alpha_{conv}$')
cb=plt.colorbar()
# plt.clim(0,3)
cb.set_label('event rate')


plt.figure()
plt.scatter(results['alpha_div_hat'], results['alpha_chain_hat'], c=results['event_rate'],s=400)
plt.xlabel(r'$\alpha_{div}$')
plt.ylabel(r'$\alpha_{chain}$')
cb=plt.colorbar()
# plt.clim(0,3)
cb.set_label('event rate')

# plt.figure(figsize=())

plt.show()