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
data_dir = "matrices/N3000_LL70_LR0_ff_alphas_all_rand/"
# data_dir = "matrices/N3000_LL70_LR0_ff_alpha_conv_div_rand/"
# data_dir = "matrices/N3000_LL70_LR70_sym_alpha_div_rand/"
# style = "FF_L70"
style = "CnsClean"
# style = "ttLongClean"

reload=False

summary_filename = "{0}Summary_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,style) 


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
        ('event_rate', np.zeros(n_indices)), ('IEI excess_kurtosis', np.zeros(n_indices)),]) 
        #('j', np.zeros(n_indices)),('ext_rate', np.zeros(n_indices)),('ext_mag', np.zeros(n_indices))])
    
    for w_index in range(start_index, end_index+1):
        
        print("w_index = {0}".format(w_index))

        samp_results = ar.load_results(N, p_AVG, w_index, style, data_dir)
        
        # results_filename = "{0}Results_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p_AVG,w_index) 

        # with open(results_filename, 'rb') as sf:
        #     try:
        #         stats = pickle.load(sf) # load in the stats for the W matrix (L, p_hat, alpha values, alpha_hat values)
        #     except (EOFError):
        #         print("unpickling error")
        
        ind = w_index-start_index
        
        for key in results:
            results[key][ind]=samp_results[key]

    # save results (pickle new stats dictionary)
    with open(summary_filename, "wb") as rf:
        pickle.dump(results, rf)

kurtosis_mean = np.nanmean(results['IEI excess_kurtosis'])
print("mean kurtosis: {0}".format(kurtosis_mean))

# results['event_rate'] = np.maximum(0,results['event_rate'])

#inds = inds=results['L_left'] > 200
#for key in results:
#    results[key] = results[key][inds]




plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=100)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure()
# plt.subplot(221)
plt.plot(results['alpha_chain_hat'], results['event_rate'], 'o', markersize=30)
plt.xlabel(r'$\alpha_{chain}$')
plt.ylabel('event rate')
# plt.tight_layout()



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
plt.rc('font', family='serif', size=100)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure()
# plt.subplot(223)
plt.plot(results['alpha_conv_hat'], results['event_rate'], 'o', markersize=30)
plt.xlabel(r'$\alpha_{conv}$')
plt.ylabel('event rate')
plt.xticks(np.arange(0,0.5,0.1))
# plt.tight_layout()




plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=100)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure()
# plt.subplot(224)
plt.plot(results['alpha_div_hat'], results['event_rate'], 'o', markersize=30)
plt.xlabel(r'$\alpha_{div}$')
plt.ylabel('event rate')
plt.xticks(np.arange(0,0.5,0.1))
# plt.tight_layout()



plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=100)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure()
plt.scatter(results['alpha_conv_hat'], results['alpha_chain_hat'], c=results['event_rate'],s=400)
plt.xlabel(r'$\alpha_{conv}$')
plt.ylabel(r'$\alpha_{chain}$')
cb=plt.colorbar()
plt.clim(0,25)
cb.set_label('event rate')


plt.figure()
plt.scatter(results['alpha_div_hat'], results['alpha_conv_hat'], c=results['event_rate'],s=400)
plt.xlabel(r'$\alpha_{div}$')
plt.ylabel(r'$\alpha_{conv}$')
cb=plt.colorbar()
plt.clim(0,25)
cb.set_label('event rate')


plt.figure()
plt.scatter(results['alpha_div_hat'], results['alpha_chain_hat'], c=results['event_rate'],s=400)
plt.xlabel(r'$\alpha_{div}$')
plt.ylabel(r'$\alpha_{chain}$')
cb=plt.colorbar()
plt.clim(0,25)
cb.set_label('event rate')

plt.show()