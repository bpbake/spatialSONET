## -*- coding: utf-8 -*-
"""
Created on Thurs May 16 15:54 2019

@author: Brittany
"""

try:
   import cPickle as pickle ## used to store python data types
except:
   import pickle


try:
    del input ## brian overwrites this, we want to reset it to the python default
except:
    pass
input_orig = input ## rename the python default for input (brian will overwrite it when imported)

import numpy as np
np.set_printoptions(threshold=np.nan)

import matplotlib.pyplot as plt

import analyze_results as ar

N = 3000 ## Number of excitatory neurons
p =50/N ## average probability of connectivity between neurons

# data_dir = 'matrices/N3000_LL100_LR0_ff_alphas_all_zero/'
# data_dir = 'matrices/N3000_LL50_LR50_recurr_alphas_all_zero/'
data_dir = 'matrices/N3000_Linf_homogeneous_alphas_all_zero/'
print("data_dir: "+data_dir)

style = "Regular5s_Clean_"
# style = "Irregular100s_Clean_"
print("style: "+style)

reload=False
# reload=True

iei_filename = "{0}IEI_Summary_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p,style) 


if reload:
    with open(iei_filename, "rb") as rf:
        results = pickle.load(rf)

else:
    ## Load in matrices one at a time
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
        ('average firing rate', np.zeros(n_indices)), ('event_rate', np.zeros(n_indices)), 
        ('IEI excess_kurtosis', np.zeros(n_indices)), ('IEI skew', np.zeros(n_indices)),#])
        # ('IEIs', np.zeros(n_indices))])#,
        ('IEI mean', np.zeros(n_indices)),('IEI std', np.zeros(n_indices)), ('IEI coeff of variation', np.zeros(n_indices))]) 
        #('j', np.zeros(n_indices)),('ext_rate', np.zeros(n_indices)),('ext_mag', np.zeros(n_indices))])
    actual_index = -1

    skipped = 0
    num_no_events = 0
    for w_index in range(start_index, end_index+1):
        if w_index%30 == 0:    
            print("\nw_index = {0}".format(w_index))

        try:
            samp_results = ar.load_results(N, p, w_index, style, data_dir)
        except:
            print("couldn't load {0}Results_W_N{1}_p{2}_{3}{4}.pickle".format(data_dir,N,p,style,w_index))
            continue

        try:
            if samp_results['average firing rate'] > 100:
                skipped += 1
                print("skipped index {0} because average firing rate {1} > 100".format(w_index, samp_results['average firing rate']))
                continue
        except:
            if samp_results['saturated']:
                skipped += 1
                print("skipped index {0} because saturated".format(w_index))
                continue
        
        actual_index += 1
        
        for key in results:
          try:
            results[key][actual_index]=samp_results[key]
          except:
            continue

        # if not samp_results['IEIs']:
        #   num_no_events += 1
        #   # imean = samp_results['simulation_time']
        #   imean = float('nan')
        #   istd = float('nan')
        #   # istd = float('inf')
        #   # print('no events in index {0}. IEI stats not computed.'.format(w_index))
        # else:
        #   imean=np.mean(samp_results["IEIs"])
        #   istd=np.std(samp_results["IEIs"])
        # icoeffvar = np.true_divide(istd,imean)
        # results["IEI mean"][actual_index]=imean
        # results["IEI std"][actual_index]=istd
        # results["IEI coeff of variation"][actual_index]=icoeffvar
        # # print('IEI coeff of variation for index {0}: {1}'.format(w_index, icoeffvar))   
        # print('IEI excess kurtosis for index {0}: {1}'.format(w_index, samp_results['IEI excess_kurtosis']))  
        # print('IEI mean: {0}'.format(imean))    
        print("IEI skew: {0}".format(samp_results['IEI skew']))


    for key in results:
        results[key].resize(actual_index+1, refcheck=False)

    results["skipped"] = skipped
    results['num_no_events'] = num_no_events
    ## save results (pickle new stats summary dictionary for future plots)
    with open(iei_filename, "wb") as rf:
        pickle.dump(results, rf)

print("\n{0} networks were skipped due to saturation or high average firing rate".format(results['skipped']))
print("there were {0} simulations with no events".format(results['num_no_events']))

print("\ndata_dir: {0}".format(data_dir))
print("Style: {0}".format(style))
# print("average IEI mean: {0}".format(np.nanmean(results['IEI mean'])))
print("IEI excess kutosis: {0}".format(results['IEI excess_kurtosis']))
# print("average IEI excess kurtosis: {0}".format(np.nanmean(results['IEI excess_kurtosis'])))
print("IEI coefficients of variation: {0}".format(results['IEI coeff of variation']))
# print("mean IEI coefficient of variation: {1}".format(style, np.nanmean(results["IEI coeff of variation"])))
print("\nIEI skew: {0}".format(results['IEI skew']))
print("mean IEI skew: {0}".format(np.nanmean(results['IEI skew'])))

print("\nevent rates: {0}".format(results['event_rate']))
print("mean event rate: {0}".format(np.nanmean(results['event_rate'])))