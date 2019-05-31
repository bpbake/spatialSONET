## -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:06:45 2018

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
from scipy import stats

import matplotlib.pyplot as plt

import analyze_results as ar

N = 3000 ## Number of excitatory neurons
p =50/N ## average probability of connectivity between neurons

# data_dir = 'matrices/N3000_LL100_LR0_ff_alphas_all_rand/'
# data_dir = 'matrices/N3000_LL50_LR50_recurr_alphas_all_rand/'
data_dir = 'matrices/N3000_Linf_homogeneous_alphas_all_rand/'
print("data_dir: "+data_dir)

# style = "Regular5s_Clean_"
style = "Irregular50s_Clean_"
print("style: "+style)

reload=False
# reload=True

summary_filename = "{0}Summary_W_N{1}_p{2}_{3}.pickle".format(data_dir,N,p,style) 


if reload:
    with open(summary_filename, "rb") as rf:
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
        ('alpha_chain_hat', np.zeros(n_indices)), 
        ('average firing rate', np.zeros(n_indices)), ('event_rate', np.zeros(n_indices)), 
        ('IEI excess_kurtosis', np.zeros(n_indices)), ('IEI skew', np.zeros(n_indices))]) 
    actual_index = -1

    skipped = 0
    for w_index in range(start_index, end_index+1):
        if w_index%30 == 0:    
            print("w_index = {0}".format(w_index))

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
            results[key][actual_index]=samp_results[key]

        if samp_results["alpha_chain_hat"]>=0.05 and samp_results["event_rate"]<=5 and samp_results["alpha_conv_hat"]>=0.3:
            print("\nindex {0} has \nlow event_rate {1} \nand high alpha_chain {2}\n".format(
                w_index, samp_results["event_rate"], samp_results["alpha_chain_hat"]))

        # if samp_results["alpha_conv_hat"]>=0.4 and samp_results["event_rate"]<=150:
        #     print("\nhigh-heel situation! \nindex {0} has \nlow event_rate {1} \nand high alpha_conv {2}\n".format(
        #         w_index, samp_results["event_rate"], samp_results["alpha_conv_hat"]))


    for key in results:
        results[key].resize(actual_index+1, refcheck=False)

    results["skipped"] = skipped
    ## save results (pickle new stats summary dictionary for future plots)
    with open(summary_filename, "wb") as rf:
        pickle.dump(results, rf)
    


print("\n{0} networks were skipped due to saturation or high average firing rate".format(results['skipped']))

kurtosis_mean = np.nanmean(results['IEI excess_kurtosis'])
print("\nmean kurtosis: {0}".format(kurtosis_mean))



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

# print("stats are: slope, intercept, r_value, p_value, std_err")


## Plot event rate vs. alpha_div
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=80)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure(figsize=(12,10))
# plt.suptitle('style: {0}, data_dir: {1}'.format(style, data_dir))
plt.plot(results['alpha_div_hat'], results['event_rate'], 'o', markersize=30)
plt.xlabel(r'$\hat\alpha_{\mathrm{div}}$')
plt.ylabel('event rate')
plt.xticks(np.arange(0.1,0.6,0.1))
plt.xlim(left=0.02)
# plt.yticks(np.arange(0,140,30))
# plt.ylim(-10,140)
# plt.yticks(np.arange(0,1.25,.5))
# plt.ylim(-0.1,1.4)
plt.tight_layout()
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
print("event rate vs. alpha_div stats: {0}".format(stats.linregress(results['alpha_div_hat'],results['event_rate'])))


## Plot event rate vs. alpha_conv
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=80)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure(figsize=(12,10))
# plt.suptitle('style: {0}, data_dir: {1}'.format(style, data_dir))
plt.plot(results['alpha_conv_hat'], results['event_rate'], 'o', markersize=30)
plt.xlabel(r'$\hat\alpha_{\mathrm{conv}}$')
plt.ylabel('event rate')
plt.xticks(np.arange(0.1,0.6,0.1))
plt.xlim(left=0.02)
# plt.yticks(np.arange(0,140,30))
# plt.ylim(-10,140)
# plt.yticks(np.arange(0,170,50))
# plt.ylim(-15,180)
# plt.yticks(np.arange(0,10,2))
# plt.ylim(-0.5,10)
plt.tight_layout()
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
print("event rate vs. alpha_conv stats: {0}".format(stats.linregress(results['alpha_conv_hat'],results['event_rate'])))

## Plot event rate vs alpha_chain
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=80)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure(figsize=(12,10))
# plt.suptitle('style: {0}, data_dir: {1}'.format(style, data_dir))
plt.plot(results['alpha_chain_hat'], results['event_rate'], 'o', markersize=30)
plt.xlabel(r'$\hat\alpha_{\mathrm{chain}}$')
plt.ylabel('event rate')
plt.xticks(np.arange(-0.4,0.4,0.2))
plt.xlim(left=-0.45)
# plt.yticks(np.arange(0,140,30))
# plt.ylim(-10,140)
# plt.yticks(np.arange(0,170,50))
# plt.ylim(-15,180)
plt.tight_layout()
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
print("event rate vs. alpha_chain stats: {0}".format(stats.linregress(results['alpha_chain_hat'],results['event_rate'])))




## Plot alpha_conv vs. alpha_div with event rate color
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=80)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure(figsize=(12,10))
# plt.suptitle('style: {0}, data_dir: {1}'.format(style, data_dir))
plt.scatter(results['alpha_div_hat'], results['alpha_conv_hat'], c=results['event_rate'],s=400)
plt.xlabel(r'$\hat\alpha_{\mathrm{div}}$')
plt.xticks(np.arange(0.1,0.6,0.1))
# plt.xlim(left=0.02)
plt.ylabel(r'$\hat\alpha_{\mathrm{conv}}$')
plt.yticks(np.arange(0.1,0.6,0.1))
# plt.ylim(bottom=0.02)
# plt.clim(0,90)
cb=plt.colorbar()
cb.set_label('event rate')
plt.tight_layout()
mng = plt.get_current_fig_manager()
mng.window.showMaximized()



## Plot alpha_chain vs. alpha_conv with event rate color
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=80)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure(figsize=(12,10))
# plt.suptitle('style: {0}, data_dir: {1}'.format(style, data_dir))
plt.scatter(results['alpha_conv_hat'], results['alpha_chain_hat'], c=results['event_rate'],s=500)
plt.xlabel(r'$\hat\alpha_{\mathrm{conv}}$')
plt.xticks(np.arange(0.1,0.6,0.1))
# plt.xlim(left=0.02)
plt.ylabel(r'$\hat\alpha_{\mathrm{chain}}$')
plt.yticks(np.arange(-0.4,0.6,0.2))
# plt.ylim(bottom=-0.45)
# plt.clim(0,90)
cb=plt.colorbar()
cb.set_label('event rate')
plt.tight_layout()
mng = plt.get_current_fig_manager()
mng.window.showMaximized()




## Plot alpha_chain vs. alpha_div with event rate color
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=80)
plt.rc('xtick', labelsize=50)
plt.rc('ytick', labelsize=50)
plt.figure(figsize=(12,10))
# plt.suptitle('style: {0}, data_dir: {1}'.format(style, data_dir))
plt.scatter(results['alpha_div_hat'], results['alpha_chain_hat'], c=results['event_rate'],s=400)
plt.xlabel(r'$\hat\alpha_{\mathrm{div}}$')
plt.xticks(np.arange(0.1,0.6,0.1))
# plt.xlim(left=0.02)
plt.ylabel(r'$\hat\alpha_{\mathrm{chain}}$')
plt.yticks(np.arange(-0.4,0.6,0.2))
# plt.ylim(bottom=-0.45)
# plt.clim(0,90)
cb=plt.colorbar()
cb.set_label('event rate')
plt.tight_layout()
mng = plt.get_current_fig_manager()
mng.window.showMaximized()

# plt.figure(figsize=())

plt.show()