# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 00:02:45 2017

@author: Brittany

# to reimport a module, use importlib.reload(module), you must first import importlib
"""

def results(N, p, start_index, end_index, style, data_dir='matrices/'):
    try:
       import cPickle as pickle # used to store python data types
    except:
       import pickle
       
    import numpy as np
    from scipy import stats

    import matplotlib.pyplot as plt
    import matplotlib

    # try:
    #     del input # brian overwrites this, we want to reset it to the python default
    # except:
    #     pass
    # input_orig = input

    # data_dir = 'matrices/N10000_LL70_LR0_ff_randadiv/'

    # N = 10000
    # p = 50/N

    results_header = ['w_index', 'N', 'L_left', 'L_right', #'largest eigenvalue', 
    'alpha_chain', 'alpha_chain_hat', 'alpha_conv', 'alpha_conv_hat', 'alpha_div', 'alpha_div_hat', 
    'alpha_recip', 'alpha_recip_hat', 'p_AVG', 'p_hat', 
    'simulation_time', 'event_rate', 'event_mag', 'IEI excess_kurtosis', 'IEI skew']

    results_type = [('w_index', int), ('alpha_div_hat', float), ('event_rate', float), 
    ('event_mag', float), ('IEI excess_kurtosis', float), ('IEI skew', float)]
    #('N', int), ('L_left', float), ('L_right', float), #('largest eigenvalue', complex list), 
    #('alpha_chain', float), ('alpha_chain_hat', float), ('alpha_conv', float), ('alpha_conv_hat', float), 
    #('alpha_div', float), , ('alpha_recip', float), ('alpha_recip_hat', float), 
    #('p_AVG', float), ('p_hat', float), ('simulation_time', float), #('IEIs', list),


    #start_index = int(input_orig("enter a starting index: "))
    #end_index = int(input_orig("enter end index: "))

    results_matrix = np.zeros((end_index-start_index+1, len(results_header)))
    result_list = []

    alpha_divs = []
    alpha_chains = []
    # alpha_recips = []
    # alpha_convs = []
    event_rates = []
    event_mags = []
    IEIkurtosis = []
    IEIskews = []
                                                      
    for w_index in range(start_index, end_index+1):
       
        # results_filename = "{0}Results_W_N{1}_p{2}_slower{3}.pickle".format(data_dir,N,p,w_index)
        results_filename = "{0}Results_W_N{1}_p{2}_{3}{4}.pickle".format(data_dir, N, p, style, w_index)
        with open(results_filename, 'rb') as rf:
            try:
                results = pickle.load(rf) 
            except (EOFError):
                print("unpickling error")
              
        results_matrix[w_index-start_index, 0] = w_index
        for k in range(1, len(results_header)):    
            results_matrix[w_index-start_index, k] = results[results_header[k]]

        result = (w_index, results['alpha_div_hat'], results['event_rate'], results['event_mag'], 
            results['IEI excess_kurtosis'], results['IEI skew'])
        result_list.append(result)

        alpha_divs.append(results['alpha_div_hat'])
        alpha_chains.append(results['alpha_chain_hat'])
        # alpha_recips.append(results['alpha_recip_hat'])
        # alpha_convs.append(results['alpha_conv_hat'])
        event_rates.append(results['event_rate'])
        event_mags.append(results['event_mag'])
        IEIkurtosis.append(results['IEI excess_kurtosis'])
        IEIskews.append(results['IEI skew'])

    result_summary = np.array(result_list, results_type)

    # r_filename = "{0}Result_Matrices_{1}_{2}.pickle".format(data_dir, start_index, end_index) 
    # with open(r_filename, "wb") as rf:
    #     pickle.dump(np.array(result_list, results_type), rf)

    matrix_filename = '{0}Result_Matrices_{1}_{2}_{3}'.format(data_dir, start_index, end_index, style)
    np.savetxt('{0}.csv'.format(matrix_filename), results_matrix, delimiter=',', 
        header=str(results_header)) 

    # configure the plots
    matplotlib.rcParams.update({'font.size': 20})
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    plt.suptitle(matrix_filename)

    xp = np.linspace(0, 0.5, 100)

    slope, intercept, r_value, p_value, std_err = stats.linregress(alpha_divs, event_rates)
    z = np.poly1d(np.array([slope, intercept]))
    print('Best fit for alpha_div vs. event_rate: {0}'.format([slope, intercept]))
    print('r-squared:{0}'.format(r_value**2))
    plt.subplot(221)
    plt.plot(alpha_divs, event_rates, 'o', xp, z(xp), '-')
    plt.ylabel('event_rate')
    plt.ylim(0,0.02)
    plt.xlabel('alpha_div_hat')
    plt.grid(True)

    slope, intercept, r_value, p_value, std_err = stats.linregress(alpha_divs, event_mags)
    z = np.poly1d(np.array([slope, intercept]))
    print('Best fit for alpha_div vs. event_mag: {0}'.format([slope, intercept]))
    print('r-squared:{0}'.format(r_value**2))
    plt.subplot(223)
    plt.plot(alpha_divs, event_mags, 'o', xp, z(xp), '-')
    plt.ylabel('event_mag')
    plt.ylim(0,0.0125)
    plt.xlabel('alpha_div_hat')
    plt.grid(True)

    plt.subplot(222)
    plt.plot(alpha_chains, event_rates, 'o')#, xp, pcr(xp), '-')
    plt.ylabel('event_rate')
    plt.ylim(0,0.02)
    plt.xlabel('alpha_chain_hat')
    plt.grid(True)

    plt.subplot(224)
    plt.plot(alpha_divs, IEIkurtosis, 'o')
    plt.ylabel('IEI excess kurtosis')
    plt.ylim(-4, 20)
    plt.xlabel('alpha_div_hat')
    plt.grid(True)

    plt.show()

    # plt.plot(IEIkurtosis, IEIskews, 'o')
    # plt.xlabel('IEI excess kurtosis')
    # plt.ylabel('IEI skew')
    # plt.grid(True)

    # plt.show()

    return result_summary