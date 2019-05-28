## -*- coding: utf-8 -*-
"""
Created on Tues Feb  26 14:25 2019

@author: Brittany
"""
# data_dir = 'matrices/N1000_Linf_recurr_alpha_div_rand/'
# data_dir = 'matrices/N3000_LL50_LR50_recurr_alpha_div_rand/'
data_dir = 'matrices/N3000_LL100_LR0_ff_alpha_conv_div_rand/'
# data_dir = 'matrices/N3000_erdos_renyi/'
# res_dir = '/var/tmp/N3000_LL70_LR0_ff_alphas_all_rand/'
# data_dir = 'matrices/N1000_LL50_LR50_recurr_alpha_div_rand/'
# data_dir = 'matrices/N1000_LL100_LR0_ff_alpha_conv_div_rand/'
res_dir = data_dir
print("data_dir: "+data_dir)
# print("results_dir: "+res_dir)
style1 = "Regular5s_"
# style2 = "Irregular50s_"
print("styles: "+style1)#+", "+style2)

import sys

try:
   import cPickle as pickle ## used to store python data types
except:
   import pickle
#import dill ## pickle works fine

#from analyze import analyze_autocor ## used to analyze synchrony of networks
## we use our synchronous event detection algorithm instead:
import analyze_results as ar

try:
    del input ## brian overwrites this, we want to reset it to the python default
except:
    pass
input_orig = input ## rename the python default for input (brian will overwrite it when imported)



N = 3000 ## Number of excitatory neurons
p = 50/N ## average probability of connectivity between neurons

if len(sys.argv) >= 3:
   start_index = int(sys.argv[1])
   end_index = int(sys.argv[2])
else:
   start_index = int(input_orig("enter starting index: "))
   end_index = int(input_orig("enter end index: "))

## now, let's update the chosen simulations
for w_index in range(start_index, end_index+1):
	print("\n\nupdating network {0}".format(w_index))
	print("data_dir {0}\n".format(data_dir))

	try:
		ar.update_results(N, p, w_index, style1, data_dir)
		# ar.update_results(N, p, w_index, style2, data_dir)
	except Exception as e:
		print("Error: %s" % e)