try:
   import cPickle as pickle # used to store python data types
except:
   import pickle
#import dill #pickle works fine

import analyze_results as ar


N = 3000 # Number of excitatory neurons
p_AVG =50/N # average probability of connectivity between neurons

# data_dir = "matrices/N3000_LL70_LR70_sym_alpha_conv_div_rand/"
data_dir = "matrices/N3000_LL70_LR0_ff_alpha_conv_rand/"

print("data_dir: "+data_dir)

# Styles:  Cns = Irregular, CnsL = Regular
# NewStyle = "Irregular50s_"
# OldStyle = "Cns" 
NewStyle = "Regular5s_"
OldStyle = "CnsL"
print("Old Style: "+OldStyle)
print("New Style: "+NewStyle)

start_index = int(input("enter starting index: "))
end_index = int(input("enter end index: "))

for w_index in range(start_index, end_index+1): #i=start_index,start_index+1,...,start_index+num_matrices-1
	if (w_index%10)==0:
		print("w_index = {0}".format(w_index))

	if w_index not in [23,27,39,42,71,77,95]: # for Regular sym_alphas_all_rand
	# if w_index not in [23,27,39]: #for Irregular sym_alphas_all_rand
		Results = ar.load_results(N, p_AVG, w_index, OldStyle, data_dir)
		ar.save_results(N, p_AVG, w_index, Results, NewStyle, data_dir)

		CleanResults = ar.load_results(N, p_AVG, w_index, OldStyle+"Clean", data_dir)
		ar.save_results(N, p_AVG, w_index, CleanResults, NewStyle+"Clean_", data_dir)

	# Results = ar.load_results(N, p_AVG, w_index, OldStyle, data_dir)
	# ar.save_results(N, p_AVG, w_index, Results, NewStyle, data_dir)

	# CleanResults = ar.load_results(N, p_AVG, w_index, OldStyle+"Clean_", data_dir)
	# ar.save_results(N, p_AVG, w_index, CleanResults, NewStyle+"Clean_", data_dir)
