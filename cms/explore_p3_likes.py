#this weird effect would be explained if in the old likelihood tables there is the feature that no bin has greater than 20-fold difference 
# 11.02.16
import numpy as np
import os
def read_new_likefile(infilename):
	vals = []
	openfile = open(infilename, 'r')
	for line in openfile:
		entries = line.split()
		vals.append(float(entries[-1]))
	openfile.close()
	return vals

def get_min_bf(scores_miss, scores_hit):
	allbfs = []
	for ibin in range(len(scores_miss)):

		thishit, thismiss = scores_hit[ibin], scores_miss[ibin]
		if thishit != 1e-10 and thismiss !=1e-10:
			bf = (thishit/thismiss)
			allbfs.append(bf)
	minbf = min(allbfs)
	return minbf
for model in ['nulldefault_constantsize', 'default_112115_825am', 'nulldefault','default_default_101715_12pm', 'gradient_101915_treebase_6_best',]:
	for pop in [1, 2, 3, 4]:
		for freq in ['allfreq']:#'low', 'mid', 'hi']:
			print('\n')
			clr_sum = 0 
			for score in ['ihs', 'delihh', 'deldaf', 'fst', 'xpehh']:
				if score in ['ihs', 'delihh', 'nsl']:
					misslikesfilename = '/idi/sabeti-scratch/jvitti/likes_111516_b/'  + model + '/' + score + '/likes_sel' + str(pop) + '_' + str(freq) + '_neut.txt'
					hitslikesfilename = '/idi/sabeti-scratch/jvitti/likes_111516_b/'  + model + '/' + score + '/likes_sel' + str(pop) + '_' + str(freq) + '_causal.txt'
				else:
					misslikesfilename = '/idi/sabeti-scratch/jvitti/likes_111516_b/'  + model + '/' + score + '/likes_sel' + str(pop) + '_choose_' + str(freq) + '_neut.txt'
					hitslikesfilename = '/idi/sabeti-scratch/jvitti/likes_111516_b/'  + model + '/' + score + '/likes_sel' + str(pop) + '_choose_' + str(freq) + '_causal.txt'
					
				#print(misslikesfilename)
				#print(hitslikesfilename)
				if os.path.isfile(hitslikesfilename) and os.path.isfile(misslikesfilename):
					miss = read_new_likefile(misslikesfilename)
					hit = read_new_likefile(hitslikesfilename)
					#print(miss)
					#print(hit)
					min_bf=	get_min_bf(miss, hit)
					print(np.log(min_bf))
					clr_sum += np.log(min_bf)
			print("\t" + str(clr_sum))


