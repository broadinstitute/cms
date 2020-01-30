#this weird effect would be explained if in the old likelihood tables there is the feature that no bin has greater than 20-fold difference (?)
# 11.02.16

import numpy as np

def read_old_like_tsv(infilename):
	ihs, stdif, meanfst, freqdif, maxxpop = [], [], [], [], []
	openfile = open(infilename, 'r')
	openfile.readline() #header
	for line in openfile:
		entries = line.split()
		ihs.append(float(entries[0]))
		stdif.append(float(entries[1]))
		meanfst.append(float(entries[2]))	
		freqdif.append(float(entries[3]))
		maxxpop.append(float(entries[4]))
	openfile.close()
	return ihs, stdif, meanfst, freqdif, maxxpop

def get_min_bf(scores_miss, scores_hit):
	allbfs = []
	for ibin in range(len(scores_miss)):

		thishit, thismiss = scores_hit[ibin], scores_miss[ibin]
		#if True:
		if thishit != 1e-10 and thismiss != 1e-10:
			bf = (thishit/thismiss)
			allbfs.append(bf)
	minbf = min(allbfs)
	return minbf

for pop in [1, 4, 5]:
	#for freq in ['lo', 'hi']:
	#	misslikesfilename = '/idi/sabeti-data/ilya/gsvn/Data/Ilya_Data/sim/newagesfx/likes/missLikes_nov13_10ky_' + str(pop) + '_' + freq + '.tsv'
	#	hitslikesfilename = '/idi/sabeti-data/ilya/gsvn/Data/Ilya_Data/sim/newagesfx/likes/hitsLikes_nov13_10ky_' + str(pop) + '_' + freq + '.tsv'

	if True:

		hitslikesfilename = "/idi/sabeti-data/ilya/sabetilab/ilya/new/nsvn/Data/Ilya_Data/sim/newages/likes/hitsLikes_toneutFixed_" + str(pop) + ".tsv"
		misslikesfilename = "/idi/sabeti-data/ilya/sabetilab/ilya/new/nsvn/Data/Ilya_Data/sim/newages/likes/missLikes_toneutFixed_" + str(pop) + ".tsv"
		ihs_miss, stdif_miss, meanfst_miss, freqdif_miss, maxxpop_miss = read_old_like_tsv(misslikesfilename)
		ihs_hit, stdif_hit, meanfst_hit, freqdif_hit, maxxpop_hit = read_old_like_tsv(hitslikesfilename)

		min_bf_ihs = (get_min_bf(ihs_miss, ihs_hit))
		min_bf_stdif = (get_min_bf(stdif_miss, stdif_hit))
		min_bf_fst = (get_min_bf(meanfst_miss, meanfst_hit))
		min_bf_freqdir = (get_min_bf(freqdif_miss, freqdif_hit))		
		min_bf_xp = (get_min_bf(maxxpop_miss, maxxpop_hit))			

		comp = 1
		for min_bf in [min_bf_ihs, min_bf_stdif, min_bf_fst, min_bf_freqdir, min_bf_xp]:
			print(np.log(min_bf))
			#comp *= min_bf
		print(np.log(comp))
