##	functions for analyzing empirical/simulated CMS output
##	last updated 04.13.2017		vitti@broadinstitute.org

import matplotlib as mp 
mp.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import math

###################
## DEFINE SCORES ##
###################
def write_master_likesfile(writefilename, model, selpop, freq,basedir,  miss = "neut",):
	'''adapted from run_likes_func.py'''
	writefile = open(writefilename, 'w')
	for score in ['ihs', 'nsl', 'delihh']: 
		hitlikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_" + str(freq) + "_causal.txt"#_smoothed.txt"
		misslikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_" + str(freq) + "_" + miss + ".txt"#"_smoothed.txt"
		#assert(os.path.isfile(hitlikesfilename) and os.path.isfile(misslikesfilename))
		writefile.write(hitlikesfilename + "\n" + misslikesfilename + "\n")
	for score in ['xpehh', 'fst', 'deldaf']:
		hitlikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_choose_" + str(freq) + "_causal.txt"#_smoothed.txt"
		misslikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_choose_" + str(freq) + "_" + miss + ".txt"#"_smoothed.txt"
		#assert(os.path.isfile(hitlikesfilename) and os.path.isfile(misslikesfilename))
		writefile.write(hitlikesfilename + "\n" + misslikesfilename + "\n")
	writefile.close()
	print("wrote to: " + writefilename)
	return
'''def normalize_global(rawscore, mean, sd):
	#cms_gw 
	rawscore, mean, sd = float(rawscore), float(mean), float(sd)
	normalizedvalue = (rawscore - mean) #/ sd
	return normalizedvalue
'''
def normalize_local(values, physpos): #called by regionviz - replace?

	clean_vals, clean_phys = [], []
	for i in range(len(values)):
		item = values[i]
		phys = physpos[i]
		if not np.isnan(item) and not np.isinf(item):
			clean_vals.append(item)
			clean_phys.append(phys)


	####
	#### NORMALIZE AS Z SCORE
	####
	#mean = np.mean(clean_vals)
	#var = np.var(clean_vals)
	#stddev = var**.5
	#normalized = [(item-mean)/stddev for item in clean_vals]
	#print(values[-1])
	#print(normalized[-1])


	print('this should be replaced by other normalize function.')
	print('currently testing - 5.4.17')
	print("for amended CMS, this should be gratuitous. The values already make sense (qua probabilities), so there's no need to scale the y-axis.")

	####
	#### NORMALIZE BY SCALING TO [0, 1]
	####
	maxVal, minVal = max(clean_vals), min(clean_vals)
	normalized = [(unnormed - minVal)/maxVal for unnormed in clean_vals]
	
	#if len(normalized) !=len(clean_phys):
	#	print(str(len(normalized)))
	#	print(str(len(clean_phys)))
	#assert len(normalized) == len(clean_phys)
	return normalized, clean_phys
###############
## REGION ID ##
###############
def get_window(istart, physpos, scores, windowlen = 100000):
	window_scores = [scores[istart]]
	startpos = physpos[istart]
	pos = startpos
	iscore = istart
	while pos < (startpos + windowlen):
		iscore += 1
		if iscore >= len(scores):
			break
		window_scores.append(scores[iscore])
		pos = physpos[iscore]
	#print(str(pos) + " " + str(startpos))
	return window_scores
def check_outliers(scorelist, cutoff = 3):
	numscores = len(scorelist)
	outliers = [item for item in scorelist if item > cutoff]
	numoutliers = len(outliers)
	percentage = (float(numoutliers) / float(numscores)) * 100.
	return percentage
def check_rep_windows(physpos, scores, windowlen = 100000, cutoff = 3, totalchrlen=1000000):
	'''
	previous implementation: !!!! this is going to result in false positives whenever I have a small uptick right near the edge of the replicate
	'''
	#check window defined by each snp as starting point
	rep_percentages = []

	numSnps = len(physpos)
	numWindows = 0
	#get exhaustive windows and stop at chrom edge
	for isnp in range(numSnps):
		if physpos[isnp] + windowlen < totalchrlen:
			numWindows +=1
		else:
			#print(str(physpos[isnp]) + "\t")
			break		
	for iPos in range(numWindows): 
		window_scores = get_window(iPos, physpos, scores, windowlen)
		percentage = check_outliers(window_scores, cutoff)
		rep_percentages.append(percentage)
	return rep_percentages
def merge_windows(chrom_signif, windowlen, maxGap = 100000):
	starts, ends = [], []
	contig = False
	for i_start in range(len(chrom_signif) - 1):
		if not contig:
			starts.append(chrom_signif[i_start])
		if ((chrom_signif[i_start] + windowlen) > chrom_signif[i_start + 1]): #contiguous
			contig = True
		#or, could also be contiguous in the situation where the next snp is not within this window because there doesn't exist such a snp
		elif chrom_signif[i_start +1] >=(chrom_signif[i_start] + windowlen)  and chrom_signif[i_start +1] < (chrom_signif[i_start] + maxGap):
			contig = True
		else:
			contig = False
		if not contig:
			windowend = chrom_signif[i_start] + windowlen
			ends.append(windowend)
	return starts, ends

##########################
## POWER & SIGNIFICANCE ##
##########################
def calc_pr(all_percentages, threshhold):
	numNeutReps_exceedThresh = 0
	totalnumNeutReps = len(all_percentages)
	for irep in range(totalnumNeutReps):
		if len(all_percentages[irep]) != 0:
			if max(all_percentages[irep]) > threshhold:
				numNeutReps_exceedThresh +=1
	numNeutReps_exceedThresh, totalnumNeutReps = float(numNeutReps_exceedThresh), float(totalnumNeutReps)
	if totalnumNeutReps != 0:
		pr = numNeutReps_exceedThresh / totalnumNeutReps 
	else:
		pr = 0
		print('ERROR; empty set')
	return pr
def get_causal_rank(values, causal_val):
	if np.isnan(causal_val):
		return(float('nan'))
	assert(causal_val in values)
	cleanvals = []
	for item in values:
		if not np.isnan(item) and not np.isinf(item):
			cleanvals.append(item)
	values = cleanvals

	values.sort()
	values.reverse()
	causal_rank = values.index(causal_val)
	return causal_rank
def get_cdf_from_causal_ranks(causal_ranks):
	numbins = max(causal_ranks) #? heuristic
	counts, bins = np.histogram(causal_ranks, bins=numbins, normed = True) #doublecheck
	cdf = np.cumsum(counts)
	return bins, cdf
def get_pval(all_simscores, thisScore):
	r = np.searchsorted(all_simscores,thisScore)
	n = len(all_simscores)
	pval = 1. - ((r + 1.) / (n + 1.)) 
	if pval > 0:
		#pval *= nSnps #Bonferroni
		return pval
	else:
		#print("r: "  +str(r) + " , n: " +  str(n))
		pval =  1. - (r/(n+1))
		#pval *= nSnps #Bonferroni
		return pval

###############
## VISUALIZE ##
###############
def quick_plot(ax, pos, val, ylabel,causal_index=-1):
	
	ax.scatter(pos, val, s=.8)
	if causal_index != -1:
		ax.scatter(pos[causal_index], val[causal_index], color='r', s=4)
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize('6')
	ax.set_ylabel(ylabel, fontsize='6')
	#ax.set_xlim([0, 1500000]) #make flexible?
	ax.yaxis.set_label_position('right')
	#ax.set_ylim([min(val), max(val)])
	return ax
def plot_dist(allvals, savefilename= "/web/personal/vitti/test.png", numBins=1000):
	#print(allvals)

	#get rid of nans and infs
	#cleanvals = [item for item in allvals if not np.isnan(item)]
	#allvals = cleanvals
	allvals = np.array(allvals)
	allvals = allvals[~np.isnan(allvals)]
	allvals = allvals[~np.isinf(allvals)]
	#allvals = list(allvals)
	#print(allvals)
	if len(allvals) > 0:
		f, ax = plt.subplots(1)
		ax.hist(allvals, bins=numBins)
		plt.savefig(savefilename)
		print('plotted to ' + savefilename)
	return
def plotManhattan(ax, neut_rep_scores, emp_scores, chrom_pos, nSnps, maxSkipVal = 0, zscores = True):
	#neut_rep_scores.sort()
	#print('sorted neutral scores...')
	lastpos = 0
	for chrom in range(1,23):
		ichrom = chrom-1
		if ichrom%2 == 0: 
			plotcolor = "darkblue"
		else:
			plotcolor = "lightblue"

		if zscores == True:
			#http://stackoverflow.com/questions/3496656/convert-z-score-z-value-standard-score-to-p-value-for-normal-distribution-in?rq=1 
			#Z SCORE cf SG email 103116
			#pvals = [get_pval(neut_rep_scores, item) for item in emp_scores[ichrom]]
			pvalues = []
			for item in emp_scores[ichrom]:
				if item < maxSkipVal:  #speed up this process by ignoring anything obviously insignificant
					pval = 1
				
				else:
					#print('scipy')
					#sys.exit()
					pval = scipy.stats.norm.sf(abs(item))
				pvalues.append(pval)
				#else:
				#	pval = get_pval(neut_rep_scores, item)
				#pvalues.append(pval)

			print("calculated pvalues for chrom " + str(chrom))
			chrom_pos = range(lastpos, lastpos + len(pvalues))

			logtenpvals = [(-1. * math.log10(pval)) for pval in pvalues]
			ax.scatter(chrom_pos, logtenpvals, color =plotcolor, s=.5)
			lastpos = chrom_pos[-1]
		else:

			chrom_pos = range(lastpos, lastpos + len(emp_scores[ichrom]))
			ax.scatter(chrom_pos, emp_scores[ichrom], color=plotcolor, s=.5)
			lastpos = chrom_pos[-1]
	return ax
def plotManhattan_extended(ax, emp_scores, chrom_pos, chrom):
	''' makes a figure more like in Karlsson 2013 instead of Grossman 2013'''
	ax.plot(chrom_pos, emp_scores, linestyle='None', marker=".", markersize=.3, color="black")
	ax.set_ylabel('chr' + str(chrom), fontsize=6, rotation='horizontal')
	labels = ax.get_yticklabels()
	ax.set_yticklabels(labels, fontsize=6)
	ax.set_axis_bgcolor('LightGray')
	return ax
