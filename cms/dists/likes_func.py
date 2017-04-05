## helper functions for generating visualizing component score likelihoods
## last updated: 04.02.2017 vitti@broadinstitute.org

from scipy.stats.kde import gaussian_kde
from random import choice
import numpy as np
import os

######################################
## GENERATE LIKELIHOODS FROM SCORES ##
######################################
def get_plot_pdf_params(score):
	''' defines likelihood score ranges etc. '''
	"""For within-region localization, the value ranges used to define the bins were as follows: 
	iHS, [-3,3]; XP-EHH, [-3,3]; iHH difference, [-3,3]; FST, [-2,2]; and DDAF, [-1,1]; 
	for CMSGW, the value ranges were as follows: 
	iHS, [-6,6]; XP-EHH, [-3,8]; iHH difference, [-3,5]; FST, [-1,6]; and DDAF, [-1,1]. 
	Values outside the range were binned into the nearest bin at the end of the range.
	"""

	if score == "ihs":
		scorerange = [0,6.] #folding for region-detection: the ancestral allele at this SNP may tag a nearby causal derived variant, so polarity isn't useful to us #[-6., 6.]
	elif score == "delihh":
		scorerange = [-3., 5.]
	elif score == "fst":
		scorerange = [0, 1.]
	elif score == "deldaf" or score =="fst_deldaf":
		scorerange = [-1., 1.]
	elif score == "xp" or score =="xpehh":
		scorerange = [-3., 8.]
	elif score=="nsl":
		scorerange = [0,5.]#[-5., 5.]
	minVal, maxVal = scorerange
	nProbBins = 50 
	annotate = False
	return minVal, maxVal, nProbBins, annotate
def plot_pdf_comparison_from_scores(ax, neutvals, causalvals, linkedvals, minVal, maxVal, nProbBins, ax_ylabel, savebase, annotate = False, saveFiles=False, heuristic=False):
	''' inputs arrays of sorted scores and outputs estimated probability density function for component score (as image and datafile)'''
	colorDict = {1:'#FFB933', 2:'#0EBFF0', 3:'#ADCD00', 4:'#8B08B0'}
	################### 
	### CURATE DATA ###
	###################
	neutvals = np.array(neutvals)
	neutvals = neutvals[~np.isnan(neutvals)]
	neutvals = np.clip(neutvals, minVal, maxVal)
	causalvals = np.array(causalvals)
	causalvals = causalvals[~np.isnan(causalvals)]
	causalvals = np.clip(causalvals, minVal, maxVal)
	linkedvals = np.array(linkedvals)
	linkedvals = linkedvals[~np.isnan(linkedvals)]
	linkedvals = np.clip(linkedvals, minVal, maxVal)

	##################### 
	### CALCULATE PDF ###
	#####################
	if heuristic: #plot both smoothed function and coarse-grained histogram
		dist_space = np.linspace(minVal, maxVal, nProbBins)
		kde_neut = gaussian_kde( neutvals )
		kde_causal = gaussian_kde( causalvals )
		kde_linked = gaussian_kde( linkedvals )
		neut_dist = kde_neut(dist_space)
		causal_dist = kde_causal(dist_space)
		linked_dist = kde_linked(dist_space)
		ax.plot( dist_space,  neut_dist, color='blue',)
		ax.plot( dist_space, causal_dist , color='red')
		ax.plot( dist_space, linked_dist , color='green')
		ax.hist(neutvals, color="blue", alpha=.5, normed=True)
		ax.hist(causalvals, color='red', alpha=.5, normed=True)
		ax.hist(linkedvals, color='green', alpha=.5, normed=True)

	else: #plot fine-grained histogram
		#get weights to plot pdf from hist
		weights_neut = np.ones_like(neutvals)/len(neutvals) 
		weights_causal = np.ones_like(causalvals)/len(causalvals)
		weights_linked = np.ones_like(linkedvals)/len(linkedvals)

		score_range = float(maxVal) - float(minVal)
		interval = score_range / float(nProbBins)
		givenBins =[float(minVal) + ibin * interval for ibin in range(nProbBins+1)] #starts and ends
		# minVal maxVal, nPro
		starts = givenBins[:-1]
		ends = givenBins[1:]
		#for i_bin in range(len(starts)):
		#	print(str(starts[i_bin]) + "\t" + str(ends[i_bin]))
		#print(str(givenBins))

		xlims = [minVal, maxVal]
		n_causal, bins_causal = np.histogram(causalvals, range=xlims, bins=givenBins, weights = weights_causal)
		n_linked, bins_linked = np.histogram(linkedvals, range=xlims,  bins=givenBins, weights = weights_linked)
		n_neut, bins_neut = np.histogram(neutvals,range=xlims, bins=givenBins, weights = weights_neut)

		#midpoints = [float(minVal) + (interval/2) + (ibin * interval) for ibin in range(nProbBins)]
		midpoints = []
		for i in range(nProbBins):
			midpoint = (givenBins[i] + givenBins[i+1])/2.
			midpoints.append(midpoint)
		#print(midpoints)
		bar_width = interval
		#print(str(bar_width))
		#print(str(n_neut))
		ax.bar(midpoints, n_neut, bar_width, alpha=.5, color="blue")
		ax.bar(midpoints, n_linked, bar_width, alpha=.5, color="green")
		ax.bar(midpoints, n_causal, bar_width, alpha=.5, color="red")

		#midpoints = [(starts[i] + ends[i])/2 for i in range(len(starts))]
		#bar_width = (starts[0] - ends[0])/5
		#if color == "blue":
		#	offset = 0
		#elif color == "red":
		#	offset = bar_width
		#elif color == "green":
		#	offset = (-1 * bar_width)
		#plot_xvals = [midpoints[i] + offset for i in range(len(midpoints))]
		#ax.scatter(midpoints, vals, color=color)
		#ax.bar(plot_xvals, vals, bar_width, color=color, edgecolor=color)
		#ax.set_xlim(xlims)
		#ax.set_ylim(ylims)
		linked_dist, causal_dist, neut_dist = n_linked, n_causal, n_neut #unsmoothed
		
	ax.set_ylabel(ax_ylabel, color=colorDict[int(ax_ylabel[0])], fontsize=8)
	ax.yaxis.set_label_position("right")
	max_bound = max(max(linked_dist), max(causal_dist), max(neut_dist))
	this_bound = .1#.5 #pad out y-limits to accomodate maximum 
	while this_bound < max_bound:
		this_bound += .01
	ax.set_ylim([0,this_bound])
	ax.tick_params(labelsize=6)
	#ax.set_yticks(fontsize=6)
	ax.grid()
	if annotate:
		mean = np.mean(neutvalues)
		var = np.var(neutvalues)
		sd = np.sqrt(var)
		annotationstring = "NEUT max: " + str(max(values)) + "\nmin: " + str(min(values)) + "\nmean: " + str(np.mean(values)) + "\nvar: " + str(np.var(values))
		xlims = ax.get_xlim()
		ylims = ax.get_ylim()
		fullxrange = xlims[1] - xlims[0]
		fullyrange = ylims[1] - ylims[0]
		anno_x = xlims[1] - (fullxrange/3)
		anno_y = ylims[1] - (fullyrange/2)
		ax.annotate(annotationstring, xy=(anno_x, anno_y), fontsize=6)

	####
	#### Should build in a check here to enforce
	#### monotonicity/smoothness, bounds etc.

	#################
	### SAVE DIST ### (this becomes its own functionality?)
	#################
	if saveFiles:
		neutwritefilename = savebase + "likes_neut.txt"
		causalwritefilename = savebase + "likes_causal.txt"
		linkedwritefilename = savebase + "likes_linked.txt"
		neutfile = open(neutwritefilename, 'w')
		causalfile = open(causalwritefilename, 'w')
		linkedfile = open(linkedwritefilename, 'w')
		#for item in dist_space:
		for ibin in range(len(starts)): #convert midpoints to bin bounds?
			writeprefix = str(starts[ibin])  + "\t" + str(ends[ibin]) + "\t" #+ str(kde_neut(dist_space)[ibin] + "\n")
			#writefile.write(str(writestring))
			neutline = writeprefix + str(neut_dist[ibin]) + "\n"
			causalline = writeprefix + str(causal_dist[ibin]) + "\n"
			linkedline = writeprefix + str(linked_dist[ibin]) + "\n"				
			neutfile.write(neutline)
			causalfile.write(causalline)
			linkedfile.write(linkedline)
		neutfile.close()
		causalfile.close()
		linkedfile.close()
		print('wrote to e.g.: ' + neutwritefilename)
	return

##########################
## VALIDATE LIKELIHOODS ##
##########################
def quick_cl_from_likes(causallikes, linkedlikes, nSnp = 1000, takeLn = True):
	assert len(causallikes) == len(linkedlikes)
	prior = 1/nSnp
	composite_likes = []
	for ibin in range(len(causallikes)):
		numerator = causallikes[ibin] * prior
		denominator = (causallikes[ibin] * prior) + (linkedlikes[ibin] * (1-prior))
		if denominator != 0:
			cl = numerator / denominator
		else:
			cl = float('nan')
		#if cl < 1e-5:
		#	pass
		if takeLn:
			cl = np.log(cl)
		composite_likes.append(cl)
	return composite_likes
def quick_clr_from_likes(causallikes, neutlikes, takeLn = True):
	assert len(causallikes) == len(neutlikes)
	composite_like_ratios = []
	for ibin in range(len(causallikes)):
		if neutlikes[ibin] != 0:
			clr = causallikes[ibin] / neutlikes[ibin]
		else:
			clr = float('nan')
		#if clr < 1e-5:
		#	pass
		if takeLn:
			clr = np.log(clr)
		composite_like_ratios.append(clr)
	return composite_like_ratios
def quick_load_likes(infilename, prob_index = -1):
	likes = []
	#print(infilename)
	openfile = open(infilename, 'r')
	for line in openfile:
		entries = line.split()
		prob = float(entries[prob_index])
		likes.append(prob)
	openfile.close()
	return likes



"""
def write_hists_to_files(writePrefix, givenBins, n_causal, n_linked, n_neut):
	assert len(givenBins) == (len(n_causal) + 1)
	for status in ['causal', 'linked', 'neut']:
		writefilename = writePrefix + "_" + status + ".txt"
		writefile = open(writefilename, 'w')
		n_scores = eval('n_' + status)
		for index in range(len(n_scores)):
			num_in_bin = n_scores[index]
			if (num_in_bin <= 1):
				writeprob = 1e-10
			else:
				if (index < (len(n_scores) - 1) and index > 0): #check neighbors
					if (n_scores[index+1] == n_scores[index-1]) and (n_scores[index+1] <=1):
						writeprob = 1e-10
					else:
						writeprob = float(n_scores[index])/(sum(n_scores)) #
				else:
					writeprob = float(n_scores[index])/(sum(n_scores)) #
			towritestring =  str(givenBins[index]) + "\t" + str(givenBins[index+1]) + "\t" + str(writeprob)+ "\n"
			writefile.write(towritestring)
		writefile.close()
	return
def read_likes_file(likesfilename):
	'''parses a file from e.g. write_hists_to_files'''
	starts, ends, vals = [], [], []
	openfile = open(likesfilename, 'r')
	for line in openfile:
		entries = line.split()
		start, end, val = [float(x) for x in entries]
		for item in ['start', 'end', 'val']:
			eval(item + "s").append(eval(item))
	openfile.close()
	return starts, ends, vals
"""