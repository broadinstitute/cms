## helper functions for generating visualizing component score likelihoods
## last updated: 03.28.2017 vitti@broadinstitute.org

from scipy.stats.kde import gaussian_kde
from random import choice
import numpy as np
import os

######################################
## GENERATE LIKELIHOODS FROM SCORES ##
######################################
def get_plot_pdf_params(score):
	''' defines likelihood score ranges etc. '''
	if score == "ihs":
		scorerange = [-6., 6.]
	elif score == "delihh":
		scorerange = [-3., 5.]
	elif score == "fst":
		scorerange = [0, 1.]
	elif score == "deldaf" or score =="fst_deldaf":
		scorerange = [-1., 1.]
	elif score == "xp" or score =="xpehh":
		scorerange = [-3., 8.]
	elif score=="nsl":
		scorerange = [-5., 5.]
	minVal, maxVal = scorerange
	nProbBins = 100 
	annotate = False
	return minVal, maxVal, nProbBins, annotate
def plot_pdf_comparison_from_scores(ax, neutvals, causalvals, linkedvals, minVal, maxVal, nProbBins, ax_ylabel, savebase, annotate = False):
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
	dist_space = np.linspace(minVal, maxVal, nProbBins)
	kde_neut = gaussian_kde( neutvals )
	kde_causal = gaussian_kde( causalvals )
	kde_linked = gaussian_kde( linkedvals )
	ax.set_ylabel(ax_ylabel, color=colorDict[int(ax_ylabel[0])], fontsize=8)
	ax.yaxis.set_label_position("right")
	neut_dist = kde_neut(dist_space)
	causal_dist = kde_causal(dist_space)
	linked_dist = kde_linked(dist_space)
	ax.plot( dist_space,  neut_dist, color='blue',)
	ax.plot( dist_space, causal_dist , color='red')
	ax.plot( dist_space, linked_dist , color='green')
	ax.hist(neutvals, color="blue", alpha=.5, normed=True)
	ax.hist(causalvals, color='red', alpha=.5, normed=True)
	ax.hist(linkedvals, color='green', alpha=.5, normed=True)

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
	max_bound = max(max(linked_dist), max(causal_dist), max(neut_dist))
	this_bound = .5 #pad out y-limits to accomodate maximum 
	while this_bound < max_bound:
		this_bound += .1
	ax.set_ylim([0,this_bound])
	ax.tick_params(labelsize=6)
	#ax.set_yticks(fontsize=6)
	ax.grid()

	#################
	### SAVE DIST ### (this becomes its own functionality?)
	#################
	neutwritefilename = savebase + "likes_neut.txt"
	causalwritefilename = savebase + "likes_causal.txt"
	linkedwritefilename = savebase + "likes_linked.txt"
	neutfile = open(neutwritefilename, 'w')
	causalfile = open(causalwritefilename, 'w')
	linkedfile = open(linkedwritefilename, 'w')
	#for item in dist_space:
	for ibin in range(len(dist_space)): #convert midpoints to bin bounds?
		writeprefix = str(dist_space[ibin])  + "\t" #+ str(kde_neut(dist_space)[ibin] + "\n")
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
		cl = numerator / denominator
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
		clr = causallikes[ibin] / neutlikes[ibin]
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
def calc_hist_from_scores(causal_scores, linked_scores, neut_scores, xlims, givenBins, thinToSize = False):
	if thinToSize:
		limiting = min(len(causal_scores), len(linked_scores), len(neut_scores))
		print("Thinning data to " + str(limiting) + " SNPs...")
		short_causal, short_linked, short_neut = [], [], []

		for parentlist in ['causal', 'linked', 'neut']:
			shortlist = eval('short_' + parentlist)
			takenindices = []
			while len(takenindices) < limiting:
				chosenIndex = randint(0, limiting-1)
				if chosenIndex not in takenindices:
					shortlist.append(eval(parentlist + "_scores")[chosenIndex])
					takenindices.append(chosenIndex)

		causal_scores, linked_scores, neut_scores = short_causal, short_linked, short_neut

	#print(causal_scores)

	#check for nans
	causal_scores = np.array(causal_scores)
	linked_scores = np.array(linked_scores)
	neut_scores = np.array(neut_scores)
	causal_scores = causal_scores[~np.isnan(causal_scores)]
	linked_scores = linked_scores[~np.isnan(linked_scores)]
	neut_scores = neut_scores[~np.isnan(neut_scores)]

	#get weights to plot pdf from hist
	#weights_causal = np.ones_like(causal_scores)/len(causal_scores)
	#weights_linked = np.ones_like(linked_scores)/len(linked_scores)
	#weights_neut = np.ones_like(neut_scores)/len(neut_scores)

	causal_scores = np.clip(causal_scores, xlims[0], xlims[1]) #np.clip
	linked_scores = np.clip(linked_scores, xlims[0], xlims[1])
	neut_scores = np.clip(neut_scores, xlims[0], xlims[1])

	n_causal, bins_causal = np.histogram(causal_scores, range=xlims, bins=givenBins)#, weights = weights_causal)
	n_linked, bins_linked = np.histogram(linked_scores, range=xlims,  bins=givenBins)#, weights = weights_linked)
	n_neut, bins_neut = np.histogram(neut_scores,range=xlims, bins=givenBins)#, weights = weights_neut)

	#debug_array = [n_causal, n_linked, n_neut, bins_causal, bins_linked, bins_neut]
	#print(debug_array)

	#totalNsnps = len(causal_scores) + len(linked_scores) + len(neut_scores)

	#for ibin in range(len(n_causal)):
	#	if n_causal[ibin] <= 1:
	#		n_causal[ibin] = 1e-10
	#	if n_linked[ibin] <= 1:
	#		n_linked[ibin] = 1e-10
	#	if n_neut[ibin] <=1:
	#		n_neut[ibin] = 1e-10

	return n_causal, n_linked, n_neut, bins_causal, bins_linked, bins_neut
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