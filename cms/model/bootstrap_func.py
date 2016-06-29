## {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS}
## this script contains functions called by script: get_neutral_targetstats_from_bootstrap.py
## last updated: 06.14.16 	vitti@broadinstitute.org

import numpy as np
import random
import os

#######################
## VERSATILE / PARSE ##
#######################

def flattenList(list3d):
	if type(list3d[0]) == list:
		flat = []
		for list1d in list3d:
			for item in list1d:
				flat.append(item)
		return flat
	else:
		return list3d
def checkFileExists(filename):
	if os.path.isfile(filename) and os.path.getsize(filename) > 0:
		return True
	else:
		return False
def readFreqsFile(filename):
	"""returns parallel lists of lists, separated by region"""
	#print "reading from: " + filename
	nregions, nsnps, totallen = 0, 0, 0
	allpi, allnderiv, allnanc = [], [], []
	seqlens = []
	openfile = open(filename, 'r')
	pi, nderiv, nanc = [], [], [] # per region
	for line in openfile:
		if line != "\n":
			entries = line.split()
			if int(entries[0]) <= 22 and len(entries) == 5: #chrom
				thispi, thisnderiv, thisnanc = float(entries[2]), int(entries[3]), int(entries[4])
				pi.append(thispi)
				nderiv.append(thisnderiv)
				nanc.append(thisnanc)
			else:
				seqlen = int(entries[0])
				
		else:
			assert len(pi) == len(nderiv) and len(nanc) == len(pi) 
			#only take non-empty regions
			if pi:
				allpi.append(pi)
				allnderiv.append(nderiv)
				allnanc.append(nanc)
				nregions +=1
				nsnps += len(pi)
				totallen += seqlen
				seqlens.append(seqlen)
			pi, nderiv, nanc = [], [], []
	openfile.close()
	assert len(seqlens) == nregions
	#print "\tlogged frequency values for " + str(nsnps) + " SNPS across " + str(nregions) + " regions covering " + str(totallen) + " bp."
	return allpi, allnderiv, allnanc, nregions, seqlens
def readLDFile(filename, dprimecutoff = .2):
	"""returns parallel lists of lists, separated by region
	D' filters on MAF > .2 for both loci. (r2 previously filtered sings/doubs)"""
	#print "reading from: " + filename
	nregions, npairs = 0, 0
	alldists, allr2 = [], []
	allgendists, alldprime = [], []
	openfile = open(filename, 'r')
	dist, gendist, r2, dprime = [], [], [], []# per region
	r2_snppairs, dprime_snppairs = 0, 0
	r2_reg, dprime_reg = 0, 0
	for line in openfile:
		if line != "\n":
			entries = line.split()
			if len(entries) == 8: #quick fix for truncated entries at end of file.
				start, end = int(entries[1]), int(entries[2])
				thisdist = end - start
				thisgendist = float(entries[3])
				freq1, freq2 = float(entries[4]), float(entries[5])
				thisr2, thisdprime = float(entries[6]), float(entries[7])

				dist.append(thisdist)
				r2.append(thisr2)

				if freq1 > dprimecutoff and freq2 > dprimecutoff:
					gendist.append(thisgendist)
					dprime.append(thisdprime)
		else:
			assert len(r2) == len(dist) and len(dprime) == len(gendist)
			#only take non-empty regions
			if r2:
				alldists.append(dist)
				allr2.append(r2)
				r2_reg +=1
				r2_snppairs += len(r2)
			if dprime:
				allgendists.append(gendist)
				alldprime.append(dprime)
				dprime_reg +=1
				dprime_snppairs += len(dprime)

			dist, gendist, r2, dprime = [], [], [], []

	alldists.append(dist)
	allr2.append(r2)
	allgendists.append(gendist)
	alldprime.append(dprime)
	openfile.close()

	#print "\tlogged r2 values for " + str(r2_snppairs) + " SNP pairs across " + str(r2_reg) + " regions."
	#print "\tlogged dprime values for " + str(dprime_snppairs) + " SNP pairs across " + str(dprime_reg) + " regions, using MAF cutoff of " + str(dprimecutoff) +"."
	return alldists, allr2, allgendists, alldprime, r2_reg, dprime_reg
def readFstFile(filename):
	"""returns parallel lists of lists, separated by regions"""
	#print "reading from: " + filename
	nregions, nsnps = 0, 0
	allfst = []
	fst = [] #per region
	openfile = open(filename, 'r')
	for line in openfile:
		if line != "\n":
			entries = line.split()
			pos, thisfst = int(entries[1]), float(entries[2])
			fst.append(thisfst)
		else:
			if fst:
				allfst.append(fst)
				nregions +=1
				nsnps += len(fst)
			fst = []
	openfile.close()
	#print "\tlogged fst values for " + str(nsnps) + " snps across " + str(nregions) + " regions."
	return allfst, nregions

####################################
## CALCULATE BOOTSTRAP ESTIMATES ###
####################################

def estimateFstByBootstrap(allRegionValues, nrep = 100):
	numRegions = len(allRegionValues)
	allValues = flattenList(allRegionValues)
	target_mean = np.mean(allValues)
	estimates = []
	for i in range(nrep):
		num_bs_regions = 0
		this_rep_values = []
		while num_bs_regions < numRegions:
			thisregion = random.choice(allRegionValues)
			this_rep_values.extend(thisregion)
			num_bs_regions +=1
		estimates.append(np.mean(this_rep_values))
	target_se = np.std(estimates)
	return target_mean, target_se
def estimateFstByBootstrap_bysnp(allRegionValues, nrep = 100):
	numRegions = len(allRegionValues)
	allValues = flattenList(allRegionValues)
	target_mean = np.mean(allValues)
	estimates = []
	for i in range(nrep):
		this_rep_values = []
		for w in range(len(allValues)):
			thissnp = random.choice(allValues)
			this_rep_values.append(thissnp)
		estimates.append(np.mean(this_rep_values))
	target_se = np.std(estimates)
	return target_mean, target_se
def estimateFreqSpectrum(allnderiv, allnanc, nhist = 6):
	"""translated from calc_pop_stats_1kg_regions.c"""
	starts = [-1, 1e-9, .1, .2, .3, .4]
	ends = [-1, .1, .2, .3, .4, .5]
	mafhist = [0, 0, 0, 0, 0, 0]
	anchist = [0, 0, 0, 0, 0, 0]
	if type(allnderiv[0]) == list:
		allnderiv = flattenList(allnderiv)
	if type(allnanc[0]) == list:
		allnanc = flattenList(allnanc)
	for i in range(len(allnderiv)):
		nderiv = allnderiv[i]
		nanc = allnanc[i]
		nminor = min(nanc, nderiv)
		nmajor = max(nanc, nderiv)
		freq = float(nminor) / nmajor
		#print str(nminor) + "\t" + str(nmajor)
		if nminor > 0:
			if nderiv == 1 or nanc == 1:
				mafhist[0] += 1
			else:
				for ibin in range(1, nhist):
					if freq >= starts[ibin] and freq < ends[ibin]:
						mafhist[ibin] += 1
			if nanc < nderiv: #if minor allele is ancestral, update anchist
				if nanc == 1:
					anchist[0] += 1
				else:
					for ibin in range(1, nhist):
						if freq >= starts[ibin] and freq < ends[ibin]:
							anchist[ibin] += 1		
	#print anchist
	return mafhist, anchist
def estimatePi(allpivals, seqlens):
	pi_sum, length_sum = 0, 0
	assert len(allpivals) == len(seqlens)
	for i in range(len(allpivals)):
		length_sum += seqlens[i]
		for j in range(len(allpivals[i])):
			pi_sum += allpivals[i][j]
	return (pi_sum/length_sum)
def estimater2decay(allr2, allphysdist):
	nDistHist = 14
	starts = [x * 5000 for x in range(nDistHist)] #0, 5000...65000
	ends = [x*5000 for x in range(1,nDistHist+1)] #5000...70000
	r2sums = [0 for x in range(nDistHist)]
	physDistHist = [0 for x in range(nDistHist)]
	allr2 = flattenList(allr2)
	allphysdist = flattenList(allphysdist)
	for i in range(len(allr2)):
		dist = allphysdist[i]
		for ibin in range(0, nDistHist):
			if dist >= starts[ibin] and dist < ends[ibin]:
				#print str(dist) + "\t" + str(ibin) + "\t" + str(allr2[i]) + "\n"
				r2sums[ibin] += allr2[i]
				physDistHist[ibin] += 1
	return r2sums, physDistHist
def estimatedprimedecay(alldprime, allgendists):
	assert len(alldprime) == len(allgendists)
	nGenDistHist = 17
	starts = [x * .001 for x in range(nGenDistHist)]
	ends = [x * .001 for x in range(1, nGenDistHist+1)]
	compLDhist = [0 for x in range(nGenDistHist)]
	genDistHist = [0 for x in range(nGenDistHist)]
	alldprime = flattenList(alldprime)
	allgendists = flattenList(allgendists)
	for i in range(len(alldprime)):
		gendist = allgendists[i]
		dprime = alldprime[i]
		#print "DEBUG: gendist: " + str(gendist) + " dprime: " + str(dprime)
		for ibin in range(0, nGenDistHist):
			if gendist >= starts[ibin] and gendist < ends[ibin]:
				genDistHist[ibin] += 1
				if dprime == 1.0:
					compLDhist[ibin] += 1
	return compLDhist, genDistHist