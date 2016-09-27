## functions for transforming component score calculations as part of composite.py
## last updated 09.27.16 	vitti@broadinstitute.org

from math import fabs
import sys
import os

def write_delIHH_file(readfilename, writefilename):
	"""taken from JV func_scores.py. given a selscan iHS file, parses it and writes delihh file"""
	readfile = open(readfilename, 'r')
	writefile = open(writefilename, 'w')
	for line in readfile:
		entries = line.split()
		#handle input with/without ihh decomp
		if len(entries) == 8:
			locus, phys, freq_1, ihh_1, ihh_0, ihs_unnormed, ihs_normed, lastcol = entries
		elif len(entries) == 11:
			locus, phys, freq_1, ihh_1, ihh_0, ihs, der_ihh_l, der_ihh_r, anc_ihh_l, anc_ihh_r, manually_normed = entries

					#ancestral - derived
		unstand_delIHH = fabs(float(ihh_1) - float(ihh_0))
		writeline = locus + "\t" + phys + "\t" + freq_1 + "\t" + str(ihs_unnormed) + "\t" + str(unstand_delIHH) +"\t" + str(unstand_delIHH) +  "\n" #6 columns for selscan norm
		writefile.write(writeline)
	writefile.close()
	print("wrote to " + writefilename)
	readfile.close()
def interpolate_haps(starts, ends, scores, ihh1s, ihh0s, thisPos, freqs):
	'''given input windows from winIhs, assigns a value to a position with missing ihs'''
	for istart in range(len(starts)):
		if thisPos < starts[0]: #handle first few snps
			interpolated_ihs = scores[0]
			interpolated_ihh1 = ihh1s[0]
			interpolated_ihh0 = ihh0s[0]
			interpol_freq = freqs[0]
			return interpolated_ihs, interpolated_ihh1, interpolated_ihh0, interpol_freq
		if thisPos >= starts[istart] and thisPos < ends[istart]: #in window
			interpolated_ihs = scores[istart]
			interpolated_ihh1 = ihh1s[istart]
			interpolated_ihh0 = ihh0s[istart]
			interpol_freq = freqs[istart]
			return interpolated_ihs, interpolated_ihh1, interpolated_ihh0, interpol_freq
		elif thisPos >= ends[istart]:
			if thisPos < starts[istart + 1]: #in next gap 
				interpolated_ihs = (scores[istart] + scores[istart+1])/2.
				interpolated_ihh1 = (ihh1s[istart] + ihh1s[istart+1])/2.
				interpolated_ihh0 = (ihh0s[istart] + ihh0s[istart+1])/2.
				interpol_freq = (freqs[istart] + freqs[istart+1])/2.
				return interpolated_ihs, interpolated_ihh1, interpolated_ihh0, interpol_freq
			else: #all the way at end. use final window value.
				interpolated_ihs = scores[-1]
				interpolated_ihh1 = ihh1s[-1]
				interpolated_ihh0 = ihh0s[-1]
				interpol_freq = freqs[-1]
				return interpolated_ihs, interpolated_ihh1, interpolated_ihh0, interpol_freq
		else:
			pass #advance window.
	return float('nan'), float('nan'), float('nan'), float('nan')
def calc_windowscore(score_array):
	"""returns mean of absolute value of all scores in window as part of winiHS (see Schlebusch et al 2012)"""
	absolutes = [fabs(x) for x in score_array]
	total = sum(absolutes)
	return total/len(score_array)
def windows(windowsize, jumplen, infilename, writefilename):
	'''from winiHS.py (/windelihh?)'''
	posindex, scoreindex = 1, 6
	ihh1_index, ihh0_index, freq1_index = 3,4,2

	positions, scores = [], []
	ihh1s, ihh0s, freq1s = [], [], []
	infile = open(infilename, 'r')
	for line in infile:
		entries = line.split()
		position, score = int(entries[posindex]), float(entries[scoreindex])
		positions.append(position)
		scores.append(score)

		ihh1, ihh0, freq1 = float(entries[ihh1_index]), float(entries[ihh0_index]), float(entries[freq1_index])
		ihh1s.append(ihh1)
		ihh0s.append(ihh0)
		freq1s.append(freq1) #not sure that this even counts for anything. 

	infile.close()

	###################
	## DEFINE WINDOWS #
	#this is what I interpret 
	#Schlebusch to be doing (?)

	numSnps = len(scores)
	print(str(numSnps) + " scores available")
	all_phys_windows, all_score_windows = [], []
	all_ihh1_windows, all_ihh0_windows, all_freq1_windows = [], [], []
	iSnp = 0
	while iSnp < numSnps:
		window_phys, window_scores = [], []
		window_ihh1, window_ihh0, window_freq1 = [], [], []
		while len(window_scores) < windowsize:
			if iSnp < numSnps: #?
				window_phys.append(positions[iSnp])
				window_scores.append(scores[iSnp])
				window_ihh1.append(ihh1s[iSnp])
				window_ihh0.append(ihh0s[iSnp])
				window_freq1.append(freq1s[iSnp])
			else:
				break
			iSnp +=1
		all_phys_windows.append(window_phys)
		all_score_windows.append(window_scores)
		all_ihh1_windows.append(window_ihh1)
		all_ihh0_windows.append(window_ihh0)
		all_freq1_windows.append(window_freq1)
		iSnp += jumplen
		#print(str(len(all_score_windows)))
	numWindows = len(all_score_windows)
	print("chunked " + str(numSnps) + " SNPs with available haplotype scores into " + str(numWindows) + " windows of size " + str(windowsize) + " SNPs with jump length " + str(jumplen) + ".")

	################################
	## CALC AVE AND WRITE TO FILE  #
	################################

	writefile = open(writefilename, 'w')
	for iWindow in range(numWindows):
		positions = all_phys_windows[iWindow]
		startPos, endPos = positions[0], positions[1]
		scores = all_score_windows[iWindow]
		winiHS = calc_winIhs(scores)

		ihh1_scores = all_ihh1_windows[iWindow] #since these are all positive,
		ihh0_scores = all_ihh0_windows[iWindow] #we can just reuse calc_winIhs func
		freq1_vals = all_freq1_windows[iWindow] #to get averages.

		winihh1 = calc_winIhs(ihh1_scores)
		winihh0 = calc_winIhs(ihh0_scores)
		freq_ave = calc_winIhs(freq1_vals)
		
		#writeline = str(startPos) + "\t" + str(endPos) + "\t" + str(winiHS) +'\n'
		writeline = str(startPos) + "\t" + str(endPos) + "\t" + str(winiHS) +'\t' + str(winihh1) + "\t" + str(winihh0) + "\t" + str(freq_ave) + "\n"
		writefile.write(writeline)
	writefile.close()
	print("wrote to file: " + writefilename)
def interpolate_from_windows(inputTpedFilename, inputIhsFilename, inputWinihsFilename, outputFilename):
	'''from interpolate.py'''
	all_snps = []
	openfile = open(inputTpedFilename, 'r')
	for line in openfile:
		entries = line.split()
		pos = int(entries[1])
		all_snps.append(pos)
	openfile.close()
	nsnps = len(all_snps)
	print('loaded ' + str(nsnps) + ' snps from' + inputTpedFilename)

	# load scores and windows from winIhs 
	openfile = open(inputWinihsFilename, 'r')
	starts, ends, scores = [], [], []
	ihh1s, ihh0s = [], []
	freqs = []

	for line in openfile:
		entries = line.split()
		startBp, endBp, winIhs = int(entries[0]), int(entries[1]), float(entries[2])
		
		ihh1, ihh0 = float(entries[3]), float(entries[4])
		freq_win = float(entries[5])

		starts.append(startBp)
		ends.append(endBp)
		scores.append(winIhs)
		ihh1s.append(ihh1)
		ihh0s.append(ihh0)
		freqs.append(freq_win)
	openfile.close()

	readfile = open(inputIhsFilename, 'r')
	writefile = open(outputFilename, 'w')
	readline = readfile.readline()
	entries = readline.split()
	thisScorePos = int(entries[0])

	for iSnp in range(nsnps):
		thisPos = all_snps[iSnp]
		if thisPos < thisScorePos: #interpolate until we advance to the first one for which we have a score
			interpolated_ihs, interpolated_ihh1, interpolated_ihh0, interpol_freq = interpolate_haps(starts, ends, scores, ihh1s, ihh0s, thisPos, freqs)

			unnormed_ihs= 0 #dummy to pass to selscan so norm will function
			writeline = str(thisPos) + "\t" + str(thisPos) + "\t" + str(interpol_freq) + "\t" + str(interpolated_ihh1) + "\t" + str(interpolated_ihh0) + "\t" + str(unnormed_ihs) + "\t" + str(interpolated_ihs) + "\t9\n" #index interpolated in outputfile, in case it's not already obvious?
			writefile.write(writeline)
		elif thisPos == thisScorePos: #found a match
			writefile.write(readline) #propagate to writefile
			readline = readfile.readline() #advance to next SNP for which we have a full iHS calc
			entries = readline.split()
			if len(entries) < 1:
				break
			else:
				thisScorePos = int(entries[0])
		else: #thisPos > thisScorePos
			print("ERROR: does your input iHS score file contain SNPs missing from your input TPED? Shame.")

	readfile.close()
	writefile.close()
	print('wrote to ' + outputFilename)
def calc_winIhs(score_array):
	"""returns mean of absolute value of all scores in window"""
	absolutes = [fabs(x) for x in score_array]
	total = sum(absolutes)
	return total/len(score_array)
