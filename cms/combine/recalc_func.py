##
## 07.23.16 	vitti@broadinstitute.org
from math import fabs
import sys
import os
def write_delIHH_file(readfilename, writefilename):
	"""taken from JV func_scores.py. given a selscan iHS file, parses it and writes delihh file"""
	readfile = open(readfilename, 'r')
	writefile = open(writefilename, 'w')
	for line in readfile:
		entries = line.split()
		locus, phys, freq_1, ihh_1, ihh_0, ihs_unnormed, ihs_normed, lastcol = entries
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
def calc_winIhs(score_array):
	"""returns mean of absolute value of all scores in window"""
	absolutes = [fabs(x) for x in score_array]
	total = sum(absolutes)
	return total/len(score_array)
