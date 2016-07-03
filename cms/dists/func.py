## helper functions for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 07.03.16 vitti@broadinstitute.org

import subprocess
import sys
import os

def give_command(commandstr):
	commands = commandstr.split(' ')
	subprocess.check_output(commands)
def check_make_dir(dirpath):
	"""function to handle creation of (sub)folders, avoiding duplication.
	"""
	if os.path.exists(dirpath):
		contents = os.listdir(dirpath)
		if len(contents) == 0:
			return dirpath
		else:
			alt, iAlt = False, 1
			while alt == False:			
				altpath = dirpath.strip('/') + "_alt" + str(iAlt) #best way to avoid recursive issues?
				if os.path.exists(altpath):
					iAlt +=1
				else:
					alt = True
			mkdircommand = "mkdir " + altpath
			give_command(mkdircommand)
			return altpath
	else:
		mkdircommand = "mkdir " + dirpath
		give_command(mkdircommand)
	return dirpath
def check_file_len(filename):
	'''counts number of lines in file'''
	if os.path.isfile(filename) and os.path.getsize(filename) > 0:
		openfile = open(filename, 'r')
		iline = 0
		for line in openfile:
			iline +=1
		openfile.close()
		return iline
	else:
		return 0
def get_bin_strings(bin_medians):
	'''given input bin medians, returns minimal len strs needed to create unique bin labels'''
	stringlen = 3 #minimum is '0.x' -- not necessarily the most descriptive way to round these, but it's functional.
	unique = False
	while unique == False:
		labels = [str(x)[:stringlen] for x in bin_medians]
		if len(labels) == len(set(labels)):
			unique = True
		else:
			stringlen +=1
	return labels
def get_bins(freqRange=".05-.95", numbins=9):
	fullrange = [float (x) for x in freqRange.split('-')]
	binlen = (fullrange[1] - fullrange[0]) / float(numbins)
	bin_starts = [fullrange[0] + binlen*i for i in range(numbins)]
	bin_ends = [fullrange[0] + binlen*i for i in range(1,numbins+1)]
	bin_medians = [(float(bin_starts[i]) + float(bin_ends[i]))/2. for i in range(len(bin_starts))]
	bin_medians_str = get_bin_strings(bin_medians)
	return fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str	
def check_bin_filled(directory, numsims):
	'''counts all non-empty files in directory, returns True if >= numsims'''
	filenames = os.listdir(directory)
	nonempty = [x for x in filenames if check_file_len(x) !=0]
	if len(filenames) < numsims:
		return False
	else:
		return True