## helper functions for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 07.26.16 vitti@broadinstitute.org

import os, subprocess

def run_sel_trajs_snakemake(outputdir, cosibuild, paramfile, numSims, maxAttempts=100000):
	if outputdir [-1] != "/":
		outputdir += "/"
	writefilename = outputdir + "Snakefile"
	writefile = open(writefilename, 'w')	
	writefile.write('SAMPLES = range(1,' + str(numSims+1) + ')\n\n')
	#top-level condition
	writefile.write('rule all:\n\tinput:\n\t\texpand(\'' + outputdir + 'rep{sample}.txt\', sample=SAMPLES)\n\n')
	#run replicate
	writefile.write('rule make_traj:\n\toutput:\n\t\t\'' + outputdir + 'rep{sample}.txt\'\n')
	writefile.write('\tshell:\n\t\t\'python run_traj.py {output} ' + cosibuild + ' ' + paramfile + ' ' + str(maxAttempts)+ '\'\n')
	writefile.close()
	snakemake_command = "snakemake -s " + writefilename
	print(snakemake_command)
	#subprocess.check_output(snakemake_command.split())
	return
def run_sel_sims_snakemake(trajDir, cosibuild, paramfilename, nSimsPerBin, runDir, genmapRandomRegions=False, dropSings=None):
	#use cosi->tped direct option
	#dispatch as sel_trajs. (parallelize? relaunch if fail)?
	if outputdir [-1] != "/":
		outputdir += "/"
	#################
	## RUN ONE SIM ##
	argumentstring = "-p " + paramfilename + " --output-gen-map "
	if genmapRandomRegions:
		argumentstring += " --genmapRandomRegions"
	if dropSings is not None:
		argumentstring += " --drop-singletons " + str(dropSings)
	argumentstring += " --tped -o " + runDir + "rep"
	cmdstring = cosibuild + argumentstring
	print(cmdstring)
	#subprocess.check_output(cmdstring.split())	
	return
def write_bin_paramfile(readfilename, writefilename, bounds):
	'''given an inclusive cosi sweep parameter file, writes another with stricter final frequency rejection criteria'''
	readfile = open(readfilename, 'r')
	writefile = open(writefilename, 'w')
	for line in readfile:
		if "sweep_mult" not in line:
			writefile.write(line)
		else:
			writestring = ""
			entries = line.split()
			for ientry in range(len(entries) - 1):
				writestring += str(entries[ientry]) + " "
			writestring += str(bounds[0]) + '-'+str(bounds[1]) +'\n'
			writefile.write(writestring)
	readfile.close()
	writefile.close()
	return

###################
### SELFREQ BINS ##
###################
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

##################
### BOOKKEEPING ##
##################
def check_make_dir(dirpath):
	"""function to handle creation of (sub)folders, avoiding duplication.
	"""
	if os.path.exists(dirpath):
		contents = os.listdir(dirpath)
		if len(contents) == 0:
			return dirpath
		else:
			print dirpath + " ALREADY EXISTS; creating alt folder..."
			alt, iAlt = False, 1
			while alt == False:			
				altpath = dirpath.strip('/') + "_alt" + str(iAlt) #best way to avoid recursive issues?
				if os.path.exists(altpath):
					iAlt +=1
				else:
					alt = True
			mkdircommand = "mkdir " + altpath
			subprocess.check_output(mkdircommand.split())
			return altpath
	else:
		mkdircommand = "mkdir " + dirpath
		subprocess.check_output(mkdircommand.split())
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
