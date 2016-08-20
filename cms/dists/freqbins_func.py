## helper functions for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 08.16.16 vitti@broadinstitute.org

import os, subprocess

def write_slurm_array(writefilename = "slurm_array.sh", arrayCmd = "./my_process $SLURM_ARRAY_TASK_ID", numTasks=1, dispatch=False):
	writefile = open(writefilename, 'w')
	writefile.write('#!/bin/bash\n')
	writefile.write('#SBATCH --array=1-' + str(numTasks) + "\n")
	writefile.write(arrayCmd)
	writefile.write('\n')
	writefile.close()
	subcommand = "sbatch --array=1-" + str(numTasks) + " " + writefilename
	print(subcommand)
	if dispatch:
		subprocess.check_output(subcommand.split())
	#print('wrote to: ' + writefilename)
	return
def run_seltraj_arrays(trajdir, cosibuild, paramfile, numSims, scriptDir = "/home/users/vitti/cms/cms/", taskIndexStr = "$SLURM_ARRAY_TASK_ID", maxAttempts = 1000):
	if trajdir [-1] != "/":
		trajdir += "/"
	traj_outputname = trajdir + 'rep' + taskIndexStr + '.txt'
	traj_cmd = "python " + scriptDir + "run_traj.py " + traj_outputname + " " + cosibuild + " " + paramfile + " " + str(maxAttempts) +'\n'#
	write_slurm_array(trajdir + "run_sel_traj.sh", traj_cmd, numSims, dispatch=True)
	return
def run_selsim_arrays(trajdir, cosibuild, paramfile, numSims, runDir, scriptDir = "/home/users/vitti/cms/cms/", taskIndexStr = "$SLURM_ARRAY_TASK_ID", maxAttempts = 1000, genmapRandomRegions=False, dropSings=None):
	'''I think I may be over-complicating things with Snakemake.'''
	#traj_outputname = trajdir + 'rep' + taskIndexStr + '.txt'
	if runDir [-1] != "/":
		runDir += "/"
	sim_cmd = "env COSI_NEWSIM=1 env COSI_LOAD_TRAJ=" + trajdir + "/rep" + taskIndexStr + ".txt " + cosibuild + " -p " + paramfile + " --output-gen-map " 
	if genmapRandomRegions:
		sim_cmd += " --genmapRandomRegions"
	if dropSings is not None:
		sim_cmd += " --drop-singletons " + str(dropSings)
	sim_cmd += " --tped " + runDir + 'rep' + taskIndexStr
	write_slurm_array(runDir + "run_sel_sims.sh", sim_cmd, numSims, dispatch=True)
	return
def run_sel_trajs_snakemake(outputdir, cosibuild, paramfile, numSims, cluster="", jobs=1, maxAttempts=10000, waitSec = 1000, scriptDir = "/home/users/vitti/cms/cms/"):
	'''this happens in each folder to be populated'''
	
	#####################
	## WRITE SNAKEFILE ##
	#####################
	if outputdir [-1] != "/":
		outputdir += "/"
	writefilename = outputdir + "Snakefile"
	writefile = open(writefilename, 'w')	
	writefile.write('SAMPLES = range(1,' + str(numSims+1) + ')\n\n')
	#top-level condition
	writefile.write('rule all:\n\tinput:\n\t\texpand(\'' + outputdir + 'rep{sample}.txt\', sample=SAMPLES)\n\n')
	#run replicate
	writefile.write('rule make_traj:\n\toutput:\n\t\t\'' + outputdir + 'rep{sample}.txt\'\n')
	#make_traj_cmd = '\tshell:\n\t\t\'env COSI_NEWSIM=1 COSI_MAXATTEMPTS=' + str(maxAttempts) + " COSI_SAVE_TRAJ={output} " + cosibuild + " -p " + paramfile + "'\n"
	#writefile.write(make_traj_cmd)
	writefile.write('\tshell:\n\t\t\'python ' + scriptDir + 'run_traj.py {output} ' + cosibuild + ' ' + paramfile + ' ' + str(maxAttempts)+ '\'\n')
	writefile.close()
	
	#############################
	## WRITE SNAKEMAKE COMMAND ##
	#############################
	snakemake_command1 = "snakemake -s " + writefilename + " --directory " + outputdir + " --unlock\n"	
	subprocess.check_output(snakemake_command1.split())

	snakemake_command2 = "snakemake -s " + writefilename 
	if cluster is not None:
		snakemake_command2 += " --cluster " + str(cluster) + " -j " + str(jobs)
	snakemake_command2 += " --output-wait " + str(waitSec)  #DEBUGGING
	#cores_string = "-j 999"
	#+ " " + cores_string + " --output-wait " + str(waitSec) + " --rerun-incomplete "  --cluster qsub"
	snakemake_command2 += " --rerun-incomplete"
	print(snakemake_command2)
	subprocess.check_output(snakemake_command2.split())
	return
def run_sel_sims_snakemake(trajDir, cosibuild, paramfilename, nSimsPerBin, runDir, genmapRandomRegions=False, dropSings=None, waitSec = 5000):
	cores_string = "-j 999"

	#use cosi->tped direct option
	#dispatch as sel_trajs. (parallelize? relaunch if fail)?

	#####################
	## WRITE SNAKEFILE ##
	#####################
	if runDir [-1] != "/":
		runDir += "/"
	writefilename = runDir + "Snakefile"
	writefile = open(writefilename, 'w')	
	writefile.write('SAMPLES = range(1,' + str(nSimsPerBin+1) + ')\n\n')
	writefile.write('POPS = range(1,4)\n\n') #FOR DEBUG
	#top-level condition
	writefile.write('rule all:\n\tinput:\n\t\texpand(\'' + runDir + 'rep{sample}.pos-{pop}.tped\', sample=SAMPLES, pop=POPS)\n\n')
	#run replicate
	writefile.write('rule run_rep:\n\tinput:\n\t\t\'' + trajDir + '/rep{sample}.txt\'\n\toutput:\n\t\t\'' + runDir + 'rep{sample}\'\n')
	argumentstring = " -p " + paramfilename + " --output-gen-map"
	if genmapRandomRegions:
		argumentstring += " --genmapRandomRegions"
	if dropSings is not None:
		argumentstring += " --drop-singletons " + str(dropSings)
	argumentstring += " --tped {output}" #+ runDir + "rep{sample}"
	commandstring = "env COSI_LOAD_TRAJ={input} " + cosibuild
	writefile.write('\tshell:\n\t\t\'' + commandstring + argumentstring+"'\n")
	#writefile.write('\tshell:\n\t\t\'python ' + scriptDir + 'run_traj.py {output} ' + cosibuild + ' ' + paramfile + ' ' + str(maxAttempts)+ '\'\n')
	writefile.close()


	#############################
	## WRITE SNAKEMAKE COMMAND ##
	#############################
	snakemake_command1 = "snakemake -s " + writefilename + " --unlock\n"
	subprocess.check_output(snakemake_command1.split())
	snakemake_command2 = "snakemake -s " + writefilename + " " + cores_string + " --output-wait " + str(waitSec)# + " --unlock"#+ " --jobs " + str(numSims) + " --cluster qsub"
	print(snakemake_command2)

	#cmdstring = cosibuild + argumentstring
	#print(cmdstring)
	#subprocess.check_output(cmdstring.split())	
	return
def write_bin_paramfile(readfilename, writefilename, bounds):
	'''given an inclusive cosi sweep parameter file, writes another with stricter final frequency rejection criteria'''
	readfile = open(readfilename, 'r')
	writefile = open(writefilename, 'w')
	foundSweep = False
	for line in readfile:
		if "sweep_mult" not in line:
			writefile.write(line)
			foundSweep = True
		else:
			writestring = ""
			entries = line.split()
			for ientry in range(len(entries) - 1):
				writestring += str(entries[ientry]) + " "
			writestring += str(bounds[0]) + '-'+str(bounds[1]) +'\n'
			writefile.write(writestring)
	readfile.close()
	writefile.close()
	if foundSweep == False:
		print("ERROR: did not find sweep parameters in inputfile.")
	return

###################
### SELFREQ BINS ##
###################
def get_bin_strings(bin_medians):
	'''given input bin medians, returns minimal len strs needed to create unique bin labels'''
	stringlen = 3 #minimum is '0.x' -- not necessarily the most descriptive way to round these, but it's functional.
	unique = False
	bin_medians_rewritten = ["%.2f" % a for a in bin_medians]
	#print('foobar!')
	#print(bin_medians_rewritten)
	#while unique == False:
	#	labels = ["%.2f" % a for a in bin_medians]
	#	[str(x)[:stringlen] for x in bin_medians]
	#	if len(labels) == len(set(labels)):
	#		unique = True
	#	else:
	#		stringlen +=1
	#return labels
	return bin_medians_rewritten
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
	** COME BACK TO THIS JV 	"""
	if os.path.exists(dirpath):
		return dirpath
		#contents = os.listdir(dirpath)
		#if len(contents) == 0:
		#	return dirpath
		#else:
		#	print(dirpath + " ALREADY EXISTS; creating alt folder...")
		#	alt, iAlt = False, 1
		#	altpath = dirpath + "alt" + str(iAlt) #best way to avoid recursive issues?
		#	while alt == False:			
		#		if os.path.exists(altpath):
		#			iAlt +=1
		#			altpath = dirpath + "alt" + str(iAlt) #best way to avoid recursive issues?
		#		else:
		#			alt = True
		#	mkdircommand = "mkdir " + altpath
		#	subprocess.check_output(mkdircommand.split())
		#	return altpath
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
