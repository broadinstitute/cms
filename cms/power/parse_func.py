## last updated 1.6.17

import subprocess
import numpy as np
import os
import matplotlib as mp 
mp.use('agg')
import matplotlib.pyplot as plt

def plot_dist(allvals, savefilename= "/web/personal/vitti/test.png", numBins=10000):
	if len(allvals) > 0:
		f, ax = plt.subplots(1)
		ax.hist(allvals, bins=10000)
		plt.savefig(savefilename)
		print('plotted to ' + savefilename)
	return

##################
## LOCATE FILES ##
##################
def get_concat_files(model, pop, score, altpop = '', basedir = '/idi/sabeti-scratch/jvitti/clean/scores/'):
	if score in ['ihs', 'delihh', 'nsl']:
		concatfilebase = basedir + model + "/neut/concat_" + str(pop) + "_"
	elif score in ['xpehh', 'fst']:
		concatfilebase = basedir + model + "/neut/concat_" + str(pop) + "_" + str(altpop) + "_"
	else:
		concatfilebase = ""
	concatfilename = concatfilebase + score + ".txt"
	binfilename = concatfilebase + score + ".bins"
	return concatfilename, binfilename
def get_component_score_files(model, irep, pop, altpop, selbin = "neut", filebase = "/idi/sabeti-scratch/jvitti/clean/scores/", normed=False):
	if selbin == 'neut':
		basedir =  filebase + model + "/neut/"
	else:
		basedir = filebase + model + "/sel" + str(pop) + "/sel_" + str(selbin) + "/"
	in_ihs_file = basedir + "ihs/rep" + str(irep) + "_" + str(pop) + ".ihs.out"
	in_delihh_file =  basedir + "delihh/rep" + str(irep) + "_" + str(pop) + ".txt"
	in_nsl_file = basedir + "nsl/rep" + str(irep) + "_" + str(pop) + ".nsl.out"
	in_xp_file = basedir + "xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out"
	in_fst_deldaf_file = basedir + "fst_deldaf/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop)

	if normed:
		in_ihs_file += ".norm" 
		in_nsl_file += ".norm"
		in_delihh_file += ".norm" 
		in_xp_file += ".norm"

	for returnfile in [in_ihs_file, in_nsl_file, in_delihh_file, in_xp_file, in_fst_deldaf_file]:
		if not os.path.isfile(returnfile):
			print("Missing: " + returnfile)

	return in_ihs_file, in_nsl_file, in_delihh_file, in_xp_file, in_fst_deldaf_file 
def get_neut_repfile_name(model, irep, pop, normed = False, basedir = "/idi/sabeti-scratch/jvitti/clean/scores/", suffix = ""):
	repfilename = basedir +"composite/"+ model + "/neut/rep" + str(irep) + "_" + str(pop) + ".cms" + suffix
	#print(repfilename)
	if normed:
		repfilename += ".norm"
	#if not os.path.isfile(repfilename):
	#	return
	#else:
	return repfilename
def get_sel_repfile_name(model, irep, pop, freqbin, normed = False,  basedir = "/idi/sabeti-scratch/jvitti/scores/", suffix=""):
	repfilename = basedir +"composite/"+ model + "/sel" + str(pop) + "/sel_" + str(freqbin) + "/rep" + str(irep) + "_" + str(pop) + ".cms" + suffix
	#print(repfilename)
	if normed:
		repfilename += ".norm"
	#if not os.path.isfile(repfilename):
	#	return	
	#else:
	return repfilename
def get_pr_filesnames(key, basedir):
	regionlen, percentage, cutoff, model, pop = key
	fprfile = basedir + "/fpr/" + model + "/sel" + str(pop) + "/fpr_" + str(regionlen) + "_" + str(percentage) + "_" + str(cutoff)
	tprfile = basedir + "/tpr/" + model + "/sel" + str(pop) + "/tpr_" + str(regionlen) + "_" + str(percentage) + "_" + str(cutoff)
	for filename in [tprfile, fprfile]:
		if not os.path.isfile(filename):
			print("missing: " + filename)
	return fprfile, tprfile 
def get_emp_cms_file(selpop, model, likessuffix, chrom, normed = False, basedir = "/idi/sabeti-scratch/jvitti/scores_composite4_b/"):#"/idi/sabeti-scratch/jvitti/synth/cms_composite/"):
	filename = basedir + selpop  +"_" + str(model) + "_" + likessuffix + ".chr" + str(chrom) + ".txt"
	if normed:
		filename += ".norm"
	if not os.path.isfile(filename):
		print("MISSING empirical file : " + filename)
	return filename
def get_likesfiles(model, selpop, likessuffix, likesdir, allfreqs = True):
	if not allfreqs:
		hi_likesfile = likesdir+ model + "/master/likes_" + str(selpop) + "_hi_vs_" + likessuffix + ".txt"
		mid_likesfile = likesdir + model + "/master/likes_" + str(selpop) + "_mid_vs_" + likessuffix + ".txt"
		low_likesfile = likesdir + model + "/master/likes_" + str(selpop) + "_low_vs_" + likessuffix + ".txt"
		return hi_likesfile, mid_likesfile, low_likesfile 
	else:
		likesfile = likesdir + model + "/master/likes_" + str(selpop) + "_allfreq_vs_" + likessuffix + ".txt_2"
		return likesfile
	return
def get_neutinscorefiles(model, selpop, irep, pairbasedir = "/idi/sabeti-scratch/jvitti/clean/scores/"):
	inscorefiles = []
	pops = [1, 2, 3, 4]
	altpops = pops[:]
	altpops.remove(int(selpop))
	for altpop in altpops:
		pairfilename = pairbasedir + model + "/neut/pairs/rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop) + ".pair"
		#print(pairfilename)
		if os.path.isfile(pairfilename):
			inscorefiles.append(pairfilename)
	return inscorefiles	
def get_selinscorefiles(model, selpop, sel_freq_bin, irep, pairbasedir = "/n/regal/sabeti_lab/jvitti/clear/"):
	inscorefiles = []
	pops = [1, 2, 3, 4]
	altpops = pops[:]
	altpops.remove(int(selpop))
	for altpop in altpops:
		pairfilename = pairbasedir + "scores/" + model + "/sel" + str(selpop) + "/sel_" + str(sel_freq_bin) + "/pairs/rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop) + ".pair"
		#print(pairfilename)
		if os.path.isfile(pairfilename):
			inscorefiles.append(pairfilename)
	return inscorefiles	


##################
## FILE PARSING ##
##################
def read_cms_repfile(infilename):
	physpos, genpos, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = [], [], [], [], [], [], [], [], [], []
	if not os.path.isfile(infilename):
		return physpos, genpos, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed
	openfile = open(infilename, 'r')
	for line in openfile:
		entries = line.split()
		physpos.append(int(entries[0]))
		genpos.append(float(entries[1]))
		ihs_normed.append(float(entries[2]))
		delihh_normed.append(float(entries[3]))
		nsl_normed.append(float(entries[4]))
		xpehh_normed.append(float(entries[5]))
		fst.append(float(entries[6]))
		deldaf.append(float(entries[7]))
		cms_unnormed.append(float(entries[8]))
		if ".norm" in infilename:
			cms_normed.append(float(entries[9]))
	openfile.close()
	return physpos, genpos, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed
def load_simscores(model, pop, vsNeut = False, numRep = 500, normed =False):
	all_simscores = []
	for irep in range(1, numRep + 1):
		neutfilename = get_neut_repfile_name(model, irep, pop, vsNeut = vsNeut, normed=normed)
		if os.path.isfile(neutfilename):
			#print(neutfilename)
			openfile = open(neutfilename, 'r')
			for line in openfile:
				entries = line.split()
				val = float(entries[-1])
				all_simscores.append(val)
			openfile.close()
	infostring = 'loaded a total of ' + str(len(all_simscores)) 
	if normed:
		infostring += " normed"
	infostring += " scores."
	print(infostring)
	return all_simscores
def load_empscores(model, selpop, likessuffix, normed=False):
	chroms = range(1,23)
	scores = []
	for chrom in chroms:
		scorefile = get_emp_cms_file(selpop, model, likessuffix, chrom, normed=normed)
		print('loading from ' + scorefile)
		assert os.path.isfile(scorefile)
		openfile = open(scorefile, 'r')
		for line in openfile:
			entries=line.split()
			scores.append(float(entries[-1]))
		openfile.close()
	return scores
def read_pr(infilename):
	openfile = open(infilename, 'r')
	firstline = openfile.readline()
	pr = float(firstline)
	openfile.close()
	return pr
def readvals_lastcol(filename):
	allvals = []
	openfile = open(filename, 'r')
	for line in openfile:
		entries = line.split()
		value=float(entries[-1])#float(entries[4]) #thisihs, thisihh, thisxpehh, thisfst, thisdelDaf,
		if not np.isnan(value):
			if not np.isnan(np.log(value)):
				allvals.append(np.log(value))
	openfile.close()
	return allvals

######################
## UGER/BOOKKEEPING ##
######################
def check_create_file(filename, checkOverwrite):
	if checkOverwrite == False:
		return True #make it anyway
	else:
		if not os.path.isfile(filename) or os.path.getsize(filename) == 0:
			return True
		else:
			return False
	return False
def check_create_dir(directory):
	if not os.path.isdir(directory):
		subprocess.check_output(['mkdir', '-p', directory])
def writeJobsToTaskArray_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-P sabeti_lab'], memory="5g", discardOut=False):
	"""flexible function for generating scripts and argfiles to dispatch parallel jobs on Univa Grid Engine.
	arglist: a list of strings. each string represents one argstring, for one task on the array."""
	#write file with arguments for all tasks
	t = 0
	taskargfilename = filewd + fileprefix + "_taskargs.txt"
	taskargfile = open(taskargfilename, "w")
	for argstring in arguments:
		taskargfile.write(argstring + "\n")
		t += 1
	taskargfile.close() 
	#write script for one job that points to task-argument-file
	scriptname = filewd + fileprefix + '.sh'
	scriptfile = open(scriptname, 'w')
	for UGERoption in UGERoptions:
		scriptfile.write("#$ " + UGERoption + "\n")
	scriptfile.write("#$ -t 1-" + str(t) + "\n")
	scriptfile.write("\n")

	scriptfile.write('export PATH=/idi/sabeti-scratch/jvitti/cms/cms:$PATH\n')
	scriptfile.write('export PATH=/home/unix/vitti/miniconda3/bin:$PATH\n')
	scriptfile.write('source activate /home/unix/vitti/miniconda3/envs/cms-env3\n')
	#scriptfile.write("source activate cms-env3\n")
	scriptfile.write('source /broad/software/scripts/useuse\n')
	#scriptfile.write('use .scipy-0.15.1-python-3.4.3\n')

	scriptfile.write("MYFILE=" + str(taskargfilename) + "\n")
	scriptfile.write("SAMPLE=$(awk \"NR==$SGE_TASK_ID\" $MYFILE)\n")
	writestring = commandstring + " $SAMPLE\n"
	scriptfile.write(writestring)
	scriptfile.close()
	if discardOut:
		subcommand = "qsub -e /tmp/test_error -o /tmp/test_out -l h_vmem=" + str(memory) + " " + scriptname
	else:
		subcommand = "qsub -l h_vmem=" + str(memory) + " " + scriptname
	print(subcommand)
	return
def execute(commandstring):
	#print(commandstring)
	subprocess.check_output(commandstring.split())
	return
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


"""
def get_component_score_files_old(model, repNum, selPop, altPop, scenario = "neut", filebase = "/idi/sabeti-scratch/jvitti/clean/scores/"):
	'''PHASE OUT'''
	if scenario == "neut":
		in_ihs_file = filebase + model + "/neut/ihs/rep" + str(repNum) + "_" + str(selPop) + "_normed.txt"
		in_delihh_file = filebase + model + "/neut/delihh/rep" + str(repNum) +  "_" + str(selPop) + ".txt.norm"
		in_xp_file = filebase + model + "/neut/xpehh/rep" + str(repNum) +  "_" + str(selPop) + "_" + str(altPop) + "_normed.txt"
		in_fst_deldaf_file = filebase + model + "/neut/fst_deldaf/rep" + str(repNum) + "_" + str(selPop) + "_" + str(altPop)

	else:
		in_ihs_file = filebase + model + "/sel" + str(selPop) + "/ihs/sel_" + str(scenario) + "/rep" + str(repNum) + "_" + str(selPop) + "_truncOk.ihs.out.norm"
		in_delihh_file = filebase + model + "/sel" + str(selPop) + '/delihh/sel_' + str(scenario) + '/rep' + str(repNum) + "_" + str(selPop) + ".txt.norm"
		in_xp_file = filebase + model + "/sel" + str(selPop) + '/xpehh/sel_' + str(scenario) +'/rep' + str(repNum) + "_" + str(selPop) + "_" + str(altPop) + ".xpehh.out.norm"
		in_fst_deldaf_file = filebase + model + "/sel" + str(selPop) + '/fst_deldaf/sel_' + str(scenario) +'/rep' + str(repNum) + "_" + str(selPop) + "_" + str(altPop)

	for returnfile in [in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file]:
		if not os.path.isfile(returnfile):
			print("Missing: " + returnfile)

	return in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file 
def getSelFiles_asList(model, score, pop, listfilename, selbins, binclustername, altpop = "", overwrite = True, numPerBin = 500):
	if not overwrite:
		if os.path.isfile(listfilename):
			print(listfilename + " exists.")
			return listfilename
	listfile = open(listfilename, 'w')
	nFile = 0
	if score in ['ihs', 'delihh', 'nsl']:
		assert(altpop == "")
	if score in ['xpehh', 'deldaf', 'fst']:
		assert(altpop != "")
	if score in ['fst', 'deldaf']:
		score = 'fst_deldaf'

	for selbin in selbins:
		selfilebase = '/idi/sabeti-scratch/jvitti/scores/' + model + "/sel" + str(pop) + '/' + score + '/sel_' + str(selbin) + '/'
		for irep in range(1, numPerBin+1):
			if score == "ihs":
				repfilename = selfilebase + "rep" + str(irep) + "_" + str(pop) + "_truncOk.ihs.out.norm"
			elif score == "delihh":
				repfilename = selfilebase + "rep" + str(irep) + "_" + str(pop) + ".txt.norm"
			elif score == "nsl":
				repfilename = selfilebase + "rep" + str(irep) + "_" + str(pop) + "." + score + ".out.norm"
			elif score == "xpehh":
				repfilename = selfilebase + "rep"+ str(irep) + "_" + str(pop) +"_" + str(altpop) + ".xpehh.out.norm"
			elif score == "fst_deldaf":
				repfilename = selfilebase + "rep"+ str(irep) + "_" + str(pop) +"_" + str(altpop)
			if not os.path.isfile(repfilename):
				print("missing: " + repfilename)
			if os.path.isfile(repfilename):	 
				listfile.write(repfilename + "\n")
				nFile +=1
	listfile.close()
	print('wrote ' + str(nFile) + ' files to ' + listfilename)
	return listfilename
def getNeutFiles_asList(model, score, pop, listfilename, altpop = "", overwrite = True, numPerBin = 500):
	if not overwrite:
		if os.path.isfile(listfilename):
			print(listfilename + " exists.")
			return listfilename
	listfile = open(listfilename, 'w')
	nFile = 0
	if score in ['ihs', 'delihh', 'nsl']:
		assert(altpop == "")
	if score in ['xpehh', 'deldaf', 'fst']:
		assert(altpop != "")
	if score in ['fst', 'deldaf']:
		score = 'fst_deldaf'
	neutfilebase = '/idi/sabeti-scratch/jvitti/scores/' + model + "/neut/" + score + '/'
	for irep in range(1, numPerBin+1):
		if score == "ihs":
			repfilename = neutfilebase + "rep" + str(irep) + "_" + str(pop) + "_normed.txt"
		elif score == "delihh":
			repfilename = neutfilebase + "rep" + str(irep) + "_" + str(pop) + ".txt.norm"
		elif score == "nsl":
			repfilename = neutfilebase + "rep" + str(irep) + "_" + str(pop) + "_normed.txt"#"." + score +# ".out.norm"
		elif score == "xpehh":
			repfilename = neutfilebase + "rep"+ str(irep) + "_" + str(pop) +"_" + str(altpop) + "_normed.txt"
		elif score == "fst_deldaf":
			repfilename = neutfilebase + "rep"+ str(irep) + "_" + str(pop) +"_" + str(altpop)
		if not os.path.isfile(repfilename):
			print("missing: " + repfilename)
		if os.path.isfile(repfilename):	 
			listfile.write(repfilename + "\n")
			nFile +=1
	listfile.close()
	print('wrote ' + str(nFile) + ' files to ' + listfilename)
	return listfilename
def getNeutCompFiles(model, pop,  writedir, writeprefix, overwrite=False, scoredir = "/idi/sabeti-scratch/jvitti/scores_composite2/"):
	if writedir[-1] != "/":
		writedir += "/"
	writefilename = writedir + writeprefix + "_sel" + str(pop) +"_comp_neut.list"
	if os.path.isfile(writefilename) and overwrite == False:
		print("already exists: " + writefilename)
		return writefilename
	else:
		modeldir = scoredir + model + "/neut/"
		#print(modeldir)
		modelcontents = os.listdir(modeldir)
		selidstring = "_" + str(pop) +".cms.out"
		filecontents = [item for item in modelcontents if selidstring in item and ".norm" not in item and ".decompose" not in item]
		print('found ' + str(len(filecontents)) + ' files for pop...')
		writefile = open(writefilename, 'w')
		for filename in filecontents:
			writefile.write(modeldir + filename + "\n")
		writefile.close()
		print('wrote to: ' + writefilename)
	return writefilename
def getSelCompFiles(model, pop,  writedir, writeprefix, selbins, binclustername, overwrite=False, scoredir = "/idi/sabeti-scratch/jvitti/scores_composite2/"):
	if writedir[-1] != "/":
		writedir += "/"
	writefilename = writedir + writeprefix + "_sel" + str(pop) +"_comp_" + binclustername + ".list"
	if os.path.isfile(writefilename) and overwrite == False:
		print("already exists: " + writefilename)
		return writefilename
	else:
		allcontents = []
		modeldir = scoredir + model + "/sel" + str(pop) + "/"
		for selbin in selbins:
			bindir = modeldir + "sel_" + str(selbin) + '/'
			#print(bindir)
			bincontents = os.listdir(bindir)
			selidstring = "_" + str(pop) +"_vsneut.cms.out"
			#print(selidstring)
			filecontents = [bindir + item for item in bincontents if selidstring in item and ".norm" not in item and ".decompose" not in item]
			allcontents.extend(filecontents)
		print('found ' + str(len(allcontents)) + ' files for pop...')
		writefile = open(writefilename, 'w')
		for filename in allcontents:
			writefile.write(filename + "\n")
		writefile.close()
		print('wrote to: ' + writefilename)
	return writefilename
"""

