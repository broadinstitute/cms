## /idi/sabeti-scratch/jvitti 10.25.16 11.4 choose comp
import subprocess
import os
import sys
import random

###########################
## HANDLE PER-REP FILES ##
###########################
def check_file_concordance(in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file, ncheck = 30):
	'''in-place to make sure pulled score files match'''

	for checkfile in [in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file]:
		if not os.path.isfile(checkfile):
			print("missing: " +checkfile )
			return False

	ihs_snp_nos, ihs_positions, ihs_freqs_1 = read_rep_info(in_ihs_file)
	delihh_snp_nos, delihh_positions, delihh_freqs_1 = read_rep_info(in_delihh_file)
	xp_snp_nos, xp_positions, xp_freqs_1 = read_rep_info(in_xp_file, stripHeader=True)	
	freqs_positions = read_fst_deldaf_sites(in_fst_deldaf_file)

	#pick central snp
	len_ihs_snp = len(ihs_positions)
	icheck = 0
	while icheck<ncheck:
		check_index = int(len_ihs_snp/2)
		barrier = int(len_ihs_snp/4)
		check_index += random.randint(-1*barrier, barrier) #random near-central snp
		if check_index < len(ihs_positions)-1:
			check_snp_no, check_position, check_freq_1 = ihs_snp_nos[check_index], ihs_positions[check_index], ihs_freqs_1[check_index]
			if check_position in freqs_positions and check_position in delihh_positions and check_position in xp_positions:
				delihh_index = delihh_positions.index(check_position)
				xp_index = xp_positions.index(check_position)
				if (check_snp_no == delihh_snp_nos[delihh_index]) and (check_snp_no == xp_snp_nos[xp_index]) and (check_freq_1 == delihh_freqs_1[delihh_index]) and (check_freq_1 == xp_freqs_1[xp_index]):
					return True
				else:
					icheck +=1
		print('no dice on '  + in_ihs_file + " after " + str(ncheck) + " checks")
		return False

def get_score_files(model, repNum, selPop, altPop, scenario = "neut", filebase = "/idi/sabeti-scratch/jvitti/scores/"):
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

def getSelFiles_asList_poolalt(model, score, pop, listfilename, selbins, binclustername, overwrite = True, numPerBin = 500):
	if not overwrite:
		if os.path.isfile(listfilename):
			print(listfilename + " exists.")
			return listfilename
	listfile = open(listfilename, 'w')
	nFile = 0
	if score in ['fst', 'deldaf']:
		score = 'fst_deldaf'
	altpops = [1, 2, 3, 4]
	altpops.remove(pop)
	for selbin in selbins:
		selfilebase = '/idi/sabeti-scratch/jvitti/scores/' + model + "/sel" + str(pop) + '/' + score + '/sel_' + str(selbin) + '/'
		for irep in range(1, numPerBin+1):
			for altpop in altpops:
				if score == "xpehh":
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

def getNeutFiles_asList_poolalt(model, score, pop, listfilename, altpop = "", overwrite = True, numPerBin = 500):
	if not overwrite:
		if os.path.isfile(listfilename):
			print(listfilename + " exists.")
			return listfilename
	listfile = open(listfilename, 'w')
	nFile = 0
	if score in ['fst', 'deldaf']:
		score = 'fst_deldaf'
	neutfilebase = '/idi/sabeti-scratch/jvitti/scores/' + model + "/neut/" + score + '/'
	altpops = [1, 2, 3, 4]
	altpops.remove(pop)
	for irep in range(1, numPerBin+1):
		for altpop in altpops:
			if score == "xpehh":
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

def getNeutFile_concat(model, score, pop, altpop = "", numPerBin = 500):
	"""I'm not sure I even want to use this anymore."""
	if score == 'ihs' or score == 'nsl':
		neutfilename = "/idi/sabeti-scratch/jvitti/scores/" + model + "/neut/" + score + "/n500_" + str(pop) + "_concat." + score + ".100bins.norm"
	elif score == 'delihh':
		neutfilename = "/idi/sabeti-scratch/jvitti/scores/" + model + "/neut/" + score + '/delihh_n500_' + str(pop) + ".100bins.norm"
	elif score == "xpehh":
		neutfilename = "/idi/sabeti-scratch/jvitti/scores/" + model + "/neut/xpehh/n500_" + str(pop) + "_" + str(altpop) + "_concat.xpehh.norm"
	elif score == "fst" or score == 'deldaf':
		neutfilename = "/idi/sabeti-scratch/jvitti/scores/" + model + "/neut/fst_deldaf/fst_deldaf_"+ str(pop) + "_" + str(altpop) + ".concat"#n500_" + str(pop) + "_" + str(altpop) + ".concat"
		if os.path.isfile(neutfilename):
			return neutfilename
		else:
			print(neutfilename)
			writefile = open(neutfilename, 'w')
			for irep in range(1, numPerBin+1):
				repfilename = "/idi/sabeti-scratch/jvitti/scores/" + model + "/neut/fst_deldaf/rep" + str(irep) +"_" + str(pop) + "_" + str(altpop)
				repfile = open(repfilename, 'r')
				for line in repfile:
					writefile.write(line)
				repfile.close()
			writefile.close()
			print("wrote to: " + neutfilename)
	if not os.path.isfile(neutfilename):
		print("missing: " + neutfilename)
	return neutfilename

###########
## ETC ##
##########

def write_master_likesfile(writefilename, model, selpop, freq,basedir,  miss = "neut", alt = "pool"):
	writefile = open(writefilename, 'w')
	for score in ['ihs', 'delihh']: 
		hitlikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_" + str(freq) + "_causal.txt"
		misslikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_" + str(freq) + "_" + miss + ".txt"
		writefile.write(hitlikesfilename + "\n" + misslikesfilename + "\n")
	for score in ['xpehh', 'fst', 'deldaf']:
		if alt != "pool":
			print("must config")
			sys.exit(0)
		else:
			hitlikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_choose_" + str(freq) + "_causal.txt"
			misslikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_choose_" + str(freq) + "_" + miss + ".txt"
			writefile.write(hitlikesfilename + "\n" + misslikesfilename + "\n")
	writefile.close()
	print("wrote to: " + writefilename)
	return
def read_rep_info(filename, stripHeader = False):
	''' pulls SNP info from scorefile to ensure concordance of simulated replicate '''
	openfile = open(filename, 'r')
	snp_nos, positions, freqs_1 = [], [], [] 
	if stripHeader:
		openfile.readline() 
	for line in openfile:
		entries = line.split()
		this_snp, this_pos, this_1freq = int(entries[0]), int(entries[1]), float(entries[2])
		#catch XP columns
		if stripHeader:
			this_1freq = float(entries[3])
		snp_nos.append(this_snp)
		positions.append(this_pos)
		freqs_1.append(this_1freq)
	openfile.close()
	return snp_nos, positions, freqs_1 
def read_fst_deldaf_sites(filename):
	openfile = open(filename, 'r')
	openfile.readline() #strip header
	positions = []
	for line in openfile:
		entries = line.split()
		if entries[0] != "chrom":
			thisPos = int(entries[0])
			positions.append(thisPos)
	openfile.close()
	return positions
def check_create_dir(directory):
	if not os.path.isdir(directory):
		subprocess.check_output(['mkdir', directory])
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

