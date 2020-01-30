##script for making likelihood tables for per-population statistics
## last updated /idi/sabeti-scratch/jvitti 11.15.16

runPrefix = "likes_111516"
scorelistdir = "/idi/sabeti-scratch/jvitti/" + runPrefix + "/"
writedir = "/idi/sabeti-scratch/jvitti/" + runPrefix + "_b/"
discardOut, overwriteList = True,False

import subprocess
import os

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
def check_create_dir(directory):
	if not os.path.isdir(directory):
		subprocess.check_output(['mkdir', directory])


selbin_clusters = [["0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90"]] #["0.10", "0.20", "0.30"], ["0.40", "0.50", "0.60",], ["0.70", "0.80", "0.90",],
cluster_names = [ "allfreq"]#"low", "mid", "hi",]
mem = "3g"

scores = ['ihs','nsl', 'delihh',]
numPerBin, selPos, numLikesBins = 500, 500000, 60 #for histograms
cmd = "python /idi/sabeti-scratch/jvitti/cms/cms/likes_from_model.py likes_from_scores"


models = ['nulldefault_constantsize', 'nulldefault','default_112115_825am','default_default_101715_12pm', 'gradient_101915_treebase_6_best',]
pops = [1, 2, 3, 4] 

def main():
	arguments = []
	check_create_dir(writedir)
	for model in models:
		model_write_dir = writedir + model + "/"
		check_create_dir(model_write_dir)
		for score in scores:
			

			model_write_likes_dir = model_write_dir + score + "/"
			check_create_dir(model_write_likes_dir)
			for pop in pops:
				neutlistfilename = scorelistdir + model +"/" + score + "/repfiles_neut_" + str(pop) + ".list"	
				neutfilename = getNeutFiles_asList(model, score, pop, neutlistfilename, overwrite=overwriteList, numPerBin=500)
				for icluster in range(len(cluster_names)):
					selbins = selbin_clusters[icluster]
					binclustername = cluster_names[icluster]
					listfilename = scorelistdir +model +"/" + score + "/repfiles_sel" + str(pop) + "_" + binclustername + ".list"	
					selfilename = getSelFiles_asList(model, score, pop, listfilename, selbins, binclustername, overwrite=False)

					outprefix = model_write_likes_dir + "likes_sel" + str(pop) + "_" + binclustername 		

					outfilenames = [outprefix +"_causal.txt", outprefix + "_neut.txt"]
					if not (os.path.isfile(outfilenames[0]) and os.path.isfile(outfilenames[1])):
					#if True:
						argstring = "--" + score + " " + neutfilename + " " + selfilename + " " + str(selPos) +  " " + outprefix + " --nLikesBins " + str(numLikesBins)
						arguments.append(argstring)
			
	writeJobsToTaskArray_fromArglist(cmd, arguments, runPrefix, discardOut=discardOut, memory=mem)
main()



