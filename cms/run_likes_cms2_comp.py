##script for making likelihood tables for population comparison-based statistics
## /idi/sabeti-scratch/jvitti 11.4.16: just pass .cms file

runPrefix = "likes_111516"
scorelistdir = "/idi/sabeti-scratch/jvitti/" + runPrefix + "/"
writedir = "/idi/sabeti-scratch/jvitti/" + runPrefix + "_b/"
discardOut, overwriteList = False, False
import os
import subprocess

def check_create_dir(directory):
	if not os.path.isdir(directory):
		subprocess.check_output(['mkdir', directory])
def writeJobsToTaskArray_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-P sabeti_lab', '-q long'], memory="5g", discardOut=False):
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

mem = "10g"

scores = ['deldaf', 'fst',  'xpehh']
cmd = "python /idi/sabeti-scratch/jvitti/cms/cms/likes_from_model.py likes_from_scores"
models = ['default_default_101715_12pm', 'nulldefault_constantsize', 'nulldefault', 'default_112115_825am', 'gradient_101915_treebase_6_best',]
pops = [1, 2, 3, 4] 
numPerBin, selPos, numLikesBins = 500, 500000, 60 #for histograms
selbin_clusters = [["0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90"]] #["0.10", "0.20", "0.30"], ["0.40", "0.50", "0.60",], ["0.70", "0.80", "0.90",]
cluster_names = ["allfreq"]#"low", "mid", "hi", 

def main():
	arguments = []
	check_create_dir(writedir)
	for model in models:
		model_likes_dir = scorelistdir + model + "/"
		model_write_dir = writedir + model + "/"
		check_create_dir(model_likes_dir)	
		check_create_dir(model_write_dir)
		for score in scores:
			model_score_likes_dir = model_likes_dir + score + "/"
			check_create_dir(model_score_likes_dir)
			write_score_dir = model_write_dir + score + "/"
			check_create_dir(write_score_dir)
			for pop in pops:
				
				neutfilename = getNeutCompFiles(model, pop, model_score_likes_dir, runPrefix, overwrite=overwriteList)

				for icluster in range(len(cluster_names)):
					selbins = selbin_clusters[icluster]
					binclustername = cluster_names[icluster]

					selfilename = getSelCompFiles(model, pop, model_score_likes_dir, runPrefix, selbins, binclustername, overwrite=overwriteList)
				
					outprefix = write_score_dir  + "likes_sel" + str(pop) + "_choose_" + binclustername

					#CHECK IF COMPLETED
					#if True:
					if not os.path.isfile(outprefix + "_causal.txt") or not os.path.isfile(outprefix + "_neut.txt"):
						argstring = "--" + score + " " + neutfilename + " " + selfilename + " " + str(selPos) +  " " + outprefix + " --nLikesBins " + str(numLikesBins)
						arguments.append(argstring)

	writeJobsToTaskArray_fromArglist(cmd, arguments, runPrefix + "_comp_pops", discardOut=discardOut, memory=mem)
main()



