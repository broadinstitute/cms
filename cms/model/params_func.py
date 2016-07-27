## !!-----this is where the user specifies the INPUT PARAMETERS for the demographic model to be fit-------!!
## By way of example, currently configured to fit African populations from 1000 Genomes. (YRI, LWK, MSL, ESN)
## User should duplicate/adapt this file to reflect the populations being modeled. 
## See http://broad-cms.readthedocs.io/en/latest/workflow.html 		last updated: 07.05.16		vitti@broadinstitute.org

def read_lines(openfile, numlines):
	for i in range(numlines):
		line = openfile.readline()
	return line

##########################
## DEFINE TARGET VALUES ##
##########################
def get_target_values(bootstrap_targetval_filename):
	#print "MUST CONFIG FOR AWS: currently just sfs"
	popDict = {'LWK': 1, 'GWD': 2, 'MSL': 3, 'YRI':4, 'ESN':5}
	stats = {}
	openfile = open(bootstrap_targetval_filename, 'r')
	for ipop in range(1,5):

		popline = read_lines(openfile, 1)
		popline = popline.strip('\n')
		ipop = popDict[popline]
		
		piline = read_lines(openfile, 1)
		entries = piline.split()
		pi_mean, pi_se = float(entries[0]), float(entries[1])
		stats[('pi', ipop)] = [pi_mean]
		stats[('pi_var', ipop)] = [pi_se ** 2]

		sfs_mean_line = read_lines(openfile, 1)
		entries = sfs_mean_line.split()
		sfs_mean = eval(sfs_mean_line)
		stats[('sfs', ipop)] = sfs_mean
		sfs_se_line = read_lines(openfile, 1)
		sfs_se = eval(sfs_se_line)
		sfs_var = [float(x)**2 for x in sfs_se]
		stats[('sfs_var', ipop)] = sfs_var

		anc_mean_line = read_lines(openfile, 1)
		anc_mean = eval(anc_mean_line)
		stats[('anc', ipop)] = anc_mean
		anc_se_line = read_lines(openfile, 1)
		anc_se = eval(anc_se_line)
		anc_var = [float(x)**2 for x in anc_se]
		stats[('anc_var', ipop)] = anc_var


	for ipop in range(1, 5):
		r2_mean_line = read_lines(openfile, 1)
		r2_mean = eval(r2_mean_line)
		stats[('r2', ipop)] = r2_mean
		r2_se_line = read_lines(openfile, 1)
		r2_se = eval(r2_se_line)
		r2_var = [float(x)**2 for x in r2_se]
		stats[('r2_var', ipop)] = r2_var

		dprime_mean_line = read_lines(openfile, 1)
		dprime_mean = eval(dprime_mean_line)
		stats[('dprime', ipop)] = dprime_mean
		dprime_se_line = read_lines(openfile, 1)
		dprime_se = eval(dprime_se_line)
		dprime_var = [float(x)**2 for x in dprime_se]
		stats[('dprime_var', ipop)] = dprime_var

	popPairs = []
	for ipop in range(1, 5):
		for jpop in range(ipop, 5):
			if ipop != jpop: 
				popPairs.append((ipop, jpop))

	for popPair in popPairs:
		fstline = read_lines(openfile, 1)
		entries = fstline.split()
		fst, fst_se = float(entries[2]), float(entries[3])
		stats[("fst", popPair)] = [fst]
		stats[("fst_var", popPair)] =[fst_se ** 2]
	
	return stats

##########################
## DEFINE TREE TOPOLOGY ##
##########################
def generate_params():
	############################
	##DEFINE GLOBAL PARAMETERS##
	############################

	paramDict = {'chromlength':100000, 'mutation_rate':1.25e-08, 'num_indivs_per_sample':170, 'gene_conv_rel_rate':2.3, 'singrate':.25,
	'recomfilename':'results_inclusive_test_filter.recom'}

	paramDict['numPops'] = 4
	paramDict['labels'] = {1:'LWK', 2:'GWD', 3:'MSL', 4:'YRI'}	

	###################
	##DEFINE BRANCHES##
	###################

	paramDict[('split', 2)] = [175]
	paramDict[('split', 3)] = [105]
	paramDict[('split', 4)] = [75]
									
									#Indexing: At THIS TIME, arrive at THIS SIZE via THIS KIND of change (forward not coalescent perspective)
									#0			#1			#2			#3					#4					#5				#6
									#present 	#midexp		#pre exp 	minimum				#transitional		#anc			#base
	paramDict[('Ne', 1)] = 			[75000000,	40000,		8000,		18400,				15666,				30700,			15000]
	paramDict[('Ne_times', 1)] = 	[1,			60,			390,		400,				2200,				3300,			15000]
	paramDict[('Ne_change', 1)] = 	[1,			1,			1,			1,					0,					0,					1] #code = {1: 'exp', 0:'normal'}

	paramDict[('Ne', 2)] = 			[300000,	40000,		8000]
	paramDict[('Ne_times', 2)] =	[1, 		50,			paramDict[('split', 2)][0]-1]
	paramDict[('Ne_change', 2)] = 	[1, 		1,			1]

	paramDict[('Ne', 3)] = 			[300000,	40000,		8000]
	paramDict[('Ne_times', 3)] = 	[1, 		75, 		paramDict[('split', 3)][0]-1]
	paramDict[('Ne_change', 3)] = 	[1, 		1,			1]

	paramDict[('Ne', 4)] = 			[300000,	40000,		8000]
	paramDict[('Ne_times', 4)] = 	[1,			50, 		paramDict[('split', 4)][0]-1]
	paramDict[('Ne_change', 4)] = 	[1, 		1, 			1]

	##############
	##POP EVENTS##
	##############
 
	paramDict[('bn', 1)] = [1e-10, 1e-10, 0.001]
	paramDict[('bn_times', 1)] = [50, 3000, 5000]
	paramDict[('bn', 2)] = [1e-10, 1e-10,]
	paramDict[('bn_times', 2)] = [paramDict[('split', 2)][0] - 1, 50]
	paramDict[('bn', 3)] = [1e-10, 1e-10,]
	paramDict[('bn_times', 3)] = [paramDict[('split', 3)][0] - 1, 50]
	paramDict[('bn', 4)] = [1e-10, 1e-10,]
	paramDict[('bn_times', 4)] = [paramDict[('split', 4)][0] - 1, 50] 

	#############
	##GENE FLOW##
	#############

	paramDict[('mig', '1->2')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '1->2')] = [0, 10, 100]
	paramDict[('mig', '2->1')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '2->1')] = [0, 10, 100]

	paramDict[('mig', '1->3')] =  [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '1->3')] = [0, 10, 100]
	paramDict[('mig', '3->1')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '3->1')] = [0, 10, 100]

	paramDict[('mig', '1->4')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '1->4')] = [0, 10, 100]
	paramDict[('mig', '4->1')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '4->1')] = [0, 10, 100]

	paramDict[('mig', '2->3')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '2->3')] = [0, 10, 100]
	paramDict[('mig', '3->2')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '3->2')] = [0, 10, 100]

	paramDict[('mig', '2->4')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '2->4')] = [0, 10, 100]
	paramDict[('mig', '4->2')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '4->2')] = [0, 10, 100]

	paramDict[('mig', '3->4')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '3->4')] = [0, 10, 100]
	paramDict[('mig', '4->3')] = [0, 0, 0]#[1e-10, 1e-10, 1e-10]
	paramDict[('mig_time', '4->3')] = [0, 10, 100]

	return paramDict
def get_ranges():
	paramDict = generateParams()
	split2time = paramDict[('split', 2)][0]
	split3time = paramDict[('split', 3)][0]
	split4time = paramDict[('split', 4)][0]

	rangeDict = {}
	rangeDict[('split', 2)] = [[175, 500]]
	rangeDict[('split', 3)] = [[100, 174]]
	rangeDict[('split', 4)] = [[10, 99]]
	rangeDict[('Ne', 1)] = [[50000, 10000000], [5000, 500000],  [1000, 100000], [1000, 100000], [1000, 1000000], [1000, 1000000], [100, 100000]]
	rangeDict[('Ne_times', 1)] = [[1, 5], [6, 199], [200, 1000], [500, 8000], [2000, 12000], [2500, 25000], [12000, 50000]]
	rangeDict[('Ne', 2)] = [[50000, 10000000], [5000, 500000],  [1000, 100000]]
	rangeDict[('Ne_times', 2)] = [[1, 5], [6, 200], [200, split2time]]
	rangeDict[('Ne', 3)] = [[50000, 10000000], [5000, 500000],  [1000, 100000]]
	rangeDict[('Ne_times', 3)] = [[1, 5], [6, split3time], [10, split3time]]
	rangeDict[('Ne', 4)] = [[50000, 10000000], [5000, 500000],  [1000, 100000]]
	rangeDict[('Ne_times', 4)] = [[1, 5], [6, split4time], [10, split4time]]
	rangeDict[('bn', 1)] = [[1e-6, .5], [1e-6, .5],  [1e-6, .5], [1e-6, .5]]
	rangeDict[('bn_times', 1)] = [[0,50000],[0,50000],[0,50000],[0,50000]]
	rangeDict[('bn', 2)] = [[1e-6, .5], [1e-6, .5], [1e-6, .5], [1e-6, .5]]
	rangeDict[('bn_times', 2)] = [[0,split2time],[0,split2time],[0,split2time],[0,split2time]]
	rangeDict[('bn', 3)] = [[1e-6, .5], [1e-6, .5], [1e-6, .5], [1e-6, .5]]
	rangeDict[('bn_times', 3)] = [[0,split3time],[0,split3time],[0,split3time],[0,split3time]]
	rangeDict[('bn', 4)] = [[1e-6, .5], [1e-6, .5], [1e-6, .5], [1e-6, .5]]
	rangeDict[('bn_times', 4)] = [[0,split4time],[0,split4time],[0,split4time],[0,split4time]]
	rangeDict[('mig', '1->2')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '1->3')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '1->4')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '2->1')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '2->3')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '2->4')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '3->1')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '3->2')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '3->4')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '4->1')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '4->2')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig', '4->3')] = [[0,.1],[0,.1], [0,.1]]
	rangeDict[('mig_time', '1->2')] = [[0,split2time],[0,split2time],[0,split2time]]
	rangeDict[('mig_time', '1->3')] = [[0,split3time],[0,split3time],[0,split3time]]
	rangeDict[('mig_time', '1->4')] = [[0,split4time],[0,split4time],[0,split4time]]
	rangeDict[('mig_time', '2->1')] = [[0,split2time],[0,split2time],[0,split2time]]
	rangeDict[('mig_time', '2->3')] = [[0,split3time],[0,split3time],[0,split3time]]
	rangeDict[('mig_time', '2->4')] = [[0,split4time],[0,split4time],[0,split4time]]	
	rangeDict[('mig_time', '3->1')] = [[0,split3time],[0,split3time],[0,split3time]]
	rangeDict[('mig_time', '3->2')] = [[0,split3time],[0,split3time],[0,split3time]]
	rangeDict[('mig_time', '3->4')] = [[0,split4time],[0,split4time],[0,split4time]]
	rangeDict[('mig_time', '4->1')] = [[0,split4time],[0,split4time],[0,split4time]]
	rangeDict[('mig_time', '4->2')] = [[0,split4time],[0,split4time],[0,split4time]]
	rangeDict[('mig_time', '4->3')] = [[0,split4time],[0,split4time],[0,split4time]]
	rangeDict[('Ne_change', 1)] = [[0,1] for i in rangeDict[('Ne', 1)]]
	rangeDict[('Ne_change', 2)] = [[0,1] for i in rangeDict[('Ne', 2)]]
	rangeDict[('Ne_change', 3)] = [[0,1] for i in rangeDict[('Ne', 3)]]
	rangeDict[('Ne_change', 4)] = [[0,1] for i in rangeDict[('Ne', 4)]]
	return rangeDict

############################################
## INPUTFILE <-model params-> PYTHON DICT ##
############################################
def write_paramfile(paramfilename, paramDict):
	#####################
	### Global params ###
	#####################
	openfile = open(paramfilename, 'w')
	openfile.write('length ' + str(paramDict['chromlength']) + '\n')
	openfile.write('mutation_rate ' + str(paramDict['mutation_rate']) + '\n')
	openfile.write('recomb_file ' + str(paramDict['recomfilename']) + '\n')
	openfile.write('gene_conversion_relative_rate ' + str(paramDict['gene_conv_rel_rate']) + '\n')
	openfile.write('\n')

	####################
	### Define pops ####
	####################

	numPops = paramDict['numPops']
	popPairs = []
	for i in range(1, numPops+1):
		for j in range(i, numPops+1):
			if i != j:
				popPair = (i, j)
				popPairs.append(popPair)

	for i in range(1, numPops+1):
		defineline = 'pop_define ' + str(i) + " " + paramDict['labels'][i] + "\n" + "sample_size " + str(i) + " " + str(paramDict['num_indivs_per_sample']) + "\n"
		defineline += 'pop_size ' + str(i) + " " +str(paramDict[('Ne', i)][0]) + "\n" #+ str(paramDict['presentSizes'][i-1]) + "\n"
		openfile.write(defineline)

	openfile.write('\n')

	###########################
	### Define demography ####
	###########################

	towriteLines = []
	towriteAges = []

	# divergence events
	split1time = 1e10
	split2time = paramDict[('split', 2)][0]
	split3time = paramDict[('split', 3)][0]
	split4time = paramDict[('split', 4)][0]

	for i in range(2, numPops+1):
		label = '\"' + paramDict['labels'][i] + ' split\"'
		parentPop = i -1
		timeSplit = paramDict[('split', i)][0]
		splitline = 'pop_event split ' + label + " " + str(parentPop) + " " + str(i) + " " + str(timeSplit) + "\n"
		towriteLines.append(splitline)
		towriteAges.append(timeSplit)

	# changes in pop size

	for i in range(1, numPops+1):
		Ne = paramDict[('Ne', i)]
		Ne_times = paramDict[('Ne_times', i)]
		Ne_change_coded = paramDict[('Ne_change', i)]
		code = {1: 'exp', 0:'normal'}
		Ne_change = [code[int(round(item))] for item in Ne_change_coded]
		assert len(Ne) == len(Ne_times)
		lastgrowthtype = ""
		for j in range(len(Ne)):
			growthtype = Ne_change[j]
			if growthtype == "exp":
				label = '\"' + paramDict['labels'][i] + ' exp Ne change\"'
				if j != len(Ne) - 1:
					changesizeline = "pop_event exp_change_size2 " + label + " " + str(i) + " " + str(Ne_times[j]) +  " " + str(Ne_times[j+1]) + " " + str(Ne[j]) + " " + str(Ne[j +1]) + "\n"
				else: #if the last entry in Ne_change is exp, just treat it as if it were normal
					if i == 1: 
						label = '\"' + paramDict['labels'][i] + ' Ne change\"'
						changesizeline = "pop_event change_size " + label + " " + str(i) + " " + str(Ne_times[j]) + " " + str(Ne[j]) + "\n"
				#	changesizeline = "pop_event exp_change_size2 " + label + " " + str(i) + " " + str(Ne_times[j]) +  " " + str(Ne_times[j] + 50) + " " + str(Ne[j]) + " " + str(Ne[j]) + "\n"
			elif growthtype == "normal":
				label = '\"' + paramDict['labels'][i] + ' Ne change\"'
				changesizeline = "pop_event change_size " + label + " " + str(i) + " " + str(Ne_times[j]) + " " + str(Ne[j]) + "\n"
				if lastgrowthtype == "exp":
					#catch redundant information and comment it out
					changesizeline2 = "#" + changesizeline
					changesizeline = changesizeline2
			if (Ne_times[j] > eval('split' + str(i) + "time")):
				#catch redundant information and comment it out
				changesizeline2 = "#" + changesizeline
				changesizeline = changesizeline2
			lastgrowthtype = growthtype
			towriteLines.append(changesizeline)
			towriteAges.append(Ne_times[j])

		#bottlenecks

		bn = paramDict[('bn', i)]
		bn_times = paramDict[('bn_times', i)]
		assert len(bn) == len(bn_times)
		for j in range(len(bn)):
			label = '\"' + paramDict['labels'][i] + ' bottleneck\"'
			bnline = "pop_event bottleneck " + label + " " + str(i) + " " + str(bn_times[j]) + " " + str(bn[j]) + "\n"
			if (bn_times[j] > eval('split' + str(i) + "time")):
				#catch redundant information and comment it out
				bnline2 = "#" + bnline
				bnline = bnline2
			towriteLines.append(bnline)
			towriteAges.append(bn_times[j])


	#migration

	for popPair in popPairs:
		oneDir = str(popPair[0]) + "->" + str(popPair[1])
		otherDir = str(popPair[1]) + "->" + str(popPair[0])
		timelimit = min(eval("split" + str(popPair[0]) + "time"), eval("split" + str(popPair[1]) + "time"))

		mig = paramDict['mig', oneDir]
		mig_time = paramDict['mig_time', oneDir]
		assert len(mig) == len(mig_time)
		for i in range(len(mig)):
			if mig[i] != 0: #ignore null mig lines, we'll get them at the end
				label = '\"mig ' + paramDict['labels'][popPair[0]] + "->" + paramDict['labels'][popPair[1]] + '\"'
				migline = "pop_event migration_rate " + label + " " + str(popPair[0]) + " " + str(popPair[1]) + " " + str(mig_time[i]) + " " + str(mig[i]) + "\n"

				if (mig_time[i] > timelimit):
					#catch redundant information and comment it out
					migline2 = "#" + migline
					migline = migline2

				towriteLines.append(migline)
				towriteAges.append(mig_time[i])

		mig = paramDict['mig', otherDir]
		mig_time = paramDict['mig_time', otherDir]				
		assert len(mig) == len(mig_time)
		for i in range(len(mig)):
			if mig[i] != 0: #ignore null mig lines, we'll get them at the end
				label = '\"mig ' + paramDict['labels'][popPair[1]] + "->" + paramDict['labels'][popPair[0]]+ '\"'
				migline = "pop_event migration_rate " + label + " " + str(popPair[1]) + " " + str(popPair[0]) + " " + str(mig_time[i]) + " " + str(mig[i]) + "\n"

				if (mig_time[i] > timelimit):
					#catch redundant information and comment it out
					migline2 = "#" + migline
					migline = migline2

				towriteLines.append(migline)
				towriteAges.append(mig_time[i])

	###################################
	### Sort and record demography ####
	###################################

	#####PYTHON 2:
	#from itertools import izip
	#sorted_lines = sorted(izip(towriteAges, towriteLines))
	sorted_lines = sorted(zip(towriteAges, towriteLines))
	for line in sorted_lines:
		openfile.write(line[1])
	openfile.write('\n')

	#migration automatically is set to zero when populations coalesce

	nomiglines = ["pop_event migration_rate \"no mig YRI->GWD\" 1 2 " + str(split2time -1) + " 0\n",
	"pop_event migration_rate \"no mig GWD->YRI\" 2 1 " + str(split2time -1) + " 0\n",
	"pop_event migration_rate \"no mig YRI->MSL\" 1 3 " + str(split2time -1) + " 0\n",
	"pop_event migration_rate \"no mig MSL->YRI\" 3 1 " + str(split2time -1) + " 0\n",
	"pop_event migration_rate \"no mig GWD->MSL\" 2 3 " + str(split2time -1) + " 0\n",
	"pop_event migration_rate \"no mig MSL->GWD\" 3 2 " + str(split2time -1) + " 0\n",
	"pop_event migration_rate \"no mig YRI->YRI\" 1 4 " + str(split4time -1) + " 0\n",
	"pop_event migration_rate \"no mig YRI->YRI\" 4 1 " + str(split4time -1) + " 0\n",
	"pop_event migration_rate \"no mig GWD->YRI\" 2 4 " + str(split4time -1) + " 0\n",
	"pop_event migration_rate \"no mig YRI->GWD\" 4 2 " + str(split4time -1) + " 0\n",
	"pop_event migration_rate \"no mig MSL->YRI\" 3 4 " + str(split4time -1) + " 0\n",
	"pop_event migration_rate \"no mig YRI->MSL\" 4 3 " + str(split4time -1) + " 0\n"]
	for line in nomiglines:
		openfile.write(line)
	openfile.write('\n')
	openfile.close()
	return
def get_dict_from_paramfile(paramfilename):
	#initialize empty lists so we can append from paramfile as we parse
	paramDict = {'numPops':4,'labels':{1:'LWK', 2:'GWD', 3:'MSL', 4:'YRI'}, 'singrate':.5, 'presentSizes':[], 'num_indivs_per_sample':200}
	for pop in range(1,5):
		for parameter in ['Ne', 'Ne_times', 'Ne_change', 'bn', 'bn_times']:
			key = (parameter, pop)
			paramDict[key] = []
		for jpop in range(1,5):
			if pop != jpop:
				migstring = str(pop) + "->" + str(jpop)
				key1, key2 = ('mig', migstring), ('mig_time', migstring)
				paramDict[key1] = []
				paramDict[key2] = []
	openfile = open(paramfilename, 'r')
	for line in openfile:
		entries = line.split()
		if "length" in entries:
			chromlength = int(entries[-1])
			paramDict['chromlength'] = chromlength
		elif "mutation_rate" in entries:
			mutation_rate = float(entries[-1])
			paramDict['mutation_rate'] = mutation_rate
		elif 'recomb_file' in entries:
			recomb_file = entries[-1]
			paramDict['recomfilename'] = recomb_file
		elif "gene_conversion_relative_rate" in entries:
			gene_conversion_relative_rate = float(entries[-1])
			paramDict['gene_conv_rel_rate'] = gene_conversion_relative_rate
		elif "pop_size" in entries:
			currentsize = int(entries[2])
			paramDict['presentSizes'].append(currentsize)
		elif "migration_rate" in entries:
			sourcepop, targetpop, time, rate = int(entries[-4]), int(entries[-3]), float(entries[-2]), float(entries[-1])
			migstring = str(sourcepop) + "->" + str(targetpop)
			key1 = ('mig', migstring)
			key2 = ('mig_time', migstring)
			paramDict[key1].append(rate)
			paramDict[key2].append(time)
		elif "bottleneck" in entries:
			pop, time, rate = int(entries[-3]), float(entries[-2]), float(entries[-1])
			key1 = ('bn', pop)
			key2 = ('bn_times', pop)
			paramDict[key1].append(rate)
			paramDict[key2].append(time)
		elif "exp_change_size" in entries or "exp_change_size2" in entries:
			pop = int(entries[-5])
			finaltime, starttime = float(entries[-4]), float(entries[-3]) #forward-looking concept of 'final/start'
			finalsize, startsize = float(entries[-2]), float(entries[-1])
			key1 = ('Ne', pop)
			key2 = ('Ne_times', pop)
			key3 = ('Ne_change', pop)
			paramDict[key1].append(finalsize)
			paramDict[key2].append(finaltime)
			paramDict[key3].append(1)
		elif "change_size" in entries:
			pop, time, size = int(entries[-3]), float(entries[-2]), float(entries[-1])
			key1 = ('Ne', pop)
			key2 = ('Ne_times', pop)
			key3 = ('Ne_change', pop)
			paramDict[key1].append(size)
			paramDict[key2].append(time)
			paramDict[key3].append(0)
		elif "split" in entries:
			parentpop, childpop, time = int(entries[-3]), int(entries[-2]), float(entries[-1])
			assert parentpop == childpop - 1
			key = ('split', childpop)
			paramDict[key] = [time]
		elif "admix" in entries:
			targetpop, sourcepop, time, rate = int(entries[-4]), int(entries[-3]), float(entries[-2]), float(entries[-1])
			paramDict['admix_rate'] = [rate]
			paramDict['admix_time'] = [time]
	openfile.close()
	#print paramDict
	return paramDict
def update_params(paramDict, keys, indices, values):
	"""paramDict: gets modified and returned.
	keys: indicate which items in paramDict to be modified.
	indices: parallel to keys; indicate which element in paramDict[key] to be modified.
	values: parallel to above; indicate new value to substitute."""
	assert len(keys) == len(indices) and len(keys) == len(values)
	for i in range(len(keys)):
		thisKey, thisIndex, thisValue = keys[i], indices[i], values[i]
		paramDict[thisKey][thisIndex] = thisValue
		#manually handle dependencies (better way to do this?)
		#if thisKey == 'presentSizes':
		#	pop = thisIndex + 1
		#	paramDict[('Ne', pop)][0] = thisValue
		if 'split' in thisKey:
			bnkey = ('bn_times', thisKey[1])
			nekey = ('Ne_times', thisKey[1])
			paramDict[bnkey][0] = thisValue-1
			paramDict[nekey][-1] = thisValue - 1
	return paramDict
