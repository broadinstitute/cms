##	12.17.2018
bootstrap_targetval_filename = "/idi/sabeti-scratch/jvitti/remodel/fit/targetstats.txt"
#"/idi/sabeti-scratch/jvitti/PEL_model/target_test_ldfixed2.txt"
#"/idi/sabeti-scratch/jvitti/PEL_model/target_test.txt"

##########################
## DEFINE TARGET VALUES ##
##########################

def read_lines(openfile, numlines):
	for i in range(numlines):
		line = openfile.readline()
	return line

def get_target_values():
	#print('getting TARGET VALS from ' + bootstrap_targetval_filename) #load VARIANCE as it says!
	popDict = {'YRI': 1, 'CEU': 2, 'CHB': 3, 'BEB':4}
	stats = {}
	openfile = open(bootstrap_targetval_filename, 'r')
	for ipop in range(1,5):

		popline = read_lines(openfile, 1)
		popline = popline.strip('\n')
		ipop = popDict[popline]
		
		piline = read_lines(openfile, 1)
		entries = piline.split()
		#pi_mean, pi_se = float(entries[0]), float(entries[1])
		pi_mean = float(entries[0])
		stats[('pi', ipop)] = [pi_mean]
		#stats[('pi_var', ipop)] = [pi_se ** 2]

		sfs_mean_line = read_lines(openfile, 1)
		entries = sfs_mean_line.split()
		sfs_mean = eval(sfs_mean_line)
		stats[('sfs', ipop)] = sfs_mean
		#sfs_se_line = read_lines(openfile, 1)
		#sfs_se = eval(sfs_se_line)
		#sfs_var = [float(x)**2 for x in sfs_se]
		#stats[('sfs_var', ipop)] = sfs_var

		anc_mean_line = read_lines(openfile, 1)
		anc_mean = eval(anc_mean_line)
		stats[('anc', ipop)] = anc_mean
		#anc_se_line = read_lines(openfile, 1)
		#anc_se = eval(anc_se_line)
		#anc_var = [float(x)**2 for x in anc_se]
		#stats[('anc_var', ipop)] = anc_var


	#for ipop in range(1, 5):
		r2_mean_line = read_lines(openfile, 1)
		r2_mean = eval(r2_mean_line)
		stats[('r2', ipop)] = r2_mean
		#r2_se_line = read_lines(openfile, 1)
		#r2_se = eval(r2_se_line)
		#r2_var = [float(x)**2 for x in r2_se]
		#stats[('r2_var', ipop)] = r2_var

		dprime_mean_line = read_lines(openfile, 1)
		dprime_mean = eval(dprime_mean_line)
		stats[('dprime', ipop)] = dprime_mean
		#dprime_se_line = read_lines(openfile, 1)
		#dprime_se = eval(dprime_se_line)
		#dprime_var = [float(x)**2 for x in dprime_se]
		#stats[('dprime_var', ipop)] = dprime_var


		pi_se_line = read_lines(openfile, 1)
		entries = pi_se_line.split()
		pi_se = float(entries[0])
		stats[('pi_var', ipop)] = [pi_se ** 2]

		sfs_se_line = read_lines(openfile, 1)
		sfs_se = eval(sfs_se_line) #WAIT RIGHT? IS THIS HOW I RECORDED IT?
		sfs_var = [float(x)**2 for x in sfs_se]
		stats[('sfs_var', ipop)] = sfs_var

		anc_se_line = read_lines(openfile, 1)
		anc_se = eval(anc_se_line)
		anc_var = [float(x)**2 for x in anc_se]
		stats[('anc_var', ipop)] = anc_var
		
		r2_se_line = read_lines(openfile, 1)
		r2_se = eval(r2_se_line)
		r2_var = [float(x)**2 for x in r2_se]
		stats[('r2_var', ipop)] = r2_var

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
	############################ #FAST FOR FITTING
	paramDict = {'chromlength':100000, 'mutation_rate':1.25e-08, 'num_indivs_per_sample':170, 'gene_conv_rel_rate':2.3, 'singrate':.25,
	#paramDict = {'chromlength':1000000, 'mutation_rate':1.25e-08, 'num_indivs_per_sample':170, 'gene_conv_rel_rate':2.3, 'singrate':.25,
	'recomfilename':'results_inclusive_test_filter.recom'}

	paramDict['numPops'] = 4
	paramDict['labels'] = {1:'YRI', 2:'CEU', 3:'CHB', 4:'PEL'}	

	###################
	##DEFINE BRANCHES##
	###################

	paramDict[('split', 2)] = [2250]#[2025]
	paramDict[('split', 3)] = [1400]#[1800]
	paramDict[('split', 4)] = [400]#500]#[700]#[1600]#[1200]  #this wants to be more like 400
									
									#Indexing: At THIS TIME, arrive at THIS SIZE via THIS KIND of change (forward not coalescent perspective)
									#0			#1			#2			#3							#5				#6
									#present 	#midexp		#pre exp 	minimum						#anc			#base
	paramDict[('Ne', 1)] = 			[500000,	10000,		6000,		25000,				25000,			17500]
	paramDict[('Ne_times', 1)] = 	[1,			80,			85,			600,				4000,			10000]
	paramDict[('Ne_change', 1)] = 	[1,			1,			1,			1,					0,					0] #code = {1: 'exp', 0:'normal'}

	paramDict[('Ne', 2)] = 			[500000,	2500,		4000,	 1250]
	paramDict[('Ne_times', 2)] =	[1, 		50,			300,	2000]
	paramDict[('Ne_change', 2)] = 	[1, 		1,			1,		1]

	paramDict[('Ne', 3)] = 			[750000,	15000,		2000, 	1000]
	paramDict[('Ne_times', 3)] = 	[1, 		75, 		300, 	900]
	paramDict[('Ne_change', 3)] = 	[1, 		1,			0, 		1]

	paramDict[('Ne', 4)] = 			[100000,		5000,		1000]
	paramDict[('Ne_times', 4)] = 	[1,			25, 		100]#paramDict[('split', 4)][0]-1]
	paramDict[('Ne_change', 4)] = 	[1, 		1, 			1]

	##############
	##POP EVENTS##
	##############
 
	paramDict[('bn', 1)] = [.005, .001, 1e-10]
	paramDict[('bn_times', 1)] = [100, 500, 5000]

	paramDict[('bn', 2)] = [0.01, 0.0005]
	paramDict[('bn_times', 2)] = [100, 2000]

	paramDict[('bn', 3)] = [0.035, 0.1]
	paramDict[('bn_times', 3)] = [100, 1500]

	paramDict[('bn', 4)] = [.05, .05]#[1.00000052029e-06, 1.00000031898e-06]
	paramDict[('bn_times', 4)] = [10, 100]#[paramDict[('split', 4)][0] - 1, 25] 

	#############
	##GENE FLOW##
	#############

	#YRI,CEU
	paramDict[('mig', '1->2')] = [.005, .0001, .1e-10]
	paramDict[('mig_time', '1->2')] = [0, 10, 149]
	paramDict[('mig', '2->1')] = [.01, .0001, 1e-10]
	paramDict[('mig_time', '2->1')] = [0, 10, 149]

	#YRI,CHB
	paramDict[('mig', '1->3')] =  [.005, .0001, 1e-10]
	paramDict[('mig_time', '1->3')] = [0, 10, 149]
	paramDict[('mig', '3->1')] = [.00975, .0001, 1e-10]
	paramDict[('mig_time', '3->1')] = [0, 10, 149]

	#YRI,PEL
	paramDict[('mig', '1->4')] = [.005, .0001, 0]
	paramDict[('mig_time', '1->4')] = [0, 10, 149]
	paramDict[('mig', '4->1')] = [.016, .0001, 0]
	paramDict[('mig_time', '4->1')] = [0, 10, 149]

	#CEU,CHB
	paramDict[('mig', '2->3')] = [.001, .0001, 0]
	paramDict[('mig_time', '2->3')] = [0, 10, 149]
	paramDict[('mig', '3->2')] = [.001, .0001, .00025]
	paramDict[('mig_time', '3->2')] = [0, 10, 149]

	#CEU,PEL (ALSO: admix)
	paramDict[('mig', '2->4')] = [.001, 1e-09, 0] #LESS?
	paramDict[('mig_time', '2->4')] = [0, 10, 149]
	paramDict[('mig', '4->2')] = [.00005, 1e-09, .0005]
	paramDict[('mig_time', '4->2')] = [0, 10, 149]

	#CHB,PEL
	paramDict[('mig', '3->4')] = [.001, 0, 0] #LESS?
	paramDict[('mig_time', '3->4')] = [0, 10, 149]
	paramDict[('mig', '4->3')] = [.001, 0, 0]
	paramDict[('mig_time', '4->3')] = [0, 10, 149]

	paramDict['admix_rate'] = [.3]#[1e-10]#[.267]
	paramDict['admix_time'] = [10]

	return paramDict
def get_ranges():
	paramDict = generate_params()
	split2time = paramDict[('split', 2)][0]
	split3time = paramDict[('split', 3)][0]
	split4time = paramDict[('split', 4)][0]

	rangeDict = {}
	rangeDict[('split', 2)] = [[2000, 5000]]
	rangeDict[('split', 3)] = [[1400, 2000]]
	rangeDict[('split', 4)] = [[400, 1400]]
	rangeDict[('Ne', 1)] = [[50000, 10000000], [5000, 500000],  [1000, 100000], [1000, 100000], [1000, 1000000], [1000, 1000000]]#, [100, 100000]]
	rangeDict[('Ne_times', 1)] = [[1, 5], [6, 84], [85, 399], [400, 1999], [2000, 11999], [12000, 25000]]#, [25000, 50000]]
	rangeDict[('Ne', 2)] = [[50000, 10000000], [5000, 500000],  [1000, 100000], [1000, 1000000]]
	rangeDict[('Ne_times', 2)] = [[1, 5], [6, 199], [200, 499], [500, split2time]]
	rangeDict[('Ne', 3)] = [[50000, 10000000], [5000, 500000],  [1000, 100000], [1000, 1000000],]
	rangeDict[('Ne_times', 3)] = [[1, 5], [6, 249], [250, 749], [750, split3time]] 
	rangeDict[('Ne', 4)] = [[50000, 10000000], [5000, 500000],  [1000, 100000]]
	rangeDict[('Ne_times', 4)] = [[1, 5], [6, 25], [25, split4time]]
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

	rangeDict['admix_rate'] = [[0,1]]
	rangeDict['admix_time'] = [[0, split4time]]

	return rangeDict

############################################
## INPUTFILE <-model params-> PYTHON DICT ##
############################################
def write_paramfile(paramfilename, paramDict):
	#print(paramDict[('Ne_times', 3)])
	#print(paramDict[('Ne', 3)])

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

	## ADMIX???
	admixTime = paramDict['admix_time'][0]
	#print(paramDict['admix_rate'])
	admixRate = paramDict['admix_rate'][0]
	#print(admixTime, admixRate)
	#if type(admixRate) == list:
	#	admixRate = float(admixRate[0])
	admix_line = "pop_event admix \"PEL admix\" 4 2 " + str(admixTime) + " " + str(admixRate)  + "\n"
	towriteLines.append(admix_line)
	towriteAges.append(admixTime)

	#print("PARAMLINESN=",str(len(towriteLines)))

	###################################
	### Sort and record demography ####
	###################################
	sorted_lines = sorted(zip(towriteAges, towriteLines))
	for line in sorted_lines:
		openfile.write(line[1])
	openfile.write('\n')

	#migration automatically is set to zero when populations coalesce
	nomiglines = ["pop_event migration_rate \"no mig YRI->GWD\" 1 2 " + str(split2time -1) + " 0\n",
	"pop_event migration_rate \"no mig CEU->YRI\" 2 1 " + str(split2time -1) + " 0\n",
	"pop_event migration_rate \"no mig YRI->CHB\" 1 3 " + str(split3time -1) + " 0\n",
	"pop_event migration_rate \"no mig CHB->YRI\" 3 1 " + str(split3time -1) + " 0\n",
	"pop_event migration_rate \"no mig CEU->CHB\" 2 3 " + str(split3time -1) + " 0\n",
	"pop_event migration_rate \"no mig CHB->CEU\" 3 2 " + str(split3time -1) + " 0\n",
	"pop_event migration_rate \"no mig YRI->PEL\" 1 4 " + str(split4time -1) + " 0\n",
	"pop_event migration_rate \"no mig PEL->YRI\" 4 1 " + str(split4time -1) + " 0\n",
	"pop_event migration_rate \"no mig CEU->PEL\" 2 4 " + str(split4time -1) + " 0\n",
	"pop_event migration_rate \"no mig PEL->CEU\" 4 2 " + str(split4time -1) + " 0\n",
	"pop_event migration_rate \"no mig CHB->PEL\" 3 4 " + str(split4time -1) + " 0\n",
	"pop_event migration_rate \"no mig PEL->CHB\" 4 3 " + str(split4time -1) + " 0\n"]
	for line in nomiglines:
		openfile.write(line)
	openfile.write('\n')
	openfile.close()

	#nLines = 0
	#openfile = open(paramfilename, 'r')
	#for line in openfile:
	#	nLines +=1
	#	if "exp_change_size2" in line and "CHB" in line:
	#		print(line)
	#openfile.close()
	#print('param file has ' + str(nLines) + " lines in it")
	return
def get_dict_from_paramfile(paramfilename):
	#initialize empty lists so we can append from paramfile as we parse
	paramDict = {'numPops':4,'labels':{1:'YRI', 2:'CEU', 3:'CHB', 4:'PEL'}, 'singrate':.25, 'presentSizes':[], 'num_indivs_per_sample':200}
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
			paramDict['admix_rate'] = rate#[rate]
			paramDict['admix_time'] = time#[time]
	openfile.close()
	return paramDict
def update_params(paramDict, keys, indices, values):
	"""paramDict: gets modified and returned.
	keys: indicate which items in paramDict to be modified.
	indices: parallel to keys; indicate which element in paramDict[key] to be modified.
	values: parallel to above; indicate new value to substitute."""

	# COULD ADD A CHECK HERE TO ENFORCE BOUNDS by first calling (or passing) (get_)rangeDict

	assert len(keys) == len(indices) and len(keys) == len(values)
	for i in range(len(keys)):
		thisKey, thisIndex, thisValue = keys[i], indices[i], values[i]
		paramDict[thisKey][thisIndex] = thisValue
		print("\tupdateParams:\t", thisKey, thisIndex,thisValue)
		if 'split' in thisKey:
			bnkey = ('bn_times', thisKey[1])
			nekey = ('Ne_times', thisKey[1])
			paramDict[bnkey][0] = thisValue-1
			paramDict[nekey][-1] = thisValue - 1
	return paramDict


