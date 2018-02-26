## !!-----this is where the user specifies the ERROR FUNCTION to quantify model goodness-of-fit-------!!
## by default, we use the root mean square error statistic (cf Schaffner 2005) for all statistics and populations. 
## this can be toggled in order to fit certain parts of the model: e.g., we might explore parameter-space by letting ancient migration rates vary (and no other parameters) and focus our error function on the effects on Fst
## by default we up-weigh the importance of heterozygosity; this statistic has wide variance across simulates of a given model and can fail to inflate the error function while traversing unrealistic corners of parameter-space if not compensated for. [the 100x scale factor is pragmatic and can easily be toggled]
## last updated: 07.04.16 	vitti@broadinstitute.org

#from model.params_func import get_target_values
import math

#generated for dataset using cms_modeller.py, e.g.
bootstrap_targetval_filename = "/idi/sabeti-scratch/jvitti/targetstats_120317_quick_TEST10.txt"
#"/idi/sabeti-scratch/jvitti/targetstats_120317_quick_TEST.txt"#"/idi/sabeti-scratch/jvitti/test.txt"
#"/Users/vitti/n2stats.txt"
#"results_inclusive_targetstats_thinned_bootstrap_bysnp_2reps_120915_3pm_chr21-22.txt"


##########################
## DEFINE TARGET VALUES ##
##########################

def read_lines(openfile, numlines):
	for i in range(numlines):
		line = openfile.readline()
	return line
def get_target_values():
	popDict = {'YRI': 1, 'CEU': 2, 'CHB': 3, 'BEB':4}
	#popDict = {'LWK': 1, 'GWD': 2, 'MSL': 3, 'YRI':4, 'ESN':5}
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

def read_custom_statsfile(statfilename, numPops):
	stats = {}
	openfile = open(statfilename, 'r')
	for ipop in range(1, numPops+1): 
		openfile.readline() #pop label
		stats[('pi', ipop)] = [float(openfile.readline())]
		stats[('sfs', ipop)]  = [float(x) for x in openfile.readline().split()]
		stats[('anc', ipop)]  = [float(x) for x in openfile.readline().split()]
		stats[('r2', ipop)]  = [float(x) for x in openfile.readline().split()]
		stats[('dprime', ipop)]  = [float(x) for x in openfile.readline().split()]
		stats[('pi_var', ipop)] = [float(openfile.readline())]
		stats[('sfs_var', ipop)]  = [float(x) for x in openfile.readline().split()]
		stats[('anc_var', ipop)]  = [float(x) for x in openfile.readline().split()]
		stats[('r2_var', ipop)]  = [float(x) for x in openfile.readline().split()]
		stats[('dprime_var', ipop)]  = [float(x) for x in openfile.readline().split()]
	popPairs = []
	for ipop in range(1, numPops+1):
		for jpop in range(ipop, numPops+1):
			if ipop != jpop: 
				popPairs.append((ipop, jpop))
	for i in range(len(popPairs)):
		line = openfile.readline()
		entries = line.split('\t')
		fst, fst_var = float(entries[1]), float(entries[2])
		stats[('fst', popPairs[i])] = [fst]
		stats[('fst_var', popPairs[i])] = [fst_var]
	openfile.close()
	#print(stats)
	return stats
def root_mean_square_error(output, target, target_var, verbose = False):
	#each of these arguments is a list of the same length (num bins)
	n = int(len(output))
	sum_bins = 0
	for i in range(n):
		a = float(output[i])
		c = float(target[i])
		d = float(target_var[i])
		partialsum = ((a-c)**2)/(math.sqrt(d))
		if verbose == True:
			print("(" + str(a) + "-" + str(c) + ")^2  / " + str(d))
			print(str(partialsum))
		sum_bins += partialsum
	return (sum_bins/n)**.5
def calc_error(statfilename, stats = ['pi', 'sfs', 'anc', 'r2', 'dprime', 'fst'], pops = [1, 2, 3, 4], piScaleFactor = 100, verbose=False): 
	"""takes in an array of target values (y). extracts y' values (i.e. model
	outputs) from statfile and compares to get root mean square error."""
	targetStats = get_target_values() 
	simStats = read_custom_statsfile(statfilename, len(pops))
	if verbose:
		print("getting statfile from " + statfilename)
	popPairs = []
	for i in range(len(pops)):
		for j in range(i, len(pops)):
			if i != j:
				popPair = (pops[i], pops[j])
				popPairs.append(popPair)
	if verbose:
		print("pop pairs: ")
		print(popPairs)
	tot = 0
	counter = 0
	for stat in stats:
		if stat == 'fst':
			for popPair in popPairs:
				item = (stat, popPair)
				sim_val = simStats[item]
				target_val = targetStats[item]
				varkey = (item[0] + "_var", item[1])
				sim_var = simStats[varkey]
				target_var = targetStats[varkey]

				RMS = root_mean_square_error(sim_val, target_val, target_var)
				if verbose:
					print("RMS, " + str(item) + ": " + str(RMS) + "\t")
				tot += (RMS ** 2)
				counter +=1
		else:
			for pop in pops:
				item = (stat, pop)
				sim_val = simStats[item]
				target_val = targetStats[item]
				varkey = (item[0] + "_var", item[1])
				sim_var = simStats[varkey]
				target_var = targetStats[varkey]

				RMS = root_mean_square_error(sim_val, target_val, target_var)
				if stat == "pi": 
					RMS *= piScaleFactor
				if verbose:
					print("RMS, " + str(item) + ": " + str(RMS) + "\t")
				tot += (RMS ** 2)
				counter +=1

	ave = tot / counter
	return ave**.5
def read_error_dimensionsfile(filename):
	''' passes parameters for error function. '''
	openfile = open(filename, 'r')
	statline = openfile.readline().strip('\n')
	stats = statline.split(',')
	popline = openfile.readline().strip('\n')
	pops = openfile.split(',')
	openfile.close()
	return stats, pops
