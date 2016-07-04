## {{POP GEN DATA --> DEM MODEL}}
## this script contains functions pertaining to the error function used to evaluate the goodness-of-fit of a model to the specified dataset.
## by default, we use the root mean square error statistic (cf Schaffner 2005) for all statistics and populations. 
## this can be toggled in order to fit certain parts of the model: e.g., we might explore parameter-space by letting ancient migration rates vary (and no other parameters) and focus our error function on the effects on Fst
## by default we up-weigh the importance of heterozygosity; this statistic has wide variance across simulates of a given model and can fail to inflate the error function while traversing unrealistic corners of parameter-space if not compensated for. [the 100x scale factor is pragmatic and can easily be toggled]
## last updated: 07.04.16 	vitti@broadinstitute.org

import math
from params import get_target_values
from parse import read_custom_statsfile

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
			print "(" + str(a) + "-" + str(c) + ")^2  / " + str(d) 
			print partialsum
		sum_bins += partialsum
	return (sum_bins/n)**.5
def calc_error(statfilename, stats = ['pi', 'sfs', 'anc', 'r2', 'dprime', 'fst'], pops = [1, 2, 3, 4], piScaleFactor = 100, verbose=False): 
	"""takes in an array of target values (y). extracts y' values (i.e. model
	outputs) from statfile and compares to get root mean square error."""
	targetStats = getTargetValues() 
	simStats = readcustomstatfile(statfilename, len(pops))
	if verbose:
		print "getting statfile from " + statfilename
	popPairs = []
	for i in range(len(pops)):
		for j in range(i, len(pops)):
			if i != j:
				popPair = (pops[i], pops[j])
				popPairs.append(popPair)
	if verbose:
		print "pop pairs: "
		print popPairs
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
					print "RMS, " + str(item) + ": " + str(RMS) + "\t"
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
					print "RMS, " + str(item) + ": " + str(RMS) + "\t"
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
