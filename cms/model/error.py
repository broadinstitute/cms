## AWS xp_run: ~/COSI_STATIC/
## 12.04.15

import math
from parse import getTargetValues, readIScustomstatfile

def RMSE(output, target, target_var, DEBUG = False):
	#each of these arguments is a list of the same length (num bins)
	n = int(len(output))
	sum_bins = 0

 	for i in range(n):
		a = float(output[i])
		c = float(target[i])
		d = float(target_var[i])
		partialsum = ((a-c)**2)/(math.sqrt(d))
		if DEBUG == True:
			print "(" + str(a) + "-" + str(c) + ")^2  / " + str(d) 
			print partialsum
		sum_bins += partialsum

	return (sum_bins/n)**.5

def calcError(statfilename, stats = ['pi', 'sfs', 'anc', 'r2', 'dprime', 'fst'], pops = [1, 2, 3, 4], DEBUG=False): 
	"""takes in an array of target values (y). extracts y' values (i.e. model
	outputs) from statfile and compares to get root mean square error."""
	targetStats = getTargetValues() 
	simStats = readIScustomstatfile(statfilename, 4)
	print "getting statfile from " + statfilename
	
	popPairs = []
	for i in range(len(pops)):
		for j in range(i, len(pops)):
			if i != j:
				popPair = (pops[i], pops[j])
				popPairs.append(popPair)
	if DEBUG == True:
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

				RMS = RMSE(sim_val, target_val, target_var)
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

				RMS = RMSE(sim_val, target_val, target_var)
				if stat == "pi": 
					RMS *= 100
				print "RMS, " + str(item) + ": " + str(RMS) + "\t"
				tot += (RMS ** 2)
				counter +=1

	ave = tot / counter
	return ave**.5