## this file contains functions for cms_modeller
## last updated: 09.10.16	vitti@broadinstitute.org

from model.error_func import calc_error
from model.params_func import get_ranges, generate_params, update_params, write_paramfile

########################
## PERTURBING A MODEL ##
########################
def sample_point(nRep, keys, indices, values, writeparamfilename = "/tmp/point.par", scaleParams = False):
	'''given a perturbation of demographic model, writes paramfile and runs point. '''
	if scaleParams:
		ranges = get_ranges()
		actualValues = []
		for i in range(len(keys)):
			key, index, value = keys[i], indices[i], values[i]
			interval = ranges[key][index]
			low, high = float(interval[0]), float(interval[1])
			actual = get_real_value(value, low, high)
			actualValues.append(actual)
	else:
		actualValues = values

	paramDict = generate_params()
	paramDict_updated = update_params(paramDict, keys, indices, actualValues)
	write_paramfile(writeparamfilename, paramDict_updated)

	singrate = paramDict_updated['singrate']
	statfilename = "/tmp/point.err"
	runstatscommand = "coalescent -p " + writeparamfilename + " -n " + str(nRep) + " --drop-singletons " + str(singrate) + "  --custom-stats --genmapRandomRegions > " + statfilename
	#print(runstatscommand)
	subprocess.check_output(runstatscommand, shell=True)
	error = calc_error(statfilename, stats, pops)
	return error
def get_real_value(scaledVal, low, high):
	interval = high - low
	return low + (scaledVal * interval)
def get_scaled_value(actualVal, low, high):
	interval = high - low
	return (actualVal - low) / interval
def get_scaled_values():
	ranges = get_ranges()
	params = generate_params()
	scaledParams = {}
	for key in list(ranges.keys()):
		low, high = float(ranges[key][0]), float(ranges[key][1])
		interval = high - low
		value = params[key]
		scaled = (value - low) / interval
		scaledParams[key] = scaled
	return scaledParams

################################
## EXPLORING PARAMETER-SPACE  ##
################################
def read_dimensionsfile(filename, runType='grid'):
	''' passes parameters for grid search'''
	keys, indices, values = [], [], []
	openfile = open(filename, 'r')
	runnameline = openfile.readline().strip('\n')
	runname = str(runnameline)
	for line in openfile:
		parameterstring = line.strip('\n')
		entries = parameterstring.split('\t')
		key, index = eval(entries[0]), eval(entries[1])
		if runType == 'grid':
			value_array = eval(entries[2:])
			values.append(value_array)
		keys.append(key)
		indices.append(index)
	openfile.close()
	if runType == 'grid':
		return runname, keys, indices, values
	elif runType == 'optimize':
		return runname, keys, indices
