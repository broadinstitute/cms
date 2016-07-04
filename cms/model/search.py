## this file contains functions for cms_modeller
## last updated: 07.04.16	vitti@broadinstitute.org

#from error import read_error_dimensionsfile

def read_grid_dimensionsfile(filename):
	''' passes parameters for grid search'''
	openfile = open(filename, 'r')
	gridnameline = openfile.readline().strip('\n')
	gridname = str(gridnameline)
	keys = eval(openfile.readline())
	indices = eval(openfile.readline())
	#values = eval(openfile.readline())
	openfile.close()
	return gridname, keys, indices#, values
