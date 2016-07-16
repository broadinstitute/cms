## functions for composite
## last updated 07.15.16 	vitti@broadinstitute.org

def get_likesfiles_frommaster(masterfile):
	'''ad hoc, for now'''
	openfile = open(masterfile, 'r')
	ihs_hit_filename = openfile.readline().strip('\n') 
	ihs_miss_filename= openfile.readline().strip('\n') 
	delihh_hit_filename = openfile.readline().strip('\n') 
	delihh_miss_filename= openfile.readline().strip('\n') 
	xpehh_hit_filename= openfile.readline().strip('\n') 
	xpehh_miss_filename= openfile.readline().strip('\n') 
	fst_hit_filename= openfile.readline().strip('\n') 
	fst_miss_filename= openfile.readline().strip('\n') 
	deldaf_hit_filename= openfile.readline().strip('\n')  
	deldaf_miss_filename= openfile.readline().strip('\n') 
	openfile.close()
	return delihh_hit_filename, delihh_miss_filename, ihs_hit_filename, ihs_miss_filename, xpehh_hit_filename, xpehh_miss_filename, fst_hit_filename, fst_miss_filename, deldaf_hit_filename, deldaf_miss_filename


