## {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS}
# JV: -> should connect this also to scans.py/initial filtration. 
## script to dispatch job to thin tped files, in order to ensure consistent sample size for the purpose of calculating summary stats (target vals for dem model pipeline)
## last updated: 06.04.16 	vitti@broadinstitute.org
#also I need to centralize the function that dispatches taskarrays...

#####################
## TOGGLE THIS RUN ##
#####################

dateprefix, fileprefix = "120915_3pm_chr21-22",'results_inclusive'
fileloc = "/idi/sabeti-scratch/jvitti/dem_model/target_vals/neutral/bootstrap/thinned/"
prefix = "thin_to_MSL_size_120315"
pops = ['GWD', 'YRI', 'MSL', 'ESN', 'LWK']
chroms = range(1,23)
num_indivs = {'GWD': 113, 'ACB': 96, 'BEB': 86, 'PEL': 85, 'LWK': 99, 'MSL': 85, 'GBR': 91, 'IBS': 107, 'ASW': 61, 'TSI': 107, 'KHV': 99, 'CEU': 99, 'YRI': 108, 'CHB': 103, 'STU': 102, 'CHS': 105, 'ESN': 99, 'FIN': 99, 'GIH': 103, 'PJL': 96, 'MXL': 64, 'PUR': 104, 'ITU': 102, 'CDX': 93, 'JPT': 104, 'CLM': 94}
print "should store and load num_indivs some place centrally."

##########
## MAIN ##
##########

def main():
	limiting = min([num_indivs[x] for x in pops])
	nChromToTake = 2 * limiting #diploid 

	print "also I need to centralize the function that dispatches taskarrays..."
	#write file with arguments for all tasks
	t = 0
	taskargfilename = prefix + "_taskargs.txt"
	taskargfile = open(taskargfilename, "w")
	for chrom in chroms:
		chrom = str(chrom)
		for pop in pops:
			taskargfile.write(pop + "\t" + str(chrom) + "\t" + str(nChromToTake) + "\n")
			t += 1
	taskargfile.close() 

	#write script for one job that points to task-argument-file
	scriptname = prefix + '.sh'
	scriptfile = open(scriptname, 'w')
	for UGERoption in UGERoptions:
		scriptfile.write("#$ " + UGERoption + "\n")
	scriptfile.write("#$ -t 1-" + str(t) + "\n")
	scriptfile.write("\n")
	scriptfile.write("MYFILE=" + str(taskargfilename) + "\n")
	scriptfile.write("SAMPLE=$(awk \"NR==$SGE_TASK_ID\" $MYFILE)\n")
	commandstring = "./thin_tped $SAMPLE\n"
	scriptfile.write(commandstring)
	scriptfile.close()
	subcommand = "qsub " + scriptname

	print subcommand

main()