# quick adhoc script to rename binfiles

modelpops = [1, 2, 3, 4]
models = [ 'nulldefault', 'gradient_101915_treebase_6_best', 'default_default_101715_12pm', 'default_112115_825am','nulldefault_constantsize',]
basedir = '/idi/sabeti-scratch/jvitti/clean/scores/'

import subprocess

def execute(command):
	subprocess.check_output(command.split())
	return

def get_concat_files(model, pop, score, altpop = '', basedir = '/idi/sabeti-scratch/jvitti/clean/scores/'):
	if score in ['ihs', 'delihh', 'nsl']:
		concatfilebase = basedir + model + "/neut/concat_" + str(pop) + "_"
	elif score in ['xpehh', 'fst']:
		concatfilebase = basedir + model + "/neut/concat_" + str(pop) + "_" + str(altpop) + "_"
	else:
		concatfilebase = ""
	concatfilename = concatfilebase + score + ".txt"
	binfilename = concatfilebase + score + ".bins"
	return concatfilename, binfilename


def main():
	irep = 1
	for model in models:
		modeldir = basedir + model
		#check_create_dir(modeldir)
		neutmodeldir = modeldir + "/neut/"
		#check_create_dir(neutmodeldir)
		for pop in modelpops:
			for score in ['ihs', 'delihh', 'nsl',]:
				binfilename = "/idi/sabeti-scratch/jvitti/cms/cms/norm_componentscores_112216.sh.e4022856." + str(irep)
				concat, binpostname = get_concat_files(model, pop, score)
				movecommand = "mv " + binfilename + " " + binpostname
				#print(movecommand)
				execute(movecommand)
				irep +=1

			for score in ['fst', 'xpehh']:
				altpops = modelpops[:]
				altpops.remove(pop)
				for altpop in altpops:
					binfilename = "/idi/sabeti-scratch/jvitti/cms/cms/norm_componentscores_112216.sh.e4022856." + str(irep)
					concat, binpostname = get_concat_files(model, pop, score, altpop)
					movecommand = "mv " + binfilename + " " + binpostname
					#print(movecommand)	
					execute(movecommand)		
					irep +=1

main()
