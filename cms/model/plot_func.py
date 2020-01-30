## this file contains functions for cms_modeller
## last updated: 09.10.16	vitti@broadinstitute.org
#### ASW 11.01.2018

from model.bootstrap_func import flattenList
from model.params_func import get_target_values
from model.error_func import read_custom_statsfile
import matplotlib.pyplot as plt
from matplotlib import pylab
import numpy as np

#basedir = "/web/personal/vitti/PEL_"

def plot_comparison(statfilename, numReps, basedir, stats =['pi', 'sfs', 'anc', 'fst', 'r2', 'dprime'], pops=[1,2,3,4], poplabels=['1', '2', '3', '4'], ):
	targetcolor, modelcolor = 'purple', 'orange'
	if "/" in statfilename:
		entries = statfilename.split('/')
		scenariobase = entries[-1]
	else:
		scenariobase = statfilename
	scenariobase = scenariobase.replace('.stat', '')
	
	targetstats = get_target_values()
#	print('\n\n')

#	print(targetstats.keys())
#	print('\n\n')

	modelstats = read_custom_statsfile(statfilename, 4)
	print('loaded stats from '  + statfilename)
	#print(modelstats.keys())
	popPairs = []
	for i in range(len(pops)):
		for j in range(i+1, len(pops)):
			popPair = (pops[i], pops[j])
			popPairs.append(popPair)
	#model_popPairs = [(1,4), (2,4), (3,4)]
	#target_popPairs = [(4,1), (4,2), (4,3)]


	for stat in stats:

		target_vals, target_ses, model_vals, model_ses = [], [], [], []
		if stat != "fst":
			for pop in pops:
				item = (stat, pop)
				target_val = targetstats[item]
				target_vals.append(target_val)
				varkey = (item[0] + "_var", item[1])
				target_se = np.sqrt(targetstats[varkey]) 
				target_ses.append(target_se) 
				model_val = modelstats[item]
				model_vals.append(model_val)
				varkey = (item[0] + "_var", item[1])
				model_se = np.sqrt(modelstats[varkey])
				model_ses.append(model_se)
		elif stat == 'fst':
			for ipair in range(len(popPairs)):#len(model_popPairs)):
				model_popPair =popPairs[ipair]#model_popPairs[ipair]
				target_popPair =popPairs[ipair] #target_popPairs[ipair]
				#print(target_popPair)
			#pair in popPairs:
				target_item = (stat, target_popPair)
				target_val = targetstats[target_item]
				target_vals.append(target_val)
				varkey = (target_item[0] + "_var", target_item[1])
				target_se = np.sqrt(targetstats[varkey]) 
				target_ses.append(target_se) 

				model_item = (stat, model_popPair)
				model_val = modelstats[model_item]
				model_vals.append(model_val)
				varkey = (model_item[0] + "_var", model_item[1])
				model_se = np.sqrt(modelstats[varkey])
				model_ses.append(model_se)
			#print(target_vals)

		if stat == "pi":
			#variance vs standard err
			#num pops
			target_vals = flattenList(target_vals)
			model_vals = flattenList(model_vals)
			target_ses = flattenList(target_ses)
			model_ses = flattenList(model_ses)
			savefilename = basedir + "_pi.png"
			plt.figure()
			plt.title("Pi, target vs. model calcs, " + str(numReps) + " reps")
			ind = [x + .6 for x in np.arange(len(target_vals))]
			ind2 = [x + 1 for x in np.arange(len(target_vals))]
			plt.bar(ind, target_vals, color = targetcolor, width =.4,  label="target")#, yerr=target_ses)
			plt.bar(ind2, model_vals, color = modelcolor, width =.4,  label="model") 
			plt.legend(loc='best')
			plt.xticks(ind2, poplabels)
			plt.xlabel('pop')
			plt.ylabel('ave. per-site heterozygosity')
			plt.savefig(savefilename)
			print("saved to " + savefilename)
			#plt.show()
			plt.close()

		elif stat == 'sfs':
			xlabels = ['sing','<0.1','<0.2','<0.3', '<0.4','<0.5']
			ind = [x + .6 for x in np.arange(len(xlabels))]
			ind2 = [x + 1 for x in np.arange(len(xlabels))]
			savefilename = basedir + "_sfs.png"
		
			fig = plt.figure()
			pylab.subplots_adjust(hspace=0.250)
			#plt.ylabel('Fraction of SNPs')
			for i in range(len(target_vals)):
				#print(i, len(target_vals[i]))
				#print(target_vals[i])
				#target_vals = [float(item) for item in target_vals]
				#print(target_vals[i])
				#print(len(target_vals[i]))
				#print(len(target_ses[i]))

				#print(model_vals[i])
				#print(len(model_vals[i]))
				#print(len(model_ses[i]))

				ax = fig.add_subplot(len(target_vals), 1, i+1, ylabel = poplabels[i], xticks = ind2, xticklabels = xlabels)
				ax.bar(ind, target_vals[i],  label="target", width=.4, color=targetcolor, yerr=target_ses[i])
				ax.bar(ind2, model_vals[i], label="model", width=.4, color=modelcolor, yerr=model_ses[i])
				ax.set_ylim([0,.5])
			plt.xlabel('Minor allele frequency')
			plt.legend(loc='best')
			plt.suptitle("Site frequency spectra, target vs. model calcs, " + str(numReps) + " reps")
			plt.savefig(savefilename)
			print("saved to " + savefilename)
			#plt.show()
			plt.close()


		elif stat == 'anc':
			xlabels = ['sing','<0.1','<0.2','<0.3', '<0.4','<0.5']
			ind = [x + .6 for x in np.arange(len(xlabels))]
			ind2 = [x + 1 for x in np.arange(len(xlabels))]
			savefilename = basedir + "_panc.png"
		
			fig = plt.figure()
			pylab.subplots_adjust(hspace=0.250)
			for i in range(len(target_vals)):
				ax = fig.add_subplot(len(target_vals), 1, i+1, ylabel = poplabels[i], xticks = ind2, xticklabels = xlabels)
				ax.bar(ind, target_vals[i],  label="target", width=.4, color=targetcolor, yerr=target_ses[i])
				ax.bar(ind2, model_vals[i], label="model", width=.4, color=modelcolor, yerr=model_ses[i])
				ax.set_ylim([0,.5])
			plt.legend(loc='best')
			plt.xlabel('Allele frequency')
			plt.suptitle("Percentrage ancestral by frequency, target vs. model calcs, " + str(numReps) + " reps")
			plt.savefig(savefilename)
			print("saved to " + savefilename)
			#plt.show()
			plt.close()

		elif stat == 'r2':
			savefilename = basedir + "_r2.png"
		
			plt.figure()
			plt.title("LD, r2 vs phys dist, target vs. model calcs, " + str(numReps) + " reps")
			xlabels = ['<5kb','<10kb','<15kb','<20kb','<25kb','<30kb','<35kb','<40kb','<45kb','<50kb','<55kb','<60kb','<65kb','<70kb']
			ind = np.arange(len(xlabels))
			fig = plt.figure()
			pylab.subplots_adjust(hspace=0.250)
			for i in range(len(target_vals)):
				ax = fig.add_subplot(len(target_vals), 1, i+1, ylabel = poplabels[i], xticks = ind, xticklabels = xlabels, ylim=[0,.4])
				ax.errorbar(ind, target_vals[i], label="target", color=targetcolor, yerr=target_ses[i])
				ax.errorbar(ind, model_vals[i],  label="model", color=modelcolor, yerr=model_ses[i])
				ax.set_xticklabels(xlabels, rotation=20, fontsize =6)

			plt.suptitle("r2 vs physical distance, target vs. model calcs, " + str(numReps) + " reps")
			plt.savefig(savefilename)
			print("saved to " + savefilename)
			#plt.show()
			plt.legend(loc = 'best')
			plt.close()

		elif stat == 'dprime':
			savefilename = basedir + "_dprime.png"
		
			plt.figure()
			plt.title("LD, P(D'=1) vs gen dist, target vs. model calcs, " + str(numReps) + " reps")
			xlabels = ['<.001cM','<.002cM','<.003cM','<.004cM','<.005cM','<.006cM','<.007cM','<.008cM','<.009cM','<.010cM','<.011cM','<.012cM','<.013cM','<.014cM','<.015cM','<.016cM','<.017cM']
			ind = np.arange(len(xlabels))
			fig = plt.figure()
			pylab.subplots_adjust(hspace=0.250)
			for i in range(len(target_vals)):
				ax = fig.add_subplot(len(target_vals), 1, i+1, ylabel = poplabels[i], xticks = ind, ylim=[0,1]) #xticklabels = xlabels,
				ax.set_xticklabels(xlabels, rotation=20, fontsize =5)
				ax.errorbar(ind, target_vals[i], label="target", color=targetcolor, yerr=target_ses[i])
				ax.errorbar(ind, model_vals[i],  label="model", color=modelcolor, yerr=model_ses[i])
			plt.legend(loc='best')
			plt.suptitle("P(D'=1) vs genetic distance, target vs. model calcs, " + str(numReps) + " reps")
			plt.savefig(savefilename)
			print("saved to " + savefilename)
			#plt.show()
			plt.close()

		elif stat == 'fst':
			popPairLabels = ['(YRI CEU)', '(YRI CHB)', '(YRI PEL)', '(CEU CHB)', '(CEU PEL)', '(CHB PEL)']
			savefilename = basedir + "_fst.png"
		
			#popPairLabels = ['(YRI ASW)', '(CEU ASW)', '(CHB ASW)']
			#print(target_vals)

			target_vals = flattenList(target_vals)
			#print(target_vals)
			model_vals = flattenList(model_vals)
			target_ses = flattenList(target_ses)
			model_ses = flattenList(model_ses)
			#popPairs = target_popPairs
			plt.figure()
			plt.title("Fst, target vs. model calcs, " + str(numReps) + " reps")
			ind = [x + .6 for x in np.arange(len(popPairs))]
			ind2 = [x + 1 for x in np.arange(len(popPairs))]
			plt.bar(ind, target_vals, color = targetcolor, width =.4,  label="target", yerr=target_ses)
			plt.bar(ind2, model_vals, color = modelcolor, width =.4, yerr=model_ses, label="model")
			plt.legend(loc='best')
			plt.xticks(ind2, popPairLabels, fontsize=10)
			plt.xlabel('pop pair')
			plt.ylabel('Fst')
			plt.savefig(savefilename)
			print("saved to " + savefilename)
			#plt.show()
			plt.close()
	return

