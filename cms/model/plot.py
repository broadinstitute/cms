## {{POP GEN DATA --> DEM MODEL}}
## this script contains functions for visualizing (changes in) goodness-of-fit during model-fitting
## currently configured to fit 1KG P3 populations: {YRI, CEU, CHB, BEB}
## last updated: 06.04.16	vitti@broadinstitute.org

from parse import readIScustomstatfile, getTargetValues
import matplotlib as mp
mp.use('agg')
import matplotlib.pyplot as plt
import pylab
import numpy 
import itertools
#plt.style.use('ggplot')

###################
## HELER METHODS ##
###################

def flattenList(list3d):
	flat = []
	for list1d in list3d:
		for item in list1d:
			flat.append(item)
	return flat

#######################################
## VISUALIZE MODEL OUTPUT STATISTICS ##
#######################################

def plot_before_and_after_comparison(statfilenamepre, statfilenamepost, savestring):
	stats = ['pi', 'sfs', 'anc', 'r2', 'dprime', 'fst']
	numReps = 1000
	pops = [1,2,3,4]#['CEU', 'BEB', 'CHB', 'YRI'] #etc
	
	poplabels = ['YRI', 'CEU', 'CHB', 'BEB']
	popPairs = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
	popPairLabels = ['(YRI CEU)', '(YRI CHB)', '(YRI BEB)', '(CEU CHB)', '(CEU BEB)', '(CHB BEB)']
	tosavedir = "/home/unix/vitti/public_html/dem_model/NRE_110915/prepost/"

	targetcolor, modelcolor = 'purple', 'orange'
	modelcolorpost = 'green'

	if "/" in statfilenamepre:
		entries = statfilenamepre.split('/')
		scenariobase = entries[-1]
	else:
		scenariobase = statfilenamepre

	scenariobase = statfilenamepre.replace('_pre.stat', '')

	print "Plotting figures at: http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/prepost"
	print "Saving files as: " + savestring
	import numpy as np

	targetstats = getTargetValues()
	modelstatspre = readIScustomstatfile(statfilenamepre, 4)
	modelstatspost = readIScustomstatfile(statfilenamepost, 4)

	for stat in stats:
		tosavefilename = tosavedir + savestring + "_" + stat + "_prepost.png"
		target_vals, target_ses = [], []
		model_vals_pre, model_ses_pre = [], []
		model_vals_post, model_ses_post = [], []

		if stat != "fst":
			for pop in pops:
				item = (stat, pop)

				target_val = targetstats[item]
				target_vals.append(target_val)
				varkey = (item[0] + "_var", item[1])
				target_se = np.sqrt(targetstats[varkey]) 
				target_ses.append(target_se) 
				
				model_val_pre = modelstatspre[item]
				model_val_post = modelstatspost[item]
				model_vals_pre.append(model_val_pre)
				model_vals_post.append(model_val_post)
				varkey = (item[0] + "_var", item[1])
				model_se_pre = np.sqrt(modelstatspre[varkey])
				model_se_post = np.sqrt(modelstatspost[varkey])
				model_ses_pre.append(model_se_pre)
				model_ses_post.append(model_se_post)

		elif stat == 'fst':
			for pair in popPairs:
				item = (stat, pair)
				target_val = targetstats[item]
				target_vals.append(target_val)
				varkey = (item[0] + "_var", item[1])
				target_se = np.sqrt(targetstats[varkey]) 
				target_ses.append(target_se)

				model_val_pre = modelstatspre[item]
				model_val_post = modelstatspost[item]
				model_vals_pre.append(model_val_pre)
				model_vals_post.append(model_val_post)
				varkey = (item[0] + "_var", item[1])
				model_se_pre = np.sqrt(modelstatspre[varkey])
				model_se_post = np.sqrt(modelstatspost[varkey])
				model_ses_pre.append(model_se_pre)
				model_ses_post.append(model_se_post)

		if stat == "pi":
			#variance vs standard err
			#num pops
			target_vals = flattenList(target_vals)
			model_vals_pre = flattenList(model_vals_pre)
			model_vals_post = flattenList(model_vals_post)
			target_ses = flattenList(target_ses)
			model_ses_post = flattenList(model_ses_post)
			model_ses_pre = flattenList(model_ses_pre)

			plt.figure()
			plt.title("Pi, target vs. model calcs, " + str(numReps) + " reps")
			ind = [x + .3 for x in np.arange(len(target_vals))]
			ind2 = [x + .6 for x in np.arange(len(target_vals))]
			ind3 = [x + .9 for x in np.arange(len(target_vals))]
			plt.bar(ind, model_vals_pre, color=modelcolor, width =.3, label = "model, pre")
			plt.bar(ind2, target_vals, color = targetcolor, width =.3,  label="target")#, yerr=target_ses)
			plt.bar(ind3, model_vals_post, color = modelcolorpost, width =.3,  label="model, post") 
			plt.legend(loc='best')
			plt.xticks(ind2, poplabels)
			plt.xlabel('pop')
			plt.ylabel('ave. per-site heterozygosity')
			plt.savefig(tosavefilename)
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/prepost/" + savestring + "_" + stat + "_prepost.png"
			plt.close()

		elif stat == 'sfs':
			xlabels = ['sing','<0.1','<0.2','<0.3', '<0.4','<0.5']
			ind = [x + .3 for x in np.arange(len(xlabels))]
			ind2 = [x + .6 for x in np.arange(len(xlabels))]
			ind3 = [x + .9 for x in np.arange(len(xlabels))]


			fig = plt.figure()
			pylab.subplots_adjust(hspace=0.250)
			#plt.ylabel('Fraction of SNPs')
			for i in range(len(target_vals)):
				ax = fig.add_subplot(len(target_vals), 1, i+1, ylabel = poplabels[i], xticks = ind2, xticklabels = xlabels)
				ax.bar(ind, model_vals_pre[i], label="model, pre", width=.3, color=modelcolor, yerr=model_ses_pre[i])
				ax.bar(ind2, target_vals[i],  label="target", width=.3, color=targetcolor, yerr=target_ses[i])
				ax.bar(ind3, model_vals_post[i], label="model, post", width=.3, color=modelcolorpost, yerr=model_ses_post[i])
				ax.set_ylim([0,.5])
				#ax.legend()
			plt.xlabel('Minor allele frequency')

			#plt.legend(loc='best')

			#plt.suptitle("Site frequency spectra, " + plot_savestring)
			plt.suptitle("Site frequency spectra, target vs. model calcs, " + str(numReps) + " reps")
			plt.savefig(tosavefilename)
			plt.close()
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/prepost/" + savestring + "_" + stat + "_prepost.png"

		elif stat == 'anc':
			xlabels = ['sing','<0.1','<0.2','<0.3', '<0.4','<0.5']
			ind = [x + .3 for x in np.arange(len(xlabels))]
			ind2 = [x + .6 for x in np.arange(len(xlabels))]
			ind3 = [x + .9 for x in np.arange(len(xlabels))]
			fig = plt.figure()
			#plt.title(item[0] + ", JV vs. IS calcs, " + str(numReps) + " reps")
			pylab.subplots_adjust(hspace=0.250)
			for i in range(len(target_vals)):
				ax = fig.add_subplot(len(target_vals), 1, i+1, ylabel = poplabels[i], xticks = ind2, xticklabels = xlabels)
				ax.bar(ind, model_vals_pre[i], label="model, pre", width=.3, color=modelcolor, yerr=model_ses_pre[i])
				ax.bar(ind2, target_vals[i],  label="target", width=.3, color=targetcolor, yerr=target_ses[i])
				ax.bar(ind3, model_vals_post[i], label="model, post", width=.3, color=modelcolorpost, yerr=model_ses_post[i])
				ax.set_ylim([0,.5])
				#ax.legend()
			#plt.legend(loc='best')
			plt.xlabel('Allele frequency')
			plt.ylabel('P(ancestral|frequency')
			plt.suptitle("Percentrage ancestral by frequency, target vs. model calcs, " + str(numReps) + " reps")
			plt.savefig(tosavefilename)
			plt.close()
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/prepost/" + savestring + "_" + stat + "_prepost.png"

		elif stat == 'r2':
			plt.figure()
			plt.title("LD, r2 vs phys dist, target vs. model calcs, " + str(numReps) + " reps")
			xlabels = ['<5kb','<10kb','<15kb','<20kb','<25kb','<30kb','<35kb','<40kb','<45kb','<50kb','<55kb','<60kb','<65kb','<70kb']
			ind = np.arange(len(xlabels))
			fig = plt.figure()
			pylab.subplots_adjust(hspace=0.250)
			for i in range(len(target_vals)):
				ax = fig.add_subplot(len(target_vals), 1, i+1, ylabel = poplabels[i], xticks = ind, xticklabels = xlabels, ylim=[0,.4])
				ax.errorbar(ind, model_vals_pre[i],  label="model, pre", color=modelcolor, yerr=model_ses_pre[i])
				ax.errorbar(ind, target_vals[i], label="target", color=targetcolor, yerr=target_ses[i])
				ax.errorbar(ind, model_vals_post[i],  label="model, post", color=modelcolorpost, yerr=model_ses_post[i])

			plt.suptitle("r2 vs physical distance, target vs. model calcs, " + str(numReps) + " reps")
			plt.savefig(tosavefilename)
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/prepost/" + savestring + "_" + stat + "_prepost.png"
			#plt.legend(loc = 'best')
			plt.close()

		elif stat == 'dprime':
			plt.figure()
			plt.title("LD, P(D'=1) vs gen dist, target vs. model calcs, " + str(numReps) + " reps")
			xlabels = ['<.001cM','<.002cM','<.003cM','<.004cM','<.005cM','<.006cM','<.007cM','<.008cM','<.009cM','<.010cM','<.011cM','<.012cM','<.013cM','<.014cM','<.015cM','<.016cM','<.017cM']
			ind = np.arange(len(xlabels))
			fig = plt.figure()
			pylab.subplots_adjust(hspace=0.250)
			for i in range(len(target_vals)):
				ax = fig.add_subplot(len(target_vals), 1, i+1, ylabel = poplabels[i], xticks = ind, ylim=[0,1]) #xticklabels = xlabels,
				ax.set_xticklabels(xlabels, rotation=20, fontsize =5)
				ax.errorbar(ind, model_vals_pre[i],  label="model, pre", color=modelcolor, yerr=model_ses_pre[i])
				ax.errorbar(ind, target_vals[i], label="target", color=targetcolor, yerr=target_ses[i])
				ax.errorbar(ind, model_vals_post[i],  label="model, post", color=modelcolorpost, yerr=model_ses_post[i])
			plt.legend(loc='best')
			plt.suptitle("P(D'=1) vs genetic distance, target vs. model calcs, " + str(numReps) + " reps")
			plt.savefig(tosavefilename)
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/prepost/" + savestring + "_" + stat + "_prepost.png"
			plt.close()

		elif stat == 'fst':
			target_vals = flattenList(target_vals)
			model_vals_pre = flattenList(model_vals_pre)
			model_vals_post = flattenList(model_vals_post)
			target_ses = flattenList(target_ses)
			model_ses_post = flattenList(model_ses_post)
			model_ses_pre = flattenList(model_ses_pre)


			plt.figure()
			plt.title("Fst, target vs. model calcs, " + str(numReps) + " reps")
			ind = [x + .3 for x in np.arange(len(popPairs))]
			ind2 = [x + .6 for x in np.arange(len(popPairs))]
			ind3 = [x + .9 for x in np.arange(len(popPairs))]
			plt.bar(ind, model_vals_pre, color = modelcolor, width =.3, yerr=model_ses_pre, label="model, pre")
			plt.bar(ind2, target_vals, color = targetcolor, width =.3,  label="target", yerr=target_ses)
			plt.bar(ind3, model_vals_post, color = modelcolorpost, width =.3, yerr=model_ses_post, label="model, post")
			plt.legend(loc='best')
			plt.xticks(ind2, popPairs)
			plt.xlabel('pop pair')
			plt.ylabel('Fst')
			plt.savefig(tosavefilename)
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/prepost/" + savestring + "_" + stat + "_prepost.png"
			plt.close()
def plot_comparison(statfilename, savestring, numReps):
	stats = ['pi', 'sfs', 'anc', 'fst', 'r2', 'dprime']
	pops = [1, 2, 3, 4]
	poplabels = ['YRI', 'CEU', 'CHB', 'BEB']
	popPairs = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
	popPairLabels = ['(YRI CEU)', '(YRI CHB)', '(YRI BEB)', '(CEU CHB)', '(CEU BEB)', '(CHB BEB)']
	tosavedir = "/home/unix/vitti/public_html/dem_model/NRE_110915/points"

	targetcolor, modelcolor = 'purple', 'orange'
	if "/" in statfilename:
		entries = statfilename.split('/')
		scenariobase = entries[-1]
 	else:
 		scenariobase = statfilename
 	scenariobase = scenariobase.replace('.stat', '')

	print "Plotting figures at: http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/points/"
	print "Saving files as: " + scenariobase
	import numpy as np

	targetstats = getTargetValues()
	modelstats = readIScustomstatfile(statfilename, 4)

	for stat in stats:
		tosavefilename = tosavedir + "/"  + savestring + "_" + stat + ".png"
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
			for pair in popPairs:
				item = (stat, pair)
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

		if stat == "pi":
			#variance vs standard err
			#num pops
			target_vals = flattenList(target_vals)
			model_vals = flattenList(model_vals)
			target_ses = flattenList(target_ses)
			model_ses = flattenList(model_ses)

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
			plt.savefig(tosavefilename)
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/points/" + scenariobase + "_" + stat + ".png"
			plt.close()

		elif stat == 'sfs':
			xlabels = ['sing','<0.1','<0.2','<0.3', '<0.4','<0.5']
			ind = [x + .6 for x in np.arange(len(xlabels))]
			ind2 = [x + 1 for x in np.arange(len(xlabels))]
			fig = plt.figure()
			pylab.subplots_adjust(hspace=0.250)
			#plt.ylabel('Fraction of SNPs')
			for i in range(len(target_vals)):
				ax = fig.add_subplot(len(target_vals), 1, i+1, ylabel = poplabels[i], xticks = ind2, xticklabels = xlabels)
				ax.bar(ind, target_vals[i],  label="target", width=.4, color=targetcolor, yerr=target_ses[i])
				ax.bar(ind2, model_vals[i], label="model", width=.4, color=modelcolor, yerr=model_ses[i])
				ax.set_ylim([0,.5])
				#ax.legend()
			plt.xlabel('Minor allele frequency')

			plt.legend(loc='best')

			#plt.suptitle("Site frequency spectra, " + plot_savestring)
			plt.suptitle("Site frequency spectra, target vs. model calcs, " + str(numReps) + " reps")
			plt.savefig(tosavefilename)
			plt.close()
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/points/" + scenariobase + "_" + stat + ".png"

		elif stat == 'anc':
			xlabels = ['sing','<0.1','<0.2','<0.3', '<0.4','<0.5']
			ind = [x + .6 for x in np.arange(len(xlabels))]
			ind2 = [x + 1 for x in np.arange(len(xlabels))]
			fig = plt.figure()
			#plt.title(item[0] + ", JV vs. IS calcs, " + str(numReps) + " reps")
			pylab.subplots_adjust(hspace=0.250)
			for i in range(len(target_vals)):
				ax = fig.add_subplot(len(target_vals), 1, i+1, ylabel = poplabels[i], xticks = ind2, xticklabels = xlabels)
				ax.bar(ind, target_vals[i],  label="target", width=.4, color=targetcolor, yerr=target_ses[i])
				ax.bar(ind2, model_vals[i], label="model", width=.4, color=modelcolor, yerr=model_ses[i])
				ax.set_ylim([0,.5])
				#ax.legend()
			plt.legend(loc='best')
			plt.xlabel('Allele frequency')
			#plt.ylabel('P(ancestral|frequency')
			plt.suptitle("Percentrage ancestral by frequency, target vs. model calcs, " + str(numReps) + " reps")
			plt.savefig(tosavefilename)
			plt.close()
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/points/" + scenariobase + "_" + stat + ".png"

		elif stat == 'r2':
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
			plt.savefig(tosavefilename)
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/points/" + scenariobase + "_" + stat + ".png"
			plt.legend(loc = 'best')
			plt.close()

		elif stat == 'dprime':
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
			plt.savefig(tosavefilename)
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/points/" + scenariobase + "_" + stat + ".png"
			plt.close()

		elif stat == 'fst':
			target_vals = flattenList(target_vals)
			model_vals = flattenList(model_vals)
			target_ses = flattenList(target_ses)
			model_ses = flattenList(model_ses)

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
			plt.savefig(tosavefilename)
			print "Saved file to http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/points/" + scenariobase + "_" + stat + ".png"
			plt.close()
	return

#######################################
## VISUALIZE SCIPY.OPTIMIZE PIPELINE ##
#######################################

def plotPipeline(errors, errorDict, pointParamss, keys, indices, savestring, descentmethod, stepsize):
	tosavedir = "/home/unix/vitti/public_html/dem_model/NRE_110915/pipelines/" 
	pubdir = "http://www.broadinstitute.org/~vitti/dem_model/NRE_110915/pipelines/"

	if True:
	#if len(errors) < 10: 
	#	pass
	#else:
		numParams = len(pointParamss[0])

		savename = tosavedir + savestring + "_pipeline.pdf"
		f, axarr = plt.subplots(3, sharex=True)
		f.suptitle('scipy.optimize.minimize output for pipeline, ' + savestring + ",\nstepsize: " + str(stepsize)[:5] + "; method: " + descentmethod)
		ax1, ax2, ax3 = axarr[0], axarr[1], axarr[2]

		#subplot 1: error
		ax1.plot(errors, color ='red')
		ax1.set_ylabel('error function', fontsize = 8)
		
		#subplot 2: error components
		errorKeys = errorDict.keys()
		for i in range(len(errorKeys)):
			key = errorKeys[i]
			stat, pop = key[0], key[1]
			ax2.plot(errorDict[key], label = key)
		ax2.set_ylabel('error components', fontsize = 8)
		#ax2.legend(prop={'size':6}, loc='best') 
		paramTrajectories = x = [[] for i in range(len(pointParamss[0]))]

		#print paramTrajectories
		#print pointParamss
		#subplot 3: parameter values
		for i in range(numParams):
			for j in range(len(pointParamss[i])):
				paramTrajectories[j].append(pointParamss[i][j])

		params = itertools.izip(keys, indices)
		for i in range(len(paramTrajectories)):
			key, index = keys[i], indices[i]
			plotlabel = str(key) + " [" + str(index) + "]" #can make this more informative
			ax3.plot(paramTrajectories[i], label = plotlabel)

		#ax3.legend(prop={'size':6}, loc='best') 
		ax3.set_ylabel('parameter trajectories', fontsize = 8)
		ax3.set_xlabel('sequential points in parameter space')

		f.subplots_adjust(hspace=.250)
		f.savefig(savename)
		plt.close()
		print "plotted to " + pubdir + savestring + "_pipeline.pdf"
	return

