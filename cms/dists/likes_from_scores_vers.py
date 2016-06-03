## {{DEM MODEL --> COMPONENT SCORE PROBABILITY DISTS}} -- script iii.
## last updated 05.31.16 	vitti@broadinstitute.org

#########################
### CONFIGURE THIS RUN ##
#########################

scriptName = "dryrun_053116" 
neutScenarios = ['default_default_101715_12pm']#, 'nulldefault', 'nulldefault_constantsize','gradient_101915_treebase_6_best', 'default_112115_825am', 'DEBUG_default_default_101715_12pm']
numBins, selPops = 50, [1, 2, 3, 4]
verbose = True
scores = ['ihs']#, 'deldaf', 'delihh', 'xp', 'fst',]#'xp']
scriptName = "debug_default_default"
plotFig, transparency = True, .5
freqGroups = ['all']#['lo', 'mid', 'hi']

##########
## MAIN ##
##########

from func_clean import *
import os

def main():

	########################
	### MAIN -- LOAD DATA ##
	########################
	data = {}
	if True:

		# DEFINE BINS
		fullrange, binlen = [.05, .95], .1 #define selfreq bins
		numbins = int((fullrange[1] - fullrange[0])/binlen)
		bin_starts = [fullrange[0] + binlen * i for i in range(numbins+1)]
		bin_ends = [fullrange[0] + binlen * i for i in range(1, numbins+2)]
		bin_medians = [fullrange[0] + (binlen/2) + i*.1 for i in range(numbins+1)]
		bin_medians_str = [str(x)[:3] for x in bin_medians]

		if verbose:
			print "LOADING DATA..\n"
		for iScore in range(len(scores)):
			################################
			## LOOP OVER COMPONENT SCORES ##
			################################
			score = scores[iScore]
			if score in ['xp', 'fst', 'deldaf']: #per pop pair
				comparison  = True
			else:
				comparison = False
			if score in ['xp']: #Toggle which scoretypes need to have a header struck
				stripHeader = True
			else:
				stripHeader = False
			if verbose:
				print "score: " + str(score) 

			for iScenario in range(len(neutScenarios)):	
				##################################
				## LOOP OVER DEMOGRAPHIC MODELS ##
				##################################
				dem_scenario = neutScenarios[iScenario]
				if verbose:
					print "\tmodel: " + dem_scenario 

				for iPop in range(len(selPops)):
					################################
					## LOOP OVER PUTATIVE SELPOPS ##
					################################
					selpop = selPops[iPop]
					if verbose:
						print "\t\tselpop: " + str(selpop)

					if comparison == False:
						altpops = ['']

					else:
						altpops = selPops[:]
						altpops.remove(selpop)
						#print altpops

					for altpop in altpops:

						#############################
						## GET DATA FROM NEUT SIMS ##
						#############################
						neutral_filename = getNeutralFilename(score, dem_scenario, selpop, altpop)
						if verbose:
							print "\t\tloading from: " + str(neutral_filename)
						assert os.path.isfile(neutral_filename)
						numCols, takeindices = getIndices(score, "neut")
						val_array = loadValsFromFile(neutral_filename, numCols, takeindices, stripHeader, verbose)	
						assert len(val_array) == 2
						neut_positions, neut_score_final = val_array[0], val_array[1]

						#neut_positions, neut_freq_anc, neut_score_unnormed, neut_score_final = val_array[0], val_array[1], val_array[2], val_array[3]
						#for debugitem in [neut_positions, neut_anc_freqs, neut_values, neut_freq_anq]:
						#		print str(len(debugitem))
						#		print debugitem[0:5]

						##############################
						## GET DATA FROM SWEPT SIMS ##
						##############################
						for ibin in range(len(bin_medians_str)):
							selfreq = bin_medians_str[ibin]
							lowerBound, upperBound = bin_starts[ibin], bin_ends[ibin]
							selscenario = "sel" + str(selpop) + "_" + str(selfreq)
							sel_filename = getSelFilename(score, dem_scenario, selpop, selfreq, altpop)
							if verbose:
								print "\t\tloading from: " + str(sel_filename)

							assert os.path.isfile(sel_filename)
							numCols, takeindices = getIndices(score, selscenario)
							val_array = loadValsFromFile(sel_filename, numCols, takeindices, stripHeader, verbose)	
							assert len(val_array) == 2
							sel_positions, sel_score_final = val_array[0], val_array[1]				

							#sel_positions, sel_freq_anc, sel_score_unnormed, sel_score_final = val_array[0], val_array[1], val_array[2], val_array[3]
							#for debugitem in [sel_positions, sel_freq_anc, sel_ihs_unnormed, sel_ihs_normed]:
							#	print str(len(debugitem))
							#	print debugitem[0:5]

							#### PARSE LINKED VS CAUSAL VARS
							causal_indices = [i for i, x in enumerate(sel_positions) if x == 500000]
							linked_indices = [i for i, x in enumerate(sel_positions) if x != 500000] 

							causal_positions = [sel_positions[variant] for variant in causal_indices]
							#causal_freq_anc = [sel_freq_anc[variant] for variant in causal_indices]
							#causal_score_unnormed = [sel_score_unnormed[variant] for variant in causal_indices]
							causal_score_final = [sel_score_final[variant] for variant in causal_indices]

							linked_positions = [sel_positions[variant] for variant in linked_indices]
							#linked_freq_anc = [sel_freq_anc[variant] for variant in linked_indices]
							#linked_score_unnormed = [sel_score_unnormed[variant] for variant in linked_indices]
							linked_score_final = [sel_score_final[variant] for variant in linked_indices]

							datatypes = []
							for status in ['causal', 'linked', 'neut']:
								for dpoint in ['positions', 'score_final']: #'score_unnormed', freq_anc $?!
									datatypes.append(status + "_" + dpoint)
						
							for datatype in datatypes:
								key = (score, dem_scenario, selpop, altpop, datatype, selfreq)
								data[key] = eval(datatype)

	####################
	### MAIN -- PLOT  ##
	####################
	if True:
		if verbose:
			print "PLOTTING..\n"
		for iScore in range(len(scores)):

			for freqGroup in freqGroups:
				if freqGroup == "lo":
					takebins = [0, 1, 2] #indices of bin_medians_str
				elif freqGroup == "mid":
					takebins = [3, 4, 5]
				elif freqGroup == 'hi':
					takebins = [6, 7, 8]
				else:
					takebins = [0, 1, 2, 3, 4, 5, 6, 7, 8]

				score = scores[iScore]
				if score not in ['xp', 'fst', 'deldaf']: #per pop scores - plot 4x5; pop by model
					f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7 ,ax8), (ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16), (ax17, ax18,ax19, ax20)) = plt.subplots(5, 4, sharex='col', sharey='row')
					f.suptitle("p( component score | demographic model + putative selpop )\n" + score + " " + freqGroup)
					iAxis = 1
					#titlestring = score + "_" + "likes"
					#plt.xlabel('putative selpop') 
					#plt.ylabel('demographic model')
					for iScenario in range(len(neutScenarios)):	
						dem_scenario = neutScenarios[iScenario]

						for iPop in range(len(selPops)):
							selpop = selPops[iPop]

							bins, scorerange, ylims = getBins(score, numBins)
							causalScoresToPlot, linkedScoresToPlot, neutScoresToPlot = [], [], []
							

							for ibin in takebins:#range(len(bin_medians_str)): #0,1,2, 3, 4, 5]:
								selfreq = bin_medians_str[ibin]
								causalScoresToPlot.extend(data[(score, dem_scenario, selpop, altpop, 'causal_score_final', selfreq)])
								linkedScoresToPlot.extend(data[(score, dem_scenario, selpop, altpop, 'linked_score_final', selfreq)])
								neutScoresToPlot.extend(data[(score, dem_scenario, selpop, altpop, 'neut_score_final', selfreq)])
							ax = eval('ax' + str(iAxis))
							titlestring = score + "_" + dem_scenario + "_" + str(selpop) + "_likes_" + freqGroup
							plotFig(ax, bins, causalScoresToPlot, linkedScoresToPlot, neutScoresToPlot, dem_scenario, str(selpop), scorerange, ylims, transparency, thinToSize=False, show=False, titlestring = titlestring)
							iAxis +=1
					if plotFig:
						plt.show()
				
				else: #per pop pair -- plot 6x5; poppair by model #NO
				# per putative selpop!
					for iPop in range(len(selPops)):
						selpop = selPops[iPop]
						altpops = selPops[:]
						altpops.remove(selpop)

						f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12), (ax13, ax14,ax15)) = plt.subplots(5, 3, sharex='col', sharey='row')
						f.suptitle("p( component score | demographic model + putative selpop/outgroup )\n" + score + " " + str(selpop))
						iAxis = 1

						for iScenario in range(len(neutScenarios)):	
							dem_scenario = neutScenarios[iScenario]
							bins, scorerange, ylims = getBins(score, numBins)
							for altpop in altpops:

								causalScoresToPlot, linkedScoresToPlot, neutScoresToPlot = [], [], []
								for ibin in takebins:#range(len(bin_medians_str)): #0,1,2, 3, 4, 5]:
									selfreq = bin_medians_str[ibin]
									causalScoresToPlot.extend(data[(score, dem_scenario, selpop, altpop, 'causal_score_final', selfreq)])
									linkedScoresToPlot.extend(data[(score, dem_scenario, selpop, altpop, 'linked_score_final', selfreq)])
									neutScoresToPlot.extend(data[(score, dem_scenario, selpop, altpop, 'neut_score_final', selfreq)])
								ax = eval('ax' + str(iAxis))
								titlestring = score + "_" + dem_scenario + "_" + str(selpop) + "_" +str(altpop) + "_likes_" + freqGroup
								plotFig(ax, bins, causalScoresToPlot, linkedScoresToPlot, neutScoresToPlot, dem_scenario, str(selpop) + "-" + str(altpop), scorerange, ylims, transparency, thinToSize=False, show=False, titlestring = titlestring)
								iAxis +=1
					if plotFig:
						plt.show()
	
main()