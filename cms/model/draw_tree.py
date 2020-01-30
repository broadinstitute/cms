## {{POP GEN DATA --> DEM MODEL}}
## this script contains functions to facilitate visualization of demographic models 
## currently configured to fit 1KG P3 African populations: {PEL, ESN, CHB, CEU, YRI}
## last updated: 06.04.16 	vitti@broadinstitute.org

colorDict = {'YRI':'#CC9933', 'CEU':'#FFB900', 'CHB':'#E1B919', 'PEL':'#FFB933'}
colorDict = {'YRI':'gold', 'CEU':'dodgerblue', 'CHB':'springgreen', 'PEL':'mediumorchid'}

#from params import generateParams
#from params_func import generate_params
import matplotlib as mp 
mp.use('agg')
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import sys

####################
## HELPER METHODS ##
####################

def calculateExpGrowthRate(ne0, ne1, t0, t1):
	ne0, ne1, t0, t1 = float(ne0), float(ne1), float(t0), float(t1)
	return (np.log(ne1/ne0))/(np.absolute(t1-t0))
def drawRectangle(centerx, starty, width, length, color, rectlabel, limits, printvals):
	Path = mpath.Path
	startx = centerx - (width/2)
	endx = centerx + (width/2)
	endy = starty + length
	path_data = [
	    (Path.MOVETO, (startx, starty)), 
	    (Path.LINETO, (endx, starty)), #x axis is Ne
	    (Path.LINETO, (endx, endy)), ##y axis is T
	    (Path.LINETO, (startx, endy)),
	    (Path.CLOSEPOLY, (startx, starty)),
    ]
	codes, verts = zip(*path_data)
	path = mpath.Path(verts, codes)
	patch = mpatches.PathPatch(path, facecolor=color, alpha=1, label=rectlabel)
	if printvals == True:
		writePopSize(width, centerx + 15000, starty + (length/5), limits)
	return patch 
def drawTrapezoid(centerx, starty, startwidth, length, endwidth, color, limits, rate, printvals):
	Path = mpath.Path
	startx1 = centerx - (startwidth/2)
	startx2 = centerx + (startwidth/2)
	endx1 = centerx - (endwidth/2)
	endx2 = centerx + (endwidth/2)
	endy = starty + length
	p1 = (startx1, starty)
	p2 = (endx1, endy)
	p3 = (endx2, endy)
	p4 = (startx2, starty)
	path_data = [
		(Path.MOVETO, p1),
		(Path.LINETO, p2),
		(Path.LINETO, p3),
		(Path.LINETO, p4),
		(Path.CLOSEPOLY, p1),
	]
	codes, verts = zip(*path_data)
	path = mpath.Path(verts, codes)
	patch = mpatches.PathPatch(path, facecolor=color, alpha=1, label="")
	if printvals == True:
		writePopGrowth(endwidth, startwidth, centerx-15000, starty+20, limits, rate)
	return patch
def drawCoalescent(coaltime, originpopaxis, divergentpopaxis, startpopsize, travellen, color):
	Path = mpath.Path
	p1 = (originpopaxis, coaltime)
	#p2 = (originpopaxis + startpopsize, coaltime)
	p3 = (divergentpopaxis + (startpopsize/2), coaltime-travellen)
	p4 = (divergentpopaxis - (startpopsize/2), coaltime-travellen)
	path_data = [
		(Path.MOVETO, p1),
		#(Path.LINETO, p2),
		(Path.LINETO, p3),
		(Path.LINETO, p4),
		(Path.CLOSEPOLY, p1),
	]
	codes, verts = zip(*path_data)
	path = mpath.Path(verts, codes)
	patch = mpatches.PathPatch(path, facecolor=color, alpha=1, label="")
	return patch
def drawBottleneck(bntime, popaxis, bnsize, limits, width=10000, height=500, printvals=False):
	Path = mpath.Path
	center = (popaxis, bntime)
	p1 = (popaxis + width, bntime + height)
	p2 = (popaxis + width, bntime - height)
	p3 = (popaxis - width, bntime + height)
	p4 = (popaxis - width, bntime - height)
	path_data = [
		(Path.MOVETO, center),
		(Path.MOVETO, p1),
		(Path.MOVETO, p2),
		(Path.MOVETO, center),
		(Path.MOVETO, p3),
		(Path.MOVETO, p4),
		(Path.CLOSEPOLY, center),
	]
	codes, verts = zip(*path_data)
	path = mpath.Path(verts, codes)
	patch = mpatches.PathPatch(path, facecolor = 'red', alpha=.5, label="")
	if printvals == True:
		writebnSize(bnsize, popaxis, bntime+50, limits)
	return patch
def drawAdmixture(startx, endx, rate, time, color, height=50, printvals = False):
	Path = mpath.Path
	starty = time
	endy = time + (rate*10)
	p1 = (startx, starty)
	p2 = (endx, starty)
	p3 = (endx, endy)
	p4 = (startx, endy)
	path_data = [
		(Path.MOVETO, p1),
		(Path.LINETO, p2),
		(Path.LINETO, p3),
		(Path.LINETO, p4),
		(Path.CLOSEPOLY, p1),
	]
	codes, verts = zip(*path_data)
	path = mpath.Path(verts, codes)
	patch = mpatches.PathPatch(path, facecolor=color, alpha=rate, label="")
	if printvals == True:
		plt.text(startx + 10000, starty + 20, "admix: " + str(rate)[:8] + " > > ", fontsize=8)
	return patch
def writebnSize(bn, x, y, limits):
	if y > limits[0] and y < limits[1]:
		plt.text(x, y, "bn: " + str(bn)[:8], fontsize =8, horizontalalignment='center',
	        verticalalignment='center')
def writePopGrowth(startsize, finishsize, x, y, limits, rate):
	if y > limits[0] and y < limits[1]:
	#	plt.text(x, y, "Ne: " + str(int(startsize)) + "->" + str(int(finishsize)),fontsize =8, horizontalalignment='center',
       # verticalalignment='center')
		plt.text(x, y, "Ne: " + str(int(startsize)) + "->" + str(int(finishsize)) +"\nrate: " + str(rate)[:6],fontsize =6, horizontalalignment='center',
        verticalalignment='center')
def writePopSize(Ne, x, y, limits):
	if y > limits[0] and y < limits[1]:
		plt.text(x, y, "Ne:" + str(int(Ne)), fontsize =8, horizontalalignment='center',
    	    verticalalignment='center')

###############################
#### DRAW TREE, MIGRATIION ####
###############################

def draw_tree(paramDict, fullsaveloc, recent = False, veryrecent = False, printvals = True):
	##################################
	### Set up figure and get data ###
	##################################
	
	nbuffer = 5000 #how much space to leave between branches
	basepadding = 2000 #how much of the tree trunk to show
	leaf_widths = [paramDict[('Ne', i)][0] for i in range(1, 5)]
	#paramDict['presentSizes']
	YRI_axis = 0 + leaf_widths[0]/2 
	CEU_axis = (nbuffer + leaf_widths[0]) + (leaf_widths[1]/2)
	PEL_axis = (2*nbuffer + leaf_widths[0] + leaf_widths[1]) + (leaf_widths[3]/2)
	CHB_axis = (3*nbuffer + leaf_widths[0] + leaf_widths[1] + leaf_widths[3]) + (leaf_widths[2]/2)
	xlims = [0 - nbuffer, CHB_axis + leaf_widths[2]/2 + nbuffer]

	ax2 = plt.subplot2grid((5,10), (0, 0), colspan = 10, rowspan =10)

	YRI_Ne = paramDict[('Ne', 1)] 
	YRI_Ne_times = paramDict[('Ne_times', 1)]
	YRI_Ne_change = paramDict[('Ne_change', 1)]
	CEU_Ne = paramDict[('Ne', 2)]
	CEU_Ne_times = paramDict[('Ne_times', 2)]
	CEU_Ne_change = paramDict[('Ne_change', 2)]	
	CHB_Ne = paramDict[('Ne', 3)]
	CHB_Ne_times = paramDict[('Ne_times', 3)]
	CHB_Ne_change = paramDict[('Ne_change', 3)]	
	PEL_Ne = paramDict[('Ne', 4)] 
	PEL_Ne_times = paramDict[('Ne_times', 4)]
	PEL_Ne_change = paramDict[('Ne_change', 4)]

	YRI_bn = paramDict[('bn', 1)] 
	YRI_bn_times = paramDict[('bn_times', 1)]
	CEU_bn = paramDict[('bn', 2)]
	CEU_bn_times = paramDict[('bn_times', 2)]
	CHB_bn = paramDict[('bn', 3)]
	CHB_bn_times = paramDict[('bn_times', 3)]
	PEL_bn = paramDict[('bn', 4)] 
	PEL_bn_times = paramDict[('bn_times', 4)]

	#configure what text gets printed depending on what's in the frame		
	if veryrecent == True:
		limits = [0, 300]
	elif recent == True:
		limits = [0, 5000]
	else:
		limits = [500, 50000]	

	###################
	### Draw leaves ###
	###################

	if paramDict[('Ne_change', 1)][0] != 'exp':
		rect1 = drawRectangle(YRI_axis, 0, leaf_widths[0], YRI_Ne_times[0], colorDict['YRI'], 'YRI', limits, False)
		ax2.add_patch(rect1)
	if paramDict[('Ne_change', 2)][0] != 'exp':	
		rect2 = drawRectangle(CEU_axis, 0, leaf_widths[1], CEU_Ne_times[0], colorDict['CEU'], 'CEU', limits, False)
		ax2.add_patch(rect2)
	if paramDict[('Ne_change', 3)][0] != 'exp':		
		rect3 = drawRectangle(CHB_axis, 0, leaf_widths[2], CHB_Ne_times[0], colorDict['CHB'], 'CHB', limits, False)
		ax2.add_patch(rect3)
	if paramDict[('Ne_change', 4)][0] != 'exp':	
		rect4 = drawRectangle(PEL_axis, 0, leaf_widths[3], PEL_Ne_times[0], colorDict['PEL'], 'PEL', limits, False)
		ax2.add_patch(rect4)

	#######################
	### Draw branchings ###
	#######################

	coalpops = ['YRI', 'CEU', 'CHB', 'PEL']
	travellendict = {'CEU':1, 'CHB':1, 'PEL':1}
	for i in range(1, len(coalpops)):
		coaltime = paramDict[('split', i + 1)][0]
		divpop = coalpops[i]
		parentpop = coalpops[i-1]
		coal = drawCoalescent(coaltime, eval(parentpop + "_axis"), eval(divpop + "_axis"), eval(divpop + "_Ne")[-1], travellendict[divpop], colorDict[divpop])
		ax2.add_patch(coal)

	for i in range(1, len(coalpops)):
		coaltime = paramDict[('split', i + 1)][0]
		divpop = coalpops[i]
		parentpop = coalpops[i-1]
		rect = drawRectangle(eval(divpop + "_axis"), eval(divpop + "_Ne_times")[-1], eval(divpop + "_Ne")[-1],  coaltime - travellendict[divpop] - eval(divpop + "_Ne_times")[-1], colorDict[divpop], '', limits, printvals)
		ax2.add_patch(rect)

	##########################
	### Draw branch widths ###
	##########################

	for pop in ['YRI', 'CEU', 'PEL', 'CHB']:
		Nesizes = eval(pop + "_Ne")
		Netimes = eval(pop + "_Ne_times")
		Nechange = eval(pop + "_Ne_change")
		lastgrowthtype = ""
		if pop == 'YRI':
			splittime = 10000000
		if pop == 'CEU':
			splittime = paramDict[('split', 2)][0]
		if pop == 'CHB':
			splittime = paramDict[('split', 3)][0]
		if pop == 'PEL':
			splittime = paramDict[('split', 4)][0]
		

		for i in range(len(Nesizes)-1):
			popsize = Nesizes[i]
			changetype = Nechange[i]
			if Netimes[i] <= splittime:
				if i != len(Nesizes)-1:
					branchlen = int(Netimes[i+1]) - int(Netimes[i])
				else:
					branchlen = basepadding #arbitrary
				if changetype == "normal" or int(changetype) == 0 or changetype == "0":
					rect = drawRectangle(eval(pop + "_axis"), Netimes[i] , popsize, branchlen, colorDict[pop], '', limits, printvals)
					ax2.add_patch(rect)
				if changetype == "exp" or int(changetype) == 1 or changetype == "1":
					pastsize, pasttime = Nesizes[i+1], Netimes[i+1]
					rate = calculateExpGrowthRate(pastsize, popsize, pasttime, Netimes[i])
					trapezoid = drawTrapezoid(eval(pop + "_axis"), Netimes[i], popsize, (pasttime - Netimes[i]), pastsize, colorDict[pop], limits, rate, printvals)
					ax2.add_patch(trapezoid)
			lastgrowthtype = changetype

		if pop == 'YRI':
			i = len(Nesizes)-1
			popsize = Nesizes[i]
			changetype = Nechange[i]
			rect = drawRectangle(eval(pop + "_axis"), Netimes[i] , popsize, branchlen, colorDict[pop], '', limits, printvals)
			ax2.add_patch(rect)

	########################
	### Draw bottlenecks ###
	########################

	for pop in ['YRI', 'CEU', 'PEL', 'CHB']:
		bnsizes = eval(pop + "_bn")
		bntimes = eval(pop + "_bn_times")
		if pop == 'YRI':
			splittime = 10000000
		if pop == 'CEU':
			splittime = paramDict[('split', 2)][0]
		if pop == 'CHB':
			splittime = paramDict[('split', 3)][0]
		if pop == 'PEL':
			splittime = paramDict[('split', 4)][0]
		for i in range(len(bnsizes)):
			bntime = bntimes[i]
			bnsize = bnsizes[i]
			popaxis = eval(pop + "_axis")
			if bntime < splittime: 
				bn = drawBottleneck(bntime, popaxis, bnsize, limits, printvals)
				ax2.add_patch(bn)
			#I'm not sure why this isn't working? but it suffices to have the text for now

	##########################################
	### Draw tree base and configure zoom ###
	##########################################

	ancrect = drawRectangle(YRI_axis,YRI_Ne_times[-1], YRI_Ne[-1], 100000, colorDict['YRI'], 'anc', limits, printvals)
	ancexprect = drawRectangle(YRI_axis, YRI_Ne_times[-2], YRI_Ne[-2], YRI_Ne_times[-1] - YRI_Ne_times[-2], colorDict['YRI'], '', limits, printvals)
	PEL_split = paramDict[('split', 4)][0]
	CEU_split = paramDict[('split', 2)][0]
	if veryrecent == True:
		ylims = [0,PEL_split+100]
	elif recent == True:
		ylims = [0,CEU_split+500]
	else:
		ylims = [0, paramDict['Ne_times', 1][-1] + basepadding]
	genlength = 25
	yearlims = [i*genlength/1000 for i in ylims]
	ax3 = ax2.twinx()
	ax3.set_ylabel("time bp (kya; gen = " + str(genlength) + "y)")
	ax3.set_ylim(yearlims)
	ax2.spines['top'].set_visible(False)
	ax2.xaxis.tick_bottom()
	ax2.grid()
	ax2.set_ylim(ylims)
	ax2.set_xlim(xlims)
	
	########################
	### Finalize and save ##
	########################
	ax2.set_title("Demographic model")
	plt.xticks([YRI_axis, CEU_axis, PEL_axis, CHB_axis], ['YRI', 'CEU', 'PEL', 'CHB'])
	ax2.set_xlabel('Width corresponds to effective population size, Ne')
	ax2.set_ylabel('time bp (gens)')
	plt.savefig(fullsaveloc)
	savename = fullsaveloc.strip('/home/unix/vitti/public_html/')
	ploturl = "http://www.broadinstitute.org/~vitti/" + savename
	plt.close()
	print("plotted: " + ploturl)
	return
def draw_tree_ax(paramDict, ax2, recent = False, veryrecent = False, printvals = True):
	##################################
	### Set up figure and get data ###
	##################################
	nbuffer = 5000 #how much space to leave between branches
	basepadding = 2000 #how much of the tree trunk to show
	leaf_widths = [paramDict[('Ne', i)][0] for i in range(1, 5)]
	#paramDict['presentSizes']
	YRI_axis = 0 + leaf_widths[0]/2 
	CEU_axis = (nbuffer + leaf_widths[0]) + (leaf_widths[1]/2)
	PEL_axis = (2*nbuffer + leaf_widths[0] + leaf_widths[1]) + (leaf_widths[3]/2)
	CHB_axis = (3*nbuffer + leaf_widths[0] + leaf_widths[1] + leaf_widths[3]) + (leaf_widths[2]/2)
	xlims = [0 - nbuffer, CHB_axis + leaf_widths[2]/2 + nbuffer]

	#ax2 = plt.subplot2grid((5,10), (0, 0), colspan = 10, rowspan =10)

	YRI_Ne = paramDict[('Ne', 1)] 
	YRI_Ne_times = paramDict[('Ne_times', 1)]
	YRI_Ne_change = paramDict[('Ne_change', 1)]
	CEU_Ne = paramDict[('Ne', 2)]
	CEU_Ne_times = paramDict[('Ne_times', 2)]
	CEU_Ne_change = paramDict[('Ne_change', 2)]	
	CHB_Ne = paramDict[('Ne', 3)]
	CHB_Ne_times = paramDict[('Ne_times', 3)]
	CHB_Ne_change = paramDict[('Ne_change', 3)]	
	PEL_Ne = paramDict[('Ne', 4)] 
	PEL_Ne_times = paramDict[('Ne_times', 4)]
	PEL_Ne_change = paramDict[('Ne_change', 4)]

	YRI_bn = paramDict[('bn', 1)] 
	YRI_bn_times = paramDict[('bn_times', 1)]
	CEU_bn = paramDict[('bn', 2)]
	CEU_bn_times = paramDict[('bn_times', 2)]
	CHB_bn = paramDict[('bn', 3)]
	CHB_bn_times = paramDict[('bn_times', 3)]
	PEL_bn = paramDict[('bn', 4)] 
	PEL_bn_times = paramDict[('bn_times', 4)]

	#configure what text gets printed depending on what's in the frame		
	if veryrecent == True:
		limits = [0, 300]
	elif recent == True:
		limits = [200, 5000]
	else:
		limits = [500, 50000]	

	###################
	### Draw leaves ###
	###################

	if paramDict[('Ne_change', 1)][0] != 'exp':
		rect1 = drawRectangle(YRI_axis, 0, leaf_widths[0], YRI_Ne_times[0], colorDict['YRI'], 'YRI', limits, printvals)
		ax2.add_patch(rect1)
	if paramDict[('Ne_change', 2)][0] != 'exp':	
		rect2 = drawRectangle(CEU_axis, 0, leaf_widths[1], CEU_Ne_times[0], colorDict['CEU'], 'CEU', limits, printvals)
		ax2.add_patch(rect2)
	if paramDict[('Ne_change', 3)][0] != 'exp':		
		rect3 = drawRectangle(CHB_axis, 0, leaf_widths[2], CHB_Ne_times[0], colorDict['CHB'], 'CHB', limits, printvals)
		ax2.add_patch(rect3)
	if paramDict[('Ne_change', 4)][0] != 'exp':	
		rect4 = drawRectangle(PEL_axis, 0, leaf_widths[3], PEL_Ne_times[0], colorDict['PEL'], 'PEL', limits, printvals)
		ax2.add_patch(rect4)

	#######################
	### Draw admixture  ###
	#######################

	admix_rate = paramDict['admix_rate'][0] 
	admix_time = paramDict['admix_time'][0]
	admix = drawAdmixture(CEU_axis, PEL_axis, admix_rate, admix_time, colorDict['PEL'], printvals) #hard-coded to assume CEU->ASI/PEL
	ax2.add_patch(admix)

	#######################
	### Draw branchings ###
	#######################

	coalpops = ['YRI', 'CEU', 'CHB', 'PEL']
	travellendict = {'CEU':75, 'CHB':75, 'PEL':25}
	for i in range(1, len(coalpops)):
		coaltime = paramDict[('split', i + 1)][0]
		divpop = coalpops[i]
		parentpop = coalpops[i-1]
		coal = drawCoalescent(coaltime, eval(parentpop + "_axis"), eval(divpop + "_axis"), eval(divpop + "_Ne")[-1], travellendict[divpop], colorDict[divpop])
		ax2.add_patch(coal)

	for i in range(1, len(coalpops)):
		coaltime = paramDict[('split', i + 1)][0]
		divpop = coalpops[i]
		parentpop = coalpops[i-1]
		rect = drawRectangle(eval(divpop + "_axis"), eval(divpop + "_Ne_times")[-1], eval(divpop + "_Ne")[-1],  coaltime - travellendict[divpop] - eval(divpop + "_Ne_times")[-1], colorDict[divpop], '', limits, printvals)
		ax2.add_patch(rect)

	##########################
	### Draw branch widths ###
	##########################

	for pop in ['YRI', 'CEU', 'PEL', 'CHB']:
		Nesizes = eval(pop + "_Ne")
		Netimes = eval(pop + "_Ne_times")
		Nechange = eval(pop + "_Ne_change")
		lastgrowthtype = ""
		if pop == 'YRI':
			splittime = 10000000
		if pop == 'CEU':
			splittime = paramDict[('split', 2)][0]
		if pop == 'CHB':
			splittime = paramDict[('split', 3)][0]
		if pop == 'PEL':
			splittime = paramDict[('split', 4)][0]
		

		for i in range(len(Nesizes)-1):
			popsize = Nesizes[i]
			changetype = Nechange[i]
			if Netimes[i] <= splittime:
				if i != len(Nesizes)-1:
					branchlen = int(Netimes[i+1]) - int(Netimes[i])
				else:
					branchlen = basepadding #arbitrary
				if changetype == "normal" or int(changetype) == 0 or changetype == "0":
					rect = drawRectangle(eval(pop + "_axis"), Netimes[i] , popsize, branchlen, colorDict[pop], '', limits, printvals)
					ax2.add_patch(rect)
				if changetype == "exp" or int(changetype) == 1 or changetype == "1":
					pastsize, pasttime = Nesizes[i+1], Netimes[i+1]
					rate = calculateExpGrowthRate(pastsize, popsize, pasttime, Netimes[i])
					trapezoid = drawTrapezoid(eval(pop + "_axis"), Netimes[i], popsize, (pasttime - Netimes[i]), pastsize, colorDict[pop], limits, rate, printvals)
					ax2.add_patch(trapezoid)
			lastgrowthtype = changetype

		if pop == 'YRI':
			i = len(Nesizes)-1
			popsize = Nesizes[i]
			changetype = Nechange[i]
			rect = drawRectangle(eval(pop + "_axis"), Netimes[i] , popsize, branchlen, colorDict[pop], '', limits, printvals)
			ax2.add_patch(rect)

	########################
	### Draw bottlenecks ###
	########################

	for pop in ['YRI', 'CEU', 'PEL', 'CHB']:
		bnsizes = eval(pop + "_bn")
		bntimes = eval(pop + "_bn_times")
		if pop == 'YRI':
			splittime = 10000000
		if pop == 'CEU':
			splittime = paramDict[('split', 2)][0]
		if pop == 'CHB':
			splittime = paramDict[('split', 3)][0]
		if pop == 'PEL':
			splittime = paramDict[('split', 4)][0]
		for i in range(len(bnsizes)):
			bntime = bntimes[i]
			bnsize = bnsizes[i]
			popaxis = eval(pop + "_axis")
			if bntime < splittime: 
				bn = drawBottleneck(bntime, popaxis, bnsize, limits, printvals)
				ax2.add_patch(bn)
			#I'm not sure why this isn't working? but it suffices to have the text for now

	##########################################
	### Draw tree base and configure zoom ###
	##########################################

	ancrect = drawRectangle(YRI_axis,YRI_Ne_times[-1], YRI_Ne[-1], 100000, colorDict['YRI'], 'anc', limits, printvals)
	ancexprect = drawRectangle(YRI_axis, YRI_Ne_times[-2], YRI_Ne[-2], YRI_Ne_times[-1] - YRI_Ne_times[-2], colorDict['YRI'], '', limits, printvals)
	PEL_split = paramDict[('split', 4)][0]
	CEU_split = paramDict[('split', 2)][0]
	if veryrecent == True:
		ylims = [0,PEL_split+100]
	elif recent == True:
		ylims = [0,CEU_split+500]
	else:
		ylims = [0, paramDict['Ne_times', 1][-1] + basepadding]
	genlength = 25
	yearlims = [i*genlength/1000 for i in ylims]
	#ax3 = ax2.twinx()
	#ax3.set_ylabel("time bp (kya; gen = " + str(genlength) + "y)")
	#ax3.set_ylim(yearlims)
	#ax2.spines['top'].set_visible(False)
	#ax2.xaxis.tick_bottom()
	ax2.grid()
	ax2.set_ylim(ylims)
	ax2.set_xlim(xlims)
	
	########################
	### Finalize and save ##
	########################
	ax2.set_title("Demographic model")
	#ax2.set_xlabel('Width corresponds to effective population size, Ne', fontsize=8)
	ax2.set_ylabel('time bp (gens)')
	ax2.set_xticks([YRI_axis, CEU_axis, PEL_axis, CHB_axis])
	ax2.set_xticklabels(['YRI', 'CEU', 'PEL', 'CHB'])

	plt.close()

	return
def draw_mig_matrix(paramDict, fullsaveloc):
	"""could visualize change over time as e.g. gif or sequence.
	currently: takes time-weighted average rate."""
	pops = range(1, paramDict['numPops'] + 1)
	poplabels = [paramDict['labels'][pop] for pop in pops]
	grid = []
	for sourcepop in pops:
		row = []
		for targetpop in pops:
			if sourcepop != targetpop:
				keystring = str(sourcepop) + "->" + str(targetpop)
				ratekey = ('mig', keystring)
				timekey = ('mig_time', keystring)

				#must find point where pops coalesce
				if targetpop != 1:
					coal_a = paramDict['split', targetpop][0]
				else:
					coal_a = 5000
				if sourcepop != 1:
					coal_b = paramDict['split', sourcepop][0]
				else:
					coal_b = 5000

				migstarttime = min(coal_a, coal_b)
				rates = paramDict[ratekey]
				times = paramDict[timekey]
				#catch the possibility that the array has an inefficacious mig time value that is more ancestral than coaleseence
				for i in range(len(times)):
					if times[i] > migstarttime:
						times[i] = migstarttime
				times.append(migstarttime)
				total_mig_gens = 0
				for i in range(len(rates)):
					rate = rates[i]
					if times[i+1] <= migstarttime:
						time = times[i+1] - times[i]
						total_mig_gens += (rate * time)

				#print keystring
				#print str(rates)
				#print str(times)
				averate = (total_mig_gens/migstarttime)
				row.append(averate)
			else:
				row.append(0)
		grid.append(row)

	plt.imshow(grid, interpolation='none')
	plt.title('Model migration rates\n(dynamic; time-weighted ave is pictured)')
	plt.ylabel('source pop')
	plt.xlabel('target pop')
	plt.colorbar()
	plt.xticks(range(4), poplabels)
	plt.yticks(range(4), poplabels)
	plt.savefig(fullsaveloc)
	savename = fullsaveloc.strip('/home/unix/vitti/public_html/')
	ploturl = "http://www.broadinstitute.org/~vitti/" + savename
	plt.close()
	print("plotted: " + ploturl)
	return
