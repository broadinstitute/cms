import numpy as np
import matplotlib as mp
mp.use('agg')
import matplotlib.pyplot as pp
import pandas as pd
import os, sys
#from Operations.MiscUtil import SlurpFileLines
from Operations.Shari_Operations.localize.clean_hapmap2_region import clean_hapmap2_region
#from Operations.Shari_Operations.localize.subs import normalize
from Classes.DotData import DotData

def normalize(rawVals, ind = itertools.repeat( True ) ):
    notNan = []
    for val, useVal in itertools.izip( rawVals, ind ):
        if useVal and not isnan(val) and not isneginf(val):
            notNan.append(val)
#    print len(notNan)
    theMean = mean(notNan)
    theStd = std(notNan)
    
    normVals = []

    for val, useVal in itertools.izip( rawVals, ind ):
        normVals.append( (val - theMean)/theStd if useVal else nan )

    return normVals, theMean, theStd

z = DotData( '/idi/sabeti-data/ilya/sabetilab/ilya/new/nsvn/Data/Shari_Data/Shari_Data/localizeHapMap/LCT/LCT_likesByStat.tsv' )
z = z[ z.gdPos < 151.6 ]
cln = clean_hapmap2_region( z )
print(cln)
pp.plot( cln.gdPos, normalize(cln.lik)[0], 'k.' )
ax = pp.axes()
yvals = ax.get_yticks()
ax.set_yticklabels(['{:3.1f}'.format(y) for y in yvals])
xvals = ax.get_xticks()
ax.set_xticklabels(['{:3.1f}'.format(x) for x in xvals])
pp.title( 'LCT')
pp.xlabel('pos (cM)' )
pp.ylabel('CMS')
pp.savefig( '/web/personal/jvitti/LCT.svg' )


