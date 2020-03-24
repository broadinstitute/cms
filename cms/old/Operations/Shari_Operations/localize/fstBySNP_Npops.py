
from Operations.MiscUtil import dbg
import numpy as np
from numpy import *

def fst(pop2name, sampleSize, pop2ancFreqs, pop2sampleSize = {}):
	"""Given one SNP's ancestral allele frequencies in three populations,
	compute estimated Fst for the three population pairs.

	Returns:

	   a map which maps the string keys 'EastAsian_European',
	   'European_WestAfrican' and 'EastAsian_WestAfrican' to the
	   Fst value of the corresponding pair, showing how much that SNP's
	   ancestral allele frequency differs between the two populations.
	   
	"""

	pops = sorted( pop2name.values() )

       	n = {}
	Fst = {}
	f1 = dict( ( pop2name[ popNum ], ancFreqs ) for popNum, ancFreqs in list(pop2ancFreqs.items()) )
	f2 = {}
	
	nanc = {}
	for pop in pops:
		nanc[pop] = zeros(len(f1[pop]))
		#print nanc[pop], f1[pop]
		f2[pop] = 1 - f1[pop]
		sampleSizeThisPop = pop2sampleSize[ pop ] if pop in pop2sampleSize else sampleSize
		for i in range(len(f1[pop])):
			nanc[pop][i] = int(f1[pop][i]*sampleSizeThisPop + 0.5)
		n[pop] = sampleSizeThisPop

	# Use Weir-Hill estimator for Fst
		
	for ipop in pops:
		for jpop in pops:
			if jpop + '_' + ipop in list(Fst.keys()) or ipop == jpop:
				continue
			Fst[ipop + '_' + jpop] = zeros(len(f1[ipop]))
			for i in range(len(f1[ipop])):			
				pmean = (nanc[ipop][i] + nanc[jpop][i]) / (n[ipop] + n[jpop])
			
				nic = n[ipop] - (n[ipop]*n[ipop])/(n[ipop]+n[jpop])
				njc = n[jpop] - (n[jpop]*n[jpop])/(n[ipop]+n[jpop])
				#assert nic == njc == 60
				nc = nic + njc
				#print nc, f1[ipop][i], f1[jpop][i]
				#print f2[ipop][i], f2[jpop][i]
				msp = n[ipop] * (f1[ipop][i] - pmean) * (f1[ipop][i] - pmean) \
				    + n[jpop] * (f1[jpop][i] - pmean) * (f1[jpop][i] - pmean)
				msg = ((n[ipop] * f1[ipop][i]* f2[ipop][i]) + (n[jpop] * f1[jpop][i] * f2[jpop][i])) \
				    / (n[ipop] - 1 +n[jpop] - 1)
				num = msp - msg
				denom = msp + (msg*(nc - 1))
				#print msp, msg
				#print num, denom
				an_fst = nan
				if denom != 0:
					an_fst = num / denom
				
				Fst[ipop + '_' + jpop][i] = an_fst
			
	return Fst


def fst_onePopPair(ancFreqs, sampleSizes):
        """Compute fst between two pops, for each SNP, given the ancestral freq in each pop and the sample sizes.
        """

        dbg( 'ancFreqs sampleSizes' )

        n = sampleSizes
        n_tot = n[0] + n[1]
#        nanc =  np.ceil([ ancFreqs[0] * sampleSizes[0], ancFreqs[1] * sampleSizes[1] ])
        nanc =  ( np.array([ ancFreqs[0] * sampleSizes[0], ancFreqs[1] * sampleSizes[1] ]) + .5 ).astype( int )

        dbg( 'nanc' )

        f1 = ancFreqs
        f2 = 1 - ancFreqs

        dbg( 'f1 f2 f1.shape f2.shape' )

        # Use Weir-Hill estimator for Fst

        pmean = (nanc[0] + nanc[1]) / n_tot

        nic = n[0] - (n[0]*n[0])/n_tot
        njc = n[1] - (n[1]*n[1])/n_tot
        dbg( 'pmean.shape nic.shape njc.shape nic njc' )
        #assert nic == njc == 60
        nc = nic + njc
        #print nc, f1[0], f1[1]
        #print f2[0], f2[1]
        msp = n[0] * (f1[0] - pmean) * (f1[0] - pmean) \
            + n[1] * (f1[1] - pmean) * (f1[1] - pmean)
        msg = ((n[0] * f1[0]* f2[0]) + (n[1] * f1[1] * f2[1])) \
            / (n[0] - 1 + n[1] - 1)
        dbg( 'msp.shape msg.shape' )
        num = msp - msg
        denom = msp + (msg*(nc - 1))
        #print msp, msg
        #print num, denom
        dbg( 'num.shape denom.shape' )
        denom [ denom == 0 ] = nan
        an_fst = num / denom

        dbg( 'an_fst.shape' )

        ipop = 0
        jpop = 1
        dbg( '"KKKKKKKKKKKKK" pmean nanc[ipop] nanc[jpop] n[ipop] n[jpop] nic njc nc msp msg num denom an_fst' )
        
        return an_fst

# end: def fst_onePopPair(ancFreqs, samplesSizes)
        

def fst_oneSnp(pop2name, sampleSize, pop2ancFreq, pop2sampleSize = {}):
	"""Given one SNP's ancestral allele frequencies in three populations,
	compute estimated Fst for the three population pairs.

	Returns:

	   a map which maps the string keys 'EastAsian_European',
	   'European_WestAfrican' and 'EastAsian_WestAfrican' to the
	   Fst value of the corresponding pair, showing how much that SNP's
	   ancestral allele frequency differs between the two populations.
	   
	"""

#        dbg( '"QQQQQQQQQQQ" pop2name sampleSize pop2ancFreq pop2sampleSize' )

	pops = sorted( pop2name.values() )

       	n = {}
	Fst = {}
	f1 = dict( ( pop2name[ popNum ], ancFreq ) for popNum, ancFreq in list(pop2ancFreq.items()) )
	f2 = {}
	
	nanc = {}
	for pop in pops:
		f2[pop] = 1 - f1[pop]
		sampleSizeThisPop = pop2sampleSize[ pop ] if pop in pop2sampleSize else sampleSize
		nanc[pop] = int(f1[pop]*sampleSizeThisPop + 0.5)
		n[pop] = sampleSizeThisPop
#                dbg( '"NNNNNNNNNNN" pop f1[pop] f2[pop] sampleSizeThisPop nanc[pop] n[pop]' )

	# Use Weir-Hill estimator for Fst
		
	for ipop in pops:
		for jpop in pops:
			if jpop + '_' + ipop in list(Fst.keys()) or ipop == jpop:
				continue
			Fst[ipop + '_' + jpop] = 0
			pmean = (nanc[ipop] + nanc[jpop]) / (n[ipop] + n[jpop])

			nic = n[ipop] - (n[ipop]*n[ipop])/(n[ipop]+n[jpop])
			njc = n[jpop] - (n[jpop]*n[jpop])/(n[ipop]+n[jpop])

			nc = nic + njc

			msp = n[ipop] * (f1[ipop] - pmean) * (f1[ipop] - pmean) \
			    + n[jpop] * (f1[jpop] - pmean) * (f1[jpop] - pmean)
			msg = ((n[ipop] * f1[ipop]* f2[ipop]) + (n[jpop] * f1[jpop] * f2[jpop])) \
			    / (n[ipop] - 1 +n[jpop] - 1)
			num = msp - msg
			denom = msp + (msg*(nc - 1))


			an_fst = ( num / denom ) if denom != 0 else nan
#                        dbg( '"MMMMMMMMMMM" ipop jpop pmean nanc[ipop] nanc[jpop] n[ipop] n[jpop] nic njc nc msp msg num denom an_fst' )
                        
			Fst[ipop + '_' + jpop] = an_fst
			
	return Fst

