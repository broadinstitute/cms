import numpy, logging
from sys import exit
from Classes.DotData import DotData
from Operations.Shari_Operations.localize.xpopMerge import xpopMerge
from Operations.Shari_Operations.localize.Scenario import GetSelectionScenarios, GetScenarios
from Operations.MiscUtil import MakeAlphaNum, Dict, Sfx, progress, AddFileSfx
from functools import reduce

def mergeSims( scenario, Ddata = '../Data/Shari_Data/sim/', simsOut = 'simsOut3', nreplicas = 5,
	       thinExt = '.thin', thinSfx = '',
	       selpop = None, getio = None ):
	"""Gathers per-SNP information, for all replicas of a given scenario, and outputs it in a single DotData where each line
	gives info for one SNP.

	Specifically, reads simulation and Sweep output, collects columns needed for composite likehood test (chrom, base pair position, genetic
	distance, anc frequencies for 3 populations, xpop for each pair, and ihs, iHH_A and iHH_D for selected population)

	Input params:

	   scenario - an object of class Scenario, indicating the simulation scenario (either neutral or a selection scenario)
	       from which all replicas were simulated.
	   nreplicas - the number of replicas simulated under this scenario.
	      Each replica represents a chromosome region, with a set of SNPs on it.
	   
	   Ddata - the directory under which the simulations and the Sweep analysis results live.
	     Under this directory we expect to find:
	         iHS analysis results, under power_ihs/
		 XP-EHH analysis results, under power_xpop
		 simulation output giving SNP positions

	   thinExt - the extension appended to simulation files that describe the SNPs in the simulated replica.
	      Sometimes we create simulations and then thin them under different thinning models (to simulate SNP ascertainment
	      by the various stages of HapMap; these differently thinned versions of the same simulations might be stored in
	      simulation files with different extensions.

	   thinSfx - the suffix appended to the power_ihs and power_xpop directory names, telling where to find iHS and XP-EHH
	      analyses of the simulations.   When we analyze the same simulations after applying different thinning scenarios,
	      the iHS and XP-EHH analyses for each thinning scenario go into a separate set of directories.

        Output params:

	    Ddata - under Ddata writes a DotData named merged_scenName.data, where each line gives info
	        for one SNP, with the following columns (type of data is float unless stated otherwise):

	        CHROM_POS 1 - physical (basepair) position of the SNP within its replica.
	           Note that one merged file contains SNPs for a set of replicas (all for the same scenario),
		   so there could be multiple SNPs with the same position.  The replica number
		   is given in the Chrom column.
		FREQ1 1 - derived allele frequency in pop 1 ( European )
		FREQ1 4 - derived allele frequency in pop 4 ( EastAsian )
		FREQ1 5 - derived allele frequency in pop 5 ( WestAfrican )

		R AllEHH logratio Deviation European_WestAfrican - XP-EHH score to the right of the SNP,
		   between European and WestAfrican pops, normalized to the neutral background.
		   Analogously for the next five columns:
		L AllEHH logratio Deviation European_WestAfrican
		R AllEHH logratio Deviation EastAsian_European
		L AllEHH logratio Deviation EastAsian_European
		R AllEHH logratio Deviation EastAsian_WestAfrican
		L AllEHH logratio Deviation EastAsian_WestAfrican

		SNP pos (cM) European_WestAfrican - genetic map position of this SNP, within its replica.
		   (the European_WestAfrican suffix is irrelevant).
		SNP pos (bases) European_WestAfrican - physical (basepair) position of this SNP within its replica.
		   (the European_WestAfrican suffix is irrelevant).
		Chrom European_WestAfrican - the replica from which this SNP comes; can be nan.
		   (the European_WestAfrican suffix is irrelevant)
		Chrom - the replica from which this SNP comes; can be nan
		SNP pos (bases) - physical (basepair) position of this SNP within its replica.
		SNP pos (cM) - genetic map position of this SNP within its replica
		Both iHH_A - sum of iHH_A for both directions from this SNP
		Both iHH_D - sum of iHH_D for both directions from this SNP
		Both iHS - the value in 'Both Unstandardised iHS' (below), but binned by derived allele frequency
		   and normalized within the bin.
		Left iHH_D - iHH_D to the left of the SNP (the raw integral value).  analogously for the next three.
		Right iHH_D
		Left iHH_A
		Right iHH_A
		Both Unstandardised iHS - log( (iHH_A_left + iHH_A_right) / ( iHH_D_left + iHH_D_right ) )
		   ( see also 'Both iHS' column for the standardized iHS score )
	
	"""

	assert selpop == None  or  scenario.is_neutral()

	DataDir = Ddata + '/'
	SimDir = DataDir + simsOut + thinSfx + '/'

	if not scenario.is_neutral():
		scenName = 'sel%d_%d' % ( scenario.mutFreq, scenario.mutPop )
		scenDir  = str( scenario.mutAge ) + 'ky/' + scenName
	else:
		scenName = 'neutral'
		scenDir = 'neutral'
	
	popName = {1:'European',4:'EastAsian',5:'WestAfrican'}
	
	ihsSignifTsv = DataDir + 'power_ihs' + thinSfx + '/' + scenDir + '/ihs_sig_' + \
	    popName[ scenario.mutPop if not scenario.is_neutral() else ( selpop if selpop != None else 1 ) ] + '.tsv'
	xpopSignifTsv = [ DataDir + 'power_xpop' + thinSfx + '/' + scenDir + '/xpop_significance_' + popPair + '.tsv'
			  for popPair in ( 'EastAsian_WestAfrican', 'EastAsian_European', 'European_WestAfrican' ) ]
	posFiles = [ SimDir + scenDir + '/' + str(ichrom) + '_' + scenName + '.pos-%d%s' % ( pop, thinExt )
		     for ichrom in range( nreplicas ) for pop in ( 1, 4, 5 ) ]
	
	ageSfx = '%dky' % ( scenario.mutAge if not scenario.isNeutral() else 10 )
	mergedDotData  = AddFileSfx( Ddata + 'merged.data/', ageSfx, scenario.scenName(), selpop, thinSfx )

	fileDescrs = \
	{ mergedDotData :
		  ( 'Various per-snp statistics for SNPs in scenario $scenario, replicas 0-$nreplicas.',
		    ( ( 'CHROM_POS 1', 'physical (basepair) position of the SNP within its replica. '
			'Note that one merged file contains SNPs for a set of replicas (all for the same scenario), '
			'so there could be multiple SNPs with the same position.  The replica number '
			'is given in the Chrom column. ' ), 
		      ( 'FREQ1 1', 'derived allele frequency in pop 1 ( European )' ),
		      ( 'R AllEHH logratio Deviation European_WestAfrican', 'XP-EHH score to the R of the SNP, '
			'between European and WestAfrican pops, normalized to the neutral background.' ),
		      ( 'SNP pos (cM) European_WestAfrican', 'genetic map SNP position' ),
		      ( 'SNP pos (bases) European_WestAfrican', 'physical SNP position' ),
		      ( 'Chrom European_WestAfrican', 'chromosome (or replica number)' ),
		      ( 'Chrom', 'chromosome (or replica number)' ),
		      ( 'SNP pos (bases)', 'physical SNP position' ),
		      ( 'SNP pos (cM)', 'genetic map SNP position' ),
		      ( 'Both iHH_A', 'sum of iHH_A scores for both sides' ),
		      ( 'Both iHH_D', 'sum of iHH_D scores for both sides' ),
		      ( 'Both iHS', 'sum of iHS scores for both sides' ),
		      ( ' Left iHH_D', 'iHH_D score to the left of the SNP' ),
		      ( 'Right iHH_D', 'iHH_D score to the right of the SNP' ),
		      ( 'Left iHH_A', 'iHH_A score to the left of the SNP' ),
		      ( 'Right iHH_A', 'iHH_A score to the right of the SNP' ), 
		      ( 'Both Unstandardised iHS', 'sum of unstandardized iHS scores for both sides' ) ) ) }

	if getio: return dict( depends_on = posFiles + [ ihsSignifTsv ] + xpopSignifTsv, creates = mergedDotData,
			       mediumRuleNameSfx = scenario.scenDir(),
			       fileDescrs = fileDescrs )

	ncausal = 0

	dashFixer = lambda v: v if v != '-' else numpy.nan

	# Load iHS of selected pop
	ihsAll = DotData(SVPath = ihsSignifTsv,ToLoad=['Chrom','SNP pos (bases)','SNP pos (cM)','Both iHH_A','Both iHH_D','Both iHS','Left iHH_D','Right iHH_D','Left iHH_A','Right iHH_A','Both Unstandardised iHS'], SVValueFixer = dashFixer)
	ihsAllChrom = ihsAll.Chrom
	
	# Load xpop values
	xpopAll = xpopMerge( *xpopSignifTsv )
	logging.info( 'done with xpopMerge' )
	
	xpopAll = xpopAll[['R AllEHH logratio Deviation European_WestAfrican','L AllEHH logratio Deviation European_WestAfrican','R AllEHH logratio Deviation EastAsian_European','L AllEHH logratio Deviation EastAsian_European','R AllEHH logratio Deviation EastAsian_WestAfrican',
			   'L AllEHH logratio Deviation EastAsian_WestAfrican','SNP pos (cM) European_WestAfrican','SNP pos (bases) European_WestAfrican','Chrom European_WestAfrican']]
	xpopAllChrom = xpopAll['Chrom European_WestAfrican']
	
	replicates = []

	xpopIdx = 0
	ihsIdx = 0
	
	for ichrom in range(nreplicas):

		progress( 'Merging replicas', ichrom, nreplicas, freq = 1 )
	
		logging.info( 'looking at replica %d of %d' % ( ichrom, nreplicas ) )
		# Load in pos files for this replica.
		# They give, for each SNP in the replica, its physical (basepair) position within the replica,
		# and the frequency of the derived and the ancestral alleles.
		pos1, pos4, pos5 = [ DotData(SVPath=SimDir + scenDir + '/' + str(ichrom) + '_' + scenName + '.pos-%d%s' % ( pop, thinExt),
					     SVSkipFirstLines = 1, SVHeader = False,
					     names = ['SNP','CHROM', 'CHROM_POS', 'ALLELE1', 'FREQ1', 'ALLELE2', 'FREQ2' ]) for pop in ( 1, 4, 5 ) ]
		assert pos1.numCols() == pos4.numCols() == pos5.numCols()
		posBlank = ((numpy.nan,)*pos1.numCols(),)*3
		logging.info( 'Loaded pos files for chrom ' + str( ichrom ) +  ': ' + str( len(pos1) ) + 'snps' )

		assert set(pos1.CHROM_POS) == set(pos4.CHROM_POS) == set(pos5.CHROM_POS)

		logging.info( 'pos file sizes are: %d, %d, %d' % ( len( pos1 ), len( pos4 ), len( pos5 ) ) )
		logging.info( 'Merging on position...' )
		posAll = DotData.mergeOnKeyCols((pos1,pos4,pos5),('CHROM_POS',)*3,posBlank, suffixes = (' 1',' 4',' 5'))

		logging.info( 'Done merging.' )
		logging.info( 'type(posAll) is ' + str( type( posAll ) ) )
		print(len(posAll))
		chrom = numpy.ones(len(posAll))*ichrom
		newChrom = DotData(Columns = [chrom,],names=['newChrom',])
		print(newChrom)
		posAll = posAll[['CHROM_POS 1','FREQ1 1','FREQ1 4','FREQ1 5']]
		posAll.hstack(newChrom)

		logging.info( 'added replica number column' )
		
		print(posAll)
		posAllBlank = (numpy.nan,)*posAll.numCols()
		
		# 10-16-08 ADDED CHROM TO MERGED OUTPT  ( not now used -- can be removed? )

		#
		# From the xpop and ihs significance results, get just the rows for SNPs in the 
		# current replica
		#

		#while xpopIdx < len( xpopAllChrom ) and xpopAllChrom[ xpopIdx ] == ichrom: xpopIdx += 1
		#xpop = xpopAll[ :xpopIdx ]
		xpop = xpopAll[ xpopAllChrom == ichrom ]
		logging.info( 'selected xpop for replica %d' % ichrom )
		xpopBlank = (numpy.nan,)*xpop.numCols()

		#while ihsIdx < len( ihsAllChrom ) and ihsAllChrom[ ihsIdx ] == ichrom: ihsIdx += 1
		#ihs = ihsAll[ :ihsIdx ]
		ihs = ihsAll[ ihsAllChrom == ichrom ]
		logging.info( 'selected ihs for replica %d' % ichrom )
		ihsBlank = (numpy.nan,)*ihs.numCols()

#		if not set( ihs[  'SNP pos (bases)' ] ).issubset( set( posAll['CHROM_POS 1'] ) ):
#			print 'bad positions: ', set( posAll['CHROM_POS 1'] ) - set( ihs[  'SNP pos (bases)' ] )
#		assert set( ihs[  'SNP pos (bases)' ] ).issubset( set( posAll['CHROM_POS 1'] ) ), "bad iHS file " + ihsSignifTsv

		logging.info( 'merging replica %d' % ichrom )
		Data = DotData.mergeOnKeyCols((posAll,xpop,ihs),('CHROM_POS 1','SNP pos (bases) European_WestAfrican','SNP pos (bases)'),
					      blanks = (posAllBlank,xpopBlank,ihsBlank), suffixes = ('pos',' xpop',' ihs'),
					      verbose = True )
		logging.info( 'done merging replica %d; now have %d records' % ( ichrom, len( Data ) ) )
		
		Data = Data[ numpy.invert( numpy.isnan( Data[ 'CHROM_POS 1' ] ) ) ]
		logging.info( 'done removing snp info for SNPs not in all .pos files for replica %d; now have %d records'
			      % ( ichrom, len( Data ) ) )
		
		replicates.append(Data)

		logging.info( 'now have ' + str( len( replicates ) ) + ' replicates.' )

	# endloop: for each replica

	logging.info( 'Stacking replicates...' )
	allData = reduce( lambda x, y: x.vstack(y), replicates)
	logging.info( 'Saving merged SNP info to ' + mergedDotData )
	allData.save( mergedDotData )
	
	logging.info( 'Finished mergeSims()' )

#        print scen + ' ncausal: ' + str(ncausal)

def DefineRulesTo_MergeSims( pr, mutAges, mutPops, mutFreqs, noNeutral, nreplicas,
                             Ddata, simsOut, thinExt = '.thin', thinSfx = '' ):
    """Pipeline generator: for each scenario, create a rule to merge SNP info for all SNPs in each replica within that scenario,
    into a single table.
    """

    for scenario in ( GetSelectionScenarios if noNeutral else GetScenarios)( mutAges, mutPops, mutFreqs ):
        print('generating rule for scenario ', scenario)
        pr.addInvokeRule( invokeFn = mergeSims,
                          invokeArgs = Dict( 'scenario nreplicas Ddata simsOut thinExt thinSfx' ) )


