import logging, itertools
import numpy as np
import os
from Operations.Shari_Operations.localize.Scenario import GetSelectionScenarios, GetScenarios
from Operations.MiscUtil import MakeAlphaNum, Dict, Sfx, progress, dbg, AddFileSfx
from Operations.Shari_Operations.localize.PopConsts import pn_European, pn_EastAsian, pn_WestAfrican, pop2name, \
    AllAges, AllPops, AllFreqs
from Operations.IDotData import IDotData
from Operations.Ilya_Operations.PipeRun.python.PipeRun import GetCreates


def mergePosFilesOneSim( scenario, Ddata, replicaNum, putativeMutPop = None, simsOut = 'simsOut', thinExt = '', thinSfx = '',
			 statsSfx = '', ihsSfx = '', pop2name = pop2name,
			 getio = None ):
	"""Merge per-SNP position files for one sim.
	"""

	if not Ddata.endswith('/'): Ddata += '/'

	popNums = sorted( pop2name.keys() )
	minPopNum = popNums[ 0 ]

	SimDir = os.path.join( Ddata, simsOut + thinSfx )

	scenName = scenario.scenName()
	scenDir = scenario.scenDir()

	if putativeMutPop == None: putativeMutPop = scenario.mutPop
	
	posFiles = [ os.path.join( SimDir, scenDir, '%d_%s.pos-%d%s' % ( replicaNum, scenName, pop, thinExt ) )
		      for pop in popNums ]

	snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
	mergedPosFN = os.path.join( snpStatsDir, 'mergedPos', AddFileSfx( 'mergedPos.tsv', statsSfx, putativeMutPop, ihsSfx, replicaNum ) )

	if getio: return dict( depends_on = posFiles, creates = mergedPosFN,
			       splitByCols = dict([ ( posFile, dict( keyCols = 'CHROM_POS',
								     tableReadOpts = dict( headingSep = 'whitespace' ) ) ) for posFile in posFiles ]),
			       mediumRuleNameSfx = scenario.scenDir() )

	replicaIDotDatas = []

	posCols = ['CHROM_POS %d' % minPopNum ] + ['FREQ1 %d' % popNum for popNum in popNums]

	# Load in pos files for this replica.
	# They give, for each SNP in the replica, its physical (basepair) position within the replica,
	# and the frequency of the derived and the ancestral alleles.
	posIndiv = [ IDotData( os.path.join( SimDir, scenDir, str(replicaNum) + '_' + scenName +
					     '.pos-%d%s' % ( pop, thinExt) ), skipFirstLines = 1,
			       headings = ('SNP','CHROM', 'CHROM_POS', 'ALLELE1',
					   'FREQ1', 'ALLELE2', 'FREQ2' )) 
		     for pop in popNums ]

	posAll = IDotData.merge( iDotDatas = posIndiv,
				 cols = ('CHROM_POS',) * len( posIndiv ),
				 blanks = (np.nan,) * len( posIndiv ),
				 suffixes = [ ' %s' % pop for pop in popNums ] )

	posAll = posAll[posCols]

	IDotData.repeat( heading = 'replicaNum', value = replicaNum ).hstack( posAll ).save( mergedPosFN )


def mergeSims( scenario, Ddata, posFileFN = None, simsOut = 'simsOut', nreplicas = 100, thinExt = '', thinSfx = '',
	       putativeMutPop = None, outFile = None,
	       pop2name = pop2name, statsSfx = '', ihsSfx = '',
               limitToPop = None,
	       getio = None ):
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


	   putativeMutPop - the population in which, we think, selection is occurring, aka "putatively selected population".
	      In practice, when localizing a given region, we will usually suspect that selection has occurred in a particular
	      population.   When doing a genome-wide scan, we can do several scans assuming each population in turn to be
	      the selected population, and find regions selected in that population.

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

	if not Ddata.endswith('/'): Ddata += '/'

	assert nreplicas > 0
	dbg( 'pop2name' )

	SimDir = os.path.join( Ddata, simsOut + thinSfx )

	scenName = scenario.scenName()
	scenDir = scenario.scenDir()

	if putativeMutPop == None: putativeMutPop = scenario.mutPop
	
	ihsSignifFN = os.path.join( Ddata, 'power_ihs' + thinSfx, scenDir,
				    'ihs_sig_' + pop2name[ putativeMutPop ] + ihsSfx + '.tsv' )

	popNames = sorted( pop2name.values() )
	popNums = sorted( pop2name.keys() )
	minPopNum = popNums[ 0 ]

	posFileKeyCols = ( 'replicaNum', 'CHROM_POS %d' % minPopNum )
	xpopIhsKeyCols = ('Chrom', 'SNP pos (bases)')
	
	popPairs = [ '%s_%s' % ( popNames[ pop1idx ], popNames[ pop2idx ] )
		     for pop1idx in range( len( popNames ) ) for pop2idx in range( pop1idx+1, len( popNames ) )
                     if limitToPop is None or limitToPop in ( popNames[ pop1idx ], popNames[ pop2idx ] )  ]
	
	xpopSignifFNs = [ os.path.join( Ddata, 'power_xpop' + thinSfx, scenDir, 'xpop_significance_' + popPair + '.tsv' )
			  for popPair in popPairs ]

	snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    
	mergedData = outFile if outFile else os.path.join( snpStatsDir, AddFileSfx( 'merged.tsv', statsSfx, putativeMutPop, ihsSfx ) )

	fileDescrs = \
	{ mergedData :
		  ( 'Various per-snp statistics for SNPs in scenario $scenario, replicas 0-$nreplicas, '
		    'assuming selection in ' + pop2name[ putativeMutPop ],
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

	if posFileFN is None: posFileFN = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir(),
							AddFileSfx( 'mergedPosStacked.tsv', statsSfx, putativeMutPop, ihsSfx ) )
	
	if getio: return dict( depends_on = [ posFileFN, ihsSignifFN ] + xpopSignifFNs, creates = mergedData,
			       splitByCols = dict([ ( posFileFN, dict( keyCols = posFileKeyCols ) ) ]
						  + [ ( signifFN, dict( keyCols = xpopIhsKeyCols ) ) for signifFN in [ ihsSignifFN ] + xpopSignifFNs ] ),
			       mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ),
			       fileDescrs = fileDescrs,
                               attrs = Dict( 'putativeMutPop nreplicas pop2name' ) )

	dashFixer = lambda v: v if v != '-' else np.nan

	ihsAll = IDotData(ihsSignifFN, valueFixer = dashFixer)
	ihsAll = ihsAll[('Chrom','SNP pos (bases)','SNP pos (cM)','Both iHH_A','Both iHH_D','Both iHS',
			 'Left iHH_D','Right iHH_D','Left iHH_A','Right iHH_A','Both Unstandardised iHS')]
	def chkReplica( r, n = nreplicas ): return r.Chrom < n 
		
	ihsAll = ihsAll.takewhile( chkReplica )
	
	xpopCols = ('Chrom','SNP pos (bases)','SNP pos (cM)','L AllEHH logratio Deviation','R AllEHH logratio Deviation')

	xpopSignif = tuple( [ IDotData(xpopSignifFN, valueFixer = dashFixer)[ xpopCols].takewhile( chkReplica )
			      for xpopSignifFN in xpopSignifFNs ] )
	
	posCols = ['CHROM_POS %d' % minPopNum ] + ['FREQ1 %d' % popNum for popNum in popNums]

	result = IDotData.merge( iDotDatas =  ( IDotData( posFileFN ), ) + xpopSignif + ( ihsAll, ),
				 cols = (posFileKeyCols,) +
					 (xpopIhsKeyCols,) * ( len( popPairs ) + 1 ),
				 blanks = (None,) + (np.nan,) * ( len( popPairs ) + 1 ),
				 suffixes = ['pos'] + [ ' %s' % popPair for popPair in popPairs ] + [ '' ] )

	aPopPair = 'European_WestAfrican' if 'European_WestAfrican' in popPairs else popPairs[0]
	useCols = [ 'replicaNum' ] + posCols + \
	    [ '%s AllEHH logratio Deviation %s' % ( side, popPair ) for popPair in popPairs for side in ( 'L', 'R' ) ] + \
	    [ 'SNP pos (cM) ' + aPopPair,
	      'SNP pos (bases) ' + aPopPair,
	      'Chrom ' + aPopPair,
	      'Chrom',
	      'SNP pos (bases)',
	      'SNP pos (cM)',
	      'Both iHH_A',
	      'Both iHH_D',
	      'Both iHS' ]

	if len( popPairs ) == 1:
		result = result.renameCols( { 'L AllEHH logratio Deviation' : 'L AllEHH logratio Deviation ' + aPopPair,
					      'R AllEHH logratio Deviation' : 'R AllEHH logratio Deviation ' + aPopPair } )
					      
	result[ useCols ].save( mergedData )
	
	logging.info( 'Finished mergeSims()' )


def DefineRulesTo_MergeSims( pr, Ddata, simsOut = 'simsOut', thinSfx = '', thinExt = '',
			     nreplicas = 100, mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs,
			     pop2name = pop2name, limitToPop = None ):
	"""Define rules to merge per-SNP data for each SNP into one file"""

	assert False
	
	if not Ddata.endswith('/'): Ddata += '/'
    
	for scenario in GetScenarios( mutAges, mutPops, mutFreqs ):
		pr.addInvokeRule( invokeFn = mergePosFilesOneSim,
				  invokeArgs = Dict( 'Ddata simsOut thinSfx thinExt '
						     'scenario nreplicas pop2name' ) )

		for putativeMutPop in ( mutPops if scenario.isNeutral() else (scenario.mutPop,) ):

			       
			pr.addInvokeRule( invokeFn = mergeSims,
					  invokeArgs = Dict( 'Ddata simsOut thinSfx thinExt '
							     'scenario nreplicas putativeMutPop pop2name limitToPop' ) )
            
