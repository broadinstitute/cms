"""CMS-related code: computing CMS component statistics, creating likelihood tables, and computing CMS for given SNPs
using given likelihood tables.
"""


from __future__ import division, with_statement

__all__ = ( 'DefineRulesTo_computeCMSstats', 'DefineRulesTo_computeCMSforSims', 'DefineRulesTo_normalizeCMS', 'DefineRulesTo_createLikesTables' ) 

from Operations.IDotData import IDotData, imeanstd
from Operations.tsvutils import DefineRulesTo_normalizeColumnsWithinGroups
from numpy import array, where, isnan, nanmax, zeros, any, all, isfinite, argsort, exp, repeat, atleast_1d, mean, min, max
import numpy as np
from scipy import interpolate
from Classes.DotData import DotData
from Operations.DotDataFunctions.ondiskops import linkDotDatas
from Operations.DotDataFunctions.datavstack import datavstack
from Operations.MiscUtil import MakeAlphaNum, Sfx, AddFileSfx, dbg, Str, Histogrammer, Dict, DictGet, tabwrite, DictGet, coerceVal, \
    MakeSeq, compile_expr, MergeDicts, DictExcept, AttrDict, ReplaceFileExt, AddSlash, SumKeeper, StatKeeper, clamp, vstack_tsvs
from Operations.Shari_Operations.localize.subs import normalize, normalizeByFreq, CoM_dist, scoreVal, likeVal, likeRatio, \
    calcRank, likeValBin 
from Operations.Shari_Operations.localize.fstBySNP_Npops import fst, fst_oneSnp
from Operations.Shari_Operations.localize.PopConsts import popName, AllPops, AllFreqs, AllAges, pn_WestAfrican, pop2name, CAUSAL_POS
from Operations.Ilya_Operations.SnpStats import DefineRulesTo_histogramSnpStatistic, identifyReplicasMeetingConds, \
    splitSnpStatsFile, joinSnpStatsFiles, DefineRulesTo_histogramReplicaStatistic, scatterPlotReplicaStatistic, \
    findReplicasMatchingConds
from Operations.Shari_Operations.localize.Scenario import GetScenarios, GetSelectionScenarios, NeutralScenario, \
    SelectionScenario, Scenario
from Operations.Shari_Operations.localize.RunAndAnalyzeSims import DefineRulesTo_RunSimsAndSweep
from Operations.Shari_Operations.localize.mergeSimsP3 import mergeSims, mergePosFilesOneSim
from Operations.Shari_Operations.localize.clean_hapmap2_region import clean_hapmap2_region
from Operations.Shari_Operations.localize.localizeSpatially import localizeSpatiallyBySplineFitting, \
    localizeSpatiallyByWindows, evalSpatialLoc, DefineRulesTo_gatherCausalSnpGdPos
from Operations.Ilya_Operations.localize.gatherCausalRanks import DefineRulesTo_gatherCausalRanks
from Operations.Ilya_Operations.sweep import DefineRulesTo_extractGeneticMapsFromSims
from Operations.Ilya_Operations.SnpStats import DefineRulesTo_histogramReplicaStatistic, \
    DefineRulesTo_histogramSnpStatistic, DefineRulesTo_gatherCausalFreqs, GraphHistograms
from Operations.Ilya_Operations.PipeRun.python.PipeRun import GetCreates
from Operations.Ilya_Operations.RBTree import RBList
from Operations.tsvutils import computeSumsWithinGroups
from operator import concat, attrgetter, itemgetter, add
import logging, os, types, itertools, sys, math, operator
import matplotlib.pyplot as pp


# Bin boundaries for CMS statistics
class CMSBinsBase(object):
    #nbins = 60

    # current code has constraints that the number of bins must be the same for all stats
    #assert set( stat_nbin.values() ) == set((nbins,))  

    # field: CMSstats - string names of the CMS statistics.
    #   the order is relied upon by some code (e.g. in computeCMSstats) so don't change it.
    CMSstats = ( 'iHS', 'StdDiff', 'meanFst', 'freqDiff', 'max_xpop' )

    nonNormedStats = ( 'freqDiff', )

    #stat_numSpecialBins = { 'max_xpop': 0, 'iHS': 2, 'meanFst': 0, 'StdDiff': 2, 'dist': 0, 'freqDiff': 0, 'DAF': 0 }
    stat_numSpecialBins = { 'max_xpop': 0, 'iHS': 0, 'meanFst': 0, 'StdDiff': 0, 'dist': 0, 'freqDiff': 0, 'DAF': 0 }
    stat_specialBinNames = { 'iHS': ( 'DAF < .05', 'DAF > .95' ), 'StdDiff': ( 'DAF < .05', 'DAF > .95' ) }
    maxSpecialBins = max( stat_numSpecialBins.values() )

class CMSBins(CMSBinsBase):
    stat_start = { 'max_xpop': -3, 'iHS': -6, 'meanFst': -1.0, 'StdDiff': -3.0, 'dist': 0, 'freqDiff': -1, 'DAF': -1.0 }
    stat_end = { 'max_xpop': 8.0, 'iHS': 6, 'meanFst': 6.0, 'StdDiff': 5.0, 'dist': 2.0 , 'freqDiff': 1, 'DAF': 1.0 }
    stat_nbin = { 'max_xpop': 60, 'iHS': 60, 'meanFst': 60, 'StdDiff': 60, 'dist': 60, 'freqDiff': 60, 'DAF': 60 }

class CMSBinsLocal(CMSBinsBase):

    stat_start = { 'max_xpop': -3., 'iHS': -3.0, 'meanFst': -2.0, 'StdDiff': -3.0, 'dist': 0,'freqDiff': -1.0, 'DAF': -1.0 }
    stat_end = { 'max_xpop': 3.0, 'iHS': 3.0, 'meanFst': 2.0, 'StdDiff': 3.0, 'dist': .01,'freqDiff': 1.0, 'DAF': 1.0 }
    stat_nbin = { 'max_xpop': 60, 'iHS': 60, 'meanFst': 60, 'StdDiff': 60, 'dist': 60, 'freqDiff': 60, 'DAF': 60 }

    
def computeCMSstats_old( Ddata, thinSfx,
                         scenario, putativeMutPop = None, sampleSize = 120,
                         pop2name = pop2name,
                         pop2sampleSize = {},
                         oldMerged = False,
                         DdataMerged = None,
                         outFile = None,
                         stats = CMSBins.CMSstats,
                         limitToPop = None,
                         getio = None ):
    """Compute CMS stats for all SNPs in all replicas of a given scenario.

    Params:

       Ddata - the root directory containing simulations and their analyses
       thinSfx - specifies which thinning version of the simulations are used.
       scen - the simulation scenario.   Here we compute stats for all replicas within that scenario, for all SNPs
          within each replica.
       mutPop - if the scenario neutral, this gives the putative selected population; that is, we will do analyses
          as if we were thinking (wrongly) that this is actually a selected region with selection in mutPop.
          This is needed because, when we localize a region, we currently assume that selection occurred in one particular
          population.  So to explore false positive rates, we need to run CMS on a neutral region while
          assuming selection in a particular population.

        oldMerged - determines the path to the merged.data/ file -- is it located under an old or a new (more uniform)
           scheme.

    See also: DefineRulesTo_computeCMSstats()      

    """

    if not Ddata.endswith( '/' ): Ddata += '/'

    if putativeMutPop == None: putativeMutPop = scenario.mutPop

    putativeMutPop = int( putativeMutPop )

    if not DdataMerged: DdataMerged = Ddata
    
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    
    snpStatsDirMerged = os.path.join( DdataMerged, 'snpStats' + thinSfx, scenario.scenDir() )

    mergedData = os.path.join( snpStatsDirMerged,  AddFileSfx( 'merged.tsv', putativeMutPop ) )

    if oldMerged:
        mergedData = os.path.join( DdataMerged,
                                   AddFileSfx( 'merged.data/',
                                               '%dky' % ( scenario.mutAge if not scenario.isNeutral() else 10 ),
                                               scenario.scenName(),
                                               putativeMutPop if scenario.isNeutral() else None ) )

    cmsStatsRawFN = outFile if outFile else os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsRaw.data/', putativeMutPop ) )

    fileDescrs = { cmsStatsRawFN:
                   ( 'CMS stats for all SNPs in all replicas of scenario $scenario, assuming selection happened in '
                     + pop2name[putativeMutPop] + '.  Each line is one SNP in one replica.',
                     ( ('Chrom', 'Chromosome or simulation replica of this SNP' ),
                       ('gdPos', 'Genetic map position of this SNP on its chromosome, in cM'),
                       ('iHS', 'Both iHS' ),
                       ( 'Pos', 'Position on chromosome, plus 1e6' ),
                       ( 'StdDiff', 'deltaiHH: ( Both_iHH_A - Both_iHH_D ), normalized by SNP frequency. '
                         'Both_iHH_A is the sum of iHH_A_left and iHH_A_right, etc.' ),
                       ( 'meanFst', 'Mean Fst comparison between the selected population and the other populations.' ),
                       ( 'derFreq', 'derived allele frequency of this SNP in the selected population' ),
                       ( 'max_xpop', 'the highest xpop comparison between the selected pop and the other pops' ),
                       ( 'meanAnc', 'the mean ancestral frequency in the non-selected pops.   deltaDAF is derived freq in '
                         'the selected pop, minus this.' ),
                       ( 'freqDiff',
                         'deltaDAF: Difference between the derived allele frequency in the selected population, '
                         'and the average of derived allele frequencies in the non-selected populations.' ) ) ) }
    
    if getio: return dict( depends_on = ( mergedData, ) +
                           ( ( pop2sampleSize, ) if isinstance( pop2sampleSize, types.StringTypes ) else () ),
                           creates = cmsStatsRawFN,
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ),
                           fileDescrs = fileDescrs )


    if isinstance( pop2sampleSize, types.StringTypes ): pop2sampleSize = dict( IDotData( pop2sampleSize ) )

    dbg( 'pop2sampleSize' )
    
    # Open neutral
    Data = DotData(SVPath = mergedData)
    
    # Merge Chrom columns
    
    popNames = sorted( pop2name.values() )
    popNums = sorted( pop2name.keys() )
    minPopNum = popNums[ 0 ]
    popPairs = [ '%s_%s' % ( popNames[ pop1idx ], popNames[ pop2idx ] )
                 for pop1idx in range( len( popNames ) ) for pop2idx in range( pop1idx+1, len( popNames ) ) ]


    # Calculate Fst
    xpopComparisonPairs = []
    xpopComparisonSigns = []

    for otherPop in popNums:
        if otherPop != putativeMutPop:
            popAname, popBname = pop2name[ putativeMutPop ], pop2name[ otherPop ]
            if popAname > popBname: popAname, popBname = popBname, popAname
            xpopComparisonPairs.append( '%s_%s' % ( popAname, popBname ) )
            xpopComparisonSigns.append( 1 if popAname == pop2name[ putativeMutPop ] else -1 )

    dbg( 'zip(xpopComparisonPairs,xpopComparisonSigns)' )

    logging.info( 'Computing Fst' )
    meanFst = zeros(len(Data))
    pop2ancFreqs = dict( ( popNum, 1 - array(Data['FREQ1 %d' % popNum ]) ) for popNum in popNums )
    Fst = fst( **Dict( 'sampleSize pop2name pop2ancFreqs pop2sampleSize' ) )
    for comp in xpopComparisonPairs:
        meanFst += Fst[comp]
    meanFst = meanFst / len( xpopComparisonPairs )

    logging.info( 'Computing max xpop' )
    # Max xpop
    xpopAll = zeros( ( 2 * len( xpopComparisonPairs ), len(Data) ) )
    j = 0
    for xpopComparisonPair, xpopComparisonSign in zip( xpopComparisonPairs, xpopComparisonSigns ):
        xpopAll[ j ] = xpopComparisonSign * Data['L AllEHH logratio Deviation ' + xpopComparisonPair]
        xpopAll[ j + 1 ] = xpopComparisonSign * Data['R AllEHH logratio Deviation ' + xpopComparisonPair]
        #print xpopName[putativeMutPop][i],xpopSign[putativeMutPop][i]
        j += 2

    max_xpop = nanmax(xpopAll,axis=0)
    #print max_xpop

    ihs = array(Data['Both iHS'])

    aPopPair = 'European_WestAfrican' if 'European_WestAfrican' in popPairs else popPairs[0]
    
    gdPos = array(Data['SNP pos (cM) ' + aPopPair])
    bpPos = array(Data['CHROM_POS %d' % minPopNum])

    # Calculate iHH difference
    iHHDiff = array(Data['Both iHH_D']) - array(Data['Both iHH_A'])
    StdDiff = normalizeByFreq(iHHDiff,Data['FREQ1 %d' % putativeMutPop])

    # Derived frequencies
    derFreq = array(Data['FREQ1 %d' % putativeMutPop])

    logging.info( 'Computing mean anc' )

    # Mean ancestral freq
    mean_anc = zeros(len(Data))
    for popNum in popNums:
        if popNum != putativeMutPop:
            mean_anc += (1 - Data['FREQ1 %d' % popNum])
    mean_anc /= ( len( popNums ) - 1 )
    
    # Freq diff
    freqDiff = derFreq - (1 - mean_anc)

    # Make new pos column
    Pos = Data['CHROM_POS %d' % minPopNum]

    Chrom = map( float, Data.replicaNum if 'replicaNum' in Data.dtype.names
                 else where( isnan( Data.Chrom ),
                             Data[ 'Chrom ' + aPopPair ],
                             Data.Chrom ) )
    
    cmsStatsRaw = DotData( Columns = [Chrom, Pos, gdPos, ihs, StdDiff, meanFst, derFreq,
                                      max_xpop, mean_anc, freqDiff,iHHDiff, Data['FREQ1 %d' % putativeMutPop] ],
                           names = ['Chrom','Pos','gdPos','iHS','StdDiff','meanFst','derFreq',
                                    'max_xpop','meanAnc','freqDiff','iHHDiff', 'normingFreqs'] )

    assert set( stats ) <= set( cmsStatsRaw.dtype.names[3:] )
    cmsStatsRaw.save( cmsStatsRawFN )

def computeCMSstats( Ddata, thinSfx,
                     scenario,
                     putativeMutPop = None, sampleSize = 120,
                     pop2name = pop2name,
                     pop2sampleSize = {},
                     oldMerged = False,
                     DdataMerged = None,
                     outFile = None,
                     statsSfx = '',
                     simsOut = 'simsOut',
                     stats = CMSBins.CMSstats,
                     ihsSfx = '',
                     ihhBackgroundArgs = None,
                     limitToPop = None,
                     getio = None ):
    """Compute CMS stats for all SNPs in all replicas of a given scenario.

    Params:

       Ddata - the root directory containing simulations and their analyses
       thinSfx - specifies which thinning version of the simulations are used.
       scen - the simulation scenario.   Here we compute stats for all replicas within that scenario, for all SNPs
          within each replica.
       mutPop - if the scenario neutral, this gives the putative selected population; that is, we will do analyses
          as if we were thinking (wrongly) that this is actually a selected region with selection in mutPop.
          This is needed because, when we localize a region, we currently assume that selection occurred in one particular
          population.  So to explore false positive rates, we need to run CMS on a neutral region while
          assuming selection in a particular population.

        oldMerged - determines the path to the merged.data/ file -- is it located under an old or a new (more uniform)
           scheme.

    See also: DefineRulesTo_computeCMSstats()      

    """

    if not Ddata.endswith( '/' ): Ddata += '/'

    if putativeMutPop == None: putativeMutPop = scenario.mutPop

    putativeMutPop = int( putativeMutPop )

    if not DdataMerged: DdataMerged = Ddata
    
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    
    snpStatsDirMerged = os.path.join( DdataMerged, 'snpStats' + thinSfx, scenario.scenDir() )

    mergedData = GetCreates( addIHHDiff, Ddata = DdataMerged, **Dict( 'scenario putativeMutPop simsOut statsSfx thinSfx pop2name ihsSfx' ) )[0]

    if oldMerged:
        mergedData = os.path.join( DdataMerged,
                                   AddFileSfx( 'merged.data/',
                                               '%dky' % ( scenario.mutAge if not scenario.isNeutral() else 10 ),
                                               scenario.scenName(),
                                               putativeMutPop if scenario.isNeutral() else None ) )

    cmsStatsRawFN = outFile if outFile else os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsRaw.tsv', statsSfx, putativeMutPop, ihsSfx ) )

    args = Dict( 'Ddata thinSfx scenario putativeMutPop simsOut statsSfx pop2name ihsSfx' )
    if ihhBackgroundArgs is not None: args = MergeDicts( ihhBackgroundArgs,
                                                         Dict( 'scenario putativeMutPop' ) )
    ihhDiff_sumsByFreqFN = getFN_ihhDiff_sumsByFreq( **args )
    ihhDiff_sumsFN = getFN_ihhDiff_sums( **args )

    fileDescrs = { cmsStatsRawFN:
                   ( 'CMS stats for all SNPs in all replicas of scenario $scenario, assuming selection happened in '
                     + pop2name[putativeMutPop] + '.  Each line is one SNP in one replica.',
                     ( ('Chrom', 'Chromosome or simulation replica of this SNP' ),
                       ('gdPos', 'Genetic map position of this SNP on its chromosome, in cM'),
                       ('iHS', 'Both iHS' ),
                       ( 'Pos', 'Position on chromosome, plus 1e6' ),
                       ( 'StdDiff', 'deltaiHH: ( Both_iHH_A - Both_iHH_D ), normalized by SNP frequency. '
                         'Both_iHH_A is the sum of iHH_A_left and iHH_A_right, etc.' ),
                       ( 'meanFst', 'Mean Fst comparison between the selected population and the other populations.' ),
                       ( 'derFreq', 'derived allele frequency of this SNP in the selected population' ),
                       ( 'max_xpop', 'the highest xpop comparison between the selected pop and the other pops' ),
                       ( 'meanAnc', 'the mean ancestral frequency in the non-selected pops.   deltaDAF is derived freq in '
                         'the selected pop, minus this.' ),
                       ( 'freqDiff',
                         'deltaDAF: Difference between the derived allele frequency in the selected population, '
                         'and the average of derived allele frequencies in the non-selected populations.' ) ) ) }
    
    if getio: return dict( depends_on = ( mergedData, ihhDiff_sumsByFreqFN, ihhDiff_sumsFN ) +
                           ( ( pop2sampleSize, ) if isinstance( pop2sampleSize, types.StringTypes ) else () ),
                           creates = cmsStatsRawFN,
                           splitByCols = { mergedData: dict( keyCols = () ) },
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ),
                           fileDescrs = fileDescrs )


    if isinstance( pop2sampleSize, types.StringTypes ): pop2sampleSize = dict( IDotData( pop2sampleSize ) )

    dbg( 'pop2sampleSize' )
    
    theData = IDotData(mergedData)

    popNames = sorted( pop2name.values() )
    popNums = sorted( pop2name.keys() )
    minPopNum = popNums[ 0 ]
    popPairs = [ '%s_%s' % ( popNames[ pop1idx ], popNames[ pop2idx ] )
                 for pop1idx in range( len( popNames ) ) for pop2idx in range( pop1idx+1, len( popNames ) )
                 if limitToPop is None or limitToPop in ( popNames[ pop1idx ], popNames[ pop2idx ] ) ]

    xpopComparisonPairs = []
    xpopComparisonSigns = []

    for otherPop in popNums:
        if otherPop != putativeMutPop:
            popAname, popBname = pop2name[ putativeMutPop ], pop2name[ otherPop ]
            if popAname > popBname: popAname, popBname = popBname, popAname
            xpopComparisonPairs.append( '%s_%s' % ( popAname, popBname ) )
            xpopComparisonSigns.append( 1 if popAname == pop2name[ putativeMutPop ] else -1 )

    dbg( 'zip(xpopComparisonPairs,xpopComparisonSigns)' )

    #
    # For normalizing iHHDiff by frequency bin, get the mean and stddev
    # within each frequency bin.
    #

    ihhDiff_statsByFreq = IDotData( ihhDiff_sumsByFreqFN ).addMeanStdCols( 'iHHDiff' ).toDotData()
    ihhDiff_stats = IDotData( ihhDiff_sumsFN )
    ihhDiff_stats = ihhDiff_stats.addMeanStdCols( 'iHHDiff' )
    
    ihhDiff_stats = ihhDiff_stats[0]

    totStd = ihhDiff_stats.iHHDiff_std
    ihhDiffMeans = dict( ( r.freqBinId, r.iHHDiff_mean ) for r in ihhDiff_statsByFreq )

    dbg( 'ihhDiffMeans ihhDiff_stats totStd' )

    #
    # ok so the next thing is to rewrite this, and you should not need
    # to make a one-line DotData.  and make sure each record we're writing is a value.
    #

    #
    # then, check that all stat values we compute are almost equal to each other.
    #
    
    with IDotData.openForWrite( cmsStatsRawFN,
                                headings =
                                'Chrom Pos gdPos iHS StdDiff meanFst derFreq max_xpop meanAnc freqDiff iHS_nanReason StdDiff_nanReason' ) \
                                as cmsStatsRaw:
        for r in theData:

            pop2ancFreq = dict( ( popNum, 1 - r['FREQ1 %d' % popNum ] ) for popNum in popNums )
            Fst = fst_oneSnp( **Dict( 'sampleSize pop2name pop2ancFreq pop2sampleSize' ) )
            meanFst = np.mean([ Fst[comp] if not np.isnan( Fst[comp] ) else 0.0 for comp in xpopComparisonPairs ])

            # Max xpop
            xpopAll = np.zeros( 2 * len( xpopComparisonPairs ) )
            j = 0
            for xpopComparisonPair, xpopComparisonSign in zip( xpopComparisonPairs, xpopComparisonSigns ):
                xpopAll[ j ] = xpopComparisonSign * r['L AllEHH logratio Deviation ' + xpopComparisonPair]
                xpopAll[ j + 1 ] = xpopComparisonSign * r['R AllEHH logratio Deviation ' + xpopComparisonPair]
                j += 2

            max_xpop = np.nanmax(xpopAll)

            derFreq = r['FREQ1 %d' % putativeMutPop]
            ancFreq = 1 - derFreq

            ihs = r['Both iHS']

            ihsNanReason = -1
            if np.isnan( ihs ):
                if derFreq < .05: ihsNanReason = 0
                if derFreq > 0.95: ihsNanReason = 1

            aPopPair = 'European_WestAfrican' if 'European_WestAfrican' in popPairs else popPairs[0]

            gdPos = r['SNP pos (cM) ' + aPopPair]
            bpPos = r['CHROM_POS %d' % minPopNum]

            # Calculate iHH difference
            iHHDiff = r['Both iHH_D'] - r['Both iHH_A']

            StdDiff = ( r.iHHDiff - ihhDiffMeans[ r.freqBinWith01Id ] ) / totStd
            
            # Mean ancestral freq
            mean_anc = 0
            for popNum in popNums:
                if popNum != putativeMutPop:
                    mean_anc += (1 - r['FREQ1 %d' % popNum])
            mean_anc /= ( len( popNums ) - 1 )

            # Freq diff
            freqDiff = derFreq - (1 - mean_anc)

            # Make new pos column
            Pos = r['CHROM_POS %d' % minPopNum]

            Chrom = float( r.replicaNum if 'replicaNum' in theData.headings
                           else ( r[ 'Chrom ' + aPopPair ] if np.isnan( r.Chrom) else r.Chrom ) )

            cmsStatsRaw.writeRecord( Chrom, Pos, gdPos, ihs, StdDiff, meanFst, derFreq,
                                     max_xpop, mean_anc, freqDiff,
                                     # nan reasons
                                     ihsNanReason, ihsNanReason )
                                     
    assert set( stats ) <= set( cmsStatsRaw.headings[3:] )

    

def GetNeutralMeanStd( Ddata, thinSfx, putativeMutPop,
                       stats = CMSBins.CMSstats, pop2name = pop2name, getio = None ):
    """Get the mean and stddev of CMS scores for neutral scenario"""

    if not Ddata.endswith( '/' ): Ddata += '/'
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, NeutralScenario().scenDir() )
    cmsStatsRawFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsRaw.data/', putativeMutPop ) )

    cmsStatsMeanStdFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsMeanStd.tsv', putativeMutPop ) )

    fileDescrs = { cmsStatsMeanStdFN :
                       ( 'Mean and stddev of raw CMS stats (' + ','.join(stats) + ') computed for all SNPs in all '
                         'replicas in neutral scenario under assumption of selection in ' +
                         pop2name[ putativeMutPop ] + '.',
                         ( ( 'stat', 'column name of statistic' ),
                           ( 'mean,stddev', 'mean,stddev of statistic across all SNPs in all neutral replicas' ) ) ) }

    if getio: return dict( depends_on = cmsStatsRawFN, creates = cmsStatsMeanStdFN,
                           mediumRuleNameSfx = putativeMutPop,
                           fileDescrs = fileDescrs )

    cmsStatsRaw = DotData( Path = cmsStatsRawFN )

    Records = []
    for stat in stats:
        statNorm, statMean, statStd = normalize( array( cmsStatsRaw[ stat ] ) ) 
        Records.append( ( stat, statMean, statStd ) )

    dbg( 'Records' )
    DotData( names = ( 'stat', 'mean', 'stddev' ), Records = Records ).saveToSV( cmsStatsMeanStdFN )
    
def normalizeStats_old( Ddata, thinSfx, scenario, putativeMutPop, stats = CMSBins.CMSstats,
                        DdataNeutral = None, getio = None ):
    """Compute normalized versions of CMS statistics computed assuming selection in $putativeMutPop,
    for all SNPs in all replicas of scenario $scenario, normalizing each statistic by its mean and stddev for all neutral
    SNPs in all replicas of the neutral scenario (computed assuming selection in $putativeMutPop)"""
    
    if not Ddata.endswith( '/' ): Ddata += '/'
    if not DdataNeutral: DdataNeutral = Ddata
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    snpStatsDirNeutral = os.path.join( DdataNeutral, 'snpStats' + thinSfx, NeutralScenario().scenDir() )

    cmsStatsRawFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsRaw.data/', putativeMutPop ) )

    cmsStatsMeanStdFN = os.path.join( snpStatsDirNeutral, AddFileSfx( 'cmsStatsMeanStd.tsv', putativeMutPop ) )
    cmsStatsNormFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsNorm.data/', putativeMutPop ) )

    fileDescrs = { cmsStatsNormFN: normalizeStats.__doc__ }

    depends_on = [ cmsStatsRawFN ]
    creates = [ cmsStatsNormFN ]
    ( creates if scenario.isNeutral() else depends_on ).append( cmsStatsMeanStdFN )
    
    if getio: return dict( depends_on = depends_on, creates = creates,
                           fileDescrs = fileDescrs, mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ),
                           name = 'normalizeStats' + ( '_neutral' if scenario.isNeutral() else '' ) )

    cmsStatsRaw = IDotData( cmsStatsRawFN )
    cmsStatsMeanStd = [] if scenario.isNeutral() else DotData( SVPath = cmsStatsMeanStdFN )

    # var: stat2row - map from name of CMS statistic to the row containing its mean and stddev
    #   in cmsStatsMeanStd
    if not scenario.isNeutral(): stat2row = dict( ( stat, row ) for row, stat in enumerate( cmsStatsMeanStd.stat ) )
    stat2normed = {}
    for stat in stats:
        if scenario.isNeutral():
            stat2normed[ stat ], statMean, statStd = normalize( array( cmsStatsRaw[ stat ] ) )
            cmsStatsMeanStd.append( ( stat, statMean, statStd ) )
        else:
            statRow = cmsStatsMeanStd[ stat2row[ stat ] ]
            statMean, statStd = statRow['mean'], statRow['stddev']
            dbg( 'type(statRow) statRow type(statMean) statMean type(statStd) statStd' )
            stat2normed[ stat ] = ( cmsStatsRaw[ stat ] - statMean ) / statStd

    dbg( 'stat2normed cmsStatsRaw' )
    DotData( names = cmsStatsRaw.dtype.names,
             Columns = [ ( stat2normed if c in stats else cmsStatsRaw )[ c ]
                         for c in cmsStatsRaw.dtype.names ] ).save( cmsStatsNormFN )

    if scenario.isNeutral():
        DotData( names = ( 'stat', 'mean', 'stddev' ), Records = cmsStatsMeanStd ).saveToSV( cmsStatsMeanStdFN )

def normalizeStats( Ddata, thinSfx, scenario, putativeMutPop,
                    stats = CMSBins.CMSstats,
                    DdataNeutral = None, outFile = None, outFileMeanStd = None, statsSfx = '',
                    pop2name = pop2name, ihsSfx = '', getio = None ):
    """Compute normalized versions of CMS statistics computed assuming selection in $putativeMutPop,
    for all SNPs in all replicas of scenario $scenario, normalizing each statistic by its mean and stddev for all neutral
    SNPs in all replicas of the neutral scenario (computed assuming selection in $putativeMutPop)"""
    
    if not Ddata.endswith( '/' ): Ddata += '/'
    if not DdataNeutral: DdataNeutral = Ddata
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    snpStatsDirNeutral = os.path.join( DdataNeutral, 'snpStats' + thinSfx, NeutralScenario().scenDir() )

    cmsStatsRawFN = GetCreates( computeCMSstats, **Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx pop2name ihsSfx' ) )[0]

    cmsStatsMeanStdFN = outFileMeanStd if outFileMeanStd else os.path.join( snpStatsDirNeutral,
                                                                            AddFileSfx( 'cmsStatsMeanStd.tsv', statsSfx, putativeMutPop, ihsSfx ) )
    cmsStatsNormFN = outFile if outFile else os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsNorm.data/', statsSfx, putativeMutPop, ihsSfx ) )

    fileDescrs = { cmsStatsNormFN: normalizeStats.__doc__ }

    depends_on = [ cmsStatsRawFN ]
    creates = [ cmsStatsNormFN ]
    ( creates if scenario.isNeutral() else depends_on ).append( cmsStatsMeanStdFN )
    
    if getio: return dict( depends_on = depends_on, creates = creates,
                           fileDescrs = fileDescrs, mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ),
                           attrs = dict( scenario = scenario.scenDir() ),
                           name = 'normalizeStats' + ( '_neutral' if scenario.isNeutral() else '' ) )

    cmsStatsRaw = IDotData( cmsStatsRawFN )

    cmsStatsMeanStd = cmsStatsRaw.getColumnStats( *stats ).renameCols( {'col':'stat','std':'stddev'} ).save( cmsStatsMeanStdFN ) \
        if scenario.isNeutral() else IDotData( cmsStatsMeanStdFN )

    cmsStatsRaw.normalizeColumns( *stats, meansStds = cmsStatsMeanStd[( 'mean', 'stddev' )] ).save( cmsStatsNormFN )

def normalizeStatsWithinReplicas( Ddata, thinSfx, scenario, putativeMutPop, stats = CMSBins.CMSstats, statsSfx = '',
                                  pop2name = pop2name, ihsSfx = '', getio = None ):
    """Compute normalized versions of CMS statistics computed assuming selection in $putativeMutPop,
    for all SNPs in all replicas of scenario $scenario, normalizing each statistic within the replica."""
    
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )

    cmsStatsRawFN = GetCreates( computeCMSstats, **Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx pop2name ihsSfx' ) )[0]
    cmsStatsNormedLocalFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsNormedLocal.data/', statsSfx, putativeMutPop, ihsSfx ) )

    if getio: return dict( depends_on = cmsStatsRawFN, creates = cmsStatsNormedLocalFN,
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ) )

    IDotData( cmsStatsRawFN ).normalizeColumnsWithinGroups( cols = stats,
                                                            groupCols = ( 'Chrom', ) ).save( cmsStatsNormedLocalFN )

def normalizeStatsWithinReplicas2( Ddata, thinSfx, scenario, putativeMutPop, nreplicas, stats = CMSBins.CMSstats,
                                   getio = None ):
    """Compute normalized versions of CMS statistics computed assuming selection in $putativeMutPop,
    for all SNPs in all replicas of scenario $scenario, normalizing each statistic within the replica."""
    
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )

    cmsStatsRawFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsRaw.data/', putativeMutPop ) )
    cmsStatsNormedLocalFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsNormedLocal2.data/', putativeMutPop ) )

    if getio: return dict( depends_on = cmsStatsRawFN, creates = cmsStatsNormedLocalFN,
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ) )

    cmsStatsRaw = DotData( Path = cmsStatsRawFN )

    normedList = []
    
    for replicaNum in range( nreplicas ):

        r = cmsStatsRaw[ cmsStatsRaw.Chrom == replicaNum ]
        r_normed = DotData( names = r.dtype.names,
                            Columns = [ r[c] if c not in stats else normalize(r[c])[0] for c in r.dtype.names ] )
        normedList.append( r_normed )

    datavstack( normedList ).save( cmsStatsNormedLocalFN )


def normalizeCMS_old( Ddata, thinSfx, scenario, putativeMutPop, complikeSfx, likesTableSfx, 
                      DdataNeutral = None, getio = None ):
    """Compute normalized versions of CMS computed assuming selection in $putativeMutPop,
    for all SNPs in all replicas of scenario $scenario, normalizing each statistic by its mean and stddev for all neutral
    SNPs in all replicas of the neutral scenario (computed assuming selection in $putativeMutPop)"""
    
    if not Ddata.endswith( '/' ): Ddata += '/'
    if not DdataNeutral: DdataNeutral = Ddata
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    snpStatsDirNeutral = os.path.join( DdataNeutral, 'snpStats' + thinSfx, NeutralScenario().scenDir() )

    cmsRawFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', putativeMutPop, complikeSfx, likesTableSfx ) )

    cmsMeanStdFN = os.path.join( snpStatsDirNeutral, AddFileSfx( 'cmsNeutralMeanStd.tsv', putativeMutPop,
                                                                 complikeSfx, likesTableSfx ) )
    cmsNormFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', 'normed', putativeMutPop, complikeSfx,
                                                       likesTableSfx ) )

    depends_on = [ cmsRawFN, cmsMeanStdFN ]
    creates = [ cmsNormFN ]
    
    if getio: return dict( depends_on = depends_on, creates = creates,
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ),
                           name = 'normalizeCMS' ) 

    cmsRaw = DotData( Path = cmsRawFN )
    cmsMeanStd = DotData( SVPath = cmsMeanStdFN )

    # var: stat2row - map from name of CMS statistic to the row containing its mean and stddev
    #   in cmsStatsMeanStd
    stat2row = dict( ( stat, row ) for row, stat in enumerate( cmsMeanStd.stat ) )
    stat2normed = {}
    stats = ( 'complike', 'complikeExp', 'complikeRatio', 'complikeRatioExp' )
    for stat in stats:
        statRow = cmsMeanStd[ stat2row[ stat ] ]
        stat2normed[ stat ] = ( cmsRaw[ stat ] - statRow[ 'mean' ] ) / statRow[ 'stddev' ]

    DotData( names = cmsRaw.dtype.names,
             Columns = [ ( stat2normed if c in stats else cmsRaw )[ c ]
                         for c in cmsRaw.dtype.names ] ).save( cmsNormFN )

def normalizeCMS( Ddata, thinSfx, scenario, putativeMutPop, complikeSfx, likesTableSfx, statsSfx = '',
                  DdataNeutral = None, pop2name = pop2name, getio = None ):
    """Compute normalized versions of CMS computed assuming selection in $putativeMutPop,
    for all SNPs in all replicas of scenario $scenario, normalizing each statistic by its mean and stddev for all neutral
    SNPs in all replicas of the neutral scenario (computed assuming selection in $putativeMutPop)"""
    
    if not Ddata.endswith( '/' ): Ddata += '/'
    if not DdataNeutral: DdataNeutral = Ddata
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    snpStatsDirNeutral = os.path.join( DdataNeutral, 'snpStats' + thinSfx, NeutralScenario().scenDir() )

    cmsRawFN = GetCreates( computeCMS, **Dict( 'Ddata thinSfx scenario putativeMutPop complikeSfx likesTableSfx statsSfx pop2name' ) ).containing( 'complike' )

    cmsMeanStdFN = GetCreates( getNeutralMeanStdForCMS, Ddata = DdataNeutral, **Dict( 'thinSfx putativeMutPop complikeSfx likesTableSfx statsSfx pop2name' ) )[0]
    cmsNormFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', 'normed', putativeMutPop, complikeSfx,
                                                       likesTableSfx ) )

    depends_on = [ cmsRawFN, cmsMeanStdFN ]
    creates = [ cmsNormFN ]
    
    if getio: return dict( depends_on = depends_on, creates = creates,
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ),
                           attrs = Dict( 'putativeMutPop', scenario = scenario.scenDir() ),
                           split = ( ( cmsRawFN, () ), ),
                           name = 'normalizeCMS' ) 

    cmsRaw = IDotData( cmsRawFN )
    cmsMeanStd = DotData( SVPath = cmsMeanStdFN )
    stat2row = dict( ( stat, row ) for row, stat in enumerate( cmsMeanStd.stat ) )

    meansStds = [ (r['mean'], r['stddev']) for r in cmsMeanStd ]
    dbg( 'meansStds' )

    cols = tuple( cmsMeanStd.stat )
    cmsRaw.normalizeColumns( *cols,
                             meansStds = [ (r['mean'], r['stddev']) for r in cmsMeanStd ] ).save( cmsNormFN )

def getNeutralMeanStdForCMS( Ddata, thinSfx, putativeMutPop, complikeSfx, likesTableSfx, statsSfx, pop2name = pop2name, getio = None ):
    """Get the mean and stddev of CMS scores for neutral scenario"""

    if not Ddata.endswith( '/' ): Ddata += '/'
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, NeutralScenario().scenDir() )
    complikeFN = GetCreates( computeCMS, scenario = NeutralScenario(), **Dict( 'Ddata pop2name thinSfx putativeMutPop complikeSfx likesTableSfx statsSfx pop2name' ) )[0]
    cmsNeutralMeanStdFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsNeutralMeanStd.tsv',
                                                                 putativeMutPop, complikeSfx, likesTableSfx ) )

    if getio: return dict( depends_on = complikeFN, creates = cmsNeutralMeanStdFN,
                           mediumRuleNameSfx = putativeMutPop )

    complike = IDotData( complikeFN )

    complikeStats = complike.getColumnStats( 'complike', 'complikeExp', 'complikeRatio', 'complikeRatioExp' )
    
    Records = []
    for r in complikeStats:
        Records.append( ( r.col, r.mean, r.std ) )

    DotData( names = ( 'stat', 'mean', 'stddev' ), Records = Records).saveToSV( cmsNeutralMeanStdFN )

def histogramCMSstats( Ddata, thinSfx, scenario, putativeMutPop, stats = CMSBins.CMSstats, getio = None ):
    """For a given scenario, for selected replicas in that scenario,
    for each statistic, histogram the distribution of the statistic
    for causal and non-causal SNPs.
    """

    print 'deprecated!'
    sys.exit(1)

    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    cmsStatsNormFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsNorm.data/', putativeMutPop ) )
    cmsStatsHistFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsHist.data/', putativeMutPop ) )

    fileDescrs = { cmsStatsHistFN:
                       ( 'Histogram of CMS statistic values for all causal and non-causal SNPs in all replicas of scenario '
                         '$scenario.  Each row represents one bin.   Bins are defined in computeCMSstats.CMSBins. ',
                         [ ( stat + which, Str( 'Number of <b>$kind</b> SNPs with <i>$stat</i> in given bin, '
                                                'among all <b>$kind</b> SNPs in all '
                                                'replicas of scenario <code>$scenario</code>' ) )
                           for stat in stats for which, kind in ( 'Hits', 'causal' ), ( 'Misses', 'non-causal' ) ] ) }
    
    if getio: return dict( depends_on = cmsStatsNormFN, creates = cmsStatsHistFN,
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ),
                           fileDescrs = fileDescrs )

    cmsStatsNorm = DotData( Path = cmsStatsNormFN )

    selpos = 500000 if not scenario.isNeutral() else None

    statCols = []
    statColsNames = []
    
    for stat in stats:

        statNorm = cmsStatsNorm[ stat ]
        statHits = zeros( CMSBins.nbins, dtype = int )
        statMisses = zeros( CMSBins.nbins, dtype = int )

        scoreVal( statNorm, cmsStatsNorm.Pos, CMSBins.nbins, CMSBins.stat_start[ stat ], CMSBins.stat_end[ stat ],
                  statHits, statMisses, selpos )

        statCols += [ statHits, statMisses ]
        statColsNames += [ stat + 'Hits', stat + 'Misses' ]

    DotData( names = statColsNames, Columns = statCols ).save( cmsStatsHistFN )

    
def addHistogramsForSelectionScenarios( inFiles, outFile, stats = CMSBins.CMSstats, getio = None ):
    """Given a number of input tables, add them up and create an output table.
    """

    fileDescrs = { outFile:
                       ( 'Histogram of CMS statistic values for all causal and non-causal SNPs in all replicas of '
                         'all selection scenarios.  Each row represents one bin.   '
                         'Bins are defined in computeCMSstats.CMSBins. ',
                         [ ( stat + which, Str( 'Number of <b>$kind</b> SNPs with <i>$stat</i> in given bin, '
                                                'among all <b>$kind</b> SNPs in all '
                                                'replicas of all selection scenarios' ) )
                           for stat in stats for which, kind in ( 'Hits', 'causal' ), ( 'Misses', 'non-causal' ) ] ) }
    
    if getio: return dict( depends_on = inFiles, creates = outFile, fileDescrs = fileDescrs )

    hists = [ DotData( Path = inFile ) for inFile in inFiles ]

    dtype = hists[0].dtype
    assert all( [ hist.dtype == dtype for hist in hists ] )

    cols = dtype.names
    DotData( names = cols,
             Columns = [ reduce( add, map( itemgetter( col ), hists ) ) for col in cols ] ).save( outFile )


def LoadBins( fname ):    
    """Read in bin information from summary file (written by analsim)"""
    stat_start = {}
    stat_end = {}
    stat_nbins = {}
    
    with open( fname ) as inf:
        firstline = inf.readline().rstrip()
        parts = firstline.split(' ')
        
        header = inf.readline()
        
        lines = inf.readlines()
        for line in lines:
            parts = line.rstrip().split('\t')
            stat_nbins[parts[0]] = int(parts[3])
            stat_start[parts[0]] = float(parts[1])
            stat_end[parts[0]] = float(parts[2])


    if 'Fst' in stat_start:
        stat_start['fst'] = stat_start['Fst']
        stat_end['fst'] = stat_end['Fst']
        stat_nbins['fst'] = stat_nbins['Fst']

    return stat_start, stat_end, stat_nbins

def LoadLikesTable( likesTable, likesTableSfx = '', mutPop = None, Ddata = None  ):
    """Load a likes table, returning a dictionary that maps likes table
    attributes to values.  For example, gives the name of hitLikes and missLikes files."""
    if isinstance( likesTable, types.StringTypes ):
        if Ddata and not os.path.dirname( likesTable ): likesTable = os.path.join( Ddata, likesTable )
        dbg( 'likesTable likesTableSfx' )
        likesTable = IDotData( AddFileSfx( likesTable, likesTableSfx, mutPop ) )
    return AttrDict( ( likesParam, likesParamVal ) for likesParam, likesParamVal in likesTable )
    

def computeCMS_old( Ddata, thinSfx, scenario, putativeMutPop, nreplicas = 100, likesTable = None,
                    likesTableSfx = '',
                    likesSfx = 'Newer', likesDir = None,
                    likesBins = None, likesPopSpecific = True,
                    stats = CMSBins.CMSstats, nonNormedStats = CMSBins.nonNormedStats,
                    nonNanStats = ( 'iHS', 'max_xpop' ),  #CMSBins.CMSstats,
                    nearCausalSfx = '', freq = 'all',
                    complikeSfx = None, rawNormSfx = '',
                    normalizeWithinReplicas = False,
                    forceNormalizeWithinReplicas = False,
                    useRegionLikes = True,
                    cleanIHS = True, cleanDer = False, cleanAnc = False, keepCausal = False,
                    outputFN = None, useDefaultNumSnps = False,
                    getio = None ):
    """Compute CMS scores for all snps in all replicas of scenario $scenario, using given likes tables.

    Params:

       complikeSfx - a suffix that identifies this version of complike scores.  This helps distinguish, e.g.,
          complike scores computed from the same data using different likes tables.
          Defaults to likesSfx.   We could just use likesSfx, but it does not capture the location (directory) of the
          likes tables -- different likes tables in different directories could have the same likes sfx.

    """

    #stats = ( 'max_xpop', )
    
    print 'in computeCMS: likesBins is ', likesBins, 'likesSfx is ', likesSfx, ' likesTable is ', likesTable

    if likesTable:
        likesTable = LoadLikesTable( **Dict( 'Ddata likesTable likesTableSfx' ) )
        likesBins = likesTable.likesBins
        normalizeWithinReplicas = coerceVal( likesTable.normalizeWithinReplicas )

    if forceNormalizeWithinReplicas: normalizeWithinReplicas = True
    
    if likesDir == None: likesDir = Ddata

    if putativeMutPop == None: putativeMutPop = scenario.mutPop

    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    cmsStatsNormFN = os.path.join( snpStatsDir,
                                   AddFileSfx( 'cmsStatsNorm' + ( 'edLocal' if normalizeWithinReplicas else '' )
                                               + '.data/', putativeMutPop, rawNormSfx, likesTableSfx ) )
    cmsStatsRawFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsRaw.data/', putativeMutPop,
                                                           rawNormSfx, likesTableSfx ) )

    hitLikesFN = likesTable.hitsLikes if likesTable else  \
        os.path.join( likesDir, 'likes',
                      AddFileSfx( 'hitsLikes.tsv', likesSfx, putativeMutPop if likesPopSpecific else None ) )
    missLikesFN = likesTable[ ('region' if useRegionLikes else 'miss') + 'Likes' ] if likesTable else \
        os.path.join( likesDir, 'likes',
                      AddFileSfx( ( 'region' if normalizeWithinReplicas and useRegionLikes else 'miss')
                                  + 'Likes.tsv',
                                  likesSfx, putativeMutPop if likesPopSpecific else None ) )

    complikeFN = outputFN or os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', putativeMutPop, complikeSfx, likesTableSfx ) )

    if likesBins and not os.path.dirname( likesBins ): likesBins = os.path.join( likesDir, likesBins )
    
    if getio: return dict( depends_on = ( hitLikesFN, missLikesFN, cmsStatsNormFN, cmsStatsRawFN ) +
                           ( (likesBins,) if likesBins else () ),
                           creates = complikeFN,
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop, likesTableSfx ) )

    hitLikes = DotData(SVPath = hitLikesFN )
    missLikes = DotData(SVPath = missLikesFN )

    cmsStatsNorm = DotData( Path = cmsStatsNormFN )
    cmsStatsRaw = DotData( Path = cmsStatsRawFN )

    cmsStatsNorm = cmsStatsNorm[ cmsStatsNorm.Chrom < nreplicas ]
    cmsStatsRaw = cmsStatsRaw[ cmsStatsRaw.Chrom < nreplicas ]

    #cmsStatsNorm = cmsStatsNorm[ cmsStatsNorm.Pos == 75103 ] 
    #cmsStatsRaw = cmsStatsRaw[ cmsStatsNorm.Pos == 75103 ] 

    dbg( 'len(cmsStatsNorm) len(cmsStatsRaw) cmsStatsNorm cmsStatsRaw' )
    
    ind = any( [ isfinite( cmsStatsNorm[ stat ] ) for stat in nonNanStats ], axis = 0 ) if len(cmsStatsNorm) > 0 else []

    statLikeRatios = []
    statLikes = []
    statLikesBins = []

    oldStat = { 'iHS': 'ihs', 'StdDiff': 'StdDiff', 'meanFst': 'fst', 'freqDiff': 'freqDiff', 'max_xpop': 'max_xpop' }

    if not likesBins:
        stat_start, stat_end, stat_nbins = CMSBins.stat_start, CMSBins.stat_end, CMSBins.stat_nbin
    else:
        stat_start, stat_end, stat_nbins = LoadBins( likesBins )

    for stat in stats:
        statNorm = ( cmsStatsNorm if stat not in nonNormedStats else cmsStatsRaw )[ stat ]

        statLikeRatio = zeros(len(statNorm))
        #statName = statNameMap[ stat ] if statNameMap.has_key( stat ) else stat
        statName = stat if stat in hitLikes.dtype.names else oldStat[ stat ]
        dbg( 'hitLikes.dtype.names statName stat_start stat_end stat_nbins' )
        dbg( 'stat statName statNorm' )
        likeRatio( statNorm, hitLikes[ statName ], missLikes[ statName ],
                   stat_start[ statName ], stat_end[ statName ],
                   stat_nbins[ statName ], statLikeRatio, ind )
        statLikeRatios.append( statLikeRatio )

        statLike = zeros(len(statNorm))
        likeVal( statNorm, hitLikes[ statName ], missLikes[ statName ],
                 stat_start[ statName ], stat_end[ statName ],
                 stat_nbins[ statName ], statLike, ind, useDefaultNumSnps = useDefaultNumSnps )
        statLikes.append( statLike )

        dbg( '"RRRRR" statName hitLikes[statName] missLikes[statName]' )

        statLikesBin = zeros(len(statNorm))
        likeValBin( statNorm, hitLikes[ statName ], missLikes[ statName ],
                    stat_start[ statName ], stat_end[ statName ],
                    stat_nbins[ statName ], statLikesBin, ind )
        statLikesBins.append( statLikesBin )

        dbg( 'stat statNorm statLike statLikeRatio statLikesBin' )

    complikeRatio = reduce( add, statLikeRatios )
    complike = reduce( add, statLikes )

    dbg( 'complike' )

    indClean = clean_hapmap2_region( Z = DotData( names = ( 'complike', 'StdDiff', 'derFreq', 'meanAnc', 'Pos' ),
                                                  Columns = ( complike, cmsStatsRaw.StdDiff, cmsStatsRaw.derFreq,
                                                              cmsStatsRaw.meanAnc, cmsStatsRaw.Pos ) ),
                                     returnInd = True, **Dict( 'cleanIHS cleanDer cleanAnc' ) )

    complikeNormedGlobal = normalize( complike, indClean )[0]
    complikeRatioNormedGlobal = normalize( complikeRatio, indClean )[0]
    
    complikeExp = exp( complike )
    complikeRatioExp = exp( complikeRatio )

    Columns = [ ind, indClean, complike, complikeRatio, complikeExp, complikeRatioExp,
                complikeNormedGlobal, complikeRatioNormedGlobal ] + cmsStatsRaw.columns() + statLikes

    result = DotData( names = ( 'ind', 'indClean', 'complike', 'complikeRatio', 'complikeExp', 'complikeRatioExp' )
                      + cmsStatsRaw.dtype.names + tuple([ stat + 'like' for stat in stats ]) + tuple([ stat + 'likeRatio' for stat in stats ]) + tuple([ stat + 'likeBin' for stat in stats ]),
                      Columns = [ ind, indClean, complike, complikeRatio, complikeExp, complikeRatioExp ] +
                      cmsStatsRaw.columns() + statLikes + statLikeRatios + statLikesBins )

    if False:
        result_idd = IDotData.fromDotData( result )
        def GetNormedBlocks():
            for chrom, r in result_idd.groupby( 'Chrom' ):

                r_dotData = DotData.fromIDotData( r )

                complikeNormedLocal = normalize( r_dotData.complike, r_dotData.indClean )[0]
                complikeRatioNormedLocal = normalize( r_dotData.complikeRatio, r_dotData.indClean )[0]

                yield IDotData( names = ( 'complikeNormedLocal', 'complikeRatioNormedLocal' ),
                                Columns = ( complikeNormedLocal, complikeRatioNormedLocal ) )

        local = IDotData.vstackFromIterable( GetNormedBlocks(),
                                             headings = ( 'complikeNormedLocal', 'complikeRatioNormedLocal' ) )

        rr = result_idd.hstack( local )
        rr.save( complikeFN )
    else: result.save( complikeFN )
    

def computeCMS_smlMem( Ddata, thinSfx, scenario, putativeMutPop, nreplicas = 100, likesTable = None,
                       likesTableSfx = '',
                       likesSfx = 'Newer', likesDir = None,
                       likesBins = None, likesPopSpecific = True,
                       stats = CMSBins.CMSstats, nonNormedStats = CMSBins.nonNormedStats,
                       nonNanStats = ( 'iHS', 'max_xpop' ),  #CMSBins.CMSstats,
                       nonNanMin = 1,
                       nearCausalSfx = '', freq = 'all',
                       complikeSfx = None, rawNormSfx = '', statsSfx = '',
                       normalizeWithinReplicas = False,
                       forceNormalizeWithinReplicas = False,
                       useRegionLikes = True,
                       cleanIHS = True, cleanDer = False, cleanAnc = False, keepCausal = False,
                       outputFN = None, useDefaultNumSnps = True,
                       pop2name = pop2name,
                       ihsSfx = '',
                       getio = None ):
    """Compute CMS scores for all snps in all replicas of scenario $scenario, using given likes tables.

    Params:

       complikeSfx - a suffix that identifies this version of complike scores.  This helps distinguish, e.g.,
          complike scores computed from the same data using different likes tables.
          Defaults to likesSfx.   We could just use likesSfx, but it does not capture the location (directory) of the
          likes tables -- different likes tables in different directories could have the same likes sfx.

    """

    #stats = ( 'max_xpop', )
    
    print 'in computeCMS: likesBins is ', likesBins, 'likesSfx is ', likesSfx, ' likesTable is ', likesTable

    if likesTable:
        likesTable = LoadLikesTable( **Dict( 'Ddata likesTable likesTableSfx' ) )
        likesBins = likesTable.likesBins
        normalizeWithinReplicas = coerceVal( likesTable.normalizeWithinReplicas )

    if forceNormalizeWithinReplicas: normalizeWithinReplicas = True
    
    if likesDir == None: likesDir = Ddata

    if putativeMutPop == None: putativeMutPop = scenario.mutPop

    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )


    cmsStatsNormFN = GetCreates( normalizeStatsWithinReplicas if normalizeWithinReplicas else normalizeStats,
                                 **Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx pop2name ihsSfx' ) ).containing( 'cmsStatsNorm' )
                                     
    
    cmsStatsRawFN = GetCreates( computeCMSstats, **Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx pop2name ihsSfx' ) )[0]

    hitLikesFN = likesTable.hitsLikes if likesTable else  \
        os.path.join( likesDir, 'likes',
                      AddFileSfx( 'hitsLikes.tsv', likesSfx, putativeMutPop if likesPopSpecific else None ) )
    missLikesFN = likesTable[ ('region' if useRegionLikes else 'miss') + 'Likes' ] if likesTable else \
        os.path.join( likesDir, 'likes',
                      AddFileSfx( ( 'region' if normalizeWithinReplicas and useRegionLikes else 'miss')
                                  + 'Likes.tsv',
                                  likesSfx, putativeMutPop if likesPopSpecific else None ) )

    complikeFN = outputFN or os.path.join( snpStatsDir,
                                           AddFileSfx( 'complike.tsv', statsSfx, putativeMutPop, complikeSfx, likesTableSfx, ihsSfx ) )

    if likesBins and not os.path.dirname( likesBins ): likesBins = os.path.join( likesDir, likesBins )
    
    if getio: return dict( depends_on = ( hitLikesFN, missLikesFN, cmsStatsNormFN, cmsStatsRawFN ) +
                           ( (likesBins,) if likesBins else () ),
                           creates = complikeFN,
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop, likesTableSfx ),
                           splitByCols = { cmsStatsRawFn: dict( keyCols = ( 'Chrom', 'Pos' ) ),
                                           cmsStatsNormFN: dict( keyCols = ( 'Chrom', 'Pos' ) ) } )

    hitLikes = DotData(SVPath = hitLikesFN )
    missLikes = DotData(SVPath = missLikesFN )

    cmsStatsNorm = IDotData( cmsStatsNormFN )
    cmsStatsRaw = IDotData( cmsStatsRawFN )


    oldStat = { 'iHS': 'ihs', 'StdDiff': 'StdDiff', 'meanFst': 'fst', 'freqDiff': 'freqDiff', 'max_xpop': 'max_xpop' }
    statNames = [ stat if stat in hitLikes.dtype.names else oldStat[ stat ] for stat in stats ]

    if not likesBins:
        stat_start, stat_end, stat_nbins = CMSBins.stat_start, CMSBins.stat_end, CMSBins.stat_nbin
    else:
        stat_start, stat_end, stat_nbins = LoadBins( likesBins )

    withSpecialBins = ( len( hitLikes[ 'iHS' ] ) > stat_nbins[ 'iHS' ] )

    # precompute the CLR for each bin for each stat

    logging.info( 'Counting rows..' )
    numSnps = 10000 if useDefaultNumSnps else cmsStatsRaw.numRows()
    dbg( 'numSnps' )
    logging.info( 'Precomputing bin stats...' )
    p_causal = 1./numSnps if numSnps else np.nan
    p_noncausal = 1 - p_causal
    dbg( 'p_causal p_noncausal' )
        
    CLR = [ None ] * len( stats )
    CLV = [ None ] * len( stats )

    stat_starts = [ stat_start[ stat ] for stat in statNames ]
    stat_ends = [ stat_end[ stat ] for stat in statNames ]
    stat_nbinss = [ stat_nbins[ stat ] for stat in statNames ]

    bin_sizes = [ float(stat_end[ stat ] - stat_start[ stat ]) / stat_nbins[ stat ] for stat in statNames ]

    for statNum, stat in enumerate( statNames ):
#        assert stat_nbins[ stat ] == len( hitLikes[stat] ) == len( missLikes[ stat ] )

        dbg( '"YYYYYY" stat hitLikes[stat] missLikes[stat]' )
    
        indNaN = hitLikes[ stat ] != 1e-10
        missingVal = np.log( np.min( hitLikes[stat][indNaN] / missLikes[stat][indNaN] ) )

        CLR[ statNum ] = [ ( np.log( hitLike / missLike ) if hitLike != 1e-10 else missingVal ) if hitLike != 0.0 else np.nan
                           for hitLike, missLike in zip( hitLikes[stat], missLikes[stat] ) ]

        CLV[ statNum ] = [ np.log( ( hitLike * p_causal ) / ( hitLike * p_causal + missLike*p_noncausal ) )
                           if numSnps and p_causal > 0.0 else np.nan
                           for hitLike, missLike in zip( hitLikes[stat], missLikes[stat] ) ]

    dbg( 'CLR CLV' )

    logging.info( 'Computing CMS...' )

    with IDotData.openForWrite( complikeFN,
                                headings = ( 'ind', 'indClean', 'complike', 'complikeRatio', 'complikeExp', 'complikeRatioExp' )
                                + cmsStatsRaw.headings + tuple([ stat + 'like' for stat in stats ]) +
                                tuple([ stat + 'likeRatio' for stat in stats ]) + tuple([ stat + 'likeBin' for stat in stats ]) ) \
                                as complikeOut:
    
        for rowNum, ( rRaw, rNorm ) in enumerate( itertools.izip( cmsStatsRaw, cmsStatsNorm ) ):

            if rRaw.Chrom >= nreplicas: break 

            statLikeRatios = []
            statLikes = []
            statLikesBins = []

            statNorms = [ ( rRaw if statName in nonNormedStats else rNorm )[ statName ] for statName in statNames ]
            nanReasons = [ rRaw[ statName + '_nanReason' ]  if CMSBins.stat_numSpecialBins[ statName ] else -1 for statName in statNames ]

            num_nonNans = np.sum( [ np.isfinite( ( rRaw if statName in nonNormedStats else rNorm )[ statName ] ) for statName in nonNanStats ] )

#            for statName, statNorm, st_start, st_binSize, st_nbins in zip( statNames, statNorms, stat_starts, bin_sizes, stat_nbinss ):
#                dbg( '"EEEEE" rowNum statName statNorm st_start st_binSize st_nbins' )
            
            bins = [ (st_nbins + nanReason) if withSpecialBins and nanReason != -1 else
                     ( clamp( int((statNorm - st_start ) / st_binSize ), 0, st_nbins-1 )
                       if np.isfinite( statNorm ) else np.nan )
                     for statNorm, nanReason, st_start, st_binSize, st_nbins in zip( statNorms, nanReasons, stat_starts,
                                                                                     bin_sizes, stat_nbinss ) ]

            CLRvals = [ CLRbinVals[ bin ] if not np.isnan( bin ) else np.nan  for CLRbinVals, bin in zip( CLR, bins ) ] 
            CLVvals = [ CLVbinVals[ bin ] if not np.isnan( bin ) else np.log( 1. / nbins )
                        for CLVbinVals, bin, nbins in zip( CLV, bins, stat_nbinss ) ]

            complikeRatio = np.sum( CLRvals )
            complike = np.sum( CLVvals )
            ind = ( num_nonNans >= nonNanMin )
            if not ind: complike = np.nan 

            indClean = ( ( ( not cleanIHS or not np.isnan( rRaw.StdDiff ) ) and
                           ( not cleanDer or rRaw.derFreq > .2 ) and
                           ( not cleanAnc or rRaw.meanAnc > .4 ) ) or 
                         ( keepCausal and rRaw.Pos == CAUSAL_POS ) )

            complikeOut.writeRecord( ind, indClean, complike, complikeRatio, exp( complike), exp( complikeRatio ),
                                     *( tuple( rRaw ) + tuple( CLVvals ) + tuple( CLRvals ) + tuple( bins ) ) )

def chooseCMS( Ddata, scenario, replicaTables, replicaCond2choice,
               complikeSfx, putativeMutPop = None, nreplicas = 100,
               likesTableSfx = '', thinSfx = '', normedLocal = False, getio = None ):
    """Construct a table of CMS scores by choosing, for SNPs in each replica, one of existing
    scores based on a replica condition.

    Params:
       replicaCond2choice - sequence of pairs of the form ( replicaCond, complikeSfx ), where for replicas matching
          replicaCond we would take the scores from complikeSfx.
    """

    if putativeMutPop is None: putativeMutPop = scenario.mutPop
    args = Dict( 'Ddata thinSfx scenario putativeMutPop' )
    complikeFunc = normalizeComplikeWithinReplicas if normedLocal else computeCMS
    complikeFNs_in = [ GetCreates( complikeFunc, complikeSfx = complikeSfx_in, **args )[0]
                       for replicaCond, complikeSfx_in in replicaCond2choice ]
    complikeFN = GetCreates( complikeFunc, complikeSfx = complikeSfx, **args )[0]

    if getio: return dict( depends_on = complikeFNs_in, creates = complikeFN,
                           args = Dict( 'replicaCond2choice scenario complikeSfx normedLocal' ) )

    replica2choice = {}
    
    for r in findReplicasMatchingConds( allScens = ( scenario, ), showVals = map( itemgetter( 0 ), replicaCond2choice ),
                                        **Dict( 'Ddata replicaTables' ) ):
        gotChoice = False
        for replicaNum, choice in enumerate( r[2:] ):
            if choice:
                replica2choice[ r.replicaNum ] = replicaNum
                gotChoice = True
                break
        assert gotChoice

    dbg( 'replica2choice' )
    iDotDatas = map( IDotData, complikeFNs_in ) 
    def makeResult():
        yield iDotDatas[0].headings
        for r in itertools.izip( *iDotDatas ):
            yield r[ replica2choice[ r[0].Chrom ] ]

    IDotData.fromFn( makeResult ).save( complikeFN )



def computeCMS( Ddata = None, thinSfx = '', scenario = None, putativeMutPop = None, nreplicas = 100, likesTable = None,
                likesTableNeut = None,
                likesTableSfx = '',
                likesSfx = 'Newer', likesDir = None,
                likesBins = None, likesPopSpecific = True,
                stats = CMSBins.CMSstats, nonNormedStats = CMSBins.nonNormedStats,
                nonNanStats = ( 'iHS', 'max_xpop' ),  #CMSBins.CMSstats,
                nearCausalSfx = '', freq = 'all',
                complikeSfx = None, rawNormSfx = '',
                normalizeWithinReplicas = False,
                forceNormalizeWithinReplicas = False,
                useRegionLikes = True,
                statsSfx = '',
                cleanIHS = True, cleanDer = False, cleanAnc = False, keepCausal = False,
                outputFN = None, useDefaultNumSnps = False, pop2name = pop2name,
                ihsSfx = '',
                nonNanMin = 1,  # currently unused
                cmsStatsRawFN = None, cmsStatsNormFN = None,
                getio = None ):
    """Compute CMS scores for all snps in all replicas of scenario $scenario, using given likes tables.

    Params:

       complikeSfx - a suffix that identifies this version of complike scores.  This helps distinguish, e.g.,
          complike scores computed from the same data using different likes tables.
          Defaults to likesSfx.   We could just use likesSfx, but it does not capture the location (directory) of the
          likes tables -- different likes tables in different directories could have the same likes sfx.

    """

    #stats = ( 'max_xpop', )
    
    print 'in computeCMS: likesBins is ', likesBins, 'likesSfx is ', likesSfx, ' likesTable is ', likesTable

    if not likesTableNeut: likesTableNeut = likesTable

    if likesTable:
        likesTable = LoadLikesTable( **Dict( 'Ddata likesTable likesTableSfx' ) )
        likesBins = likesTable.likesBins
        normalizeWithinReplicas = coerceVal( likesTable.normalizeWithinReplicas )

    if likesTableNeut:
        likesTableNeut = LoadLikesTable( **Dict( 'Ddata likesTableSfx', likesTable = likesTableNeut ) )
        likesBinsNeut = likesTable.likesBins
        normalizeWithinReplicasNeut = coerceVal( likesTableNeut.normalizeWithinReplicas )
        
        
    assert ( likesTable and cmsStatsRawFN and cmsStatsNormFN ) or ( Ddata and scenario )

    assert not likesTableNeut or ( np.all( likesBins == likesBinsNeut )
                                   and normalizeWithinReplicas == normalizeWithinReplicasNeut )
        
    if forceNormalizeWithinReplicas: normalizeWithinReplicas = True
    
    if likesDir == None: likesDir = Ddata

    if putativeMutPop == None: putativeMutPop = scenario.mutPop

    scenDir = scenario.scenDir() if scenario else 'unknown_scenDir'
    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenDir ) if not outputFN else None

    if not cmsStatsRawFN: cmsStatsRawFN = GetCreates( computeCMSstats, **Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx ihsSfx pop2name' ) )[0]
    
    if not cmsStatsNormFN: cmsStatsNormFN = AddFileSfx( cmsStatsRawFN, 'normedLocal' if normalizeWithinReplicas else 'normed' )

    hitLikesFN = likesTable.hitsLikes if likesTable else  \
        os.path.join( likesDir, 'likes',
                      AddFileSfx( 'hitsLikes.tsv', likesSfx, putativeMutPop if likesPopSpecific else None ) )
    missLikesFN = likesTableNeut[ ('region' if useRegionLikes else 'miss') + 'Likes' ] if likesTableNeut else \
        os.path.join( likesDir, 'likes',
                      AddFileSfx( ( 'region' if normalizeWithinReplicas and useRegionLikes else 'miss')
                                  + 'Likes.tsv',
                                  likesSfx, putativeMutPop if likesPopSpecific else None ) )

    complikeFN = outputFN or os.path.join( snpStatsDir, AddFileSfx( 'complike.tsv', putativeMutPop, complikeSfx, likesTableSfx, ihsSfx ) )

    if likesBins and not os.path.dirname( likesBins ): likesBins = os.path.join( likesDir, likesBins )

    if getio: return dict( depends_on = ( hitLikesFN, missLikesFN, cmsStatsNormFN, cmsStatsRawFN ) +
                           ( (likesBins,) if likesBins else () ),
                           creates = complikeFN,
                           splitByCols = { cmsStatsNormFN: dict( keyCols = ( 'Chrom', 'Pos' ) ),
                                           cmsStatsRawFN: dict( keyCols = ( 'Chrom', 'Pos' ) ) },
                           name = 'computeCMS', 
                           attrs = Dict( 'putativeMutPop', scenario = scenDir, piperun_short = True ),
                           mediumRuleNameSfx = ( scenDir, putativeMutPop, likesTableSfx ) )

    hitLikes = DotData(SVPath = hitLikesFN )
    missLikes = DotData(SVPath = missLikesFN )

    cmsStatsNorm = DotData( cmsStatsNormFN )
    cmsStatsRaw = DotData( cmsStatsRawFN )

    oldStat = { 'iHS': 'ihs', 'StdDiff': 'StdDiff', 'meanFst': 'fst', 'freqDiff': 'freqDiff', 'max_xpop': 'max_xpop' }
    statNames = [ stat if stat in hitLikes.dtype.names else oldStat[ stat ] for stat in stats ]

    if not likesBins:
        stat_start, stat_end, stat_nbins = CMSBins.stat_start, CMSBins.stat_end, CMSBins.stat_nbin
    else:
        stat_start, stat_end, stat_nbins = LoadBins( likesBins )

    withSpecialBins = ( len( hitLikes[ 'iHS' ] ) > stat_nbins[ 'iHS' ] )

    #
    # Precompute the CLR for each bin for each stat
    #

    logging.info( 'Counting rows..' )
    numSnps = 10000 if useDefaultNumSnps else cmsStatsRaw.numRows()
    dbg( 'numSnps' )
    logging.info( 'Precomputing bin stats...' )
    p_causal = 1./numSnps if numSnps else np.nan
    p_noncausal = 1 - p_causal
    dbg( 'p_causal p_noncausal' )
        
    CLR = [ None ] * len( stats )
    CLV = [ None ] * len( stats )

    stat_starts = [ stat_start[ stat ] for stat in statNames ]
    stat_ends = [ stat_end[ stat ] for stat in statNames ]
    stat_nbinss = [ stat_nbins[ stat ] for stat in statNames ]

    bin_sizes = [ float(stat_end[ stat ] - stat_start[ stat ]) / stat_nbins[ stat ] for stat in statNames ]

    for statNum, stat in enumerate( statNames ):
#        assert stat_nbins[ stat ] == len( hitLikes[stat] ) == len( missLikes[ stat ] )

        dbg( '"YYYYYY" stat hitLikes[stat] missLikes[stat]' )
    
        indNaN = hitLikes[ stat ] != 1e-10
        missingVal = np.log( np.min( hitLikes[stat][indNaN] / missLikes[stat][indNaN] ) ) if np.any( indNaN ) else -1e10

        CLR[ statNum ] = np.array([ ( np.log( hitLike / missLike ) if hitLike != 1e-10 else missingVal ) if hitLike != 0.0 else np.nan
                                    for hitLike, missLike in zip( hitLikes[stat], missLikes[stat] ) ] + [ np.nan ])

        CLV[ statNum ] = np.array([ np.log( ( hitLike * p_causal ) / ( hitLike * p_causal + missLike*p_noncausal ) )
                                    if numSnps and p_causal > 0.0 else np.nan
                                    for hitLike, missLike in zip( hitLikes[stat], missLikes[stat] ) ] +
                                  [ np.log( 1. / stat_nbinss[ statNum ] ) ])

    logging.info( 'Computing CMS...' )

    likesCols = []
    likesRatioCols = []
    likesBinsCols = []
    for stat, st_start, st_end, st_nbins, st_binSize, CLRhere, CLVhere in \
            zip( stats, stat_starts, stat_ends, stat_nbinss, bin_sizes, CLR, CLV ):
        statNorm = ( cmsStatsRaw if stat in nonNormedStats else cmsStatsNorm )[ stat ]
        nanReasons = cmsStatsRaw[ stat + '_nanReason' ] if CMSBins.stat_numSpecialBins[ stat ] else np.repeat( -1, len( cmsStatsRaw ) )

        bins = np.where( nanReasons == -1,
                         np.where( np.isfinite( statNorm ),
                                   np.clip( ( ( statNorm - st_start ) / st_binSize ).astype( np.int16 ), 0, st_nbins-1 ),
                                   len( CLVhere ) - 1 ),
                         st_nbins + nanReasons )

        dbg( '"DDDDDDDDD" stat st_start st_end st_nbins st_binSize statNorm nanReasons bins' )

        likesBinsCols.append( bins )
        likesRatioCols.append( CLRhere[ bins ] )
        likesCols.append( CLVhere[ bins ] )

    complike = reduce( operator.add, likesCols )
    complikeRatio = reduce( operator.add, likesRatioCols )
    complikeExp = np.exp( complike )
    complikeRatioExp = np.exp( complikeRatio )

    indClean = clean_hapmap2_region( Z = DotData( names = ( 'complike', 'StdDiff', 'derFreq', 'meanAnc', 'Pos' ),
                                                  Columns = ( complike, cmsStatsRaw.StdDiff, cmsStatsRaw.derFreq,
                                                              cmsStatsRaw.meanAnc, cmsStatsRaw.Pos ) ),
                                     returnInd = True, **Dict( 'cleanIHS cleanDer cleanAnc' ) )


    DotData( names = ( 'indClean', 'complike', 'complikeRatio', 'complikeExp', 'complikeRatioExp' )
             + cmsStatsRaw.headings + tuple([ stat + 'like' for stat in stats ]) +
             tuple([ stat + 'likeRatio' for stat in stats ]) + tuple([ stat + 'likeBin' for stat in stats ]),
             Columns = [ indClean, complike, complikeRatio, complikeExp, complikeRatioExp ] +
             [ cmsStatsRaw[ h ] for h in cmsStatsRaw.headings ] +
             likesCols + likesRatioCols + likesBinsCols ).saveToSV( complikeFN )
    
def DefineRulesTo_chooseCMS( pr, Ddata, replicaTables, replicaCond2choice,
                             complikeSfx, nreplicas = 100 ):
    for scenario in GetSelectionScenarios():
        for normedLocal in ( False, True ):
            pr.addInvokeRule( invokeFn = chooseCMS,
                              invokeArgs = Dict( 'Ddata scenario replicaTables replicaCond2choice complikeSfx '
                                                 'normedLocal nreplicas' ) )



def plotCMS( Ddata, thinSfx, scenario, nreplicas, putativeMutPop = None, complikeSfx = '',
             likesTableSfx = '', selpos = 500000, includeSpatialLoc = True,
             whichSpatialLoc = 'Spline', fromReplica = None, toReplica = None,
             xAxisGd = True, getio = None ):
    """Plot CMS scores across each replica, against genetic or physical distance.
    Show the location of the causal SNP for selection scenarios.  Show the results of spatial localization.
    
    Params:
    
    includeSpatialLoc - if true, the results of spatial localization will be shown on the plot.
    xAxisGd - use genetic distance for the X axis if True, physical distance if False.
    """

    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats'+ thinSfx, scenario.scenDir() )
    if putativeMutPop == None: putativeMutPop = scenario.mutPop
    sfxs = ( putativeMutPop, complikeSfx, likesTableSfx )
    complikeFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', *sfxs ) )

    xAxis = 'gd' if xAxisGd else 'bp'
    xCol = 'gdPos' if xAxisGd else 'Pos'
    xColOther = 'gdPos' if not xAxisGd else 'Pos'
    xLabel = 'genetic map position (cM)' if xAxisGd else 'position (bp)'
    
    plotFilesFNs = [ os.path.join( complikeFN[:-len('.data/')] + Sfx( xAxis + 'Plots',
                                                                      whichSpatialLoc if includeSpatialLoc else None ),
                                   'cmsPlot_%s_%d.svg' % ( xAxis, replicaNum ) ) \
                         if ( fromReplica is None or replicaNum >= fromReplica ) \
                         and ( toReplica is None or replicaNum <= toReplica ) \
                         else None
                     for replicaNum in range( nreplicas ) ]

    intervalsListFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervals' + whichSpatialLoc + 'List.tsv', *sfxs ) )
    gdMapFN = os.path.join( snpStatsDir, 'gdMap.tsv' )
    causalGdPosFN = os.path.join( replicaStatsDir, 'causalGdPos.tsv' )

    if getio:
        return dict( depends_on = ( complikeFN, gdMapFN ) +
                     ( ( causalGdPosFN, ) if not scenario.isNeutral() else () ) +
                     ( ( intervalsListFN, ) if includeSpatialLoc else () ),
                     creates = tuple( filter( None, plotFilesFNs ) ),
                     name = 'plotCMS' + Sfx( xAxis, whichSpatialLoc if includeSpatialLoc else None ),
                     mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop if scenario.isNeutral() else None,
                                           xAxis ) )

    if scenario.isNeutral(): selpos = None
    
    usedReplicas = zeros( nreplicas, dtype = bool )
    for ( replicaNum, cmsForReplica ), plotFileFN, ( replicaNum2, intervals ) in \
            itertools.izip( IDotData( complikeFN ).groupby( 'Chrom' ),
                            plotFilesFNs,
                            IDotData( intervalsListFN ).groupby( 'replicaNum' ) if includeSpatialLoc \
                                else itertools.repeat( ( None, None ) ) ):
        replicaNum = int( replicaNum )
        if includeSpatialLoc:
            replicaNum2 = int( replicaNum2 )
            assert replicaNum2 == replicaNum

        if fromReplica is not None and replicaNum < fromReplica: continue
        if toReplica is not None and replicaNum > toReplica: break
        assert plotFileFN
            
        minCMS = min( cmsForReplica.complikeExp )
        maxCMS = max( cmsForReplica.complikeExp )
            
        dbg( '"used" replicaNum' )
        pp.figure( figsize = ( 10, 10 ) )
        pp.axes( ( .15, .1, .85, .65 ) )

        # for highest-scoring points, add info as URLs


        def locrtn():
            X = cmsForReplica
            ind = X.complikeExp.argsort()
            lik = X.complikeExp[ind]
            maxlik = mean(lik[-5:])
            #maxlik = mean(lik[-5:])
            like = X.complikeExp / ( maxlik if maxlik else np.nan )

            pp.hlines( maxlik, min( cmsForReplica[ xCol ] ), max( cmsForReplica[ xCol ] ),
                       color = 'g', linestyles = 'dashed' )
            
            foundLike = 0
            margin = ( maxCMS - minCMS ) / 100
            for pos, absLikeHere, likeHere in itertools.izip( cmsForReplica[ xCol ], X.complikeExp, like ):
                if likeHere > .2:
                    pp.vlines( [ pos ], absLikeHere - margin, absLikeHere + margin, color = 'g',
                               linewidth = 3.0 ).set_urls( [ 'like:_%.2f' % likeHere ] )
                    foundLike += 1
            dbg( '"FFFF" foundLike' )
        locrtn()
        
        pp.plot( cmsForReplica[ xCol ], cmsForReplica.complikeExp, 'k.' )

        if selpos:
            # determine the genetic map position of the causal SNP

            causalLine = cmsForReplica[ cmsForReplica.Pos == selpos ]
            if causalLine:
                pp.plot( causalLine[ xCol ], causalLine.complikeExp,
                         'ro' )

                pp.hlines( causalLine.complikeExp, min( cmsForReplica[ xCol ] ), max( cmsForReplica[ xCol ] ),
                           color = 'r', linestyles = 'dashed' )
                         
            causalPos_gd = next( iter( causalLine[ xCol ] ) ) if causalLine else \
              interpolate.interp1d( cmsForReplica.Pos, cmsForReplica[ xCol ] )( selpos )

            dbg( 'replicaNum causalPos_gd' )

            pp.vlines( causalPos_gd if xAxisGd else selpos, 0, maxCMS * 1.03,
                       color = 'r', linestyles = 'dashed' ).set_urls( [ 'gd:_%.5f' % causalPos_gd ] ) 

        if includeSpatialLoc:

            gdTot = 0.0
            bpTot = 0
            
            for interval in intervals:
                aLine = pp.hlines( maxCMS * 1.03, interval[ xAxis + 'From' ], interval[ xAxis + 'To' ], colors = 'b' )
                aLine.set_urls( [ 'gd:_%.3f_bp:_%d' % ( interval.gdSize, interval.bpSize ) ] )
                gdTot += interval.gdSize
                bpTot += interval.bpSize

        pp.xlabel( xLabel )
        pp.ylabel( 'CMS' )
        title = 'replica %s of scenario %s' % ( replicaNum, scenario.scenDir() )
        if includeSpatialLoc: title += '\ntot_gd: %.3f tot_bp: %d' % ( gdTot, bpTot )
        #pp.title( title )
        pp.figtext( .5, .97, title, ha = 'center' )
        #
        # Add a second axes on top, showing the alternate position coordinates:
        # if the main x-axis (at the bottom) is in bp, the one at the top will show cM,
        # and vice versa.
        #

        minOther = min( cmsForReplica[ xColOther ] )
        maxOther = max( cmsForReplica[ xColOther ] )
        otherTickPositionInterval = ( maxOther - minOther ) / 10

        minOrig = min( cmsForReplica[ xCol ] )
        maxOrig = max( cmsForReplica[ xCol ] )
        origTickPositionInterval = ( maxOrig - minOrig ) / 20
        
        otherTickPositions = []
        otherTickLabels = []
        origTickPositions = []
        for r in cmsForReplica:
            if not otherTickPositions or r[ xColOther ] - otherTickPositions[-1] > otherTickPositionInterval:
                otherTickPositions.append( r[ xColOther ] )
                origTickPositions.append( r[ xCol ] )
                otherTickLabels.append( str( otherTickPositions[ -1 ] ) )

        ax1 = pp.gca()

        for t in ax1.get_xticklabels():
            t.set_rotation( 'vertical' )
        
        ax2 = pp.twiny()
        ax2.set_xlim( ax1.get_xlim()  )
        ax1.callbacks.connect( 'xlim_changed', lambda ax: ax2.set_xlim( ax.get_xlim() ) )

        dbg( 'zip(origTickPositions,otherTickPositions)' ) 
        ax2.set_xticks( origTickPositions )
        ax2.set_xticklabels( otherTickLabels )

        for t in ax2.get_xticklabels():
            t.set_rotation( 'vertical' )
            t.set_size( 6 )
        
        dbg( 'ax2.get_xlim()' )

        dbg( '"SAVING_TO" plotFilesFNs[replicaNum]' )

        pp.savefig( plotFileFN )
        usedReplicas[ replicaNum ] = True

    numMissingReplicas = 0
    for replicaNum, plotFileFN in zip( range( nreplicas ), plotFilesFNs ):
        if plotFileFN and not usedReplicas[ replicaNum ]:
            pp.figure()
            pp.xlabel( xLabel )
            pp.ylabel( 'CMS' )
            pp.title( 'MISSING replica %s of scenario %s' % ( replicaNum, scenario.scenDir() ) )
            pp.savefig( plotFileFN )
            numMissingReplicas += 1

    dbg( 'numMissingReplicas' )


def DefineRulesTo_plotCMS( pr, Ddata, thinSfx, mutAges, mutPops, mutFreqs, nreplicas,
                           complikeSfx = '', likesTableSfx = '', includeSpatialLoc = False,
                           fromReplica = None, toReplica = None, xAxisGd = True ):
    """Define rules to plot CMS scores for given replicas."""

    for scenario in GetScenarios( **Dict( 'mutAges mutPops mutFreqs' ) ):
        for putativeMutPop in ( mutPops if scenario.isNeutral() else (scenario.mutPop,) ):
            pr.addInvokeRule( invokeFn = plotCMS,
                              invokeArgs = Dict( 'Ddata thinSfx scenario nreplicas '
                                                 'putativeMutPop complikeSfx likesTableSfx '
                                                 'includeSpatialLoc fromReplica toReplica xAxisGd' ) )
    
def plotCMSstat( Ddata, thinSfx, scenario, nreplicas, putativeMutPop = None,
                 snpStat = 'complikeExp',  complikeSfx = '',
                 likesTableSfx = '', selpos = 500000, includeSpatialLoc = True,
                 whichSpatialLoc = 'Spline', 
                 xAxisGd = True, fromReplica = None, toReplica = None, whichReplicas = None,
                 figSize = ( 10, 10 ), nonNanStats = 'ALL', plotSfx = '', getio = None ):
    """Plot CMS scores across each replica, against genetic or physical distance, for specified complikeSfxs.
    Show the location of the causal SNP for selection scenarios.  Show the results of spatial localization.

    Params:

       includeSpatialLoc - if true, the results of spatial localization will be shown on the plot.
       xAxisGd - use genetic distance for the X axis if True, physical distance if False.
    """

    if isinstance( scenario, types.StringTypes ): scenario = Scenario.fromString( scenario )
    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats'+ thinSfx, scenario.scenDir() )
    if putativeMutPop == None: putativeMutPop = scenario.mutPop
    complikeFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', 'normedLocal',
                                                        putativeMutPop, complikeSfx, likesTableSfx,
                                                        'nonNan', *MakeSeq( nonNanStats ) ) )

    xAxis = 'gd' if xAxisGd else 'bp'
    xCol = 'gdPos' if xAxisGd else 'Pos'
    xColOther = 'gdPos' if not xAxisGd else 'Pos'
    xLabel = 'genetic map position (cM)' if xAxisGd else 'position (bp)'

    dbg( 'nreplicas fromReplica toReplica whichReplicas' )
    plotFilesFNs = dict( [( replicaNum,
                            os.path.join( complikeFN[:-len('.data/')] + Sfx( xAxis + 'Plots',
                                                                             whichSpatialLoc
                                                                             if includeSpatialLoc else None,
                                                                             plotSfx ),
                                          'cmsPlot_%s_%s_%s_%d.svg' % ( snpStat, complikeSfx, xAxis, replicaNum ) ) )
                          for replicaNum in range( nreplicas ) if ( fromReplica is None or replicaNum >= fromReplica )
                          and ( toReplica is None or replicaNum <= toReplica ) and
                          ( whichReplicas is None or replicaNum in whichReplicas )  ] )
    assert( plotFilesFNs )
                          
    intervalsListFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervals' + whichSpatialLoc + 'List.tsv',
                                                                 putativeMutPop, complikeSfx, likesTableSfx ) ) 
    gdMapFN = os.path.join( snpStatsDir, 'gdMap.tsv' )
    causalGdPosFN = os.path.join( replicaStatsDir, 'causalGdPos.tsv' )

    if getio:
        return dict( depends_on = ( complikeFN, gdMapFN ) +
                     ( ( causalGdPosFN, ) if not scenario.isNeutral() else () ) +
                     ( ( intervalsListFN, ) if includeSpatialLoc else () ),
                     creates = tuple( plotFilesFNs.values() ),
                     name = 'plotCMS' + Sfx( xAxis, whichSpatialLoc if includeSpatialLoc else None ),
                     mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop if scenario.isNeutral() else None,
                                           xAxis ) )

    if scenario.isNeutral(): selpos = None
    
    usedReplicas = set()

    if whichReplicas:
        minReplica = min( whichReplicas )
        maxReplica = max( whichReplicas )

    pp.figure( figsize = figSize )
    for ( replicaNum, cmsForReplica ), ( replicaNum2, intervals ) in \
            itertools.izip( IDotData( complikeFN ).groupby( 'Chrom' ),
                            IDotData( intervalsListFN ).groupby( 'replicaNum' ) if includeSpatialLoc \
                                else itertools.repeat( ( None, None ) ) ):
        replicaNum = int( replicaNum )
        if includeSpatialLoc:
            replicaNum2 = int( replicaNum2 )
            assert replicaNum2 == replicaNum

        if fromReplica is not None and replicaNum < fromReplica: continue
        if toReplica is not None and replicaNum > toReplica: break
        if whichReplicas is not None:
            if replicaNum < minReplica: continue
            if replicaNum > maxReplica: break
            if replicaNum not in whichReplicas: continue
            
        minCMS = min( cmsForReplica[ snpStat ] )
        maxCMS = max( cmsForReplica[ snpStat ] )
            
        dbg( '"used" replicaNum' )
        pp.axes( ( .15, .1, .85, .65 ) )

        # for highest-scoring points, add info as URLs


        def locrtn():
            X = cmsForReplica
            ind = X[ snpStat ].argsort()
            lik = X[ snpStat ][ind]
            maxlik = mean(lik[-5:])
            #maxlik = mean(lik[-5:])
            like = X[ snpStat ] / maxlik

            pp.hlines( maxlik, min( cmsForReplica[ xCol ] ), max( cmsForReplica[ xCol ] ),
                       color = 'g', linestyles = 'dashed' )
            
            foundLike = 0
            margin = ( maxCMS - minCMS ) / 100
            for pos, absLikeHere, likeHere in itertools.izip( cmsForReplica[ xCol ], X[ snpStat ], like ):
                if likeHere > .2:
                    pp.vlines( [ pos ], absLikeHere - margin, absLikeHere + margin, color = 'g',
                               linewidth = 3.0 ).set_urls( [ 'like:_%.2f' % likeHere ] )
                    foundLike += 1
            dbg( '"FFFF" foundLike' )
        locrtn()
        
        pp.plot( cmsForReplica[ xCol ], cmsForReplica[ snpStat ], 'k.' )

        if selpos:
            # determine the genetic map position of the causal SNP

            causalLine = cmsForReplica[ cmsForReplica.Pos == selpos ]
            if causalLine:
                pp.plot( causalLine[ xCol ], causalLine[ snpStat ],
                         'ro' )

                pp.hlines( causalLine[ snpStat ], min( cmsForReplica[ xCol ] ), max( cmsForReplica[ xCol ] ),
                           color = 'r', linestyles = 'dashed' )
                         
            causalPos_gd = next( iter( causalLine[ xCol ] ) ) if causalLine else \
              interpolate.interp1d( cmsForReplica.Pos, cmsForReplica[ xCol ] )( selpos )

            dbg( 'replicaNum causalPos_gd' )

            pp.vlines( causalPos_gd if xAxisGd else selpos, 0, maxCMS * 1.03,
                       color = 'r', linestyles = 'dashed' ).set_urls( [ 'gd:_%.5f' % causalPos_gd ] ) 

        if includeSpatialLoc:

            gdTot = 0.0
            bpTot = 0
            
            for interval in intervals:
                aLine = pp.hlines( maxCMS * 1.03, interval[ xAxis + 'From' ], interval[ xAxis + 'To' ], colors = 'b' )
                aLine.set_urls( [ 'gd:_%.3f_bp:_%d' % ( interval.gdSize, interval.bpSize ) ] )
                gdTot += interval.gdSize
                bpTot += interval.bpSize

        pp.xlabel( xLabel )
        pp.ylabel( snpStat )
        title = 'replica %s of scenario %s%s' % ( replicaNum, scenario.scenDir(),
                                                  '' if not complikeSfx else ' (%s)' % complikeSfx )
        if includeSpatialLoc: title += '\ntot_gd: %.3f tot_bp: %d' % ( gdTot, bpTot )
        #pp.title( title )
        pp.figtext( .5, .97, title, ha = 'center' )
        #
        # Add a second axes on top, showing the alternate position coordinates:
        # if the main x-axis (at the bottom) is in bp, the one at the top will show cM,
        # and vice versa.
        #

        minOther = min( cmsForReplica[ xColOther ] )
        maxOther = max( cmsForReplica[ xColOther ] )
        otherTickPositionInterval = ( maxOther - minOther ) / 10

        minOrig = min( cmsForReplica[ xCol ] )
        maxOrig = max( cmsForReplica[ xCol ] )
        origTickPositionInterval = ( maxOrig - minOrig ) / 20
        
        otherTickPositions = []
        otherTickLabels = []
        origTickPositions = []
        for r in cmsForReplica:
            if not otherTickPositions or r[ xColOther ] - otherTickPositions[-1] > otherTickPositionInterval:
                otherTickPositions.append( r[ xColOther ] )
                origTickPositions.append( r[ xCol ] )
                otherTickLabels.append( str( otherTickPositions[ -1 ] ) )

        ax1 = pp.gca()

        for t in ax1.get_xticklabels():
            t.set_rotation( 'vertical' )
        
        ax2 = pp.twiny()
        ax2.set_xlim( ax1.get_xlim()  )
        ax1.callbacks.connect( 'xlim_changed', lambda ax: ax2.set_xlim( ax.get_xlim() ) )

        dbg( 'zip(origTickPositions,otherTickPositions)' ) 
        ax2.set_xticks( origTickPositions )
        ax2.set_xticklabels( otherTickLabels )

        for t in ax2.get_xticklabels():
            t.set_rotation( 'vertical' )
            t.set_size( 6 )
        
        dbg( 'ax2.get_xlim()' )

        dbg( '"SAVING_TO" plotFilesFNs[replicaNum]' )

        pp.savefig( plotFilesFNs[ replicaNum ] )
        usedReplicas.add( replicaNum )

    numMissingReplicas = 0
    for replicaNum in range( nreplicas ):
        if fromReplica is not None and replicaNum < fromReplica: continue
        if toReplica is not None and replicaNum > toReplica: break
        if whichReplicas is not None:
            if replicaNum < minReplica: continue
            if replicaNum > maxReplica: break
            if replicaNum not in whichReplicas: continue
        
        if replicaNum not in usedReplicas:
            pp.figure()
            pp.xlabel( xLabel )
            pp.ylabel( 'CMS' )
            pp.title( 'MISSING replica %s of scenario %s' % ( replicaNum, scenario.scenDir() ) )
            pp.savefig( plotFilesFNs[ replicaNum ] )
            numMissingReplicas += 1

    dbg( 'numMissingReplicas' )
    



    
def getBinPlotFilesFNs( Ddata, thinSfx, scenario, putativeMutPop, complikeSfx, likesTableSfx = ''):
    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats'+ thinSfx, scenario.scenDir() )
    if putativeMutPop == None: putativeMutPop = scenario.mutPop
    sfxs = ( putativeMutPop, complikeSfx, likesTableSfx )
    
    plotFilesFNs = [ os.path.join( complikeFN[:-len('.data/')] + '_gdPlots2', 'cmsPlot_gd2_%d.svg' % replicaNum )
                     for replicaNum in range( nreplicas ) ]

    return plotFilesFNs

    
def plotCMS_gd2( Ddata, thinSfx, scenario, putativeMutPop, nreplicas, complikeSfx = '',
                 likesTableSfx = '', selpos = 500000, getio = None ):
    """Plot CMS scores across each replica, against genetic distance; also plot various
    auxiliary/debugging information."""

    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats'+ thinSfx, scenario.scenDir() )
    if putativeMutPop == None: putativeMutPop = scenario.mutPop
    sfxs = ( putativeMutPop, complikeSfx, likesTableSfx )
    complikeFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', *sfxs ) )

    plotFilesFNs = [ os.path.join( complikeFN[:-len('.data/')] + '_gdPlots2', 'cmsPlot_gd2_%d.svg' % replicaNum )
                     for replicaNum in range( nreplicas ) ]

    intervalsListFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervalsSplineList.tsv', *sfxs ) )
    intervalsStatsFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervalsSplineStats.tsv', *sfxs ) )

    posteriorSplineFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervalsSplineSpline.tsv', *sfxs ) )
    binInfoFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervalsSplineBinInfo.tsv', *sfxs ) )
    
    if getio:
        return dict( depends_on = ( complikeFN, intervalsListFN, intervalsStatsFN, posteriorSplineFN, binInfoFN ),
                     creates = plotFilesFNs,
                     mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop if scenario.isNeutral() else None ) )

    if scenario.isNeutral(): selpos = None
    
    usedReplicas = zeros( nreplicas, dtype = bool )
    for ( replicaNum, cmsForReplica ), plotFileFN, ( replicaNum2, intervals ), \
            ( replicaNum3, posteriorSpline ), ( replicaNum4, binInfo ), intervalsStats in \
            itertools.izip( IDotData( complikeFN ).groupby( 'Chrom' ),
                            plotFilesFNs,
                            IDotData( intervalsListFN ).groupby( 'replicaNum' ),
                            IDotData( posteriorSplineFN ).groupby( 'replicaNum' ),
                            IDotData( binInfoFN ).groupby( 'replicaNum' ),
                            IDotData( intervalsStatsFN ) ):
        replicaNum, replicaNum2, replicaNum3, replicaNum4 = \
            map( int, ( replicaNum, replicaNum2, replicaNum3, replicaNum4 ) )

        assert replicaNum2 == replicaNum == replicaNum3 == replicaNum4 == intervalsStats.replicaNum
            
        dbg( '"used" replicaNum' )
        pp.figure()
        pp.subplot(211)

        theLabels = []
        theHandles = []
        
        theHandles.append( pp.plot( cmsForReplica[ xCol ], cmsForReplica.complikeExp, 'k.')[0] )
        theLabels.append( 'CMS values' )
        minCMS = min( cmsForReplica.complikeExp )
        maxCMS = max( cmsForReplica.complikeExp )
        if selpos:
            # determine the genetic map position of the causal SNP

            causalLine = cmsForReplica[ cmsForReplica.Pos == selpos ]
            if causalLine:
                pp.plot( causalLine[ xCol ], causalLine.complikeExp, 'ro' )

                pp.hlines( causalLine.complikeExp, min( cmsForReplica[ xCol ] ), max( cmsForReplica[ xCol ] ),
                           color = 'r', linestyles = 'dashed' )
                         
            causalPos_gd = next( iter( causalLine.gdPos ) ) if causalLine else \
              interpolate.interp1d( cmsForReplica.Pos, cmsForReplica.gdPos )( selpos )

            dbg( 'replicaNum causalPos_gd' )

            pp.vlines( causalPos_gd, 0, maxCMS,
                       color = 'r', linestyles = 'dashed' )

        theHandles.append( pp.plot( posteriorSpline.gdPos, posteriorSpline.complikeExp, 'c-' )[0] )
        theLabels.append( 'posterior spline' )

        aLines = pp.hlines( [ maxCMS * 1.03] * len( intervals ), intervals.gdFrom, intervals.gdTo, colors = 'b', label = 'confidence intervals' )
        aLines.set_urls( [ 'gd:_%.3f_bp:_%d_nbins:_%d_area:_%.3f_binsMax:_%.3f_binsMaxFrac:_%.3f' % ( interval.gdSize, interval.bpSize, interval.binsInInterval,
                                                                      interval.binsArea, interval.binsMax, interval.binsMaxFrac ) for interval in intervals ] )
        gdTot = sum( intervals.gdSize )
        bpTot = sum( intervals.bpSize )

        pp.xlabel( 'genetic map position (cM)' )
        pp.ylabel( 'CMS' )
        title = 'replica %s of scenario %s' % ( replicaNum, scenario.scenDir() )
        title += '\ntot_gd: %.3f tot_bp: %d' % ( gdTot, bpTot )
        pp.title( title )

        fig = pp.gcf()
        ax1 = fig.add_subplot(212)
        ax1.set_xlabel( 'genetic map position (cM)' )
        theHandles.append( ax1.plot( binInfo.binCenters, binInfo.binIntegralNormed, 'bo' )[0] )
        theLabels.append( 'normed bin integral' )
        ax1.set_ylabel( 'normed integral', color = 'b' )
        for t1 in ax1.get_yticklabels(): t1.set_color('b')

        ax2 = ax1.twinx()
        ax2.set_ylabel( 'binRank', color = 'r' )
        theHandles.append( ax2.plot( binInfo.binCenters, -binInfo.binRank, 'g.' )[0] )
        theLabels.append( 'bin rank' )
        for t1 in ax2.get_yticklabels(): t1.set_color('r')

        pp.figlegend( loc = 'lower center', labels = theLabels, handles = theHandles )
        
        pp.savefig( plotFilesFNs[ replicaNum ] )
        
        usedReplicas[ replicaNum ] = True

    numMissingReplicas = 0
    for replicaNum in range( nreplicas ):
        if not usedReplicas[ replicaNum ]:
            pp.figure()
            pp.xlabel( 'genetic map position (cM)' )
            pp.ylabel( 'CMS' )
            pp.title( 'MISSING replica %s of scenario %s' % ( replicaNum, scenario.scenDir() ) )
            pp.savefig( plotFilesFNs[ replicaNum ] )
            numMissingReplicas += 1

    dbg( 'numMissingReplicas' )
    

    
            
def normalizeComplikeRatio( Ddata, thinSfx, scenario, putativeMutPop,
                            complikeSfx = None,
                            DdataNeutral = None, getio = None ):
    """Normalize the complike ratio scores."""

    if not DdataNeutral: DdataNeutral = Ddata
    
    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    complikeFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', putativeMutPop, complikeSfx ) )

    complikeWithLogCLRFN = os.path.join( snpStatsDir,
                                         AddFileSfx( 'complikeWithLogCLR.data/', putativeMutPop, complikeSfx ) )
    
    snpStatsDirNeutral = os.path.join( DdataNeutral, 'snpStats' + thinSfx, NeutralScenario().scenDir() )
    complikeRatioMeanStdFN = os.path.join( snpStatsDirNeutral, AddFileSfx( 'complikeRatioMeanStd.tsv',
                                                                           putativeMutPop, complikeSfx ) )

    depends_on = [ complikeFN ]
    creates = [ complikeWithLogCLRFN ]
    ( creates if scenario.isNeutral() else depends_on ).append( complikeRatioMeanStdFN )
    if getio: return dict( depends_on = depends_on, creates = creates, attrs = Dict( 'scenario putativeMutPop' ),
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ),
                           name = 'normalizeComplikeRatio' + ( '_neutral' if scenario.isNeutral() else '' ) )

    complike = DotData( Path = complikeFN )
    if scenario.isNeutral():
        logCLR, meanCL, stdCL = normalize( complike.complikeRatio )
        DotData( names = ( 'mean', 'std' ), Records = ( ( meanCL, stdCL ), ) ).saveToSV( complikeRatioMeanStdFN )
    else:
        meanCL, stdCL = DotData( SVPath = complikeRatioMeanStdFN )[0]
        logCLR = ( complike.complikeRatio - meanCL ) / stdCL

    complike.hstack( DotData( names = ( 'logCLR', ), Columns = ( logCLR, ) ) ).save( complikeWithLogCLRFN )


def normalizeComplikeWithinReplicas( Ddata, thinSfx, scenario, putativeMutPop, complikeSfx = None,
                                     likesTableSfx = None,
                                     nonNanStats = 'ALL', nonNanMin = 1, statsSfx = '', ihsSfx = '',
                                     pop2name = pop2name,
                                     getio = None ):
    """Normalize the complike scores within replicas."""

    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )

    complikeFN = GetCreates( computeCMS, **Dict( 'Ddata thinSfx scenario putativeMutPop complikeSfx likesTableSfx statsSfx ihsSfx pop2name' ) )[0]

    complikeNormedLocalFN = os.path.join( snpStatsDir,
                                          AddFileSfx( 'complike.data/', 'normedLocal', statsSfx, ihsSfx, putativeMutPop, complikeSfx,
                                                      'nonNan', nonNanMin, *MakeSeq( nonNanStats ) ) )
    
    if getio: return dict( depends_on = complikeFN, creates = complikeNormedLocalFN,
                           attrs = Dict( 'putativeMutPop', scenario = scenario.scenDir() ),
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop, complikeSfx ) )

    cmsScores = IDotData( complikeFN )
    if nonNanStats == 'ALL': nonNanStats = cmsScores.headings

    #cmsScores = cmsScores.filter( lambda r: np.sum( [ isfinite( r[c] ) for c in nonNanStats ] ) >= nonNanMin )

    cmsScores.normalizeColumnsWithinGroups( cols = ( 'complike', 'complikeRatio', 'complikeExp', 'complikeRatioExp' ) +
                                            CMSBinsBase.CMSstats +
                                            tuple([ c + 'like' for c in CMSBinsBase.CMSstats ]) +
                                            tuple([ c + 'likeRatio' for c in CMSBinsBase.CMSstats ]),
                                            groupCols = 'Chrom',
                                            normedSfx = 'normedLocal',
                                            meanSfx = 'mean', stdSfx = 'std',
                                            countSfx = 'count' ).save( complikeNormedLocalFN )
    
def normalizeComplikeRatioWithinReplicas( Ddata, thinSfx, scenario, putativeMutPop, complikeSfx = None,
                                          statsSfx = '', ihsSfx = '',
                                          getio = None ):
    """Normalize the complike ratio scores."""

    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    #complikeFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', putativeMutPop, complikeSfx ) )
    complikeFN = GetCreates( computeCMS, **Dict( 'Ddata thinSfx scenario putativeMutPop complikeSfx '
                                                 'likesTableSfx statsSfx ihsSfx') )[0]
    os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', putativeMutPop, complikeSfx ) )

    complikeWithLogCLRFN = os.path.join( snpStatsDir,
                                         AddFileSfx( 'complikeWithLogCLRnormedLocal.data/', putativeMutPop, complikeSfx ) )
    
    if getio: return dict( depends_on = complikeFN, creates = complikeWithLogCLRFN,
                           attrs = Dict( 'scenario putativeMutPop' ),
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ) )

    complike = IDotData( Path = complikeFN )

    def GetNormedBlocks():
        for chrom, r in complike.groupby( 'Chrom' ):
            z = DotData( names = r.headings, Records = r )
            z = clean_hapmap2_region( z, cleanIHS = True, cleanDer = False, cleanAnc = False, keepCausal = True )
            
            yield IDotData.fromIterable( headings = z.dtype.names, iterable = itertools.imap( tuple, z ) )

    IDotData.vstackFromIterable( GetNormedBlocks() ).save( complikeWithLogCLRFN )

    
def computeSNPranks( Ddata, thinSfx, scenario, putativeMutPop, complikeSfx = None,
                     likesTableSfx = None, getio = None ):
    """Within each replica, compute the rank of each SNP by CMS score.
    """

    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    complikeFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', putativeMutPop, complikeSfx, likesTableSfx ) )
    complikeWithRanksFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/',
                                                                 putativeMutPop, complikeSfx, likesTableSfx, 'ranks' ) )

    if getio: return dict( depends_on = complikeFN, creates = complikeWithRanksFN,
                           mediumRuleNameSfx = ( scenario.scenDir(), putativeMutPop ) )

    complikeWithLogCLR = IDotData( complikeFN )

    def RanksColumn( colName ):
        for chrom, r in complikeWithLogCLR.groupby( 'Chrom' ):
            replaceNans = r.mapRecords( lambda rec: -1e10 if any( map( isnan, rec ) ) else rec.complike,
                                        headings = ( colName, ) )

            dbg( 'chrom' )
            for x in replaceNans[ colName ].argsort()[::-1]:
                yield x

    complikeWithLogCLR.hstack( IDotData.fromIterable( headings = ( 'complikeRank', ),
                                                      iterable = RanksColumn( 'complike' ) ) ).save( complikeWithRanksFN )


def addIHHDiff_calc( mergedData, putativeMutPop ):
    """Add for each SNP the unstandardized iHHDiff value, and the id of the frequency bin."""
    
    Frequency = np.arange(0.05, 1.05, 0.05)
    
    freqBinSize = 0.05
    numFreqBins = int( 1.0 / freqBinSize )

    freqBinIDs = np.repeat( -1, len( mergedData ) )
    freqBinWith01IDs = np.repeat( -1, len( mergedData ) )
    for i, r in enumerate( mergedData ):
        derFreq = 1 - r[ 'FREQ1 %d' % putativeMutPop ]
        iHHDiff = r['Both iHH_D'] - r['Both iHH_A']

        for freqBin in range(len(Frequency)):
            if ((Frequency[freqBin] - derFreq) < .05 and Frequency[freqBin] - derFreq > 0):
                freqBinIDs[ i ] = freqBin
            if ((Frequency[freqBin] - derFreq) < .05 and Frequency[freqBin] - derFreq >= 0):
                freqBinWith01IDs[ i ] = freqBin

    iHHDiff = mergedData[ 'Both iHH_D' ] - mergedData[ 'Both iHH_A' ]

    return mergedData.hstack( DotData( names = ( 'iHHDiff', 'freqBinId', 'freqBinWith01Id' ),
                                       Columns = ( iHHDiff, freqBinIDs, freqBinWith01IDs ) ) )

def addIHHDiff( Ddata, thinSfx, scenario, putativeMutPop, simsOut, statsSfx, pop2name, ihsSfx, getio = None ):
    """Add for each SNP the unstandardized iHHDiff value, and the id of the frequency bin."""

    mergedDataFN = GetCreates( mergeSims, **Dict( 'Ddata scenario putativeMutPop simsOut statsSfx thinSfx pop2name ihsSfx' ) )[0]

    mergedWithIHHFN = AddFileSfx( mergedDataFN, 'withIHH' )

    if getio: return dict( depends_on = mergedDataFN, creates = mergedWithIHHFN,
                           uses = addIHHDiff_calc,
                           splitByCols = { mergedDataFN: dict(  keyCols = () ) },
                           attrs = Dict( 'putativeMutPop', scenario = scenario.scenDir() ) )

    mergedData = DotData( mergedDataFN )
    result = addIHHDiff_calc( mergedData = mergedData, putativeMutPop = putativeMutPop )
    result.save( mergedWithIHHFN )

def getFN_ihhDiff_sumsByFreq( **kwargs ):
    return AddFileSfx( GetCreates( addIHHDiff, **kwargs )[0], 'sumsByFreq' )

def getFN_ihhDiff_sums( **kwargs ):
    return AddFileSfx( GetCreates( addIHHDiff, **kwargs )[0], 'sums' )


def DefineRulesTo_computeCMSstats( pr, Ddata, simsOut = 'simsOut', mutPops = AllPops,
                                   mutAges = AllAges, mutFreqs = AllFreqs, nreplicas = 100,
                                   scen2nreplicas = {},
                                   oldMerged = False, DdataMerged = None,
                                   thinExt = '', thinSfx = '', sampleSize = 120, pop2sampleSize = {},
                                   pop2name = pop2name,
                                   DdataNeutral = None,
                                   stats = CMSBins.CMSstats,
                                   statsSfx = '', ihsSfx = '', ihhBackgroundArgs = None,
                                   normalizeToNeutral = True,
                                   normalizeWithinReplicas = True,
                                   forceIncludeNeutral = False,
                                   limitToPop = None ):
    """Define rules to compute CMS component statistics for each SNP, from simulations and their Sweep analyses.

    Params:

       scen2nreplicas - map from scenario to number of replicas in that scenario; for scenarios not in the map,
          'nreplicas' is assumed.
          
    """

    includeNeutral = forceIncludeNeutral  or  not DdataNeutral and normalizeToNeutral

    dbg( '"YYYYYYYYYYYYY" mutAges mutPops mutFreqs includeNeutral' )

    for scenario in GetScenarios( **Dict( 'mutAges mutPops mutFreqs includeNeutral' ) ):

        
        for putativeMutPop in ( mutPops if scenario.isNeutral() else (scenario.mutPop,) ):
            nreplicas = DictGet( scen2nreplicas, scenario, nreplicas )

            replicaPosFiles = []
            for replicaNum in range( nreplicas ):
                rule = pr.addInvokeRule( invokeFn = mergePosFilesOneSim,
                                         invokeArgs = Dict( 'Ddata simsOut thinSfx thinExt '
                                                            'scenario pop2name putativeMutPop '
                                                            'statsSfx ihsSfx replicaNum' ) )
                replicaPosFiles.append( rule.creates[0] )

            posFileFN = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir(),
                                      AddFileSfx( 'mergedPosStacked.tsv', statsSfx, putativeMutPop, ihsSfx ) )
            
            pr.addInvokeRule( invokeFn = vstack_tsvs,
                              invokeArgs = dict( inFNs = replicaPosFiles, outFN = posFileFN ),
                              name = 'merged_posfiles_vstack' )

            pr.addInvokeRule( invokeFn = mergeSims,
                              invokeArgs = Dict( 'Ddata simsOut posFileFN thinSfx thinExt '
                                                 'scenario putativeMutPop pop2name nreplicas '
                                                 'statsSfx ihsSfx limitToPop' ) )

            r = pr.addInvokeRule( invokeFn = addIHHDiff,
                                  invokeArgs = Dict( 'Ddata scenario putativeMutPop simsOut statsSfx thinSfx pop2name ihsSfx' ) )

            args = Dict( 'Ddata thinSfx scenario putativeMutPop simsOut statsSfx pop2name ihsSfx' )
            ihhDiff_sumsByFreqFN = getFN_ihhDiff_sumsByFreq( **args )
            ihhDiff_sumsFN = getFN_ihhDiff_sums( **args )
                        
            pr.addInvokeRule( invokeFn = computeSumsWithinGroups,
                              invokeArgs = dict( inFN = r.creates[0], cols = 'iHHDiff', groupCols = 'freqBinId',
                                                 groupsAreContiguous = False,
                                                 outFN = ihhDiff_sumsByFreqFN ),
                              name = 'computeSumsWithinGroups_ihhDiff_byFreq' )

            pr.addInvokeRule( invokeFn = computeSumsWithinGroups,
                              invokeArgs = dict( inFN = r.creates[0], cols = 'iHHDiff', groupCols = (),
                                                 groupsAreContiguous = False,
                                                 outFN = ihhDiff_sumsFN ), name = 'computeSumsWithinGroups_ihhDiff' )

            
            pr.addInvokeRule( invokeFn = computeCMSstats,
                              #invokeFnOld = computeCMSstats_old,
                              invokeArgs = Dict( 'Ddata thinSfx scenario putativeMutPop sampleSize pop2sampleSize '
                                                 'pop2name oldMerged DdataMerged stats statsSfx ihsSfx simsOut ihhBackgroundArgs limitToPop' ) )

            # pr.addInvokeRule( invokeFn = computeCMSstats_old,
            #                   #invokeFnOld = computeCMSstats_old,
            #                   invokeArgs = Dict( 'Ddata thinSfx scenario putativeMutPop sampleSize pop2sampleSize '
            #                                      'pop2name oldMerged DdataMerged stats',
            #                                      outFile = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir(),
            #                                                              AddFileSfx( 'cmsStatsRawOld.tsv', statsSfx, putativeMutPop,
            #                                                                          thinSfx) ) ) )
            

            if normalizeToNeutral:
                pr.addInvokeRule( invokeFn = normalizeStats,
                                  invokeArgs = Dict( 'Ddata thinSfx scenario putativeMutPop DdataNeutral statsSfx ihsSfx pop2name',
                                                     outFile =
                                                     os.path.join( snpStatsDir,
                                                                   'cmsStatsRaw_%d_normed.tsv' % putativeMutPop )
                                                     if not scenario.isNeutral() else
                                                     os.path.join( Ddata, 'snpStats', 'neutral', 'cmsStatsRaw_%d_normed.tsv' % putativeMutPop )) )

            if normalizeWithinReplicas:
                cmsFN = GetCreates( computeCMSstats, **Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx pop2name ihsSfx') )[0]
                DefineRulesTo_normalizeColumnsWithinGroups( pr, inFN = cmsFN,
                                                            cols = stats, groupCols = 'Chrom',
                                                            outFN = AddFileSfx( cmsFN, 'normedLocal' ) )
                # r = pr.addInvokeRule( invokeFn = normalizeStatsWithinReplicas,
                #                       invokeArgs = Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx ihsSfx pop2name' ) )


def DefineRulesTo_computeCMSforSims( pr, Ddata, simsOut = 'simsOut', mutPops = AllPops,
                                     mutAges = AllAges, mutFreqs = AllFreqs, nreplicas = 100, scen2nreplicas = {},
                                     thinExt = '', thinSfx = '', sampleSize = 120, pop2sampleSize = {},
                                     pop2name = pop2name,
                                     oldMerged = False,
                                     DdataMerged = None,
                                     likesTable = '../Data/Common_Data/sim/likesTable_toneut5to25full_1.tsv',
                                     likesDir = None,
                                     likesSfx = '_all_ages',
                                     likesBins = 'sel20_1_Summary.txt',
                                     likesPopSpecific = True,

                                     replicaTables = (), replicaConds = (),
                                     nonNormedStats = CMSBins.nonNormedStats,
                                     freq = 'all',
                                     complikeSfx = None,
                                     nonNanStats = ( 'iHS', 'max_xpop' ),
                                     nonNanMin = 1,
                                     stats = CMSBins.CMSstats,
                                     DdataNeutral = None,
                                     useRegionLikes = True,
                                     normalizeWithinReplicas = True,
                                     forceNormalizeWithinReplicas = False,
                                     cleanIHS = True, cleanDer = False, cleanAnc = False, keepCausal = False,
                                     excludeNeutral = False,
                                     includeRulesTo_computeCMSstats = True,
                                     forceIncludeNeutral = False,
                                     useGenMap = None,
                                     useDefaultNumSnps = True,
                                     statsSfx = '', ihsSfx = '',
                                     ihhBackgroundArgs = None,
                                     piperun_memreq = None,
                                     limitToPop = None ):
    """Define rules to compute CMS scores for all SNPs in a set of simulations (or real data in simulations format).
    """

    if likesTable:
        mainLikesTable = LoadLikesTable( **Dict( 'Ddata likesTable', likesTableSfx =  '' ) )
        normalizeWithinReplicas = coerceVal( mainLikesTable.normalizeWithinReplicas )

    if forceNormalizeWithinReplicas: normalizeWithinReplicas = True
        
    if includeRulesTo_computeCMSstats:
        dbg( '"DEFINING" normalizeWithinReplicas' )
        DefineRulesTo_computeCMSstats( normalizeToNeutral = not normalizeWithinReplicas,
                                       **Dict( 'pr Ddata simsOut mutPops mutAges mutFreqs nreplicas scen2nreplicas thinExt '
                                               'thinSfx sampleSize pop2sampleSize pop2name oldMerged DdataMerged '
                                               'DdataNeutral stats statsSfx ihsSfx ihhBackgroundArgs limitToPop' ) )

    nreplicasDefault = nreplicas
    replicaTables = MakeSeq( replicaTables )
    
    for scenario in GetScenarios( includeNeutral = forceIncludeNeutral or \
                                      not DdataNeutral and not normalizeWithinReplicas and not excludeNeutral,
                                  **Dict( 'mutAges mutPops mutFreqs' ) ):
        for putativeMutPop in ( mutPops if scenario.isNeutral() else (scenario.mutPop,) ):

            dbg( '"QQQQQQQQQQQQ" scenario putativeMutPop' )

            nreplicas = DictGet( scen2nreplicas, scenario, nreplicasDefault )

            computeCMSargs = Dict( 'Ddata likesTable likesDir likesSfx likesBins likesPopSpecific thinSfx '
                                   'scenario nreplicas '
                                   'putativeMutPop nonNormedStats nonNanStats complikeSfx '
                                   'normalizeWithinReplicas forceNormalizeWithinReplicas useRegionLikes '
                                   'cleanIHS cleanDer cleanAnc keepCausal useDefaultNumSnps statsSfx '
                                   'pop2name ihsSfx')

            dbg( 'computeCMSargs' )
            
            if not replicaConds:
                pr.addInvokeRule( invokeFn = computeCMS,
                                  attrs = Dict( 'piperun_memreq' ) if piperun_memreq else {},
                                  invokeArgs = computeCMSargs )

                pr.addInvokeRule( invokeFn = plotLikesTable,
                                  invokeArgs = dict( hitsLikesFile = mainLikesTable.hitsLikes,
                                                     missLikesFile = mainLikesTable.missLikes,
                                                     regionLikesFile = mainLikesTable.regionLikes if hasattr( mainLikesTable, 'regionLikes' ) else None,
                                                     likesTable = likesTable,
                                                     normalizeWithinReplicas = normalizeWithinReplicas,
                                                     plotFile = os.path.join( Ddata,
                                                                              AddFileSfx( 'showLikes.svg',
                                                                                          complikeSfx ) ) ) )
            else:

                #
                # We need to use different likes tables for different replicas,
                # as specified in replicaConds.
                # We do this by first splitting the raw CMS stats into separate files,
                # where stats in each file meet one condition (e.g. have a hi-freq sweep)
                # and so should be processed with one likes table.
                # We then compute CMS scores separately for each raw CMS stats file.
                # Finally, we merge the computed CMS scores into a common file that
                # now has CMS scores for all the replicas.
                #
                
                condNames = map( itemgetter(0), replicaConds )
                condsFileFN = AddFileSfx( 'replicaConds.tsv', complikeSfx, *condNames )
                pr.addInvokeRule( invokeFn = identifyReplicasMeetingConds,
                                  invokeArgs = 
                                  Dict( 'replicaTables replicaConds Ddata scenario thinSfx condsFileFN nreplicas' ) )
                
                snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )

                aLikesTable = LoadLikesTable( **Dict( 'Ddata likesTable', likesTableSfx = condNames[0] ) )

                cmsStatsRawFN = os.path.join( snpStatsDir, AddFileSfx( 'cmsStatsRaw.data/', putativeMutPop,
                                                                       ihsSfx ) )
                cmsStatsNormFN = os.path.join( snpStatsDir,
                                               AddFileSfx( 'cmsStatsNorm' +
                                                           ( 'edLocal' if aLikesTable.normalizeWithinReplicas else '' )
                                                           + '.data/', putativeMutPop, ihsSfx ) )
                
                pr.addInvokeRule( invokeFn = splitSnpStatsFile,
                                  invokeArgs = dict( inFileFN = cmsStatsRawFN, sfx = complikeSfx,
                                                     **Dict( 'condNames condsFileFN Ddata scenario thinSfx' ) ),
                                  name = 'splitCMSstatsRaw' )

                pr.addInvokeRule( invokeFn = splitSnpStatsFile,
                                  invokeArgs = dict( inFileFN = cmsStatsNormFN, sfx = complikeSfx,
                                                     **Dict( 'condNames condsFileFN Ddata scenario thinSfx' ) ),
                                  name = 'splitCMSstatsNorm' )

                mainLikesTable = LoadLikesTable( **Dict( 'Ddata likesTable', likesTableSfx = '' ) )
                dbg( '"MMMMM" mainLikesTable' )
                for condName in condNames:
                    pr.addInvokeRule( invokeFn = computeCMS,
                                      invokeFnOld = computeCMSold,
                                      invokeArgs =
                                      dict( likesTableSfx = condName, rawNormSfx = complikeSfx, **computeCMSargs ) )

                    
                    pr.addInvokeRule( invokeFn = plotLikesTable,
                                      invokeArgs = dict( hitsLikesFile = AddFileSfx( mainLikesTable.hitsLikes, condName ),
                                                         missLikesFile = AddFileSfx( mainLikesTable.missLikes, condName ),
                                                         regionLikesFile = AddFileSfx( mainLikesTable.regionLikes,
                                                                                       condName ),
                                                         likesTable = AddFileSfx( likesTable, condName ),
                                                         normalizeWithinReplicas = normalizeWithinReplicas,
                                                         condName = condName,
                                                         plotFile = os.path.join( Ddata, AddFileSfx( 'showLikes.svg',
                                                                                                     condName ) ) ) )
                                      
                complikeFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', putativeMutPop, complikeSfx, ihsSfx ) )
                pr.addInvokeRule( invokeFn = joinSnpStatsFiles,
                                  invokeArgs = dict( outFileFN = complikeFN,
                                                     **Dict( 'Ddata scenario condNames condsFileFN thinSfx' ) ) )

#            pr.addInvokeRule( invokeFn = normalizeComplikeRatio,
#                              invokeArgs =
#                              Dict( 'Ddata thinSfx scenario putativeMutPop DdataNeutral complikeSfx' ),
#                              attrs = dict( piperun_memreq = 5 ) ) 

            cmsFN = GetCreates( computeCMS, **computeCMSargs )[0]

            if useRegionLikes:
                DefineRulesTo_normalizeColumnsWithinGroups( pr = pr, inFN = cmsFN,
                                                            outFN = AddFileSfx( cmsFN, 'normedLocal' ),
                                                            groupCols = ( 'Chrom', ),
                                                            cols = ( 'complike', 'complikeRatio', 'complikeExp', 'complikeRatioExp' ),
                                                            nameSfx = 'cms' )

            else:
                DefineRulesTo_normalizeColumnsWithinGroups( pr = pr, inFN = cmsFN,
                                                            outFN = AddFileSfx( cmsFN, 'normed' ),
                                                            groupCols = (),
                                                            cols = ( 'complike', 'complikeRatio', 'complikeExp', 'complikeRatioExp' ),
                                                            nameSfx = 'cms' )
            
            
            # pr.addInvokeRule( invokeFn = normalizeComplikeWithinReplicas,
            #                   invokeArgs =
            #                   Dict( 'Ddata thinSfx scenario putativeMutPop complikeSfx statsSfx '
            #                         'nonNanStats nonNanMin ihsSfx pop2name' ),
            #                   attrs = Dict( 'complikeSfx' ) )


def callSelectionByCMS( Ddata, thinSfx, scenario, putativeMutPop, complikeSfx, likesTableSfx,
                        windowSize = 100000, signifSnpFraction = .3, cmsThreshold = 4,
                        normedSfx = 'normed', chromCol = 'Chrom', cmsStat = 'complikeExp',
                        windowCol = 'Pos',
                        getio = None ):
    "Identify replicas that would be flagged as under selection by CMS"

    if not Ddata.endswith( '/' ): Ddata += '/'
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats' + thinSfx, scenario.scenDir() )

    cmsNormFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.tsv', putativeMutPop, normedSfx, complikeSfx,
                                                       likesTableSfx ) )

    cmsCallFN = os.path.join( replicaStatsDir, AddFileSfx( 'cmsCall.tsv', putativeMutPop, complikeSfx,
                                                           likesTableSfx, normedSfx, windowSize, signifSnpFraction, cmsThreshold ) )

    if getio: return dict( depends_on = cmsNormFN, creates = cmsCallFN )

    cmsNorm = IDotData( cmsNormFN )

    with IDotData.openForWrite( cmsCallFN, headings = ( 'Chrom', 'Call' ) ) as cmsCall:
        for replicaNum, idd1, idd2 in cmsNorm.groupby( chromCol, multiPass = 2):

            windowStart = iter( idd1 )
            windowEnd = iter(idd2 )

            callHere = 0

            numSnps = 0
            numSnpsOverThreshold = 0

            rStart = next( windowStart )
            rEnd = next( windowEnd )
            while( rEnd[ windowCol ] - rStart[ windowCol ] < windowSize ):
                numSnps += 1
                if rEnd[ cmsStat ] > cmsThreshold:
                    numSnpsOverThreshold += 1
                rEnd = next( windowEnd )

            try:
                while True:
                    if numSnpsOverThreshold / numSnps > signifSnpFraction:
                        callHere = 1
                        break
                    rEnd = next( windowEnd )
                    numSnps += 1
                    if rEnd[ cmsStat ] > cmsThreshold:
                        numSnpsOverThreshold += 1

                    while( rEnd[ windowCol ] - rStart[ windowCol ] > windowSize ):
                        numSnps -= 1
                        if rStart[ cmsStat ] > cmsThreshold:
                            numSnpsOverThreshold -= 1
                        rStart = next( windowStart )

            except StopIteration:
                pass
            
            cmsCall.writeRecord( replicaNum, callHere )


def callSelectionByCMS_withSplit( Ddata, thinSfx, scenario, putativeMutPop, complikeSfx, likesTableSfx,
                                  windowSize = 100000, signifSnpFraction = .3, cmsThreshold = 4,
                                  normedSfx = 'normed', chromCol = 'Chrom', cmsStat = 'complikeExp',
                                  windowCol = 'Pos',
                                  getio = None ):
    "Identify replicas that would be flagged as under selection by CMS"

    if not Ddata.endswith( '/' ): Ddata += '/'
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats' + thinSfx, scenario.scenDir() )

    cmsNormFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.tsv', putativeMutPop, normedSfx, complikeSfx,
                                                       likesTableSfx ) )

    cmsCallFN = os.path.join( replicaStatsDir, AddFileSfx( 'cmsCallWithSplit.tsv', putativeMutPop, complikeSfx,
                                                           likesTableSfx, normedSfx, windowSize, signifSnpFraction, cmsThreshold ) )

    if getio: return dict( depends_on = cmsNormFN, creates = cmsCallFN,
                           splitByCols = { cmsNormFN: dict( keyCols = ( 'Chrom', ) ) } )

    windowHalf = windowSize / 2
    calls = []
    cmsNorm = DotData( SVPath = cmsNormFN )
    for chrom in sorted( set( cmsNorm.Chrom ) ):
        z = cmsNorm[ cmsNorm.Chrom == chrom ]
        chromEndPos = z[ len( z ) - 1 ].Pos
        highScores = z[ ( z[ cmsStat ] > cmsThreshold ) & ( z.Pos + windowSize < chromEndPos ) ]
        foundWindow = False
        for r in highScores:
            windowSNPs = z[ ( r.Pos <= z.Pos ) & ( z.Pos <= r.Pos + windowSize ) ]
            fraction = np.sum( windowSNPs[ cmsStat ] >= cmsThreshold ) / len( windowSNPs )
            if fraction >= signifSnpFraction:
                foundWindow = True
                break
        calls.append( ( chrom, foundWindow ) )
    DotData( names = ( 'Chrom', 'Call' ), Records = calls ).saveToSV( cmsCallFN )
            
            
def gatherReplicaThresholds_fixedFractions( Ddata, thinSfx, scenario, putativeMutPop, complikeSfx, likesTableSfx,
                                            windowSize = 100000, signifSnpFraction = .3, cmsThreshold = 4,
                                            normedSfx = 'normed', chromCol = 'Chrom', cmsStat = 'complikeExp',
                                            windowCol = 'Pos',
                                            getio = None ):
    """For fixed window size and SNP fraction, for each replica find the highest CMS threshold that this top fraction
of SNPs in a window has."""

    if not Ddata.endswith( '/' ): Ddata += '/'
    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats' + thinSfx, scenario.scenDir() )

    cmsNormFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', normedSfx, putativeMutPop, complikeSfx,
                                                       likesTableSfx ) )

    cmsCallFN = os.path.join( replicaStatsDir, AddFileSfx( 'cmsCall2.tsv', putativeMutPop, complikeSfx,
                                                           likesTableSfx ) )

    if getio: return dict( depends_on = cmsNormFN, creates = cmsCallFN )

    cmsNorm = IDotData( cmsNormFN )

    with IDotData.openForWrite( cmsCallFN, headings = ( 'Chrom', 'maxThresh' ) ) as cmsCall:
        for replicaNum, idd1, idd2 in cmsNorm.groupby( chromCol, multiPass = 2):

            windowStart = iter( idd1 )
            windowEnd = iter(idd2 )

            maxThresh = None
            lastWindowEnd = None

            numSnps = 0

            rblist = RBList( unique = False )

            rStart = next( windowStart )
            rEnd = next( windowEnd )
            while( rEnd[ windowCol ] - rStart[ windowCol ] < windowSize ):
                numSnps += 1
                rblist.insert( rEnd[ cmsStat ] )
                rEnd = next( windowEnd )

            try:
                while True:

                    logging.info( str( numSnps ) )

                    rEnd = next( windowEnd )
                    numSnps += 1
                    rblist.insert( rEnd[ cmsStat ] )

                    while( rEnd[ windowCol ] - rStart[ windowCol ] > windowSize ):
                        numSnps -= 1
                        rblist.remove( rStart[ cmsStat ], all = False ) 
                        rStart = next( windowStart )

            except StopIteration:
                pass

            logging.info( 'done with replica ' + str( replicaNum ) )
            #cmsCall.writeRecord( replicaNum, callHere )

def DefineRulesTo_normalizeCMS( pr, Ddata, mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs,
                                thinSfx = '', complikeSfx = '', statsSfx = '',
                                likesTableSfx = '', includeNeutral = True ):
    for putativeMutPop in mutPops:
        pr.addInvokeRule( invokeFn = getNeutralMeanStdForCMS,
                          invokeArgs = Dict( 'Ddata thinSfx putativeMutPop complikeSfx likesTableSfx statsSfx' ) )
    
    for scenario in GetScenarios( **Dict( 'mutAges mutFreqs mutPops includeNeutral' ) ):
        for putativeMutPop in ( mutPops if scenario.isNeutral() else ( scenario.mutPop, ) ): 
            pr.addInvokeRule( invokeFn = normalizeCMS,
                              invokeArgs = Dict( 'Ddata scenario thinSfx putativeMutPop complikeSfx likesTableSfx statsSfx' ) )
            
def DefineRulesTo_callSelectionByCMS( pr, Ddata, mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs,
                                      thinSfx = '', complikeSfx = '',
                                      likesTableSfx = '', signifSnpFraction = .3, cmsThreshold = 3,
                                      normedSfx = 'normed', cmsStat = 'complikeRatio' ):

    #DefineRulesTo_normalizeCMS( **Dict( 'pr Ddata mutAges mutPops mutFreqs thinSfx complikeSfx likesTableSfx' ) )
    for scenario in GetScenarios( **Dict( 'mutAges mutFreqs mutPops' ) ):
        for putativeMutPop in ( mutPops if scenario.isNeutral() else ( scenario.mutPop, ) ): 
            pr.addInvokeRule( invokeFn = callSelectionByCMS_withSplit,
                              invokeArgs = Dict( 'Ddata scenario thinSfx putativeMutPop complikeSfx likesTableSfx '
                                                 'signifSnpFraction cmsThreshold normedSfx cmsStat') )

            
def DefineRulesTo_evalCMS( pr, Ddata, mutPops = AllPops,
                           mutAges = AllAges, mutFreqs = AllFreqs, nreplicas = 100,
                           thinExt = '', thinSfx = '',
                           complikeSfx = '', likesTableSfx = '', simsOut = 'simsOut',
                           cumulativeUpTo = 0.99, nonNanStats = 'ALL',
                           DdataNeutral = None, excludeNeutral = False ):
    """Define rules to evaluate the results of CMS localization."""

    DefineRulesTo_gatherCausalRanks( **Dict( 'pr Ddata mutAges complikeSfx nonNanStats' ) )
    DefineRulesTo_gatherCausalFreqs( **Dict( 'pr Ddata mutAges simsOut thinExt' ) )


    DefineRulesTo_callSelectionByCMS( **Dict( 'pr Ddata mutAges mutPops mutFreqs thinSfx complikeSfx likesTableSfx' ) )

    stat2histFiles = {}
    dfltArgs = dict( coarsenBy = 32, binSize = 10, ticksCoarsen = 2 )
    for statTable, statCol, extraArgs in \
            ( ( 'causalRank', 'causalRank', dict( coarsenBy = 8 ) ),
              ( 'causalStats', 'complikeExp_normedLocal', {} ),
              ( 'causalStats', 'complike',  {} ), ( 'causalStats', 'complikeExp',  {} ),
              ( 'causalStats', 'complike_mean', {} ), ( 'causalStats', 'complike_std', {} ),
              ( 'causalStats', 'complike_count', {} )):
        causalRankTable = statTable + Sfx( complikeSfx, 'nonNan', *MakeSeq( nonNanStats ) )
        stat2histFiles[ ( statTable, statCol ) ] = \
            DefineRulesTo_histogramReplicaStatistic( pr = pr, Ddata = Ddata,
                                                     outFile = '%s_%s.svg' % ( statTable, statCol ),
                                                     sfx = complikeSfx,
                                                     replicaTables = ( causalRankTable, 'replicaStats' ),
                                                     replicaStat = causalRankTable + '.' + statCol,
                                                     replicaConds = ( 'True', 'replicaStats.causalAlleleFreq < .5',
                                                                      'replicaStats.causalAlleleFreq >= .5' ),
                                                     replicaCondsSfxs = ( 'all', 'lo', 'hi' ),
                                                     scenCond = 'not is_neutral',
                                                     normed = True, cumulativeUpTo = cumulativeUpTo,
                                                     **MergeDicts( dfltArgs, extraArgs ) )

    for causalStat, binSize in ( 'complike', .1 ), ( 'complikeExp', 1e-3 ):
        DefineRulesTo_histogramSnpStatistic( outFile = os.path.join( Ddata, 'snpHists', AddFileSfx( 'selVsRegion.svg',
                                                                                                    complikeSfx,
                                                                                                    causalStat ) ),

                                             replicaTables = ( 'replicaStats', ),
                                             replicaCondsSfxs = ( 'all', 'lo', 'hi' ),
                                             replicaConds = ( 'True', 'replicaStats.causalAlleleFreq < .5',
                                                              'replicaStats.causalAlleleFreq >= .5' ),
                                             scenCond = 'not is_neutral',
                                             snpTables = ( 'complike_normedLocal.data', ),
                                             scen2sfxs = lambda scen: Sfx( scen.mutPop, complikeSfx, 'nonNan', 'ALL' ),
                                             snpConds = ( 'complike_normedLocal.Pos != 500000',
                                                          'complike_normedLocal.Pos == 500000' ),
                                             snpCondsSfxs = ( 'region', 'selected' ),
                                             snpStat = 'complike_normedLocal.' + causalStat,
                                             allScens = GetSelectionScenarios(),
                                             coarsenBy = 16,
                                             subplots_adjust = dict( bottom = .3 ),
                                             **Dict( 'pr Ddata binSize nreplicas' ) )


    DefineRulesTo_extractGeneticMapsFromSims( **Dict( 'pr Ddata mutAges mutFreqs mutPops' ) )
    DefineRulesTo_gatherCausalSnpGdPos( **Dict( 'pr Ddata mutAges mutFreqs mutPops' ) )

    for scenario in GetScenarios( includeNeutral = not DdataNeutral and not excludeNeutral,
                                  **Dict( 'mutAges mutPops mutFreqs' ) ):
        for putativeMutPop in ( mutPops if scenario.isNeutral() else (scenario.mutPop,) ):

            pr.addInvokeRule( invokeFn = computeSNPranks,
                              invokeArgs =
                              Dict( 'Ddata thinSfx scenario putativeMutPop complikeSfx' ) )

            pr.addInvokeRule( invokeFn = localizeSpatiallyBySplineFitting,
                              invokeArgs = Dict( 'Ddata scenario nreplicas thinSfx putativeMutPop complikeSfx' ) )

            pr.addInvokeRule( invokeFn = localizeSpatiallyByWindows,
                              invokeArgs = Dict( 'Ddata scenario nreplicas thinSfx putativeMutPop complikeSfx' ) )

            if False:
                pr.addInvokeRule( invokeFn = plotCMS,
                                  invokeArgs = Dict( 'Ddata scenario nreplicas thinSfx putativeMutPop complikeSfx',
                                                     includeSpatialLoc = False, xAxisGd = True ),
                                  name = 'plotCMS_basic_gd' )

                pr.addInvokeRule( invokeFn = plotCMS,
                                  invokeArgs = Dict( 'Ddata scenario nreplicas thinSfx putativeMutPop complikeSfx',
                                                     includeSpatialLoc = False, xAxisGd = False ),
                                  name = 'plotCMS_basic_bp' )

            for whichSpatialLoc in 'Spline', 'Windows':
                for xAxisGd in False, True:
                    pr.addInvokeRule( invokeFn = plotCMS,
                                      invokeArgs = Dict( 'Ddata scenario nreplicas thinSfx putativeMutPop complikeSfx '
                                                         'whichSpatialLoc xAxisGd' ) )

            pr.addInvokeRule( invokeFn = plotCMS_gd2,
                              invokeArgs = Dict( 'Ddata scenario nreplicas thinSfx putativeMutPop complikeSfx' ),
                              attrs =
                              dict( piperun_show_file =
                                    getBinPlotFilesFNs( **Dict( 'Ddata scenario thinSfx '
                                                                'putativeMutPop complikeSfx' )  )[:10] ) )
                              
            if not scenario.isNeutral():
                for whichSpatialLoc in 'Spline', 'Windows':
                    pr.addInvokeRule( invokeFn = evalSpatialLoc,
                                      invokeArgs = Dict( 'Ddata scenario nreplicas thinSfx putativeMutPop '
                                                         'complikeSfx whichSpatialLoc' ) )

    for whichSpatialLoc in 'Spline', 'Windows':

        spatialLocEval = 'spatialLocEval' + whichSpatialLoc
        DefineRulesTo_histogramReplicaStatistic( pr = pr, Ddata = Ddata,
                                                 outFile = 'spatialLocHist%s_causalIncluded.svg' % whichSpatialLoc,
                                                 sfx = ( complikeSfx, likesTableSfx ),
                                                 nameSfx = whichSpatialLoc,
                                                 replicaTables = ( spatialLocEval, 'replicaStats' ),
                                                 scen2sfxs =
                                                 lambda scen: { spatialLocEval : Sfx( scen.mutPop, complikeSfx, likesTableSfx ) },
                                                 replicaStat = spatialLocEval + '.causalIncluded',
                                                 replicaConds = ( 'True', 'replicaStats.causalAlleleFreq < .5',
                                                                  'replicaStats.causalAlleleFreq >= .5' ),
                                                 replicaCondsSfxs = ( 'all', 'lo', 'hi' ),
                                                 scenCond = 'not is_neutral',
                                                 binSize = .5, cumulative = False,
                                                 normed = True )

    return stat2histFiles
                    

def detectSelectionByCMS( Ddata, thinSfx, scenario, nreplicas, fraction = .3, threshold = 4, windowSize = 100000 ):
    """Detect replicas in which selection occurred, using the CMS scores.   A replica is flagged as having selection
    if in at least one window of size 'windowSize', the fraction 'fraction' of SNPs have CMS scores over 'threshold'.
    """

    


def DefineRulesTo_copyToUniformFileNames( pr, Ddata, simsOut = 'simsOut', mutPops = AllPops,
                                          mutAges = AllAges, mutFreqs = AllFreqs, nreplicas = 100 ):
    """Copy or link simulation analysis files to more uniform file names, so that each scenario directory under
    snpStats has a file with the same name.
    """

    for scen in GetSelectionScenarios( **Dict( 'mutAges mutFreqs mutPops' ) ):
        for whichFile in 'merged', 'cmsStatsRaw', 'cmsStatsNorm', 'cmsStatsNormedLocal':
            oldDir = os.path.join( Ddata, 'snpStats', scen.scenDir(), '%s_%d.data/' % ( whichFile, scen.mutPop  ) )
            newDir = os.path.join( Ddata, 'snpStats', scen.scenDir(), '%s.data/' % whichFile  )
            dbg( 'oldDir newDir' )
            pr.addInvokeRule( invokeFn = linkDotDatas, invokeArgs = dict( source = oldDir, target = newDir ) )
            
def DefineRulesTo_createLikesTables( pr, Ddata, simsOut = 'simsOut', mutPops = AllPops,
                                     mutAges = AllAges, mutFreqs = AllFreqs, nreplicas = 100,
                                     thinExt = '', thinSfx = '', sampleSize = 120, pop2sampleSize = {},
                                     pop2name = pop2name, likesSfx = 'all_ages', secondLikesSfx = '',
                                     selPos = 500000,
                                     stat_start = CMSBins.stat_start, stat_end = CMSBins.stat_end,
                                     stat_nbin = CMSBins.stat_nbin, stats = CMSBins.CMSstats,
                                     nonNormedStats = CMSBins.nonNormedStats,
                                     normalizeWithinReplicas = False, likesTable = 'likesTable.tsv',
                                     likesTableComment = None,
                                     statsSfx = '',
                                     replicaTables = (), replicaConds = ( ( '', 'True' ), ),
                                     includeRulesTo_computeCMSstats = True,
                                     copyLikesTablesTo = None ):
    """Define rules to create likes tables based on specified simset.  Can create several likes tables
    based on named subsets of the simulations, e.g. just the simulations with high-freq or with low-freq sweeps.

       Params:

            pr - the PipeRun object to which to add the rules
            Ddata - directory under which the simulations and their analyses reside

         Choosing which simulations to use:
         
            simsOut - directory under Ddata, under which the simulations reside
            mutPops, mutAges, mutFreqs - use replicas from which scenarios with selection?
               this chooses scenarios with selection in the populations 'mutPops', with
               selected allele appearing at times 'mutAges', and reaching present-day
               frequency of 'mutFreqs'.
            nreplicas - under each scenario (selection or neutral), how many replicas were simulated?
               right now we assume it's the same number of replicas for all scenarios.
            thinExt - file extension added to simulation output files (.thin if we're using
               thinned simulations)
            thinSfx - suffix 
            
            

         Output:

            likesDir - where to put the likes table.
    """

    pr.addLsfNeedDir( Ddata )

    if includeRulesTo_computeCMSstats:
        DefineRulesTo_computeCMSstats( **Dict( 'pr Ddata simsOut mutPops mutAges mutFreqs nreplicas thinExt thinSfx '
                                               'sampleSize pop2sampleSize pop2name statsSfx' ) )

    condNames = map( itemgetter(0), replicaConds )
    
    if normalizeWithinReplicas:
        stat_start = CMSBinsLocal.stat_start
        stat_end = CMSBinsLocal.stat_end
        stat_nbin = CMSBinsLocal.stat_nbin
        stats = CMSBinsLocal.CMSstats

    for condName, replicaCond in replicaConds:

        for scenario in GetScenarios( mutAges = mutAges, mutFreqs = mutFreqs, mutPops = mutPops ):
            for putativeMutPop in ( mutPops if scenario.isNeutral() else (scenario.mutPop,) ):
                analsimCreates = pr.addInvokeRule( invokeFn = analsim,
                                                   invokeArgs =
                                                   MergeDicts( Dict( 'replicaTables replicaCond condName' ) \
                                                                   if not scenario.isNeutral() else {},
                                                               Dict( 'Ddata scenario putativeMutPop thinSfx likesSfx '
                                                                     'normalizeWithinReplicas nreplicas statsSfx '
                                                                     'stat_start stat_end stat_nbin stats nonNormedStats '
                                                                     ) ) )['creates']
                likesBins = [ f for f in analsimCreates if 'Summary' in f ][0]

        for putativeMutPop in mutPops:
            pr.addInvokeRule( invokeFn = likes,
                              invokeArgs = Dict( 'Ddata mutAges mutPops mutFreqs thinSfx '
                                                 'putativeMutPop likesSfx secondLikesSfx condName statsSfx' ) )

            likesInfo = pr.addInvokeRule( invokeFn = saveLikesTable,
                                          invokeArgs = Dict( 'Ddata likesSfx secondLikesSfx condName likesBins normalizeWithinReplicas '
                                                             'likesTable putativeMutPop '
                                                             'likesTableComment pop2sampleSize statsSfx',
                                                             likesTableCreationParams =
                                                             str( Dict( 'mutAges mutFreqs mutPops nreplicas sampleSize '
                                                                        'pop2sampleSize condName replicaCond '
                                                                        'putativeMutPop '
                                                                        'stats nonNormedStats stat_start '
                                                                        'stat_end stat_nbin statsSfx' ) ) ) )

            for includeSpecialBins in ( False, True ):
                pr.addInvokeRule( invokeFn = plotLikesTable,
                                  invokeArgs = dict( hitsLikesFile = likesInfo.depends_on[ 0 ],
                                                     missLikesFile = likesInfo.depends_on[ 1 ],
                                                     regionLikesFile = likesInfo.depends_on[ 2 ],
                                                     likesBinsFile = likesInfo.depends_on[ 3 ],
                                                     likesTable = likesInfo.creates[ 0 ],
                                                     plotFile = ReplaceFileExt( AddFileSfx( likesInfo.creates[ 0 ],
                                                                                            'plot',
                                                                                            includeSpecialBins ), '.svg' ),
                                                     **Dict( 'nonNormedStats likesSfx secondLikesSfx condName putativeMutPop '
                                                             'normalizeWithinReplicas statsSfx includeSpecialBins' ) ) )

            if copyLikesTablesTo:
                for f in likesInfo.depends_on + likesInfo.creates:
                    fDir = os.path.dirname( f )
                    pr.addCpRule( comment = 'copy likes tables',  name = 'copyLikesTables',
                                  source = f,
                                  target = AddSlash( os.path.join( copyLikesTablesTo,
                                                                   'likes' if os.path.basename( fDir ) == 'likes'
                                                                   else '' ) ) )
            

def analsim_old( Ddata, scenario, putativeMutPop = None, thinSfx = '', likesSfx = '', condName = '',
                 stat_start = CMSBins.stat_start, stat_end = CMSBins.stat_end, stat_nbin = CMSBins.stat_nbin,
                 stats = CMSBins.CMSstats, nonNormedStats = CMSBins.nonNormedStats,
                 normalizeWithinReplicas = False, nreplicas = 100,
                 replicaTables = (), replicaCond = '', statsSfx = '',
                 getio = None ):

    """Compute histograms of CMS component statistics for SNPs in chosen replicas of one scenario."""

    #if not os.path.exists(Ddata + '/power_GW/'): os.makedirs(Ddata +'/' +  name + '/')

    if not putativeMutPop: putativeMutPop = scenario.mutPop

    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )

    cmsStatsNormFN = GetCreates( normalizeStatsWithinReplicas if normalizeWithinReplicas else normalizeStats,
                                 **Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx' ) ).containing( 'cmsStatsNorm' )
                                     
    
    cmsStatsRawFN = GetCreates( computeCMSstats, **Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx' ) )[0]

    ranksFN = os.path.join( Ddata, 'scenStats' + thinSfx, scenario.scenDir(),
                            AddFileSfx( 'ranks.data/', likesSfx, putativeMutPop, statsSfx, condName ) )
    summaryFN = os.path.join( Ddata, AddFileSfx( 'Summary.txt', scenario.scenDir(),
                                                 likesSfx, putativeMutPop, statsSfx, condName ) )

    replicaTables = MakeSeq( replicaTables )

    replicaTableFiles = tuple(  os.path.join( Ddata, 'replicastats', scenario.scenDir(),
                                              replicaTable + '.tsv' )
                                for replicaTable in replicaTables )

    if getio: return dict( depends_on = ( cmsStatsNormFN, cmsStatsRawFN ) + replicaTableFiles,
                           creates = ( ranksFN, summaryFN ),
                           mediumRuleNameSfx = ( scenario.scenDir(), likesSfx, putativeMutPop, condName ),
                           fileDescrs = { ranksFN: 'Histograms of CMS component stats for causal/non-causal SNPs',
                                          summaryFN: 'Bin boundaries used for CMS component histograms' } )
    
    selpos = 500000

    #rep_freqs = DotData(SVPath = '../Data/Shari_Data/sim/rep_freq.tsv', names = ['rep','low','high'])

    cmsStatsNorm = DotData( Path = cmsStatsNormFN )
    cmsStatsRaw = DotData( Path = cmsStatsRawFN )

    assert len( cmsStatsNorm ) == len( cmsStatsRaw )
    assert all( cmsStatsNorm.Pos == cmsStatsRaw.Pos )
    assert all( cmsStatsNorm.Chrom == cmsStatsRaw.Chrom )

    dbg( 'len(cmsStatsNorm)' )
    
    if replicaCond and replicaCond != 'True':

        # identify replicas matching the condition

        replicaCondExpr = compile_expr( replicaCond )
        replicaTableVals = [ DotData( SVPath = f ) for f in replicaTableFiles ]

        assert all([ len( rtv ) == nreplicas for rtv in replicaTableVals ])

        replicasToUse = array( [ eval( replicaCondExpr, globals(), dict( zip( replicaTables, replicaTableRows ) ) )
                                 for replicaTableRows in itertools.izip( *replicaTableVals ) ] ) if replicaTables \
                                 else repeat( eval( replicaCondExpr, globals(), {} ), nreplicas )
                                 
        rowsToUse = replicasToUse[ array( map( int, cmsStatsNorm.Chrom ) ) ]
        cmsStatsNorm = cmsStatsNorm[ rowsToUse ]
        cmsStatsRaw = cmsStatsRaw[ rowsToUse ]

        dbg( '"after_filtering" len(cmsStatsNorm)' )

    hitmissCols = {}

    hitmissNames = []
    for stat in stats:
        for hitmiss in 'Hits', 'Misses':
            hitmissCols[ stat + hitmiss ] = zeros( stat_nbin[ stat ], dtype = int )
            hitmissNames.append( stat + hitmiss )

        scoreVal( ( cmsStatsRaw if stat in nonNormedStats else cmsStatsNorm )[ stat ], cmsStatsNorm.Pos,
                  stat_nbin[ stat ], stat_start[ stat ], stat_end[ stat ],
                  hitmissCols[ stat+'Hits' ], hitmissCols[ stat+'Misses' ], selpos )

    DotData( names = hitmissNames, Columns = [ hitmissCols[ n ] for n in hitmissNames ] ).save( ranksFN )

    with open(summaryFN, 'w') as  outfile:
        outfile.write('N causal: 0 N missed: 0\n')  # dummy line
        tabwrite( outfile, 'Stat', 'Start', 'End', 'N_bin')
        for stat in stat_nbin.keys():
            tabwrite( outfile, stat, stat_start[stat], stat_end[stat], stat_nbin[stat] )



def analsim( Ddata, scenario, putativeMutPop = None, thinSfx = '', likesSfx = '', condName = '',
             stat_start = CMSBins.stat_start, stat_end = CMSBins.stat_end, stat_nbin = CMSBins.stat_nbin,
             stats = CMSBins.CMSstats, nonNormedStats = CMSBins.nonNormedStats,
             normalizeWithinReplicas = False, nreplicas = 100,
             replicaTables = (), replicaCond = '', statsSfx = '',
             getio = None ):

    """Compute histograms of CMS component statistics for SNPs in chosen replicas of one scenario."""

    #if not os.path.exists(Ddata + '/power_GW/'): os.makedirs(Ddata +'/' +  name + '/')

    if not putativeMutPop: putativeMutPop = scenario.mutPop

    snpStatsDir = os.path.join( Ddata, 'snpStats' + thinSfx, scenario.scenDir() )
#    cmsStatsNormFN = GetCreates( normalizeStatsWithinReplicas if normalizeWithinReplicas else normalizeStats,
#                                 **Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx' ) ).containing( 'cmsStatsNorm' )

    cmsStatsNormFN = os.path.join( snpStatsDir,
                                   'cmsStatsRaw_%d_normed%s.tsv' % ( putativeMutPop, 'Local' if normalizeWithinReplicas else '' ) )
    
    cmsStatsRawFN = GetCreates( computeCMSstats, **Dict( 'Ddata thinSfx scenario putativeMutPop statsSfx' ) )[0]

    ranksFN = os.path.join( Ddata, 'scenStats' + thinSfx, scenario.scenDir(),
                            AddFileSfx( 'ranks.data/', likesSfx, putativeMutPop, statsSfx, condName ) )
    summaryFN = os.path.join( Ddata, AddFileSfx( 'Summary.txt', scenario.scenDir(),
                                                 likesSfx, putativeMutPop, statsSfx, condName ) )

    replicaTables = MakeSeq( replicaTables )

    replicaTableFiles = tuple(  os.path.join( Ddata, 'replicastats', scenario.scenDir(),
                                              replicaTable + '.tsv' )
                                for replicaTable in replicaTables )

    if getio: return dict( depends_on = ( cmsStatsNormFN, cmsStatsRawFN ) + replicaTableFiles,
                           creates = ( ranksFN, summaryFN ),
                           mediumRuleNameSfx = ( scenario.scenDir(), likesSfx, putativeMutPop, condName ),
                           attrs = dict( scenario = scenario.scenDir() ),
                           fileDescrs = { ranksFN: 'Histograms of CMS component stats for causal/non-causal SNPs',
                                          summaryFN: 'Bin boundaries used for CMS component histograms' } )
    
    #rep_freqs = DotData(SVPath = '../Data/Shari_Data/sim/rep_freq.tsv', names = ['rep','low','high'])

    cmsStatsNorm = IDotData( Path = cmsStatsNormFN )
    cmsStatsRaw = IDotData( Path = cmsStatsRawFN )

    #assert len( cmsStatsNorm ) == len( cmsStatsRaw )
    #assert all( cmsStatsNorm.Pos == cmsStatsRaw.Pos )
    #assert all( cmsStatsNorm.Chrom == cmsStatsRaw.Chrom )

    #dbg( 'len(cmsStatsNorm)' )

    replicasToUse = None
    if replicaCond and replicaCond != 'True':

        # identify replicas matching the condition

        replicaCondExpr = compile_expr( replicaCond )
        replicaTableVals = [ DotData( SVPath = f ) for f in replicaTableFiles ]

        assert all([ len( rtv ) == nreplicas for rtv in replicaTableVals ])

        replicasToUse = array( [ eval( replicaCondExpr, globals(), dict( zip( replicaTables, replicaTableRows ) ) )
                                 for replicaTableRows in itertools.izip( *replicaTableVals ) ] ) if replicaTables \
                                 else repeat( eval( replicaCondExpr, globals(), {} ), nreplicas )

    hitmissCols = {}
    hitmissNames = []

    assert len( set( stat_nbin.values() ) ) == 1
    assert sorted( stat_start.keys() ) == sorted( stat_end.keys() ) == sorted( stat_nbin.keys() )

    for stat in stats:
        for hitmiss in 'Hits', 'Misses':
            hitmissCols[ stat + hitmiss ] = zeros( stat_nbin[ stat ] + CMSBinsBase.maxSpecialBins, dtype = int )
            hitmissNames.append( stat + hitmiss )
    
    stat_starts = [ stat_start[ stat ] for stat in stats ]
    stat_nbinss = [ stat_nbin[ stat ] for stat in stats ]
    stat_ends = [ stat_end[ stat ] for stat in stats ]
    bin_sizes = [ float(stat_end[ stat ] - stat_start[ stat ]) / stat_nbin[ stat ] for stat in stats ]

    normedOrNots = [ int( stat not in nonNormedStats ) for stat in stats ]

    hmColsHits = [ hitmissCols[ stat + 'Hits' ] for stat in stats ]
    hmColsMisses = [ hitmissCols[ stat + 'Misses' ] for stat in stats ]
            
    for rowNum, r in enumerate( itertools.izip( cmsStatsRaw, cmsStatsNorm ) ):
        if replicasToUse is not None and not replicasToUse[ int( r[0].Chrom ) ]: continue

        for stat, start, bin_size, nbin, normedOrNot, hmColsHit, hmColsMiss \
                in zip( stats, stat_starts, bin_sizes, stat_nbinss, normedOrNots,
                        hmColsHits, hmColsMisses ):
            statVal = r[ normedOrNot ][ stat ]
            nanReason = r[ 0 ][ stat + '_nanReason' ] if CMSBins.stat_numSpecialBins[ stat ] else -1

            bin = ( nbin + nanReason ) if nanReason != -1 else \
                ( -1 if not np.isfinite( statVal ) else \
                   clamp( int( (statVal - start ) / bin_size ), 0, nbin-1 ) )

            if bin >= 0: ( hmColsHit if not scenario.isNeutral() and r[0].Pos == CAUSAL_POS else hmColsMiss )[ bin ] += 1

    DotData( names = hitmissNames, Columns = [ hitmissCols[ n ] for n in hitmissNames ] ).save( ranksFN )

    with open(summaryFN, 'w') as  outfile:
        outfile.write('N causal: 0 N missed: 0\n')  # dummy line
        tabwrite( outfile, 'Stat', 'Start', 'End', 'N_bin')
        for stat in stat_nbin.keys():
            tabwrite( outfile, stat, stat_start[stat], stat_end[stat], stat_nbin[stat] )

def likes( Ddata = '../Data/Shari_Data/sim/', mutAges = AllAges, mutPops = AllPops,
           mutFreqs = AllFreqs, thinSfx = '', nearCausalSfx = '',
           putativeMutPop = pn_WestAfrican, likesSfx = '', secondLikesSfx = '', condName = '', statsSfx = '',
           stats = CMSBins.CMSstats, popSpecificLikes = False, getio = None ):
    """Creates the likelihood tables: calculates for each statistic for each bin the probability that for a given SNP that
    statistic is in that bin if the SNP is causal, and if the SNP is non-causal.  From these tables, and the statistic
    values for each SNP, the function 'complike' calculates each SNP's individual and composite likelihood of being causal.
    """

    if popSpecificLikes: mutPops = ( putativeMutPop, )
    selScenarios = GetSelectionScenarios( **Dict( 'mutFreqs mutPops mutAges' ) )
    scenarios = tuple( selScenarios ) + ( NeutralScenario(), )


    scen2ranksFile = dict([ ( scenario,
                              GetCreates( analsim, scenario = scenario,
                                          putativeMutPop = putativeMutPop if scenario.isNeutral() else scenario.mutPop,
                                          condName = '' if scenario.isNeutral() else condName,
                                          **Dict( 'Ddata statsSfx likesSfx' ) ).containing( 'ranks' ) )
                            for scenario in scenarios ])
        
    likesFiles = dict([ ( (likes, which), os.path.join( Ddata, 'likes',
                                                       AddFileSfx( which + likes + '.tsv', likesSfx, secondLikesSfx,
                                                                   putativeMutPop, statsSfx, condName ) ) )
                       for likes in ( 'Likes', '' ) for which in ( 'hits', 'miss', 'region' ) ])

    if getio: return dict( depends_on = scen2ranksFile.values(), creates = likesFiles.values(),
                           attrs = dict( piperun_short = True ),
                           mediumRuleNameSfx = ( likesSfx, putativeMutPop, condName ) )

    # var: hits - for each statistic, for each value range of that statistic, how many causal snps have that statistic
    #  in that value range; accumulated over all selection scenarios.  map from statistic name to a histogram
    #  (an array of bins, giving the number of values in each bin). 
    # var: misses - for each statistic, for each value range of that statistic, how many non-causal snps have that statistic
    #  in that value range; accumulated over all selection scenarios.
    hits = {}
    misses = {}
    region = {}

    print 'stats are ', stats
    for stat in stats:
        nbin_here = CMSBins.stat_nbin[stat] + CMSBinsBase.maxSpecialBins
        hits[stat] = zeros( nbin_here, dtype = int )
        misses[stat] = zeros( nbin_here, dtype = int )
        region[stat] = zeros( nbin_here, dtype = int )

    # var: stat_nbin - map from statistic name to number of bins
    # var: stat_start - map from statistic name to start of binned range
    # var: stat_end - map from statistic name to end of binned range.
    stat_nbin = CMSBins.stat_nbin
    stat_start = CMSBins.stat_start
    stat_end = CMSBins.stat_end
    
    neutralData = DotData( Path = scen2ranksFile[ NeutralScenario() ] )

    for stat in stats:
        misses[stat] += array(neutralData[stat + 'Misses'])

    for scenario in selScenarios:

        Data = DotData(Path = scen2ranksFile[ scenario ])
        
        for stat in stats:
            hits[stat] += array(Data[stat + 'Hits'])
            region[stat] += array(Data[stat + 'Misses'])

    # var: nhitSNPs - map from statistic to the total number of causal SNPs counted for that statistic
    #    ( i.e. total number of causal SNPs that contributed a count to one of the bins for that statistic ).
    # var: nmissSNPs - same for non-causal SNPs
    nhitSNPs = {}
    nmissSNPs = {}
    nregionSNPs = {}

    for stat in stats:
        nhitSNPs[stat] = sum(hits[stat])
        nmissSNPs[stat] = sum(misses[stat])
        nregionSNPs[stat] = sum(region[stat])

    # var: outhits - a histogram for each statistic, which for each bin gives the fraction of causal SNPs
    #     that have that statistic in that bin.  (so, this is the conditional probability that if a SNP
    #     is causal, it has that statistic in that bin).   Note that this is for SNPs in all replicas 
    #     of all scenarios.
    # var: outmisses - same for non-causal.
    outhits = {}
    outmiss = {}
    outregion = {}

    for stat in stats:
        nbins_orig = stat_nbin[stat]
        nbins = nbins_orig + CMSBins.maxSpecialBins
        dbg('stat nbins_orig nbins')
        outhits[stat] = zeros(nbins)
        outmiss[stat] = zeros(nbins)
        outregion[stat] = zeros(nbins)

        for bin in range(nbins):
            # Deal with zeros

            for vals, outvals, nvalSNPs in ( hits, outhits, nhitSNPs ), ( misses, outmiss, nmissSNPs ), \
                    ( region, outregion, nregionSNPs ):
            
                if vals[stat][bin] == 0:
                    if bin >= nbins_orig or \
                            ((bin == 0 and vals[stat][1] == 0) or (bin == nbins_orig-1 and vals[stat][nbins_orig - 2] == 0) or \
                                 0 < bin < nbins_orig-1  and  (vals[stat][bin - 1] == 0 or vals[stat][bin + 1] == 0)):
                        outvals[stat][bin] = 1e-10
                    else:
                        outvals[stat][bin] = 1. / (nvalSNPs[stat] + 1)

                elif vals[stat][bin] == 1:
                    if  bin >= nbins_orig or \
                            ( 1 < bin < nbins_orig-2 and vals[stat][bin - 1] == 0 and vals[stat][bin - 2 ] == 0 and \
                                  vals[stat][bin + 1] == 0 and vals[stat][bin + 2 ] == 0 ):
                        outvals[stat][bin] = 1e-10
                    else:
                        outvals[stat][bin] = 1. / (nvalSNPs[stat] + 1)

                else:
                    outvals[stat][bin] = float(vals[stat][bin] + 1) / (nvalSNPs[stat] + 1)
    

    DotData(Columns = [outhits[stat] for stat in stats],names=stats).saveToSV(likesFiles[('Likes','hits')])
    DotData(Columns=[outmiss[stat] for stat in stats],names=stats).saveToSV(likesFiles[('Likes','miss')])
    DotData(Columns=[outregion[stat] for stat in stats],names=stats).saveToSV(likesFiles[('Likes','region')])

    DotData(Columns = [hits[stat] for stat in stats],names=stats).saveToSV(likesFiles[('','hits')])
    DotData(Columns=[misses[stat] for stat in stats],names=stats).saveToSV(likesFiles[('','miss')])
    DotData(Columns=[region[stat] for stat in stats],names=stats).saveToSV(likesFiles[('','region')])
    



def likes2( Ddata = '../Data/Shari_Data/sim/', mutAges = AllAges, mutPops = AllPops,
            mutFreqs = AllFreqs, thinSfx = '', nearCausalSfx = '',
            putativeMutPop = pn_WestAfrican, likesSfx = '', secondLikesSfx = '', condName = '', statsSfx = '',
            stats = CMSBins.CMSstats, popSpecificLikes = False, likesRatiosFile = 'likesRatios.tsv', getio = None ):
    """Creates the likelihood tables: calculates for each statistic for each bin the probability that for a given SNP that
    statistic is in that bin if the SNP is causal, and if the SNP is non-causal.  From these tables, and the statistic
    values for each SNP, the function 'complike' calculates each SNP's individual and composite likelihood of being causal.
    """

    if popSpecificLikes: mutPops = ( putativeMutPop, )
    selScenarios = GetSelectionScenarios( **Dict( 'mutFreqs mutPops mutAges' ) )
    scenarios = tuple( selScenarios ) + ( NeutralScenario(), )


    scen2ranksFile = dict([ ( scenario,
                              GetCreates( analsim, scenario = scenario,
                                          putativeMutPop = putativeMutPop if scenario.isNeutral() else scenario.mutPop,
                                          condName = '' if scenario.isNeutral() else condName,
                                          **Dict( 'Ddata statsSfx likesSfx' ) ).containing( 'ranks' ) )
                            for scenario in scenarios ])
        
    likesFiles = dict([ ( (likes, which), os.path.join( Ddata, 'likes',
                                                       AddFileSfx( which + likes + '.tsv', likesSfx, secondLikesSfx,
                                                                   putativeMutPop, statsSfx, condName ) ) )
                       for likes in ( 'Likes', '' ) for which in ( 'hits', 'miss', 'region' ) ])

    if getio: return dict( depends_on = scen2ranksFile.values(), creates = likesFiles.values() + [ likesRatiosFile ],
                           attrs = dict( piperun_short = True ),
                           mediumRuleNameSfx = ( likesSfx, putativeMutPop, condName ) )

    # var: hits - for each statistic, for each value range of that statistic, how many causal snps have that statistic
    #  in that value range; accumulated over all selection scenarios.  map from statistic name to a histogram
    #  (an array of bins, giving the number of values in each bin). 
    # var: misses - for each statistic, for each value range of that statistic, how many non-causal snps have that statistic
    #  in that value range; accumulated over all selection scenarios.
    hits = {}
    misses = {}
    region = {}

    print 'stats are ', stats
    for stat in stats:
        nbin_here = CMSBins.stat_nbin[stat] + CMSBinsBase.maxSpecialBins
        hits[stat] = zeros( nbin_here, dtype = int )
        misses[stat] = zeros( nbin_here, dtype = int )
        region[stat] = zeros( nbin_here, dtype = int )

    # var: stat_nbin - map from statistic name to number of bins
    # var: stat_start - map from statistic name to start of binned range
    # var: stat_end - map from statistic name to end of binned range.
    stat_nbin = CMSBins.stat_nbin
    stat_start = CMSBins.stat_start
    stat_end = CMSBins.stat_end
    
    neutralData = DotData( Path = scen2ranksFile[ NeutralScenario() ] )

    for stat in stats:
        misses[stat] += array(neutralData[stat + 'Misses'])

    for scenario in selScenarios:

        Data = DotData(Path = scen2ranksFile[ scenario ])
        
        for stat in stats:
            hits[stat] += array(Data[stat + 'Hits'])
            region[stat] += array(Data[stat + 'Misses'])

    # var: nhitSNPs - map from statistic to the total number of causal SNPs counted for that statistic
    #    ( i.e. total number of causal SNPs that contributed a count to one of the bins for that statistic ).
    # var: nmissSNPs - same for non-causal SNPs
    nhitSNPs = {}
    nmissSNPs = {}
    nregionSNPs = {}

    for stat in stats:
        nhitSNPs[stat] = sum(hits[stat])
        nmissSNPs[stat] = sum(misses[stat])
        nregionSNPs[stat] = sum(region[stat])

    # var: outhits - a histogram for each statistic, which for each bin gives the fraction of causal SNPs
    #     that have that statistic in that bin.  (so, this is the conditional probability that if a SNP
    #     is causal, it has that statistic in that bin).   Note that this is for SNPs in all replicas 
    #     of all scenarios.
    # var: outmisses - same for non-causal.
    outhits = {}
    outmiss = {}
    outregion = {}

    likesRatios = []

    for stat in stats:
        nbins_orig = stat_nbin[stat]
        nbins = nbins_orig + CMSBins.maxSpecialBins
        dbg('stat nbins_orig nbins')
        outhits[stat] = zeros(nbins)
        outmiss[stat] = zeros(nbins)
        outregion[stat] = zeros(nbins)

        for bin in range(nbins):
            # Deal with zeros

            for vals, outvals, nvalSNPs, kind in ( hits, outhits, nhitSNPs, 'hits' ), ( misses, outmiss, nmissSNPs, 'misses' ), \
                    ( region, outregion, nregionSNPs, 'region' ):
            
                if vals[stat][bin] == 0:
                    if bin >= nbins_orig or \
                            ((bin == 0 and vals[stat][1] == 0) or (bin == nbins_orig-1 and vals[stat][nbins_orig - 2] == 0) or \
                                 0 < bin < nbins_orig-1  and  (vals[stat][bin - 1] == 0 or vals[stat][bin + 1] == 0)):
                        outvals[stat][bin] = 1e-10
                    else:
                        outvals[stat][bin] = 1. / (nvalSNPs[stat] + 1)

                elif vals[stat][bin] == 1:
                    if  bin >= nbins_orig or \
                            ( 1 < bin < nbins_orig-2 and vals[stat][bin - 1] == 0 and vals[stat][bin - 2 ] == 0 and \
                                  vals[stat][bin + 1] == 0 and vals[stat][bin + 2 ] == 0 ):
                        outvals[stat][bin] = 1e-10
                    else:
                        outvals[stat][bin] = 1. / (nvalSNPs[stat] + 1)

                else:
                    outvals[stat][bin] = float(vals[stat][bin] + 1) / (nvalSNPs[stat] + 1)
                    likesRatios.append( ( stat, bin, kind, vals[stat][bin] + 1, nvalSNPs[stat] + 1 ) ) 
    

    DotData(Columns = [outhits[stat] for stat in stats],names=stats).saveToSV(likesFiles[('Likes','hits')])
    DotData(Columns=[outmiss[stat] for stat in stats],names=stats).saveToSV(likesFiles[('Likes','miss')])
    DotData(Columns=[outregion[stat] for stat in stats],names=stats).saveToSV(likesFiles[('Likes','region')])

    DotData(Columns = [hits[stat] for stat in stats],names=stats).saveToSV(likesFiles[('','hits')])
    DotData(Columns=[misses[stat] for stat in stats],names=stats).saveToSV(likesFiles[('','miss')])
    DotData(Columns=[region[stat] for stat in stats],names=stats).saveToSV(likesFiles[('','region')])

    DotData( names = ( 'stat', 'bin', 'kind', 'k', 'n' ), Records = likesRatios ).saveToSV( likesRatiosFile )
    


    
def combineBinCounts( Ddata, scenario, putativeMutPop = None, thinSfx = '', likesSfx = '', condName = '',
                      getio = None ):
    """ """
    pass
    

def saveLikesTable( Ddata, likesSfx, likesBins, normalizeWithinReplicas, likesTable, putativeMutPop = None,
                    condName = '', statsSfx = '',
                    likesTableComment = None, likesTableCreationParams = None, secondLikesSfx = '',
                    pop2sampleSize = None, useRealPaths = False, getio = None ):
    """Save the likes table to a file.

    Params:

       Input:

       Output:

          likesTable - the filename of the file to which to write the likes table.
            Note that we only write references to likes files and the bins file, not their content.

       Misc:

         useRealPaths - if True, absolute paths to likes files and the bins file are recorded instead of relative paths.
            also, symlinks are resolved.

    """

    if not os.path.dirname( likesTable ): likesTable = os.path.join( Ddata, likesTable )

    likesTable = AddFileSfx( likesTable, likesSfx, secondLikesSfx, putativeMutPop, statsSfx, condName )

    headings = ( 'param', 'val' )
    recs = []
    depends_on = []
    for which in 'hits', 'miss', 'region':
        likesFileName = \
            GetCreates( likes, **Dict( 'Ddata likesSfx secondLikesSfx putativeMutPop condName statsSfx' ) ).containing( which + 'Likes' )
        if useRealPaths: likesFileName = os.path.realpath( likesFileName )
        recs.append( ( which + 'Likes', likesFileName ) )
        depends_on.append( likesFileName )
        
    if likesBins and not os.path.dirname( likesBins ): likesBins = os.path.join( Ddata, likesBins )
    depends_on.append( likesBins )
    if useRealPaths: likesBins = os.path.realpath( likesBins )
    
    recs.append( ( 'likesBins', likesBins ) )
    recs.append( ( 'normalizeWithinReplicas', normalizeWithinReplicas ) )
    if pop2sampleSize: recs.append( ( 'pop2sampleSize', pop2sampleSize ) )

    if likesTableCreationParams: recs.append( ( 'creationParams', likesTableCreationParams ) )
    if likesTableComment: recs.append( ( 'comment', likesTableComment ) )
        
    if getio: return dict( depends_on = depends_on, creates = likesTable,
                           attrs = dict( piperun_short = True ),
                           mediumRuleNameSfx = ( likesSfx, putativeMutPop, condName ) )
    
    IDotData( names = headings, Records = recs ).save( likesTable )


def plotLikesTable( hitsLikesFile, missLikesFile, plotFile, normalizeWithinReplicas, regionLikesFile = None,
                    likesBinsFile = None, statsSfx = '',
                    likesTable = None, nonNormedStats = (), likesSfx = '', secondLikesSfx = '', condName = '', putativeMutPop = None,
                    includeSpecialBins = True,
                    getio = None ):
    """Visually plot a likes table.
    """

    if getio: return dict( depends_on = filter( None, ( hitsLikesFile, missLikesFile, regionLikesFile, likesTable, likesBinsFile ) ),
                           creates = plotFile,
                           mediumRuleNameSfx = ( likesSfx, putativeMutPop, condName ),
                           attrs = Dict( 'includeSpecialBins', piperun_short = True ) )

    hitsLikes = IDotData( hitsLikesFile )
    missLikes = IDotData( missLikesFile )
    if regionLikesFile: regionLikes = IDotData( regionLikesFile )

    pp.figure( figsize = ( 16, 18 ) )

    if not likesBinsFile:
        if normalizeWithinReplicas:
            stat_start, stat_end, stat_nbins = CMSBinsLocal.stat_start, CMSBinsLocal.stat_end, CMSBinsLocal.stat_nbin
        else:
            stat_start, stat_end, stat_nbins = CMSBins.stat_start, CMSBins.stat_end, CMSBins.stat_nbin
    else:
        stat_start, stat_end, stat_nbins = LoadBins( likesBinsFile )

    regionLine = None
    for statNum, stat in enumerate( hitsLikes.headings ):

        rawStep = 1.0 / len( hitsLikes.headings ) * 0.93

        rawBottom = rawStep * statNum
        rawTop = rawBottom + rawStep

        r = ( 0.1, 0.05 + rawBottom, 0.8, rawStep * 0.6 )
        dbg( 'r' )
        pp.axes( r )

        pp.title( stat + ( ' (non-normed)' if stat in nonNormedStats else '' ) )

        assert len( hitsLikes ) == len( missLikes ) == stat_nbins[ stat ] + CMSBins.maxSpecialBins

        
        binSize = ( stat_end[stat] - stat_start[stat] ) / stat_nbins[stat]
        binStarts = [ stat_start[stat] + binSize * i
                          for i in range( stat_nbins[ stat ] + ( CMSBins.stat_numSpecialBins[ stat ] if includeSpecialBins else 0 ) )  ]

        pp.gca().set_xticks( binStarts )
        pp.gca().set_xticklabels( [ '%.2f' % b for b in binStarts[: stat_nbins[stat] ] ] +
                                  ( list( DictGet( CMSBins.stat_specialBinNames,
                                                   stat, () ) ) if includeSpecialBins else [] ),
                                  rotation = 'vertical' )
#        pp.gca().set_xticklabels( map( str, binStarts ) + [ 's%d' % i for i in range( CMSBins.stat_numSpecialBins[ stat ] ) ] )
        
        dbg( 'stat binStarts' )
        hitsLine, = pp.plot( binStarts , hitsLikes[ stat ][:len( binStarts )], 'r-' )
        missLine, = pp.plot( binStarts , missLikes[ stat ][:len( binStarts )], 'g-' )
        if regionLikesFile: regionLine, = pp.plot( binStarts, regionLikes[ stat ][:len(binStarts)], 'b-' )

    pp.figlegend( filter( None, ( hitsLine, missLine, regionLine) ),
                  ( 'selected SNPs', 'neutral SNPs in neutral regions' ) +
                  ( ( 'neutral SNPs in selected regions', ) if regionLikesFile else () ),
                  'upper center' )

    pp.savefig( plotFile )


def plotLikesTablePair( likesTableFNs,
                        plotFile, nonNormedStats = (),
                        includeSpecialBins = True,
                        getio = None ):
    """Visually plot a likes table.
    """

    if getio: return dict( depends_on = likesTableFNs,
                           creates = plotFile,
                           attrs = dict( piperun_short = True ) )

    likesTable = map( LoadLikesTable, likesTableFNs )

    hitsLikes = [ IDotData( likesTable[ i ].hitsLikes ) for i in range( 2 ) ]
    missLikes = [ IDotData( likesTable[ i ].missLikes ) for i in range( 2 ) ]
    regionLikes = [ IDotData( likesTable[ i ].regionLikes ) for i in range( 2 ) ]
    
    pp.figure( figsize = ( 16, 18 ) )

    stat_start, stat_end, stat_nbins = LoadBins( likesTable[0].likesBins )
    stat_start1, stat_end1, stat_nbins1 = LoadBins( likesTable[1].likesBins )

    assert( stat_start == stat_start1 )
    assert( stat_end == stat_end1 )
    assert( stat_nbins == stat_nbins1 )
    assert( hitsLikes[0].headings == hitsLikes[1].headings )
    assert( missLikes[0].headings == missLikes[1].headings )
    assert( regionLikes[0].headings == regionLikes[1].headings )

    regionLine = None
    for statNum, stat in enumerate( hitsLikes[0].headings ):

        rawStep = 1.0 / len( hitsLikes[0].headings ) * 0.93

        rawBottom = rawStep * statNum
        rawTop = rawBottom + rawStep

        r = ( 0.1, 0.05 + rawBottom, 0.8, rawStep * 0.6 )
        dbg( 'r' )
        pp.axes( r )

        pp.title( stat + ( ' (non-normed)' if stat in nonNormedStats else '' ) )

        assert len( hitsLikes[0] ) == len( missLikes[0] ) == stat_nbins[ stat ] + CMSBins.maxSpecialBins

        
        binSize = ( stat_end[stat] - stat_start[stat] ) / stat_nbins[stat]
        binStarts = [ stat_start[stat] + binSize * i
                          for i in range( stat_nbins[ stat ] + ( CMSBins.stat_numSpecialBins[ stat ] if includeSpecialBins else 0 ) )  ]

        pp.gca().set_xticks( binStarts )
        pp.gca().set_xticklabels( [ '%.2f' % b for b in binStarts[: stat_nbins[stat] ] ] +
                                  ( list( DictGet( CMSBins.stat_specialBinNames,
                                                   stat, () ) ) if includeSpecialBins else [] ),
                                  rotation = 'vertical' )
#        pp.gca().set_xticklabels( map( str, binStarts ) + [ 's%d' % i for i in range( CMSBins.stat_numSpecialBins[ stat ] ) ] )
        
        dbg( 'stat binStarts' )
        hitsLine = [ None, None ]
        missLine = [ None, None ]
        regionLine = [ None, None ]
        for i, style in ( ( 0, '-' ), ( 1, ':' ) ):
            hitsLine[i], = pp.plot( binStarts , hitsLikes[i][ stat ][:len( binStarts )], 'r' + style  )
            missLine[i], = pp.plot( binStarts , missLikes[i][ stat ][:len( binStarts )], 'g' + style )
            regionLine[i], = pp.plot( binStarts, regionLikes[i][ stat ][:len(binStarts)], 'b' + style )

    pp.figlegend( filter( None, ( hitsLine[0], missLine[0], regionLine[0],
                                  hitsLine[1], missLine[1], regionLine[1] ) ),
                  ( 'selected SNPs 1', 'neutral SNPs in neutral regions 1', 'region snps 1',
                    'selected SNPs 2', 'neutral SNPs in neutral regions 2', 'region snps 2',
                ),
                  'upper center' )

    pp.savefig( plotFile )

# end: def plotLikesTablePair()    


def writeBinCounts( Ddata, complikeSfx, scenario, normalizeWithinReplicas, putativeMutPop = None, getio = None ):
    """Write bin counts, as well as bin definitions, as readable tables."""

    if putativeMutPop is None and not scenario.isNeutral(): putativeMutPop = scenario.mutPop
    
    inFN = GetCreates( analsim, **Dict( 'Ddata scenario putativeMutPop', likesSfx = complikeSfx ) ).containing( 'ranks' )
    outFN = AddFileSfx( inFN, 'table' )
    dbg( 'inFN outFN' )

    if getio: return dict( depends_on = inFN, creates = outFN, attrs = dict( scenario = scenario.scenDir(),
                                                                             putativeMutPop = putativeMutPop, piperun_short = True ) )

    binCountsFile = IDotData( inFN )

    if normalizeWithinReplicas:
        stat_start, stat_end, stat_nbins = CMSBinsLocal.stat_start, CMSBinsLocal.stat_end, CMSBinsLocal.stat_nbin
    else:
        stat_start, stat_end, stat_nbins = CMSBins.stat_start, CMSBins.stat_end, CMSBins.stat_nbin


    binStartss = []
    binEndss = []
    binDescrss = []

    numBinsAll = len( binCountsFile )
        
    for stat in CMSBins.CMSstats:
        binSize = ( stat_end[stat] - stat_start[stat] ) / stat_nbins[stat]
        binStarts = [ ( stat_start[stat] + binSize * i ) if i < stat_nbins[ stat ] else np.nan
                      for i in range( numBinsAll )  ]
        binEnds = [ b + binSize for b in binStarts ]

        binStartss.append( binStarts )
        binEndss.append( binEnds )

        binDescrs = [ CMSBins.stat_specialBinNames[ stat ][ i - stat_nbins[ stat ] ] \
                          if stat in CMSBins.stat_specialBinNames and i >= stat_nbins[ stat ] \
                          and i - stat_nbins[stat] < CMSBins.stat_numSpecialBins[ stat ] 
                      else 'interval'
                      for i in range( numBinsAll )  ]
        binDescrss.append( binDescrs )

    Columns = reduce( operator.concat, [ ( binCountsFile[ stat + 'Hits' ], binCountsFile[ stat + 'Misses' ],
                                           binStarts, binEnds, binDescrs )
                                         for stat, binStarts, binEnds, binDescrs
                                         in zip( CMSBins.CMSstats,
                                                 binStartss, binEndss, binDescrss ) ] )
    dbg( 'map(len,Columns)' )
    
    DotData( names = reduce( operator.concat, [ [ stat + n for n in ( 'Hits', 'Misses', 'BinStart', 'BinEnd', 'BinDescr' ) ]
                                                for stat in CMSBins.CMSstats ] ),
             Columns = Columns ).save( outFN )


def writeLikesCounts( Ddata, complikeSfx, putativeMutPop, normalizeWithinReplicas, condName = '',
                      stats = None, getio = None ):
    """Write bin counts, as well as bin definitions, as readable tables."""

    dbg( 'likes' )
    if not stats: stats = CMSBins.CMSstats
    inFNs = GetCreates( likes, **Dict( 'Ddata condName putativeMutPop', likesSfx = complikeSfx ) )
    
    outFN = os.path.join( Ddata, 'likes', AddFileSfx( 'likesInfo.data/', complikeSfx, putativeMutPop, condName, *stats ) )
    dbg( 'inFNs outFN' )

    if getio: return dict( depends_on = inFNs, creates = outFN, attrs = Dict( 'stats complikeSfx putativeMutPop',
                                                                              piperun_short = True ) )

    likesFilesFNs = dict( ( (likesKind, which), os.path.join( Ddata, 'likes',
                                                              AddFileSfx( which + likesKind + '.tsv', complikeSfx,
                                                                      putativeMutPop, condName ) ) )
                          for likesKind in ( 'Likes', '' ) for which in ( 'hits', 'miss', 'region' ) )

    likesFiles = dict( [ ( k, IDotData( FN ) ) for k, FN in likesFilesFNs.items() ] )
    
    if normalizeWithinReplicas:
        stat_start, stat_end, stat_nbins = CMSBinsLocal.stat_start, CMSBinsLocal.stat_end, CMSBinsLocal.stat_nbin
    else:
        stat_start, stat_end, stat_nbins = CMSBins.stat_start, CMSBins.stat_end, CMSBins.stat_nbin

    Columns = []
    colNames = []

    numBinsAll =  max( stat_nbins.values() ) + CMSBins.maxSpecialBins  #len( likesFiles[ ('', 'hits' ) ] )

    for stat in stats:
        binSize = ( stat_end[stat] - stat_start[stat] ) / stat_nbins[stat]
        binStarts = [ ( stat_start[stat] + binSize * i ) if i < stat_nbins[ stat ] else np.nan
                      for i in range( numBinsAll )  ]
        binEnds = [ b + binSize for b in binStarts ]

        Columns.append( binStarts )
        colNames.append( stat + '_binStarts' )
        Columns.append( binEnds )
        colNames.append( stat + '_binEnds' )

        binDescrs = [ CMSBins.stat_specialBinNames[ stat ][ i - stat_nbins[ stat ] ] \
                          if stat in CMSBins.stat_specialBinNames and i >= stat_nbins[ stat ] \
                          and i - stat_nbins[stat] < CMSBins.stat_numSpecialBins[ stat ] 
                      else ( 'interval' if i < stat_nbins[stat] else 'None' )
                      for i in range( numBinsAll )  ]
        Columns.append( binDescrs )
        colNames.append( stat + '_binDescrs' )
        likesCols = {}
        for likesKind in ( '', 'Likes' ):
            for which in ( 'hits', 'miss', 'region' ):
                likesCol = list( likesFiles[ ( likesKind, which ) ][ stat ] )
                if len( likesCol ) < numBinsAll:
                    likesCol += [ 1e-10 if likesKind else 0 for i in range( len( likesCol ), numBinsAll ) ]
                dbg( 'likesKind which len(likesCol)' )
                likesCol = np.array( likesCol )
                likesCols[ ( likesKind, which ) ] = likesCol
                Columns.append( likesCol )
                colNames.append( stat + Sfx( which, likesKind ) )

        colNames += [ stat + '_' + which + 'Tot' for which in ( 'hits', 'miss', 'region' ) ]
        Columns += [ np.repeat( np.sum( likesFiles[ ( '', which ) ][ stat ] ), len( binStarts ) ) for which in ( 'hits', 'miss', 'region' ) ]

        ratioCol = likesCols[ ( 'Likes', 'hits' ) ] / likesCols[ ( 'Likes', 'miss' ) ]
        Columns.append( ratioCol )
        colNames.append( stat + '_ratioCausal2Neutral' )
        Columns.append( np.log( ratioCol ) )
        colNames.append( stat + '_ratioLogCausal2Neutral' )

        ratioCol = likesCols[ ( 'Likes', 'region' ) ] / likesCols[ ( 'Likes', 'miss' ) ]
        Columns.append( ratioCol )
        colNames.append( stat + '_ratioRegion2Neutral' )
        Columns.append( np.log( ratioCol ) )
        colNames.append( stat + '_ratioLogRegion2Neutral' )

    dbg( 'zip(colNames,map(len,Columns))' )
    DotData( names = colNames, Columns = Columns ).save( outFN )

    
def compareLikesTables( likesTableDefs, putativeMutPop, normalizeWithinReplicas, cmpOp = 'diff', getio = None ):
    """Plot a comparison of two or more likes tables."""

    likesDataFNs = [ GetCreates( writeLikesCounts, **Dict( 'putativeMutPop normalizeWithinReplicas', **likesTableDef ) )[0]
                     for likesTableDef in likesTableDefs ]

    outFN = os.path.join( likesTableDefs[0]['Ddata'], 'graphs',
                          AddFileSfx( 'likesCmp.svg', putativeMutPop, cmpOp, *map( operator.itemgetter( 'complikeSfx' ), likesTableDefs ) ) )

    if getio: return dict( depends_on = likesDataFNs, creates = outFN, attrs = Dict( 'putativeMutPop cmpOp', piperun_short = True ) )

    likesDataFiles = map( IDotData, likesDataFNs )

    pp.figure( figsize = ( 16, 18 ) )

    if normalizeWithinReplicas:
        stat_start, stat_end, stat_nbins = CMSBinsLocal.stat_start, CMSBinsLocal.stat_end, CMSBinsLocal.stat_nbin
    else:
        stat_start, stat_end, stat_nbins = CMSBins.stat_start, CMSBins.stat_end, CMSBins.stat_nbin

    regionLine = None
    for statNum, stat in enumerate( CMSBins.CMSstats ):

        rawStep = 1.0 / len( CMSBins.CMSstats ) * 0.93

        rawBottom = rawStep * statNum
        rawTop = rawBottom + rawStep

        r = ( 0.1, 0.05 + rawBottom, 0.8, rawStep * 0.6 )
        dbg( 'r' )
        pp.axes( r )

        pp.title( stat + ( ' (non-normed)' if stat in CMSBins.nonNormedStats else '' ) )

        #assert len( hitsLikes ) == len( missLikes ) == stat_nbins[ stat ] + CMSBins.maxSpecialBins
        
        binSize = ( stat_end[stat] - stat_start[stat] ) / stat_nbins[stat]
        binStarts = [ stat_start[stat] + binSize * i for i in range( stat_nbins[ stat ] + CMSBins.stat_numSpecialBins[ stat ] )  ]

        pp.gca().set_xticks( binStarts )
        pp.gca().set_xticklabels( [ '%.2f' % b for b in binStarts[: stat_nbins[stat] ] ] + list( DictGet( CMSBins.stat_specialBinNames,
                                                                                                          stat, () ) ),
                                  rotation = 'vertical' )
        
        dbg( 'stat binStarts' )
        r0 = likesDataFiles[0][ stat + '_ratioLogCausal2Neutral' ]
        r1 = likesDataFiles[1][ stat + '_ratioLogCausal2Neutral' ][:len(binStarts)]
        hitsLine, = pp.plot( binStarts , ( r0 - r1 ) if cmpOp == 'diff' else ( r0 / r1 ), 'r-' )

        pp.hlines( 0, min( binStarts ), max( binStarts ), linestyles = 'dotted' )

    pp.savefig( outFN )
        
def DefineRulesTo_showLikesTables( pr, Ddata, complikeSfx, normalizeWithinReplicas, secondLikesSfx = '' ):
    """Define rules to show likes tables in various ways"""
    for scenario in GetScenarios():
        for putativeMutPop in ( AllPops if scenario.isNeutral() else ( scenario.mutPop, ) ):
            pr.addInvokeRule( invokeFn = writeBinCounts, invokeArgs = Dict( 'Ddata scenario complikeSfx normalizeWithinReplicas '
                                                                            'putativeMutPop' ) )

    complikeSfx += Sfx( secondLikesSfx )

    for putativeMutPop in AllPops:
        pr.addInvokeRule( invokeFn = writeLikesCounts, invokeArgs = Dict( 'Ddata complikeSfx putativeMutPop normalizeWithinReplicas' ) )

        for stat in CMSBins.CMSstats:
            pr.addInvokeRule( invokeFn = writeLikesCounts, invokeArgs = Dict( 'Ddata complikeSfx putativeMutPop normalizeWithinReplicas',
                                                                              stats = ( stat, ) ) )

def DefineRulesTo_compareLikesTables( pr, Ddata, complikeSfxs, normalizeWithinReplicas ):
    """Define rules to compare likes tables"""
    for putativeMutPop in AllPops:
        for cmpOp in ( 'diff', 'div' ):
            pr.addInvokeRule( invokeFn = compareLikesTables,
                              invokeArgs = Dict( 'putativeMutPop normalizeWithinReplicas cmpOp',
                                                 likesTableDefs = [ Dict( 'Ddata putativeMutPop normalizeWithinReplicas',
                                                                          complikeSfx = complikeSfx )
                                                                    for complikeSfx in complikeSfxs ] ) )
    
def DefineRulesTo_evalLikesDefs( pr, Ddata, mutAges, nreplicas, simsOut, thinExt,
                                 likesDefs,
                                 thinSfx = '',
                                 nonNanStats = 'ALL',
                                 includeRulesTo_runSweep = False,
                                 includeRulesTo_computeCMS = True ):
    """Define rules to evaluate and compare several different ways of computing likes scores for the same
    set of simulations.

    Params:

       likesDefs - one or more definitions of how to compute CMS scores; each likesDef is a combination of
          a likes table and possibly subtables to be used for replicas meeting specified conditions.
          It is an IDotData with the following columns:
             - likesDefName: name of the likesDef (a string)
             - descr: a description of the likesDef
             - likesTable: path to the likesTable
             - replicaTables: tables used in the replicaConds below
             - replicaConds: sequence of conditions saying which subtables should be used for which replicas.
    """

    pr.addLsfNeedDir( Ddata )
    
    if includeRulesTo_runSweep:
        DefineRulesTo_RunSimsAndSweep( inputParamsFiles = ( Ddata + '/neutralParams.txt', ),
                                       tests = ( 'ihs', 'xpop' ), 
                                       **Dict( 'pr Ddata mutAges nreplicas simsOut thinExt' ) )


    if includeRulesTo_computeCMS:

        args = Dict( 'pr Ddata simsOut thinExt nreplicas mutAges' )
        DefineRulesTo_computeCMSstats( **args )

    for r in likesDefs:
        DefineRulesTo_computeCMSforSims( includeRulesTo_computeCMSstats = False,
                                         likesTable = r.likesTable,
                                         complikeSfx = r.likesDefName,
                                         replicaTables = r.replicaTables,
                                         replicaConds = r.replicaConds, **args )
        
        DefineRulesTo_evalCMS( excludeNeutral = True,
                               complikeSfx = r.likesDefName,
                               **args )


    for replicaCondsSfx in ( 'all', 'lo', 'hi' ):

        for table, stat, coarsenBy, ticksCoarsen in ( ( 'causalStats', 'complikeExp_normedLocal', 2, 2 ),
                                                      ( 'causalStats', 'complikeExp', 2, 2 ),
                                                      ( 'causalRank', 'causalRank', 2, 2 ) ):

            pr.addInvokeRule( invokeFn = GraphHistograms,
                              invokeArgs = dict( histFiles =
                                                 [ os.path.join( Ddata, 'hist', '_'.join((table, stat,
                                                                                          likesDef.likesDefName,
                                                                                          replicaCondsSfx )) + '.tsv' )
                                                   for likesDef in likesDefs ],
                                                 labels = [ likesDef.likesDefName for likesDef in likesDefs ],
                                                 xlabel = table + '.' + stat,
                                                 ylabel = 'fraction of replicas',
                                                 title = table + '.' + stat + ' (' + replicaCondsSfx + ')',
                                                 normed = True,
                                                 coarsenBy = coarsenBy,
                                                 ticksCoarsen = ticksCoarsen,
                                                 outFile = os.path.join( Ddata, AddFileSfx( 'hist.svg',
                                                                                            table, stat,
                                                                                            replicaCondsSfx ) )),
                              name = 'GraphReplicaHists',
                              mediumRuleNameSfx = ( table, stat, replicaCondsSfx ),
                              attrs = dict( replicaStat = table + '.' + stat,
                                            sfx = replicaCondsSfx ) )


    for likesDef1, likesDef2 in itertools.combinations( likesDefs, r = 2 ):
        sfx1 = Sfx( likesDef1.likesDefName, 'nonNan', nonNanStats )
        sfx2 = Sfx( likesDef2.likesDefName, 'nonNan', nonNanStats )

        dfltArgs = dict( binSize = .1 )

        for table, stat, extraArgs in \
                ( ( 'causalStats', 'complikeExp_normedLocal', {} ),
                  ( 'causalStats', 'complike_normedLocal', {} ),
                  ( 'causalStats', 'complikeExp', {} ),
                  ( 'causalStats', 'complike', {} ),
                  ( 'causalStats', 'complike_mean', {} ),
                  ( 'causalStats', 'complike_std', {} ),
                  ( 'causalStats', 'complike_count', {} ),
                  ( 'causalStats', 'complikeRatio', {}  ),
                  ( 'causalStats', 'complikeRatio_mean', {}  ),
                  ( 'causalStats', 'complikeRatio_std', {}  ),
                  ( 'causalStats', 'iHSlike', {} ),
                  ( 'causalStats', 'max_xpoplike', {} ),
                  ( 'causalStats', 'StdDifflike', {} ),
                  ( 'causalStats', 'freqDifflike', {} ),
                  ( 'causalStats', 'meanFstlike', {} ),
                  ( 'causalRank', 'causalRank', dict( binSize = 1, coarsenBy = 12, ticksCoarsen = 40 )  ) ):

            DefineRulesTo_histogramReplicaStatistic( pr = pr, Ddata = Ddata,
                                                     outFile = '_'.join(('cmp', likesDef1.likesDefName,
                                                                         likesDef2.likesDefName,
                                                                         table, stat)) + '.svg',
                                                     replicaTables = ( table + sfx1, table + sfx2, 'replicaStats' ),
                                                     replicaStat =
                                                     table + sfx1 + '.' + stat + '-' + table + sfx2 + '.' + stat,

                                                     replicaConds = ( 'True', 'replicaStats.causalAlleleFreq < .5',
                                                                      'replicaStats.causalAlleleFreq >= .5' ),
                                                     replicaCondsSfxs = ( 'all', 'lo', 'hi' ),
                                                     scenCond = 'not is_neutral', **MergeDicts( dfltArgs, extraArgs ) )
                                                     
                                                     
            pr.addInvokeRule( invokeFn = scatterPlotReplicaStatistic,
                              invokeArgs = Dict( 'Ddata thinSfx nreplicas',
                                                 allScens = GetSelectionScenarios( mutAges = mutAges ),
                                                 replicaTables = ( table + sfx1, table + sfx2, 'replicaStats', 'selfreqStats' ),
                                                 #replicaCond = 'replicaStats.causalAlleleFreq >= 0.2',
                                                 replicaStatX = table + sfx1 + '.' + stat,
                                                 replicaStatY = table + sfx2 + '.' + stat,
                                                 replicaColorings =
                                                 ( ( 'hi', 'replicaStats.causalAlleleFreq > .5', 'r' ),
                                                   ( 'lo', '.2 <= replicaStats.causalAlleleFreq <= .5', 'b' ),
                                                   ( 'verylow', 'replicaStats.causalAlleleFreq < .2', 'g' ) ),
                                                 replicaShow =
                                                 ( '"causalFreq"',
                                                   '"%.2f" % replicaStats.causalAlleleFreq',
                                                   '"max_xpop"', '"%.2f" % selfreqStats.max_xpop',
                                                   '"predictionMargin"',
                                                   '"%.2f" % ( ( replicaStats.causalAlleleFreq - .5 ) ' \
                                                       '* ( 1 if selfreqStats.max_xpop > 5 else -1 ) )' ),
                                                 title = '%s.%s for %s vs %s' % ( table, stat, likesDef1.likesDefName,
                                                                                  likesDef2.likesDefName ) +
                                                 '\nred = hi, blue = lo, green = under .2',
                                                 outFile =
                                                 os.path.join( Ddata,
                                                               AddFileSfx( 'scatterPlot.svg',
                                                                           table, stat, sfx1, sfx2 ) ) ),
                              attrs =
                              Dict( 'stat', likesDef1 = likesDef1.likesDefName, likesDef2 = likesDef2.likesDefName,
                                    filterUnder20 = False ) )

            pr.addInvokeRule( invokeFn = scatterPlotReplicaStatistic,
                              invokeArgs = Dict( 'Ddata thinSfx nreplicas',
                                                 allScens = GetSelectionScenarios( mutAges = mutAges ),
                                                 replicaTables = ( table + sfx1, table + sfx2, 'replicaStats' ),
                                                 replicaCond = 'replicaStats.causalAlleleFreq >= 0.2',
                                                 replicaStatX = table + sfx1 + '.' + stat,
                                                 replicaStatY = table + sfx2 + '.' + stat,
                                                 replicaColorings =
                                                 ( ( 'hi', 'replicaStats.causalAlleleFreq > .5', 'r' ),
                                                   ( 'lo', 'replicaStats.causalAlleleFreq <= .5', 'b' ), ),
                                                 replicaShow = '"%.2f" % replicaStats.causalAlleleFreq',
                                                 title = '%s.%s for %s vs %s' % ( table, stat, likesDef1.likesDefName,
                                                                                  likesDef2.likesDefName ) +
                                                 '\nred = hi, blue = lo; under .2, filtered',
                                                 outFile =
                                                 os.path.join( Ddata,
                                                               AddFileSfx( 'scatterPlot.svg',
                                                                           table, stat, sfx1, sfx2, 'filt' ) ) ),
                              attrs =
                              Dict( 'stat', likesDef1 = likesDef1.likesDefName, likesDef2 = likesDef2.likesDefName,
                                    filterUnder20 = True ) )
            

                              
    if False:
        DefineRulesTo_histogramReplicaStatistic( pr = pr, Ddata = Ddata,
                                                 outFile = 'causalComplikeExpNormedLocalHist.svg',
                                                 sfx = complikeSfx,
                                                 replicaTables = ( causalStatsTable, 'replicaStats' ),
                                                 replicaStat = causalStatsTable + '.complikeExp_normedLocal',
                                                 scenCond = 'not is_neutral',
                                                 binSize = .1, ticksCoarsen = 2, cumulative = True,
                                                 normed = True, coarsenBy = 32, cumulativeUpTo = cumulativeUpTo )


def DefineRulesTo_copyLikesTables( pr, DdataFrom, DdataTo, likesSfx, replicaCondsSfxs = ( '', ),
                                   mutPops = AllPops ):
    """Define rules to copy likes tables, including all associated files, from one location to another."""

    for replicaCondsSfx in replicaCondsSfxs:
        for mutPop in mutPops:
            likesTableFN = AddFileSfx( 'likesTable.tsv', likesSfx, mutPop, replicaCondsSfx )
            pr.addCpRule( comment = 'copy likes table files',
                          source = os.path.join( DdataFrom, likesTableFN ), 
                          target = AddSlash( DdataTo ) )
            likesTable = LoadLikesTable( Ddata = DdataFrom, likesTable = 'likesTable.tsv',
                                         likesTableSfx = likesSfx,
                                         mutPop = mutPop )
            for which in 'hits', 'miss', 'region':
                pr.addCpRule( comment = 'copy likes table component',
                              source = likesTable[ which + 'Likes' ],
                              target = os.path.join( DdataTo, 'likes/' ) )


