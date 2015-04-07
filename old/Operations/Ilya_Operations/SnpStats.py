
"""Flexible code for histogramming per-snp and per-replica statistics for selected SNPs in selected replicas in
selected scenarios and/or demographies."""

from Operations.Shari_Operations.localize.Scenario import GetScenarios, GetSelectionScenarios
from Operations.MiscUtil import Dict, compile_expr, dbg, Histogrammer, AddFileSfx, ReplaceFileExt, \
    MergeDicts, MakeSeq, SlurpFile, IsSeq, Sfx, MakeAlphaNum, DictGet, tmap, PrintDictDiff
from Classes.DotData import DotData
from Operations.Ilya_Operations.PipeRun.python.PipeRun import GetDependsOn
from Operations.Shari_Operations.localize.PopConsts import AllFreqs, AllPops, AllAges, CAUSAL_POS
from Operations.IDotData import IDotData
import operator, os, logging, contextlib, functools, collections, types, ast
from itertools import izip
import itertools, string
from UserDict import DictMixin
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as pp
import numpy as np
import math
import traceback as tb


__all__ = ( 'gatherCausalFreqs', 'DefineRulesTo_gatherCausalFreqs', 'histogramSnpStatistic', 'histogramReplicaStatistic',
            'AddUpHistograms', 'GraphHistograms', 'GraphCumulPlots', 'DefineRulesTo_histogramSnpStatistic',
            'DefineRulesTo_histogramReplicaStatistic', 'findReplicasMatchingConds', 'findSnpsMatchingConds',
            'identifyReplicasMeetingConds', 'splitSnpStatsFile',
            'DefineRulesTo_identifyReplicasMeetingCommonConds' )

def gatherCausalFreqs( scen, Ddata, simsOut, thinSfx, thinExt, nreplicas, getio = None ):
    """For all replicas within one scenario, gather some useful summary info for each replica:
    e.g. that replica's modern-day frequency of the causal allele, the genetic map position of the
    causal SNP, number of SNPs in the replica, the range of the genetic map, etc.
    """

    #hm3big/simsOutHm3big/10ky/sel100_1/
    simScenDir = Ddata + '/' + simsOut + thinSfx + '/' + scen.scenDir()
    statScenDir = Ddata + '/replicastats' + thinSfx + '/' + scen.scenDir()

    posFileNames = [ simScenDir + '/' + '%d_%s.pos-%d%s' % ( replicaNum, scen.scenName(), scen.mutPop, thinExt )
                     for replicaNum in range( nreplicas ) ] if not scen.is_neutral() else []

    replicaInfoFileName = statScenDir + '/replicaStats.tsv' 

    if getio: return dict( depends_on = posFileNames, creates = replicaInfoFileName,
                           mediumRuleNameSfx = scen )

    causalAlleleFreqs = [ ]
    replicaNums = [ ]

    selpos = 500000

    okReplicas = 0

    for replicaNum in range( nreplicas ):
        if scen.is_neutral(): causalFreq = np.nan
        else:
            posFile = DotData( SVPath = posFileNames[ replicaNum ], SVSkipFirstLines = 1, SVHeader = False,
                               names = ['SNP','CHROM', 'CHROM_POS', 'ALLELE1', 'FREQ1', 'ALLELE2', 'FREQ2' ] )
            causalLine = posFile[ posFile.CHROM_POS == selpos ]
            assert len( causalLine ) == 1
            causalFreq = causalLine[0].FREQ1
        causalAlleleFreqs.append( causalFreq )
        replicaNums.append( replicaNum )

    DotData( names = [ 'replicaNum', 'causalAlleleFreq', 'targetCausalFreq' ],
             Columns = [ replicaNums, causalAlleleFreqs,
                         (( 0 if scen.isNeutral() else scen.mutFreq),)*nreplicas ] ).saveToSV( replicaInfoFileName )

def gatherReplicaGDstats( scen, Ddata, simsOut, thinSfx, thinExt, nreplicas, getio = None ):
    """For all replicas within each scenario, gather some genetic map-related info for each replica:
    e.g. the genetic map position of the
    causal SNP, the range of the genetic map, etc.
    """

    #hm3big/simsOutHm3big/10ky/sel100_1/
    simScenDir = os.path.join( Ddata, simsOut + thinSfx, scen.scenDir() )
    statScenDir = os.path.join( Ddata,  'replicastats' + thinSfx, scen.scenDir() )

    posFileNames = [ simScenDir + '/' + '%d_%s.pos-%d%s' % ( replicaNum, scen.scenName(), scen.mutPop, thinExt )
                     for replicaNum in range( nreplicas ) ] if not scen.is_neutral() else []

    replicaInfoFileName = statScenDir + '/replicaStats.tsv' 

    if getio: return dict( depends_on = posFileNames, creates = replicaInfoFileName,
                           mediumRuleNameSfx = scen )

    causalAlleleFreqs = [ ]
    replicaNums = [ ]

    selpos = 500000



    
    
def DefineRulesTo_gatherCausalFreqs( pr, Ddata, simsOut = 'simsOut',
                                     mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs,
                                     thinSfx = '', thinExt = '', nreplicas = 100 ):
    """Define rules to gather per-replica statistics"""

    for scen in GetScenarios( mutAges, mutPops, mutFreqs ):

        pr.addInvokeRule( invokeFn = gatherReplicaStats,
                          invokeArgs = Dict( 'scen Ddata simsOut thinSfx thinExt nreplicas' ) )


# for compatibility with old code
gatherReplicaStats = gatherCausalFreqs
DefineRulesTo_gatherReplicaStats = DefineRulesTo_gatherCausalFreqs
        
        
def histogramSnpStatistic( Ddata, thinSfx, scenDir, replicaTables, replicaCond, snpTables, snpCond, snpStat,
                           outFile, nreplicas, binSize, binShift = 0.0, sfx = None, scenSfx = None, getio = None ):
    """Compute histogram of $snpStat for snps matching $snpCond in replicas matching $replicaCond in scenario $scenDir.

    Params:

       statTable - the name of the per-snp statistics table.  we assume there is a file called
          Ddata/snpstats/scenDir/statTable_pop.tsv for each scenario.

       statCol - column name to histogram.
       
    """

    replicaTables = MakeSeq( replicaTables )
    snpTables = MakeSeq( snpTables )
    
    replicaCondExpr = compile_expr( replicaCond )
    snpCondExpr = compile_expr( snpCond )
    snpStatExpr = compile_expr( snpStat )

    outFile = AddFileSfx( outFile, sfx )
    outFileStats = AddFileSfx( outFile, 'stats' )
    
    if IsSeq( scenSfx ): scenSfx = dict( scenSfx )
    
    replicaTableFiles = [ os.path.join( Ddata, 'replicastats' + thinSfx, scenDir,
                                        replicaTable + ( '.tsv' if '.' not in replicaTable else '' ) )
                          for replicaTable in replicaTables ]

    snpTableFiles = [ os.path.join( Ddata,  'snpStats' + thinSfx, scenDir,
                                    AddFileSfx( snpTable + ( '.tsv' if '.' not in snpTable else '' ),
                                                scenSfx if isinstance( scenSfx, types.StringTypes )
                                                else scenSfx[ os.path.splitext( snpTable )[0] ] ) )
                      for snpTable in snpTables ]

    #dbg('replicaTableFiles snpTableFiles')
    #dbg('"*****" replicaTableFiles+snpTableFiles')
    replicaTableFiles = [ f + '/' if f.endswith('.data') else f for f in replicaTableFiles ]
    snpTableFiles = [ f + '/' if f.endswith('.data') else f for f in snpTableFiles ]

    snpTables = [ os.path.splitext(snpTable)[0] for snpTable in snpTables ]

    if getio: return dict( depends_on = replicaTableFiles + snpTableFiles,
                           creates = ( outFile, AddFileSfx( outFile, 'stats' ) ),
                           attrs = Dict( 'scenDir snpCond replicaCond snpStat' ),
                           mediumRuleNameSfx = ( scenDir, scenSfx ) )
    
    replicaTableVals = [ DotData( SVPath = f ) for f in replicaTableFiles ]

    replicasToUse = [ eval( replicaCondExpr, globals(), dict( zip( replicaTables, replicaTableRows ) ) )
                      for replicaTableRows in izip( *replicaTableVals ) ]

    #dbg( 'sum(replicasToUse)' )

    snpTableVals = [ IDotData( SVPath = f ) for f in snpTableFiles ]

    histogramBuilder = Histogrammer( binSize = binSize, binShift = binShift )

    lastReplica = np.nan
    for snpTableRows in izip( *snpTableVals ):
        r0 = snpTableRows[ 0 ]
        assert all([ r.Chrom == r0.Chrom for r in snpTableRows ]) or all([ np.isnan( r.Chrom ) for r in snpTableRows ])
        assert all([ r.Pos == r0.Pos for r in snpTableRows ])
        
        replica = int( r0.Chrom ) if not np.isnan( r0.Chrom ) else -1
        useThisReplica = not replicaTables or replicasToUse[ replica ]
        if replica != lastReplica: dbg( 'replica useThisReplica histogramBuilder.getNumVals()' )
        if useThisReplica:
            snpDict = dict( zip( snpTables, snpTableRows ) )
            if eval( snpCondExpr, globals(), snpDict ):
                val = eval( snpStatExpr, globals(), snpDict )
                histogramBuilder.addVal( val )
        lastReplica = replica

    logging.info('saving histogram to ', outFile )
        
    histogramBuilder.save( outFile )


def histogramReplicaStatistic( Ddata, thinSfx, replicaCond, replicaStat,
                               outFile, nreplicas, binSize, scenCond = 'True',
                               replicaTables = None,
                               scen2sfxs = {}, allScens = GetScenarios(),
                               sfx = None,  replicaCondSfx = '',
                               nameSfx = '', getio = None ):
    """Compute histogram of $replicaStat for replicas matching $replicaCond in scenarios matching $scenCond.

    Saves the histogram as well as overall stats about the values of this statistic, e.g. the average.

    Params:

       statTable - the name of the per-snp statistics table.  we assume there is a file called
          Ddata/snpstats/scenDir/statTable_pop.tsv for each scenario.

       statCol - column name to histogram.

    """

    outFile = AddFileSfx( outFile, sfx, replicaCondSfx )
    outFileStats = AddFileSfx( outFile, 'stats' )

    args = Dict( 'Ddata thinSfx replicaTables scenCond replicaCond scen2sfxs allScens' )
    if getio: return dict( depends_on =
                           findReplicasMatchingConds( getio = True, **args )[ 'depends_on' ],
                           creates = ( outFile, outFileStats ),
                           mediumRuleNameSfx = sfx, attrs = dict( piperun_short = True ),
                           name = 'histogramReplicaStatistic' + Sfx( nameSfx ) )
    
    histogramBuilder = Histogrammer( binSize = binSize )
    histogramBuilder.addVals( findReplicasMatchingConds( showHeadings = 'val', showVals = replicaStat, **args ).val )
    histogramBuilder.save( outFile )


def histogramSnpStatistic2( Ddata, thinSfx, snpTables, snpCond, snpCondSfx, replicaTables, replicaCond, replicaStat,
                            outFile, nreplicas, binSize, scenCond = 'True',
                            scen2sfxs = {}, allScens = GetScenarios(),
                            sfx = None,  replicaCondSfx = '',
                            nameSfx = '', getio = None ):
    """Compute histogram of $replicaStat for replicas matching $replicaCond in scenarios matching $scenCond.

    Saves the histogram as well as overall stats about the values of this statistic, e.g. the average.

    Params:

       statTable - the name of the per-snp statistics table.  we assume there is a file called
          Ddata/snpstats/scenDir/statTable_pop.tsv for each scenario.

       statCol - column name to histogram.

    """

    outFile = AddFileSfx( outFile, sfx, replicaCondSfx, snpCondSfx )
    outFileStats = AddFileSfx( outFile, 'stats' )

    args = Dict( 'Ddata thinSfx snpTables snpCond replicaTables scenCond replicaCond scen2sfxs allScens' )
    if getio: return dict( depends_on =
                           finSnpsMatchingConds( getio = True, **args )[ 'depends_on' ],
                           creates = ( outFile, outFileStats ),
                           mediumRuleNameSfx = sfx, attrs = dict( piperun_short = True ),
                           name = 'histogramReplicaStatistic' + Sfx( nameSfx ) )
    
    histogramBuilder = Histogrammer( binSize = binSize )
    histogramBuilder.addVals( findSnpsMatchingConds( showHeadings = 'val', showVals = snpStat, **args ).val )
    histogramBuilder.save( outFile )


    
def AddUpHistograms( histFiles, outFile, getio = None ):
    """Add up histograms from separate files, write results to new file"""

    outFileStats = AddFileSfx( outFile, 'stats' )
    if getio: return dict( depends_on = histFiles, creates = ( outFile, outFileStats ),
                           attrs = dict( piperun_short = True ) )

    sumHist = reduce( operator.add, map( Histogrammer.load, histFiles ) )
    sumHist.save( outFile )
    
def GraphHistograms( histFiles, outFile = None, xlabel = '', ylabel = '', title = '',
                     labels = (), colors = 'brcmygkbrcmygkbrcmygkbrcmygk',
                     relWidth = 0.4,
                     xbound = None, ybound = None, coarsenBy = None, sfx = '',
                     ticksCoarsen = 1, log = False, normed = False,
                     cumulative = False,
                     cumulativeUpTo = None,
                     figSize = (24, 12 ),
                     subplots_adjust = {},
                     getio = None ):
    """Plot one or more histograms sharing the same bins.

    Params:

       normalizeHistograms - if true, for each histogram on the y-axis we plot not the number of
          items in a given bin, but their fraction out of the total number of items in that histogram.
          This lets us compare different histograms.

    """

    #dbg( '"at_first" labels' )

    # ( if the bins of one are strictly finer than bins of other, i.e. if they form a DAG in this
    # relationship, then we can still do the graph).

    histFiles = MakeSeq( histFiles )

    if not outFile:
        assert len( histFiles ) == 1
        outFile = ReplaceFileExt( histFiles[0], '.png' )

    outFile = AddFileSfx( outFile, sfx )

    if not labels: labels = [ os.path.splitext( os.path.basename( f ) )[0] for f in histFiles ]

    if getio: return dict( depends_on = histFiles, creates = outFile,
                           mediumRuleNameSfx = sfx,
                           attrs = dict( piperun_short = True ) )


    pp.figure(1, figsize = figSize )
    #pp.clf()
    
    pp.subplots_adjust( **MergeDicts( dict( hspace = 0.3, bottom = 0.15 ), subplots_adjust ) )

    for which, cumulative in enumerate( ( True, False ) ):

        pp.subplot( 2, 1, which + 1 )

        pp.xlabel( xlabel )
        pp.ylabel( ylabel )

        pp.hold( True )

        binSize = None
        binShift = None
        theLabels = []
        theHandles = []

        hists = map( Histogrammer.load, histFiles )
        if coarsenBy: hists = [ hist.coarsenBy( coarsenBy ) for hist in hists ]
        allBinIds = reduce( operator.concat, [ hist.bin2count.keys() for hist in hists ] )
        if not allBinIds: allBinIds = ( 0, )
        minBinId = min( allBinIds )
        maxBinId = max( allBinIds ) + 1

        if cumulativeUpTo is not None:
            maxBinId = min( maxBinId, max( [ hist.getCumulativeBinFor( cumulativeUpTo ) for hist in hists ] ) ) + 1

        for color, label, ( histFileNum, hist ) in zip( colors, labels, enumerate( hists ) ):

            # check that all histograms we're loading have the same bins
            if binSize is None: binSize = hist.binSize
            else: assert abs( hist.binSize - binSize ) < 1e-12
            if binShift is None: binShift = hist.binShift
            else: assert abs( hist.binShift - binShift ) < 1e-12

            width = binSize * relWidth / len( histFiles )

            left = np.array( hist.getAllBinLefts( minBinId = minBinId, maxBinId = maxBinId ) ) + histFileNum * width

            if histFileNum == 0: pp.xticks( [ x for i, x in enumerate( left ) if i % ticksCoarsen == 0 ] )

            height = hist.getAllBinCounts( normed = normed, cumulative = cumulative,
                                           minBinId = minBinId, maxBinId = maxBinId )

            rects = pp.bar( height = height,
                            width = width * 0.95, **Dict( 'left color log' ) )
            if rects:
                labelHere = label + ' (%d values)' % hist.getNumVals()
                if hist.getNumNaNs(): labelHere += ' (%d nans)' % hist.getNumNaNs()
                if hist.getNumInfs(): labelHere += ' (%d infs)' % hist.getNumInfs()
                rects[ 0 ].set_label( labelHere )
                theLabels.append( labelHere )
                theHandles.append( rects[0] )

        pp.title( title )
        if theLabels and theHandles:
            pp.figlegend( loc = 'lower center', labels = theLabels, handles = theHandles )
        if xbound: pp.gca().set_xbound( *xbound )
        if ybound: pp.gca().set_ybound( *ybound )

    pp.savefig( outFile )


def GraphCumulPlots( histFiles, outFile = None, xlabel = '', ylabel = '', title = '',
                     labels = (), colors = 'brcmygkbrcmygkbrcmygkbrcmygk',
                     relWidth = 0.4,
                     xbound = None, ybound = None, coarsenBy = None, sfx = '',
                     ticksCoarsen = 1, log = False, normed = True,
                     getio = None ):
    """Plot one or more cumulative plots.
    """

    # ( if the bins of one are strictly finer than bins of other, i.e. if they form a DAG in this
    # relationship, then we can still do the graph).

    histFiles = MakeSeq( histFiles )

    if not outFile:
        assert len( histFiles ) == 1
        outFile = ReplaceFileExt( histFiles[0], '.png' )

    if not labels: labels = [ os.path.splitext( os.path.basename( f ) )[0] for f in histFiles ]

    outFileTable = outFile + '.points.tsv'
    
    if getio: return dict( depends_on = histFiles, creates = ( outFile, outFileTable ),
                           mediumRuleNameSfx = sfx,
                           attrs = dict( piperun_short = True ) )


    pp.figure(1, figsize = (18,6) )
    #pp.clf()

    pp.subplots_adjust( bottom = 0.37 )

    pp.xlabel( xlabel + '\n\n\n\n' )
    pp.ylabel( ylabel )

    pp.hold( True )

    binSize = None
    theLabels = []
    theHandles = []

    for color, label, ( histFileNum, histFile ) in zip( colors, labels, enumerate( histFiles ) ):

        hist = Histogrammer.load( histFile )

        if coarsenBy: hist = hist.coarsenBy( coarsenBy )

        if not binSize: binSize = hist.binSize
        else:
            if not abs( hist.binSize - binSize ) < 1e-12:
                dbg( 'hist.binSize binSize hist.binSize-binSize' )
            assert abs( hist.binSize - binSize ) < 1e-12

        binLefts = hist.getBinLefts()
        
        if histFileNum == 0: pp.xticks( [ x for i, x in enumerate( binLefts ) if i % ticksCoarsen == 0 ] )

        binCounts = hist.getBinCounts( normed = normed, cumulative = True )
        rects = pp.plot( binLefts, binCounts, label = label, color = color )

        DotData( names = ( 'binLefts', 'binCounts' ), Columns = ( binLefts, binCounts ) ).saveToSV( outFileTable )

        if rects:
            theLabels.append( label )
            theHandles.append( rects )
        
    pp.title( title )

    if theLabels and theHandles:
        pp.figlegend( loc = 'lower center', labels = theLabels, handles = theHandles )
    if xbound: pp.gca().set_xbound( *xbound )
    if ybound: pp.gca().set_ybound( *ybound )

    pp.savefig( outFile )

def DefineRulesTo_histogramSnpStatistic( pr, Ddata,
                                         outFile, snpTables, snpStat, binSize,
                                         binShift = 0.0,
                                         scen2sfxs = lambda scen: '',
                                         scenCond = 'True',
                                         allScens = GetScenarios(),
                                         nreplicas = 100, thinSfx = '', replicaTables = (),
                                         replicaConds = 'True', replicaCondsSfxs = '',
                                         snpConds = 'True', snpCondsSfxs = '', title = '', titlePrefix = '',
                                         xlabel = '', ylabel = '',
                                         xbound = None, ybound = None, log = False, coarsenBy = None, sfx = '',
                                         ticksCoarsen = 1, cumulative = False, normed = False,
                                         colors = 'brcmygkbrcmygkbrcmygkbrcmygk',
                                         subplots_adjust = {},
                                         name = None ):
    """A generic way to plot the distribution of some per-snp statistics for some subset of SNPs.

    Params:

       statTable - the name of the per-snp statistics table.  we assume there is a file called
          Ddata/snpstats/scenDir/statTable_pop.tsv for each scenario.

       statCol - column name to histogram.

    Notes:

       - for histogramming should not need to load it all into memory.  can do a pre-pass to just get
       the range of values, define the bins, then do a second pass to count what goes in what bin.

       could also add bins as we go.  so, really just need to know bin size, and then can do all this
       with one pass.  can also, later, make this automatically parallelized.
       
    """

    if not os.path.dirname( outFile ): outFile = os.path.join( Ddata, outFile )
    
    scenCondExpr = compile_expr( scenCond )

    replicaConds = MakeSeq( replicaConds )
    replicaCondsSfxs = MakeSeq( replicaCondsSfxs )
    
    snpConds = MakeSeq( snpConds )
    snpCondsSfxs = MakeSeq( snpCondsSfxs )

    totaledHistFiles = []
    totaledLabels = []

    outFile = AddFileSfx( outFile, sfx )
    
    baseOutFile = outFile

    for replicaCond, replicaCondSfx in zip( replicaConds, replicaCondsSfxs ):
        for snpCond, snpCondSfx in zip( snpConds, snpCondsSfxs ):
    
            histFiles = []

            for scen in allScens:

                if not eval( scenCondExpr, globals(), ScenAttrs( scen ) ): continue

                scenDir = scen.scenDir()

                for scenSfx in MakeSeq( scen2sfxs( scen ) if callable( scen2sfxs ) else scen2sfxs[ scen ] ):

                    histOutFile = os.path.join( Ddata, 'hist', scenDir,
                                                AddFileSfx( ReplaceFileExt( os.path.basename( outFile ), '.tsv' ),
                                                            snpStat, 
                                                            replicaCondSfx, snpCondSfx, scenSfx, sfx ) )

                    rule = pr.addInvokeRule( invokeFn = histogramSnpStatistic,
                                             invokeArgs =
                                             dict( outFile = histOutFile,
                                                   **Dict( 'Ddata thinSfx replicaTables replicaCond snpTables snpCond '
                                                           'snpStat nreplicas binSize binShift scenDir scenSfx sfx' ) ),
                                             name = name,
                                             comment = 'Compute distribution of ' + snpStat
                                             + ' for SNPs matching ' + snpCond + ' in replicas matching ' + replicaCond )
                    histFiles.append( histOutFile )

            totaledHistFile = os.path.join( Ddata, 'hist',
                                            AddFileSfx( ReplaceFileExt( os.path.basename( outFile ), '.tsv' ),
                                                        snpCondSfx, replicaCondSfx, sfx ) )
            totaledHistFiles.append( totaledHistFile )

            totaledLabel = ''
            if replicaCondSfx:
                totaledLabel += replicaCondSfx + ' replicas' + ( (' (' + replicaCond + ') ') \
                                                                     if replicaCond != 'True' else '' )
            if snpCondSfx: totaledLabel += snpCondSfx + ' SNPs' + ( (' (' + snpCond + ') ') \
                                                                        if snpCond != 'True' else '' )
            
            totaledLabels.append( totaledLabel )
            
            pr.addInvokeRule( invokeFn = AddUpHistograms, invokeArgs = dict( histFiles = histFiles,
                                                                             outFile = totaledHistFile ),
                              mediumRuleNameSfx = ( sfx, snpStat, replicaCondSfx, snpCondSfx ), name = 'AddUpSnpHists',
                              fileDescrs = { 0:
                                                 ( 'Distribution of <b>' + snpStat + '</b> among '
                                                   + ( 'all SNPs' if snpCond == 'True'
                                                       else ' snps matching <em>' + snpCond + '</em>' )
                                                   + ' in '
                                                   + ( 'all replicas' if replicaCond == 'True' else
                                                       'replicas matching <em>' + replicaCond + '</em>' )
                                                   + ' in '
                                                   + ( 'all scenarios' if scenCond == 'True' else
                                                       'scenarios matching <em>' + scenCond + '</em>' ),
                                                   ( ( 'count', 'Number of SNPs with ' + snpStat + ' in given bin' ),) ) } )


    if not title:
        title = 'Histogram of ' + snpStat + '\n'
        if scenCond != 'True': title += ' scenCond: ' + scenCond
        if any( replicaCond != 'True' for replicaCond in replicaConds ):
            title += ' replicaConds: ' + ', '.join(replicaCondsSfxs)
        if any( snpCond != 'True' for snpCond in snpConds ): title += ' snpConds: ' + ', '.join(snpCondsSfxs)
    title = titlePrefix + title

    if not ylabel: ylabel = ('#' if not normed else 'fraction') + ' of snps'
    if not xlabel: xlabel = snpStat
            
    pr.addInvokeRule( invokeFn = GraphHistograms,
                      mediumRuleNameSfx = (snpStat,) + tuple(replicaCondsSfxs) + tuple(snpCondsSfxs),
                      name = 'GraphSnpHists',
                      invokeArgs = dict( histFiles = totaledHistFiles, labels = totaledLabels,
                                         **Dict( 'xlabel ylabel title xbound ybound coarsenBy log outFile '
                                                 'cumulative normed ticksCoarsen colors' ) ),
                      attrs = Dict( 'snpStat replicaConds snpConds scenCond subplots_adjust' ) )
    
    
def DefineRulesTo_histogramReplicaStatistic( pr, Ddata,
                                             outFile, replicaStat, binSize,
                                             scenCond = 'True',
                                             replicaTables = None,
                                             sfx = '',
                                             scen2sfxs = lambda scen: '',
                                             allScens = tuple( GetScenarios() ),
                                             nreplicas = 100, thinSfx = '', 
                                             replicaConds = 'True', replicaCondsSfxs = '',
                                             title = '', titlePrefix = '',
                                             xlabel = '', ylabel = '',
                                             xbound = None, ybound = None, log = False, coarsenBy = None,
                                             ticksCoarsen = 1, cumulative = False, normed = False,
                                             cumulativeUpTo = 0.99,
                                             subplots_adjust = {},
                                             name = None, nameSfx = '' ):
    """Define rules to plot the distribution of a specified per-replica statistic for some subsets of replicas
    in some subset of scenarios.

    Params:

       pr - the PipeRun object to which the rules should be added
       Ddata - the root folder of the genetic data in simulations format
       outFile - the filename to which the histogram plot will be written
       replicaTables - names of tables containing per-replica values.  For each such table T,
         there must be a file of the form os.path.join( Ddata, replicastats, scenario.scenDir(), T + '.tsv' )
         giving some values for each replica in the scenario.
       replicaStat - a Python expression in which the names in replicaTables may appear as variables, and refer
         to a named tuple representing the replica's row in the corresponding replicaTable.

    Notes:

       - for histogramming should not need to load it all into memory.  can do a pre-pass to just get
       the range of values, define the bins, then do a second pass to count what goes in what bin.

       could also add bins as we go.  so, really just need to know bin size, and then can do all this
       with one pass.  can also, later, make this automatically parallelized.
       
    """

    if not os.path.dirname( outFile ): outFile = os.path.join( Ddata, outFile )
    
    scenCondExpr = compile_expr( scenCond )
    ourScens = [ scen for scen in allScens if eval( scenCondExpr, globals(), ScenAttrs( scen ) ) ]

    if callable( scen2sfxs ):
        scen2sfxs = dict( ( scen, scen2sfxs( scen ) ) for scen in ourScens )

    replicaConds = MakeSeq( replicaConds )
    replicaCondsSfxs = MakeSeq( replicaCondsSfxs )

    totaledHistFiles = []
    totaledLabels = []

    for replicaCond, replicaCondSfx in zip( replicaConds, replicaCondsSfxs ):
        
        totaledHistFile = os.path.join( Ddata, 'hist',
                                        ReplaceFileExt( os.path.basename( outFile ), '.tsv' ) )

        totaledLabels.append( replicaCondSfx + ': ' + replicaCond )
            
        r = pr.addInvokeRule( invokeFn = histogramReplicaStatistic,
                              invokeArgs = Dict( 'Ddata thinSfx replicaTables replicaCond replicaStat nreplicas '
                                                 'binSize scenCond scen2sfxs allScens nameSfx sfx replicaCondSfx',
                                                 outFile = totaledHistFile ),
                              mediumRuleNameSfx = ( replicaStat, replicaCondSfx, sfx ),
                              fileDescrs = { 0:
                                                 ( 'Distribution of <b>' + replicaStat + '</b> among '
                                                   + ( 'all replicas' if replicaCond == 'True' else
                                                       'replicas matching <em>' + replicaCond + '</em>' )
                                                   + ' in '
                                                   + ( 'all scenarios' if scenCond == 'True' else
                                                       'scenarios matching <em>' + scenCond + '</em>' ),
                                                   ( ( 'count', 'Number of replicas with ' + replicaStat +
                                                       ' in given bin' ),
                                                     )) } )
        totaledHistFiles.append( r.creates[0] )

    if not title:
        if scenCond != 'True': title += ' scenCond: ' + scenCond
        if len( replicaConds ) == 1 and replicaConds[0] != 'True': title += ' replicaCond: ' + replicaConds[0]
    title = titlePrefix + title

    if not ylabel: ylabel =  ('#' if not normed else 'fraction') + ' of replicas'
    if not xlabel: xlabel = replicaStat

    pr.addInvokeRule( invokeFn = GraphHistograms,
                      invokeArgs = dict( histFiles = totaledHistFiles, labels = totaledLabels,
                                         **Dict( 'xlabel ylabel title xbound ybound coarsenBy log outFile '
                                                 'sfx ticksCoarsen cumulative normed cumulativeUpTo' ) ),
                      name = 'GraphReplicaHists' + Sfx( nameSfx ),
                      mediumRuleNameSfx = ( replicaStat, sfx ) + tuple( replicaConds ),
                      attrs = Dict( 'replicaStat sfx subplots_adjust' ) )

    return totaledHistFiles
    

def identifyReplicasMeetingConds( Ddata, scenario, replicaTables, replicaConds, condsFileFN, nreplicas,
                                  thinSfx = '', getio = None ):
    """Given a list of named replica conditions, determine for each replica which conditions it meets, and
    write out the result in an easy-to-access format.

    Input params:

      replicaConds - sequence of pairs of the form ( condName, cond ) -- for example,
         ( ( 'hi', 'replicaStats.causalAlleleFreq >= .5' ), ( 'lo', 'replicaStats.causalAlleleFreq < .5' ) )
    """

    replicaTables = MakeSeq( replicaTables )

    replicaTableFiles = [ os.path.join( Ddata, 'replicastats' + thinSfx, scenario.scenDir(),
                                        replicaTable + ( '.tsv' if not os.path.splitext( replicaTable )[1] else '' ) )
                          for replicaTable in replicaTables ]

    if not os.path.dirname( condsFileFN ): condsFileFN = os.path.join( Ddata, 'replicastats' + thinSfx, scenario.scenDir(),
                                                                       condsFileFN )

    if getio: return dict( depends_on = replicaTableFiles, creates = condsFileFN, mediumRuleNameSfx = scenario.scenDir(),
                           attrs = dict( piperun_short = True,
                                         condNames = ', '.join( map( operator.itemgetter( 0 ), replicaConds ) ) ) )
    
    replicaTableVals = [ DotData( SVPath = f ) for f in replicaTableFiles ]

    assert all([ len( replicaTableVal ) == nreplicas for replicaTableVal in replicaTableVals ])

    matchingReplicas = []
    
    for replicaCond in map( operator.itemgetter( 1 ), replicaConds ):
        replicaCondExpr = compile_expr( replicaCond )
        replicasToUse = [ int( eval( replicaCondExpr, globals(), dict( zip( replicaTables, replicaTableRows ) ) ) )
                          for replicaTableRows in izip( *replicaTableVals ) ]
        matchingReplicas.append( replicasToUse )

    Records = []
    condNames = tuple( map( operator.itemgetter( 0 ), replicaConds ) )
    for replicaNum, condResults in enumerate( izip( *matchingReplicas ) ):
        Records.append( ( replicaNum, ','.join( replicaCondName for condNum, replicaCondName
                                                in enumerate( condNames )
                                                if condResults[ condNum ]  ) )
                        + condResults )
    IDotData( names = ( 'replicaNum', 'matchingConds' ) + condNames, Records = Records ).save( condsFileFN )

def DefineRulesTo_identifyReplicasMeetingCommonConds( pr, Ddata, thinSfx = '', allScens = GetSelectionScenarios(),
                                                      nreplicas = 100 ):
    """Define rules to identify replicas meeting common conditions such as all/lo/hi freq"""

    for scenario in allScens:
        pr.addInvokeRule( invokeFn = identifyReplicasMeetingConds,
                          invokeArgs = Dict( 'Ddata scenario nreplicas thinSfx',
                                             replicaTables = ( 'replicaStats', ),
                                             replicaConds = ( ( 'all', 'True' ),
                                                              ( 'hi', 'replicaStats.causalAlleleFreq >= .5' ),
                                                              ( 'lo', 'replicaStats.causalAlleleFreq < .5' ) ),
                                             condsFileFN = 'commonReplicaConds.tsv' ) )
                          
    
def splitSnpStatsFile( Ddata, scenario, inFileFN, condsFileFN, condNames, thinSfx = '',
                       replicaColName = 'Chrom', sfx = '', getio = None ):
    """Split a file containing per-snp data for all replicas, into separate files containing the same data for each
    kind of replica."""

    if not os.path.dirname( inFileFN ): inFileFN = os.path.join( Ddata, scenario.scenDir(), inFileFN )
    if not os.path.dirname( condsFileFN ):
        condsFileFN = os.path.join( Ddata, 'replicastats' + thinSfx, scenario.scenDir(), condsFileFN )

    outFileFNs = [ AddFileSfx( inFileFN, sfx, condName ) for condName in condNames ]

    if getio: return dict( depends_on = ( inFileFN, condsFileFN ),
                           creates = outFileFNs, mediumRuleNameSfx = scenario.scenDir() )

    condsFile = IDotData( condsFileFN )
    inFile = IDotData( inFileFN )

    with contextlib.nested( *map( functools.partial( IDotData.openForWrite, headings = inFile.headings ),
                                  outFileFNs ) ) as outFiles:

        for (replica, replicaRows), condValues in izip( inFile.groupby( replicaColName, multiPass = False ), condsFile ):

            assert condValues.replicaNum == replica

            # if this replica matches more than one condition, save the replica rows so we can iterate over them more
            # than once
            if sum( condValues[ condNames ] ) > 1: replicaRows = tuple( replicaRows )

            for condName, outFile in zip( condNames, outFiles ):
                if condValues[ condName ]: outFile.writeRecords( replicaRows )


                    
def joinSnpStatsFiles( Ddata, scenario, outFileFN, condNames, condsFileFN, thinSfx = '',
                       replicaColName = 'Chrom', sfx = '', getio = None ):
    """Join several per-snp stats files, each containing data for some of the replicas,
    into a single file containing data for all the replicas.
    """

    if not os.path.dirname( outFileFN ): outFileFN = os.path.join( Ddata, scenario.scenDir(), outFileFN )
    if not os.path.dirname( condsFileFN ): condsFileFN = os.path.join( Ddata, 'replicastats' + thinSfx, scenario.scenDir(),
                                                                       condsFileFN )
                                                                       
    inFileFNs = [ AddFileSfx( outFileFN, sfx, condName ) for condName in condNames ]

    if getio: return dict( depends_on = [ condsFileFN ] + inFileFNs, creates = outFileFN,
                           mediumRuleNameSfx = scenario.scenDir() )

    inFiles = map( IDotData, inFileFNs )
    dbg( 'inFiles' )

    condsFile = IDotData( condsFileFN )
    groupIters = [ inFile.groupby( replicaColName ) for inFile in inFiles ]

    def getBlocks():

        for r in condsFile:
            for condName, groupIter in zip( condNames, groupIters ):
                if r[ condName ]:
                    replicaNum, replicaRows = next( groupIter )
                    assert replicaNum == r.replicaNum
                    yield replicaRows
                    break

    IDotData.vstackFromIterable( getBlocks() ).save( outFileFN )

def ScenAttrs( scen ):
    """Make a dictionary describing the attributes of a scenario"""
    scenAttrs = dict( scen = scen, is_neutral = scen.is_neutral(), isNeutral = scen.isNeutral() )
    if not scen.is_neutral(): scenAttrs.update( mutAge = scen.mutAge,
                                                mutPop = scen.mutPop,
                                                mutFreq = scen.mutFreq )
    return scenAttrs


    
def scatterPlotReplicaStatistic( Ddata, nreplicas, replicaStatX,
                                 replicaStatY,
                                 outFile,
                                 thinSfx = '', 
                                 scenCond = 'True',
                                 replicaTables = (), replicaCond = 'True',
                                 replicaColorings = (),
                                 replicaDefaultColor = 'b',
                                 replicaShow = None,
                                 allScens = tuple( GetScenarios() ), nameSfx = '',
                                 scen2sfxs = {},
                                 title = '', subtitle = '',
                                 highlightScen = None, highlightReplica = None,
                                 xbound = None, ybound = None,
                                 getio = None ):
    """Draw a scatter plot where for each replica we have a pair of values.

    """

    args = Dict( 'Ddata thinSfx replicaTables scenCond scen2sfxs replicaCond allScens' )
    if getio: return dict( depends_on = findReplicasMatchingConds( getio = True, **args )[ 'depends_on' ],
                           creates = outFile,
                           name = 'scatterPlotReplicaStatistic' + Sfx( nameSfx ),
                           mediumRuleNameSfx = ( replicaStatX, replicaStatY ), attrs = dict( piperun_short = True ) )

    x = []
    y = []
    urls = []
    nskipped = 0
    colors = []

    if IsSeq( replicaShow ): replicaShow = '"_".join(map(str,["%.2f" % v if isinstance(v,float) else v for v in (' + ','.join( map( str, replicaShow ) ) + ')]))'
    
    for r in findReplicasMatchingConds( showHeadings = ( 'valX', 'valY', 'valShow' ) + tmap( operator.itemgetter( 0 ),
                                                                                             replicaColorings ),
                                        showVals = ( replicaStatX, replicaStatY,
                                                     replicaShow if replicaShow is not None else '0' ) +
                                        tmap( operator.itemgetter( 1 ), replicaColorings ),
                                        **args ):
        x.append( r.valX )
        y.append( r.valY )
        urls.append( '%s_%d_x=%s_y=%s' % ( r.scenario, r.replicaNum,
                                           '%.2f' % r.valX if isinstance( r.valX, float ) else r.valX,
                                           '%.2f' % r.valY if isinstance( r.valY, float ) else r.valY ) +
                     ( ( '' if str( r.valShow).startswith('_') else '_' ) + str( r.valShow ) if replicaShow else '' ) )
        if replicaColorings:
            colorHere = None
            for name, cond, color in replicaColorings:
                if r[ name ]:
                    colorHere = color
                    break
            colors.append( colorHere if colorHere is not None else replicaDefaultColor )

    pp.scatter( **Dict( 'x y urls', c = colors if colors else 'b' ) )
    pp.axis( 'equal' )

    if xbound: pp.gca().set_xbound( *xbound )
    if ybound: pp.gca().set_ybound( *ybound )

    if not xbound and not ybound:
        start = min( min(x), min(y) )
        rng = max( max( x ) - min( x ), max(y) - min(y) )
        pp.plot( [ start, start+rng ], [ start, start+rng ], 'g--' )
        
    pp.xlabel( replicaStatX )
    pp.ylabel( replicaStatY )
    if title: pp.title( title )
    pp.savefig( outFile )

def findTableFiles( Ddata, thinSfx, whichStats, tables, scenCond, allScens, scen2sfxs ):
    """Return table files used in conditions"""

    tables = MakeSeq( tables )
    scen2sfxs = dict( scen2sfxs )
    
    scenCondExpr = compile_expr( scenCond )
    ourScens = [ scen for scen in allScens if eval( scenCondExpr, globals(), ScenAttrs( scen ) ) ]

    depends_on = []
    scen2table2file = {}

    for scen in ourScens:
        thisScenDict = {}
        for table in tables:
            # identify the scenario-specific suffix for this table
            scenSfx = DictGet( scen2sfxs, scen, '' )
            if scenSfx:
                if IsSeq( scenSfx ): scenSfx = dict( scenSfx )
                if not isinstance( scenSfx, types.StringTypes ): scenSfx = DictGet( dict( scenSfx ),
                                                                                    os.path.splitext( table )[0], '' )
            
            tableFile = os.path.join( Ddata, whichStats+ thinSfx, scen.scenDir(),
                                             AddFileSfx( table + ( '.tsv' if '.' not in table else '' ),
                                                         scenSfx )
                                             + ( '/' if table.endswith( '.data' ) else '' ) )
            depends_on.append( tableFile )
            thisScenDict[ table ] = tableFile
        scen2table2file[ scen ] = thisScenDict

    tableNames = map( operator.itemgetter( 0 ), map( os.path.splitext, tables ) )
    return tableNames, tables, ourScens, scen2table2file, depends_on
                

def FindChromCol( iDotData ):
    """Find the column representing the replica or chromosome, based on our conventions."""
    return 'replicaNum' if 'replicaNum' in iDotData.headings else ( 'Chrom' if 'Chrom' in iDotData.headings else 'chrom' )

def FindPosCol( iDotData ):
    """Find the column representing the SNP position, based on our conventions."""
    return 'Pos' if 'Pos' in iDotData.headings else 'pos'


class NameCollector(ast.NodeVisitor):
    """Gather table names used in an expression"""

    def __init__(self):
        self.names = []

    def visit_Name(self, node):
        self.names.append( node.id )

    @staticmethod
    def getNamesIn( expr ):
        nc = NameCollector()
        nc.visit( ast.parse( expr ) )
        return tuple( set( nc.names ) )

def FindTables( *exprs ):
    """Find tables referenced in specified expressions"""

    return tuple( set( reduce( operator.concat, map( NameCollector.getNamesIn, exprs ) ) ) - set( ( 'True', 'False' ) ) )
    

def findReplicasMatchingConds( Ddata, 
                               replicaTables = None, replicaCond = 'True',
                               outFile = None,
                               scenCond = 'True',
                               showHeadings = (),
                               showVals = (),
                               allScens = GetScenarios(),
                               scen2sfxs = {},
                               thinSfx = '',
                               getio = None ):
    """Make an IDotData containing specified per-replica values for replicas meeting specified conditions."""

    dbg( '"findReplicasMatchingConds" scenCond replicaTables replicaCond showHeadings showVals scen2sfxs' )

    if replicaTables is None: replicaTables = FindTables( replicaCond, *MakeSeq( showVals ) )

    replicaTables = tuple( set( MakeSeq( replicaTables ) ) )
    
    replicaTableNames, replicaTables, ourScens, scen2table2file, depends_on = \
        findTableFiles( whichStats = 'replicastats', tables = replicaTables,
                        **Dict( 'Ddata thinSfx scenCond allScens scen2sfxs' ) )
    
    if getio: return dict( depends_on = depends_on,
                           creates = outFile,
                           attrs = dict( piperun_short = True ) )

    replicaCondExpr = compile_expr( replicaCond )
    showVals = MakeSeq( showVals )
    showValsExpr = map( compile_expr, showVals )
    if not showHeadings:
        showHeadings = map( MakeAlphaNum, showVals )
        showHeadings2 = []
        for h in showHeadings:
            h_new = h
            i = 1
            while h_new in showHeadings2:
                h_new = h + Sfx( i )
                i += 1
            showHeadings2.append( h_new )
        showHeadings = showHeadings2
    
    def makeResult():
        yield ( 'scenario', 'replicaNum' ) + tuple( MakeSeq( showHeadings ) )

        numReplicasSkippedTot, numReplicasAllowedTot = 0, 0
        
        for scen in ourScens:

            logging.info( '"findReplicasMatchingConds" scen' )

            numReplicasSkipped, numReplicasAllowed = 0, 0
            
            thisScenDict = scen2table2file[ scen ]
            replicaTableVals = [ IDotData( thisScenDict[ replicaTable ] ) for replicaTable in replicaTables ]

            for replicaTableRows in \
                    IDotData.TableIterInnerJoinAuxAsTuples( tableIters = map( iter, replicaTableVals ),
                                                            cols = map( FindChromCol, replicaTableVals ),
                                                            blanks = ( None, ) * len( replicaTableVals ),
                                                            headingLens = map( IDotData.rootClass.numCols,
                                                                               replicaTableVals ) ):

                vdict = dict( zip( replicaTableNames, replicaTableRows ) )
                dbg( 'scen vdict' )
                    
                evalHere = lambda expr: eval( expr, globals(), vdict )
                if evalHere( replicaCondExpr ):
                    numReplicasAllowed += 1
                    yield [ scen.scenName(), replicaTableRows[0].replicaNum ] + map( evalHere, showValsExpr )
                else:
                    numReplicasSkipped += 1

            dbg( '"in_scenario" scen numReplicasSkipped numReplicasAllowed' )
            numReplicasSkippedTot += numReplicasSkipped
            numReplicasAllowedTot += numReplicasAllowed

        dbg( 'numReplicasSkippedTot numReplicasAllowedTot' )
                    
    r = IDotData.fromFn( makeResult )
    if outFile: r.save( outFile )
    return r


def findSnpsMatchingConds( Ddata, 
                           snpTables = (), snpCond = 'True', replicaTables = (), replicaCond = 'True',
                           outFile = None,
                           scenCond = 'True',
                           showHeadings = (),
                           showVals = (),
                           allScens = GetScenarios(),
                           scen2sfxs = {},
                           thinSfx = '',
                           getio = None ):
    """Make an IDotData containing specified per-replica values for SNPs meeting specified conditions in
    replicas meeting specified conditions."""

    snpTables = tuple( set( MakeSeq( snpTables ) ) )
    
    dbg( '"findSnpsMatchingConds" scenCond snpTables snpCond replicaTables replicaCond showHeadings showVals '
         'scen2sfxs' )

    replicaArgs = Dict( 'Ddata thinSfx scenCond allScens scen2sfxs' )
    snpTableNames, snpTables, ourScens, scen2table2file, depends_on = \
        findTableFiles( whichStats = 'snpStats', tables = snpTables, **replicaArgs )
    
    if getio: return dict( depends_on = depends_on + findTableFiles( whichStats = 'replicastats', tables = replicaTables,
                                                                     **replicaArgs )[-1],
                           creates = outFile )

    snpCondExpr = compile_expr( snpCond )
    showVals = MakeSeq( showVals )
    showValsExpr = map( compile_expr, showVals )
    if not showHeadings: showHeadings = map( MakeAlphaNum, showVals )

    numSnpsSkippedTot, numSnpsAllowedTot = 0, 0

    def makeResult():
        yield ( 'scenario', 'replicaNum', 'Pos' ) + tuple( MakeSeq( showHeadings ) )
        for scen in ourScens:

            dbg( '"findSnpsMatchingConds" scen ')

            numSnpsAllowed, numSnpsSkipped = 0, 0

            replicasHere = findReplicasMatchingConds( **MergeDicts( replicaArgs, 
                                                                    Dict( 'replicaTables replicaCond scenCond',
                                                                          allScens = ( scen, ) ) ) )

            replicasHereSet = frozenset( replicasHere.replicaNum )
            dbg( 'scen len(replicasHereSet) replicasHereSet' )

            thisScenDict = scen2table2file[ scen ]
            dbg( '#[ ( thisScenDict[ snpTable ] ) for snpTable in snpTables ]' )
            snpTableVals = [ IDotData( thisScenDict[ snpTable ] ) for snpTable in snpTables ]

            lastReplica = None
            lastReplicaResult = None
            replicaCol = FindChromCol( snpTableVals[ 0 ] )
            posCol = FindPosCol( snpTableVals[ 0 ] )

            numSnpsSkippedTot, numSnpsAllowedTot = 0, 0
            
            for snpTableRows in \
                    IDotData.TableIterInnerJoinAuxAsTuples( tableIters = map( iter, snpTableVals ),
                                                            cols = zip( map( FindChromCol, snpTableVals ),
                                                                        map( FindPosCol, snpTableVals ) ),
                                                            blanks = ( None, ) * len( snpTableVals ),
                                                            headingLens = map( IDotData.rootClass.numCols,
                                                                               snpTableVals ) ):

                thisReplica = snpTableRows[0][ replicaCol ]
                if thisReplica != lastReplica:
                    thisReplicaResult = ( thisReplica in replicasHereSet )
                    if not thisReplicaResult: dbg( '"SKIPPING_REPLICA" thisReplica' )
                    lastReplicaResult = thisReplicaResult
                    lastReplica = thisReplica
                    
                if thisReplicaResult:
                    localDict = dict( zip( snpTableNames, snpTableRows ) )
                    evalHere = lambda expr: eval( expr, globals(), localDict )
                    evalResult = evalHere( snpCondExpr )
                    if evalResult:
                        v = [ scen.scenName(), thisReplica, snpTableRows[0][ posCol ] ] \
                            + map( evalHere, showValsExpr )
                        numSnpsAllowed += 1
                        yield v
                    else: numSnpsSkipped += 1
            
            numSnpsSkippedTot += numSnpsSkipped
            numSnpsAllowedTot += numSnpsAllowed

            dbg( 'scen numSnpsSkippedTot numSnpsAllowedTot' )

        dbg( '"finalCount" numSnpsSkippedTot numSnpsAllowedTot' )

    r = IDotData.fromFn( makeResult )
    if outFile: r.save( outFile )
    return r


def gatherCausalStat( Ddata, scenario, snpStatFN, replicaCol = 'Chrom', posCol = 'Pos', getio = None ):
    """Gather a specified per-SNP statistic just for the causal SNPs, and write them out as a replicastat.
    """

    replicaStatFN = string.replace( snpStatFN, 'snpStats', 'replicastats', 1 )

    if getio: return dict( depends_on = snpStatFN, creates = replicaStatFN, attrs = dict( scenario = scenario.scenDir() ) )

    snpStatFile = IDotData( snpStatFN )
    with IDotData.openForWrite( replicaStatFN, snpStatFile.headings ) as replicaStatFile:
        for r in snpStatFile:
            if r[ posCol ] == CAUSAL_POS:
                replicaStatFile.writeRecord( r )

def DefineRulesTo_gatherCausalStat( pr, Ddata, scen2snpStatFN, posCol = 'Pos' ):
    """Define rules to gather a specified per-SNP statistic for the causal SNPs into a replica stat."""

    for scenario, snpStatFN in scen2snpStatFN.items():
        pr.addInvokeRule( invokeFn = gatherCausalStat, invokeArgs = Dict( 'Ddata scenario snpStatFN posCol' ) )
        
        
        

    
    
