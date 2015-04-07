"""Code for spatially localizing the selected SNP: identifying sub-intervals of the selected region that likely contain
the selected SNP."""

from __future__ import division
from Classes.DotData import DotData
from Operations.IDotData import IDotData
from Operations.MiscUtil import AddFileSfx, dbg, Dict
from Operations.Shari_Operations.localize.PopConsts import AllAges, AllFreqs, AllPops
from Operations.Shari_Operations.localize.Scenario import GetSelectionScenarios
import numpy as np
import os, functools, itertools, contextlib, operator
from scipy import interpolate

def localizeSpatiallyBySplineFitting( Ddata, scenario, nreplicas, thinSfx = '',
                                      putativeMutPop = None, complikeSfx = '', likesTableSfx = '',
                                      confidence = .9, minBins = 20, nbins = 200, smoothing = 0.0,
                                      getio = None ):
    """For each replica within a given scenario,
    localize the selected SNP spatially, by fitting a spline to a (smoothed version of) the CMS scores, dividing the
    region into bins, and finding the set of bins that cover 90% (or specified fraction of) area under the spline.

    Params:

       confidence - the spatially localized region will (hopefully) have this probability of containing the causal SNP;
         here, specifically, this means we'll include bins in the region that collectively cover this fraction of area
         under the posterior density curve.

       minBins - the region will include at least this many of the highest-average bins.
       
    """

    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats'+ thinSfx, scenario.scenDir() )
    if putativeMutPop == None: putativeMutPop = scenario.mutPop
    sfxs = ( putativeMutPop, complikeSfx, likesTableSfx )
    complikeFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', *sfxs ) )

    intervalsListFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervalsSplineList.tsv', *sfxs ) )
    intervalsStatsFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervalsSplineStats.tsv', *sfxs ) )

    posteriorSplineFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervalsSplineSpline.tsv', *sfxs ) )
    binInfoFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervalsSplineBinInfo.tsv', *sfxs ) )

    if getio:
        return dict( depends_on = complikeFN,
                     creates = ( intervalsListFN, intervalsStatsFN, posteriorSplineFN, binInfoFN ),
                     mediumRuleNameSfx = ( scenario.scenDir(), ) + sfxs,
                     fileDescrs = { intervalsListFN:
                                        'List of intervals in the region, one of which (hopefully) contains the causal SNP.'
                                    ' For each replica this table has one or more lines, giving intervals in that replica.',
                                    intervalsStatsFN: 'Per-replica statistic about confidence intervals in that replica',
                                    posteriorSplineFN: 'The results of interpolating posterior density as a spline',
                                    binInfoFN: 'Information about individual bins under the spline' } )
    
    complike = IDotData( complikeFN ).filter( lambda r: all( np.isfinite( ( r.iHS, r.meanFst, r.max_xpop ) ) ) )

    binInfoHeadings = 'replicaNum binNum binStart binEnd included binCenters binAvgCMS binMaxCMS binIntegral '\
        'binIntegralNormed binRank'

    nbins_orig = nbins

    with contextlib.nested( IDotData.openForWrite( posteriorSplineFN, 'replicaNum gdPos complikeExp' ),
                            IDotData.openForWrite( binInfoFN, binInfoHeadings ),
                            IDotData.openForWrite( intervalsListFN, 'replicaNum gdFrom gdTo gdSize bpFrom bpTo '
                                                   'bpSize binsInInterval binsArea binsMax binsMaxFrac binsAvg '
                                                   'binsMinRank binsMaxRank' ),
                            IDotData.openForWrite( intervalsStatsFN, 'replicaNum numSegments gdTotLen bpTotLen' ) ) \
                            as ( posteriorSplineFile, binInfoFile, intervalsListFile, intervalsStatsFile ):

        for replicaNum, complikeForReplica in complike.groupby( 'Chrom' ):

            gdMin = np.min( complikeForReplica.gdPos )
            gdMax = np.max( complikeForReplica.gdPos )
            binSize = ( gdMax - gdMin ) / nbins_orig
            binLefts = np.arange( gdMin, gdMax, binSize )
            binCenters = binLefts + binSize / 2
            binRights = binLefts + binSize
            dbg( 'len(complikeForReplica) replicaNum len(binLefts) len(binRights) len(binCenters)' )
            nbins = len( binLefts )

            #
            # Compute the mean and max CMS score in each gdPos bin
            #
            binAvgCMS = np.zeros( nbins )
            binMaxCMS = np.zeros( nbins )
            binNums = np.zeros( nbins, dtype = int )

            for bin, valsInBin in \
                    complikeForReplica.addCol( 'bin',
                                               map( functools.partial( min, nbins-1 ),
                                                    map( int,
                                                         ( complikeForReplica.gdPos - gdMin ) / binSize ))).groupby('bin'):

                binNums[ bin ] = bin
                if valsInBin:
                    binAvgCMS[ bin ] = np.mean( valsInBin.complikeExp )
                    binMaxCMS[ bin ] = max( valsInBin.complikeExp )
                else:
                    binAvgCMS[ bin ] = binAvgCMS[ bin-1 ]
                    binMaxCMS[ bin ] = binMaxCMS[ bin-1 ]

            # fit a spline to the function ( binCenter, binAvgCMS ), approximating the posterior probability of a SNP being
            # causal.
            posteriorDensitySpline = interpolate.splrep( binCenters, binAvgCMS, s = smoothing )

            splineX = np.arange( gdMin, gdMax, binSize / 2 )
            for x, y in zip( splineX, interpolate.splev( splineX, posteriorDensitySpline ) ):
                posteriorSplineFile.writeRecord( replicaNum, x, y )

            # compute the integral under the interpolated function, for each bin.
            binIntegral = np.zeros( nbins )
            for binLeft, binRight, binNum in zip( binLefts, binRights, binNums ):
                binIntegral[ binNum ] = interpolate.splint( binLeft, binRight, posteriorDensitySpline )

            # normalize the integral value above each bin so that the total posterior density of being causal
            # integrates to 1.0
            binIntegralNormed = binIntegral / sum( binIntegral )

            binsByIntegralSize = binIntegral.argsort()[::-1]

            binRank = np.zeros( nbins )
            for i in range( nbins ):
                binRank[ binsByIntegralSize[ i ] ] = i

            binsToUse = np.zeros( nbins, dtype = bool )
            binsIncluded = 0
            fractionCovered = 0.0
            for bin in binsByIntegralSize:
                if fractionCovered >= confidence and binsIncluded >= minBins: break
                binsToUse[ bin ] = True
                binsIncluded += 1
                fractionCovered += binIntegralNormed[ bin ]

            # list the confidence intervals for this replica,
            # merging adjacent intervals.

            gd2bp = interpolate.interp1d( complikeForReplica.gdPos, complikeForReplica.Pos )
            gdMax = max( complikeForReplica.gdPos )

            binInfo = IDotData( names = binInfoHeadings, 
                                Columns = ( itertools.repeat( replicaNum, nbins ),
                                            range( nbins ),
                                            binLefts, binRights, binsToUse,
                                            binCenters, binAvgCMS, binMaxCMS, binIntegral, binIntegralNormed, binRank
                                            ) )

            binInfoFile.writeRecords( binInfo )


            binsOverallMaxCMS = max( binInfo.binMaxCMS )

            gdTotLen = 0.0
            bpTotLen = 0
            numSegments = 0
            for included, bins in binInfo.groupby( 'included' ):
                if included:
                    gdFrom = min( bins.binStart )
                    gdTo = min( max( bins.binEnd ), gdMax )
                    bpFrom = int( gd2bp( gdFrom ) )
                    bpTo = int( gd2bp( gdTo ) )

                    gdTotLen += ( gdTo - gdFrom )
                    bpTotLen += ( bpTo - bpFrom )
                    numSegments += 1

                    intervalsListFile.writeRecord( replicaNum, gdFrom, gdTo, gdTo - gdFrom, bpFrom, bpTo, bpTo - bpFrom,
                                                   len( bins ), sum( bins.binIntegralNormed ),
                                                   max( bins.binMaxCMS ), max( bins.binMaxCMS ) / binsOverallMaxCMS,
                                                   np.mean( bins.binAvgCMS ),
                                                   min( bins.binRank ), max( bins.binRank ) )

            assert numSegments > 0 and bpTotLen > 0 and gdTotLen > 0.0
            intervalsStatsFile.writeRecord( replicaNum, numSegments, gdTotLen, bpTotLen )
            dbg( 'replicaNum numSegments gdTotLen bpTotLen' )


def gatherCausalSnpGdPos( Ddata, thinSfx, scenario, selpos = 500000, getio = None ):
    """For each replica in the scenario, save the genetic map position of the causal SNP
    as a replica statistic."""

    assert not scenario.isNeutral()
    
    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats'+ thinSfx, scenario.scenDir() )

    gdMapFN = os.path.join( snpStatsDir, 'gdMap.tsv' )
    causalGdPosFN = os.path.join( replicaStatsDir, 'causalGdPos.tsv' )

    if getio: return dict( depends_on = gdMapFN, creates = causalGdPosFN,
                           mediumRuleNameSfx = scenario.scenDir() )

    with IDotData.openForWrite( causalGdPosFN, 'replicaNum causalGdPos' ) as causalGdPosFile:
        gdMap = IDotData( gdMapFN )
        for r in gdMap[ gdMap.Pos == selpos ]:
            causalGdPosFile.writeRecord( r.Chrom, r.gdPos )

def DefineRulesTo_gatherCausalSnpGdPos( pr, Ddata, mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs,
                                        thinSfx = '' ):
    """Define rules to gather the causal SNP genetic position in every replica"""

    for scenario in GetSelectionScenarios( **Dict( 'mutAges mutFreqs mutPops' ) ):
        pr.addInvokeRule( invokeFn = gatherCausalSnpGdPos,
                          invokeArgs = Dict( 'Ddata scenario thinSfx' ) )
    
                
def evalSpatialLoc( Ddata, thinSfx, scenario, putativeMutPop, nreplicas, complikeSfx = '',
                    likesTableSfx = '', selpos = 500000, whichSpatialLoc = 'Spline',
                    getio = None ):
    """Evaluate spatial localization.  Compute relevant replica statistic.  For each replica,
    compute: whether the localized intervals include the causal SNP; statistics about the
    localized intervals; the position of the causal SNP relative to the localized intervals."""

    assert not scenario.isNeutral()
    
    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats'+ thinSfx, scenario.scenDir() )
    if putativeMutPop == None: putativeMutPop = scenario.mutPop
    sfxs = ( putativeMutPop, complikeSfx, likesTableSfx )

    intervalsListFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervals%sList.tsv' % whichSpatialLoc, *sfxs ) )
    causalGdPosFN = os.path.join( replicaStatsDir, 'causalGdPos.tsv' )
    spatialLocEvalFN = os.path.join( replicaStatsDir, AddFileSfx( 'spatialLocEval%s.tsv' % whichSpatialLoc, *sfxs ) )

    if getio: return dict( depends_on = ( intervalsListFN, causalGdPosFN ),
                           creates = spatialLocEvalFN,
                           mediumRuleNameSfx = scenario.scenDir(),
                           name = 'evalSpatialLoc_' + whichSpatialLoc )

    with IDotData.openForWrite( spatialLocEvalFN,
                                'replicaNum numIntervals totLenBp totLenGd causalIncluded '
                                'distanceToIntervalBoundaryBp distanceToIntervalBoundaryGd' ) as spatialLocEvalFile:

        for ( replicaNum2, replicaIntervals ), ( replicaNum3, causalGdPos, replicaNum4 ) in \
                itertools.izip( IDotData( intervalsListFN ).groupby( 'replicaNum' ),
                                IDotData.merge( iDotDatas = ( IDotData( causalGdPosFN ),
                                                              IDotData( intervalsListFN ).replicaNum.removeDups() ),
                                                cols = ( 'replicaNum', 'replicaNum' ) ) ):

            replicaNum2, replicaNum3, replicaNum4  = map( int, ( replicaNum2, replicaNum3, replicaNum4 ) )
            if not replicaNum2 == replicaNum3 == replicaNum4:
                dbg( 'replicaNum2 replicaNum3 replicaNum4 intervalsListFN complikeFN causalGdPosFN spatialLocEvalFN' )
            assert replicaNum2 == replicaNum3 == replicaNum4

            causalIncluded = False
            totLenBp = 0
            totLenGd = 0.0
            for replicaInterval in replicaIntervals:
                dbg( 'replicaInterval causalGdPos' )
                #assert bool( replicaInterval.bpFrom <= selpos <= replicaInterval.bpTo ) == bool( replicaInterval.gdFrom <= causalPos_gd <= replicaInterval.gdTo )
                assert ( replicaInterval.gdFrom <= causalGdPos <= replicaInterval.gdTo ) == ( replicaInterval.bpFrom <= selpos <= replicaInterval.bpTo )
                if replicaInterval.gdFrom <= causalGdPos <= replicaInterval.gdTo:
                    causalIncluded = True

                totLenGd += ( replicaInterval.gdTo - replicaInterval.gdFrom )
                totLenBp += ( replicaInterval.bpTo - replicaInterval.bpFrom )

            spatialLocEvalFile.writeRecord( replicaNum2, len( replicaIntervals ),
                                            totLenBp, totLenGd, int( causalIncluded ),
                                            np.min( ( np.min( np.abs( replicaIntervals.bpFrom - selpos ) ),
                                                      np.min( np.abs( replicaIntervals.bpTo - selpos ) ) ) ),
                                            np.min( ( np.min( np.abs( replicaIntervals.gdFrom - causalGdPos ) ),
                                                      np.min( np.abs( replicaIntervals.gdTo - causalGdPos ) ) ) ) )
            
def localizeSpatiallyByWindows(Ddata, scenario, nreplicas, thinSfx = '', putativeMutPop = None, complikeSfx = '',
                               likesTableSfx = '',
                               threshold = .5, numSNP = 1,
                               minGdInEachDir = .05,
                               fromReplica = None, toReplica = None, getio = None):
    """
    Spatially localize the selected variant for all replicas within a given scenario.
    The approach is to start with the highest-scoring SNP, and move left and right from it in fixed-size windows
    for as long as the windows contain at least 'numSNP' snps with score at least 'threshold'.

    Adapted from Operations.Shari_Operations.localize.hapmap_regions_0615.plotRegions() .
    """

    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats'+ thinSfx, scenario.scenDir() )
    if putativeMutPop == None: putativeMutPop = scenario.mutPop
    sfxs = ( putativeMutPop, complikeSfx, likesTableSfx )
    complikeFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', *sfxs ) )

    intervalsListFN = os.path.join( replicaStatsDir, AddFileSfx( 'intervalsWindowsList.tsv', *sfxs ) )

    if getio:
        return dict( depends_on = complikeFN,
                     creates = intervalsListFN,
                     mediumRuleNameSfx = ( scenario.scenDir(), ) + sfxs,
                     fileDescrs =
                     { intervalsListFN:
                           'List of intervals in the region, one of which (hopefully) contains the causal SNP.'
                       ' For each replica this table has one or more lines, giving intervals in that replica.' } )
    

    #complike = IDotData( complikeFN ).filter( lambda r: all( np.isfinite( ( r.iHS, r.meanFst, r.max_xpop ) ) ) )
    complike = IDotData( complikeFN )

    with IDotData.openForWrite( intervalsListFN,
                                'replicaNum gdFrom gdTo gdSize bpFrom bpTo bpSize numPositiveBins '
                                'numSnpsOver_0_2 maxSNP_Pos maxSNP_lik' ) as intervalsListFile:

        for replicaNum, complikeForReplica in complike.groupby( 'Chrom' ):


            if fromReplica is not None and replicaNum < fromReplica: continue
            if toReplica is not None and replicaNum > toReplica: break

            X = complikeForReplica.toDotData()

            minPos = np.min(X.gdPos)
            maxPos = np.max(X.gdPos)
            bins = np.arange(0,1,.01)

            ind = X.complikeExp.argsort()
            lik = X.complikeExp[ind]
            maxlik = np.mean(lik[-5:])
            #maxlik = mean(lik[-5:])
            like = X.complikeExp / maxlik

            Y = X[ind]
            maxSNP = Y.Pos[-1]

            maxScore = Y.complikeExp[-1]
            relPos = X.gdPos - Y.gdPos[-1]
            X = X.hstack(DotData(Columns = [like,relPos],names=['scaled_like','relPos']))

            topGdPos = Y.gdPos[-1]
            minPos = Y.Pos[-1]
            maxPos = Y.Pos[-1]
            minGdPos = topGdPos
            maxGdPos = topGdPos
            numPositiveBins = 0

            dbg( 'replicaNum minPos maxPos minGdPos maxGdPos' )

            for dir in -1, +1:
                for bin in bins:
                    Z = X[np.abs(X.relPos - dir*bin) <= .02]

                    dbg( 'dir bin len(Z)' )

                    if len(Z) == 0: continue
                    top = Z[Z.scaled_like > threshold]

                    dbg( 'len(top)' )

                    if len(top) <= numSNP and np.abs( topGdPos - top.gdPos ) > minGdInEachDir: break
                    if len( top ) == 0: top = Z
                    if dir == -1:
                        minPos = np.min(top.Pos)
                        minGdPos = np.min(top.gdPos)
                    else:
                        maxPos = np.max(top.Pos)
                        maxGdPos = np.max(top.gdPos)
                    numPositiveBins += 1

                    
            ind = np.all([X.Pos > minPos, X.Pos < maxPos],axis=0)
            peak = X[ind]
            
            intervalsListFile.writeRecord( replicaNum, minGdPos, maxGdPos, maxGdPos - minGdPos,
                                           minPos, maxPos, maxPos - minPos,
                                           numPositiveBins, sum(peak.scaled_like > .2),
                                           maxSNP, maxScore )
                                            
                           
    
    
    
                


