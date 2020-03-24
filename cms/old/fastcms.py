"""Fast CMS computation"""



import os, logging, functools
from functools import reduce
#import blaze as bz
import numpy as np
import pandas as pd
#from into import into
from Operations.MiscUtil import Dict, dbg, AddFileSfx, MakeSeq, StatKeeper
from Operations.tsvutils import DefineRulesTo_computeMeanStd, DefineRulesTo_normalizeOneColumn, \
    computeMeanStd_binned_tsvs, normalizeInBins_tsv
from Operations.Shari_Operations.localize import subs
from Operations.Shari_Operations.localize.fstBySNP_Npops import fst_onePopPair
from Operations.Shari_Operations.localize.CMS import CMSBins
from Operations.bioutil import genomeBuild2genMapSfx

def getFN_xpop_signif( sweepDir, chrom, pop1, pop2 ):
    """Return filename of xpop significance scores from Sweep output"""
    return os.path.join( sweepDir, 'analysis', 'chr%(chrom)s' % locals(),
                         'xpop_significance_%(pop1)s_%(pop2)s.tsv' % locals() )

def getFN_ihs_signif( sweepDir, chrom, pop ):
    """Return filename of xpop significance scores from Sweep output"""
    return os.path.join( sweepDir, 'analysis', 'chr%(chrom)s' % locals(),
                         'ihs_significance_%(pop)s.tsv' % locals() )

def gatherXPOPscores( pops, chrom, selPop, sweepDir, outFN, getio = None ):
    """Gather xpop scores into a convenient form."""

    pops = [p for p in pops if p != selPop]
    pop2FN = dict([ ( pop, getFN_xpop_signif( pop1 = selPop, pop2 = pop, sweepDir = sweepDir, chrom = chrom ) )
                    for pop in pops ])

    if getio: return dict( depends_on = list(pop2FN.values()), creates = outFN, attrs = Dict( 'chrom', pop = pops, piperun_short = True ) )

    def LoadComparison( pop ):
        """Load comparison with one pop"""

        d0 = pd.read_csv( pop2FN[ pop ], sep = '\t',
                          usecols = ( 'Pop 1', 'Pop 2', 'Chrom' ), nrows = 1 )
        dbg( 'd0' )

        assert str( d0.loc[ 0, 'Chrom' ] ) == str( chrom )
        assert ( d0.loc[ 0, 'Pop 1'] == selPop and d0.loc[0, 'Pop 2'] == pop ) or ( d0.loc[0,'Pop 1'] == pop and d0.loc[0,'Pop 2'] == selPop )

        flip = ( d0.loc[0,'Pop 1'] == pop )
        
        d = pd.read_csv( pop2FN[ pop ], sep = '\t',
                         usecols = ( 'SNP pos (bases)', 'L AllEHH logratio Deviation', 'R AllEHH logratio Deviation' ),
                         index_col = 'SNP pos (bases)',
                         na_values = ( '-', ) )
        d.info()

        if flip:
            d[ 'L AllEHH logratio Deviation' ] *= -1
            d[ 'R AllEHH logratio Deviation' ] *= -1

        return pd.DataFrame.from_dict({ pop : d.max( axis = 1 ) })

    # end: def LoadComparison( pop )

    comparisons = reduce( lambda d1, d2: d1.join( d2, how = 'inner' ),
                          list(map( LoadComparison, pops )) ).max( axis = 1,
                                                             columns = ( 'max_xpop', ) )
    comparisons.index.name = 'pos'
    comparisons.name = 'max_xpop'
#    print 'type of comparisons is', type(comparisons)
#    print comparisons
    comparisons.to_csv( outFN, sep = '\t' , header = True )

# end: def gatherXPOPscores( pops, chrom, selPop, sweepDir, outFN, getio = None ):

def gather_snp_info( pops, pop2snpInfoFN, pop2ancFreqFN, pop2sampleSizeFN, getio = None ):
    """Gather SNP freq info"""
    
    if getio: return dict( depends_on = list(pop2snpInfoFN.values()), creates = ( pop2ancFreqFN, pop2sampleSizeFN ),
                           attrs = dict( pop = pops ) )


    pop2ancFreq = pd.DataFrame( data =
                                dict([ ( pop, pd.read_csv( pop2snpInfoFN[ pop ],
                                                           sep = '\t', usecols = ( 'SNP pos (bases)', 'Ancestral Freq' ),
                                                           index_col = 'SNP pos (bases)' )[ 'Ancestral Freq' ] )
                                       for pop in pops ]) )

    pop2ancFreq.dropna( inplace = True )
    pop2ancFreq.to_csv( pop2ancFreqFN, sep = '\t', index_label = 'pos' )

    def getSampleSize( pop ):
        z = pd.read_csv( pop2snpInfoFN[ pop ],
                         sep = '\t', usecols = ( 'A0count', 'A1count' ), nrows = 1 )
        return z.at[ 0, 'A0count' ] + z.at[ 0, 'A1count' ]

    pop2sampleSize = pd.Series( dict([ ( pop, getSampleSize( pop ) ) for pop in pops ]),
                                name = 'sampleSize' )
    pop2sampleSize.to_csv( pop2sampleSizeFN, sep = '\t', header = True, index_label = 'pop' )

# end: def gather_snp_info( pops, pop2snpInfoFN, pop2ancFreqFN, pop2sampleSizeFN, getio = None )

def gather_iHS_scores( selPop, chrom, ihsFN, pop2ancFreqFN, ihsOutFN, dihhOutFN, getio = None ):
    """Gather iHS scores"""

    if getio: return dict( depends_on = ( ihsFN, pop2ancFreqFN ), creates = ( ihsOutFN, dihhOutFN ), 
                           attrs = Dict( 'chrom', pop = selPop, piperun_short = True ) )

    d0 = pd.read_csv( ihsFN, sep = '\t',
                      usecols = ( 'Population', 'Chrom' ), nrows = 1 )
    dbg( 'd0' )

    assert str( d0.loc[ 0, 'Chrom' ] ) == str( chrom )
    assert d0.loc[ 0, 'Population'] == selPop
        
    d = pd.read_csv( ihsFN, sep = '\t',
                     usecols = ( 'SNP pos (bases)', 'Ancestral Freq', 'Both iHS', 'Both iHH_D', 'Both iHH_A' ),
                     index_col = 'SNP pos (bases)',
                     na_values = ( '-', ) )
    
    d.index.name = 'pos'

    pop2ancFreq = pd.read_table( pop2ancFreqFN, index_col = 'pos', usecols = ( 'pos', selPop, ) )
#    snp2ancFreq = pd.read_table( snpInfoFN, index_col = 'SNP pos (bases)',
#                                 usecols = ( 'SNP pos (bases)', 'Ancestral Freq' ) )
#    snp2ancFreq.dropna( inplace = True )
    # dbg( 'len(pop2ancFreq) len(snp2ancFreq) pop2ancFreq.index.difference(snp2ancFreq.index)' )
    # dbg( 'len(pop2ancFreq) len(snp2ancFreq) snp2ancFreq.index.difference(pop2ancFreq.index)' )
    # dbg( 'np.all(pop2ancFreq.index.values==snp2ancFreq.index.values)' )
    # dbg( 'np.sum(pop2ancFreq.index.values==snp2ancFreq.index.values)' )
    # dbg( 'len(pop2ancFreq.index.values) len(snp2ancFreq.index.values)' )
#    pop2ancFreq.index.name = 'pos'
#    dbg( '3 pop2ancFreq selPop pop2ancFreq.columns' )
    pop2ancFreq.rename( columns = { selPop : selPop + '_ancFreq' }, inplace = True )
#    dbg( '4 pop2ancFreq' )

#    print "ii:", pop2ancFreq.info()
#    pop2ancFreq.to_csv( 'pf.tsv', sep = '\t', header = True, na_rep = 'NaN' )
#    dbg( '1 d' )
    
    d = d.join( pop2ancFreq, how = 'right', sort = True )
    
#    dbg( '2 d' )

#    af1 = d['Ancestral Freq']
    af2 = d[selPop + '_ancFreq']
#    dbg( '"GGGGGGGGGG" (af1-af2).max() (af1.isnull()==af2.isnull()).all()' )

    d_iHS = pd.DataFrame( data = dict( iHS = d[ 'Both iHS' ] ) )
    d_iHS.to_csv( ihsOutFN, sep = '\t', header = True, na_rep = 'NaN' )

#    dihh = subs.normalizeByFreq( rawVals = ( d[ 'Both iHH_D' ] - d[ 'Both iHH_A' ] ).values,
#                                 ancfreq = 1.0 - af2.values )
    d_iHH = pd.DataFrame( data = dict( iHHDiff = d[ 'Both iHH_D' ] - d[ 'Both iHH_A' ],
                                       normingFreqs = 1.0 - af2 ) )
    
    d_iHH.to_csv( dihhOutFN, sep = '\t', header = True, na_rep = 'NaN' )

# end: def gather_iHS_scores( selPop, chrom, ihsFN, snpInfoFN, ihsOutFN, dihhOutFN, getio = None ):

def computeMeanFstAndFreqDiffScores( pops, chrom, selPop, sweepDir,
                                     pop2ancFreqFN, pop2sampleSizeFN, outMeanFstFN, outFreqDiffFN, getio = None ):
    """Compute meanFst and freqDiff scores"""

    if selPop not in pops: pops = tuple( MakeSeq( pops ) ) + ( selPop, )
    cmpPops = [ pop for pop in pops if pop != selPop ]


    if getio: return dict( depends_on = ( pop2ancFreqFN, pop2sampleSizeFN ),
                           creates = ( outMeanFstFN, outFreqDiffFN ),
                           attrs = Dict( 'chrom', pop = pops ) )    

#    pop2ancFreq.to_csv( 'befdrop.tsv', sep = '\t' )
#    pop2ancFreq.fillna( value = 1.0, inplace = True )
    
#    pop2ancFreq.to_csv( 'aftdrop.tsv', sep = '\t' )

    pop2ancFreq = pd.read_table( pop2ancFreqFN, index_col = 'pos' )
    pop2sampleSize = pd.read_table( pop2sampleSizeFN, index_col = 'pop' ).sampleSize

    dbg( 'pop2sampleSize' )

    #pop2snpInfo.to_csv( 'test.tsv', sep = '\t', header = True )

    derFreq = 1.0 - pop2ancFreq[ selPop ]
    cmpAncFreqs = pop2ancFreq[ [ pop for pop in pops if pop != selPop ] ]
    meanAnc = cmpAncFreqs.mean( axis = 1 )
    freqDiff = derFreq - ( 1.0 - meanAnc )
    freqDiff.name = 'freqDiff'
    freqDiff.to_csv( outFreqDiffFN, sep = '\t', header = True )

    # compute meanFst

#    dbg( '"vvvvvvvvvvvw" selPop pop2ancFreq[selPop] pop2ancFreq["JPT+CHB"] pop2ancFreq["YRI"]' )
#    dbg( 'selPop pop2sampleSize[selPop] pop2sampleSize["JPT+CHB"] pop2sampleSize["YRI"]' )
    d = dict([ ( pop, fst_onePopPair( ancFreqs = np.array( ( pop2ancFreq[ selPop ], pop2ancFreq[ pop ] ) ),
                                      sampleSizes = ( pop2sampleSize[ selPop ], pop2sampleSize[ pop ] ) ) )
               for pop in cmpPops ])
    fstVals = pd.DataFrame( data = d, index = pop2ancFreq.index )
#    spc = fst_onePopPair( ancFreqs = np.array( ( pop2ancFreq[ 'BEB' ], pop2ancFreq[ 'ASN' ] ) ),
#                          sampleSizes = ( pop2sampleSize[ 'BEB' ], pop2sampleSize[ 'ASN' ] ) )
#    dbg( '"ddddddddddd" fstVals.loc[526736] spc' )
#    dbg( 'fstVals' )
    fstVals.fillna( value = 0.0, inplace = True )
    #fstVals.to_csv( 'fstvals.tsv', sep = '\t', header = True, na_rep = 'NaN' )
    fstMean = fstVals.mean( axis = 1 )
    dbg( 'fstVals fstMean' )
    fstMean.name = 'meanFst'
    fstMean.to_csv( outMeanFstFN, sep = '\t', header = True, na_rep = 'NaN' )
                                                
# end: def computeMeanFstAndFreqDiffScores( pops, chrom, selPop, sweepDir, outMeanFstFN, outFreqDiffFN, getio = None )    

def computeLikeRatioForStat_do( statVals, hitLikes, missLikes, bins ):
    """Compute likes ratio"""

    # Precompute the likelihood ratio corresponding to each bin
    
    dbg( 'statVals hitLikes missLikes bins' )
         
    indNaN = hitLikes != 1e-10
    missingVal = np.log( np.min( hitLikes[indNaN] / missLikes[indNaN] ) )

    
    CLR = [ ( np.log( hitLike / missLike ) if hitLike != 1e-10 else missingVal ) if hitLike != 0.0 else np.nan
            for hitLike, missLike in zip( hitLikes, missLikes ) ]
    CLR = np.array( [ CLR[0] ] + CLR + [ CLR[-1] ] )

    binIds = np.digitize( statVals.values, bins )
    st_binSize = ( bins[1] - bins[0] )
    st_nbins = len( bins ) -1
    binIds2 = np.where( np.isfinite( statVals.values ),
                        np.clip( ( ( statVals.values - bins[0] ) / st_binSize ).astype( np.int16 ), 0, st_nbins-1 ),
                        len( hitLikes ) ) + 1

    return np.where( np.isnan( statVals.values ),
                     np.repeat( np.nan, len( statVals ) ), CLR[ binIds ] ), binIds, binIds2

# end: def computeLikeRatioForStat_do( statVals, hitLikes, missLikes, bins )    

def computeLikeRatioForStat( stat, statValsFN, hitLikesFN, missLikesFN, stat_start, stat_end, stat_nbin, statLikesRatioFN, getio = None ):
    """Compute likes for one stat"""

    if getio: return dict( depends_on = ( statValsFN, hitLikesFN, missLikesFN ), creates = statLikesRatioFN,
                           uses = computeLikeRatioForStat_do )

    statVals = pd.read_table( statValsFN )
    hitLikes = pd.read_table( hitLikesFN )[ stat ]
    missLikes = pd.read_table( missLikesFN )[ stat ]

    bins = np.linspace( stat_start, stat_end, stat_nbin+1 )

    statLikeRatio, statBinIds, statBinIds2 = computeLikeRatioForStat_do( statVals = statVals[ stat ],
                                                                         **Dict( 'hitLikes missLikes bins' ) )
    statVals[ stat + 'likeRatio' ] = statLikeRatio
    statVals[ stat + 'Bin' ] = statBinIds
    statVals[ stat + 'Bin2' ] = statBinIds2

    statVals.to_csv( statLikesRatioFN, sep = '\t', columns = ( 'pos', stat, stat + 'likeRatio', stat + 'Bin',
                                                               stat + 'Bin2' ), index = False,
                     na_rep = 'NaN' )

# end: def computeLikeRatioForStat

## Bins values by frequencies and then normalizes within bins        
def normalizeByFreq_getMeanStd(rawVals,ancfreq,stdKeeper,meanKeepers):

    Frequency = np.arange(0.05, 1.05, 0.05)
    der_freq = 1 - ancfreq

    stdKeeper.addVals(rawVals)
    for i in range(len(Frequency)):
        idx = ((Frequency[i] - der_freq) < .05 ) & ( Frequency[i] - der_freq > 0)
#        dbg( 'i Frequency[i] np.sum(idx)' )
        meanKeepers[i].addVals( rawVals[  idx ] )

def normalizeByFreq_getMeanStd_tsv(iHHDiffFNs, globalStatFN, binsStatFN, getio = None):
    """Compute mean and stddev for normalizing ihhdiff within freqs"""

    if getio: return dict( depends_on = iHHDiffFNs, creates = ( globalStatFN, binsStatFN ) )

    stdKeeper = StatKeeper()
    meanKeepers = [ StatKeeper() for i in range( 20 ) ]
    for f in iHHDiffFNs:
        d = pd.read_table( f )
        normalizeByFreq_getMeanStd( d.iHHDiff.values, 1.0 - ( 1.0 - d.normingFreqs.values ), stdKeeper, meanKeepers )

#    dbg( '"ZZZZZZ" iHHDiffFNs stdKeeper.getStd map(StatKeeper.getMean,meanKeepers) map(StatKeeper.getCount,meanKeepers)' )

    pd.DataFrame( dict( std = ( stdKeeper.getStd(), ) ) ).to_csv( globalStatFN, sep = '\t', na_rep = 'NaN',
                                                                  header = True, index = False )
    pd.DataFrame( dict( mean = list(map( StatKeeper.getMean, meanKeepers )) ) ).to_csv( binsStatFN,
                                                                                  sep = '\t', na_rep = 'NaN',
                                                                                  header = True, index_label = 'binId' )

def normalizeByFreq_compute_normed(rawVals,ancfreq, StdDev, expectation):

    Frequency = np.arange(0.05, 1.05, 0.05)
    #print Frequency

    der_freq = 1 - ancfreq
    normVal = np.repeat( np.nan, len( rawVals ) )

    # Bookkeeping

    dbg( 'StdDev' )
    dbg( 'der_freq' )
    for i in range(len(Frequency)):
        idx = ((Frequency[i] - der_freq) < .05) & ( (Frequency[i] - der_freq) >= 0 ) & np.isfinite( rawVals )
        normVal[ idx ] = (rawVals[ idx ] - expectation[i])/StdDev
#        dbg( '"KKKKK" i Frequency[i] expectation[i] idx.nonzero() rawVals[idx] normVal[idx]' )
    
    return normVal

def normalizeByFreq_compute_normed_tsv(iHHDiffFN, globalStatFN, binsStatFN, StdDiffFN, getio = None):
    """Computed normed iHHDiff"""

    if getio: return dict( depends_on = ( iHHDiffFN, globalStatFN, binsStatFN ),
                           creates = StdDiffFN )

    d = pd.read_table( iHHDiffFN )
    gstat = pd.read_table( globalStatFN )
    binsStat = pd.read_table( binsStatFN )

    normVal = normalizeByFreq_compute_normed( rawVals = d.iHHDiff.values,
                                              ancfreq = d.normingFreqs.values,
                                              StdDev = gstat[ 'std' ].iloc[0],
                                              expectation = binsStat[ 'mean' ].values )

    d[ 'StdDiff' ] = normVal
    d.to_csv( StdDiffFN, sep = '\t', header = True, index = False, na_rep = 'NaN' )

def computeMeanStd( inFNs, colName, outFN, getio = None ):
    """Compute mean and std using blaze"""
    if getio: return dict( depends_on = inFNs, creates = outFN )

    filenames = inFNs
    dbg( 'inFNs' )

    sk = StatKeeper()
    for f in filenames:
        dbg( 'f' )
        d = pd.read_table( f )
        dbg( 'f len(d)' )
        sk.addVals( d[ colName ].values )

    pd.DataFrame( dict( stat = 'mean std count numNaNs'.split(),
                        val = ( sk.getMean(), sk.getStd(), sk.getCount(), sk.getNumNaNs() ) ) ).to_csv( outFN,
                                                                                                        sep = '\t',
                                                                                                        index = False,
                                                                                                        na_rep = 'NaN' )


def addLikesRatios( inFNs, colNames, outFN, getio = None ):
    """Add up likes ratios"""

    if getio: return dict( depends_on = inFNs, creates = outFN )

    result = None
    for fn, colName in zip( inFNs, colNames ):
        d = pd.read_table( fn, index_col = 0 ).dropna()[ colName ]
        d.name = 'likesRatio'
        if result is None:
            result = d
        else:
            result += d

    result.to_csv( outFN, header = True, na_rep = 'NaN', sep = '\t' )
    
#    data = bz.chunks(bz.CSV)([bz.CSV(fn) for fn in filenames])

#    d = bz.Data(data)
#    into( outFN, pd.DataFrame( dict( mean = d[ colName ].mean(), std = d[ colName ].std() ) ) )


def joinStats( snpInfoFN, statLikesFNs, likesRatioFN, outFN, getio = None ):
    """Join stats into one file"""

    if getio:
        return dict( depends_on = ( snpInfoFN, likesRatioFN ) + tuple( MakeSeq( statLikesFNs ) ),
                     creates = outFN )

    snpInfo = pd.read_table( snpInfoFN, index_col = 'SNP pos (bases)' )
    snpInfo.index.rename( 'pos', inplace = True )

    statLikes = [ pd.read_table( statLikeFN, index_col = 'pos' ) for statLikeFN in statLikesFNs ]
    likesRatio = pd.read_table( likesRatioFN, index_col = 'pos' )

    result = snpInfo.join( statLikes + [ likesRatio ], how = 'outer' )
    result.info()
    dbg( 'result.describe()' )
    
    result.to_csv( outFN, sep = '\t', na_rep = 'NaN', header = True )

# end: def joinStats( snpInfoFN, stats, statFNs, statLikesFNs, outFN, getio = None )    


def DefineRulesTo_fastCMS( pr, pops, chroms, selPop, sweepDir, cmsDir, genomeBuild = 'hg19' ):
    """Define rules to do fast CMS computation.

    Params:

       pr - the PipeRun object to which to add rules

       selPop - testing selection in which pop?
       pops - comparing selPop to which pops?
       sweepDir - the sweep directory
       cmsDir - the directory under which CMS stats go
    """

    pops = list( MakeSeq( pops ) )
    if selPop not in pops: pops.append( selPop )

    allPops = tuple( MakeSeq( pops ) )
    if selPop not in allPops: allPops += ( selPop, )
    cmpPops = [ pop for pop in allPops if pop != selPop ]

    rawScoresFN = {}

    genMapSfx = genomeBuild2genMapSfx[ genomeBuild ]
    for pop in allPops:
        for chrom in chroms:
            with pr.settingAttrs( 'pop chrom' ):
                snpInfoFN = os.path.join( sweepDir, 'analysis/chr%(chrom)s/snps_%(pop)s.tsv' % locals() )
                projDir = os.path.join( sweepDir, 'data/chr%(chrom)s' % locals() )
                ancestralImportedFN = os.path.join( projDir, 'ancestral.tsv.imported' )
                genotypesImportedFN = os.path.join( projDir, 'genotypes_chr%(chrom)s_%(pop)s_r21_nr_fwd_phased_all.imported' % locals() )
                genMapImportedFN = os.path.join( projDir, 'genetic_map_chr%(chrom)s_%(genMapSfx)s.txt.imported' % locals() )
                pr.addRule( name = 'extractSnpInfo',
                            commands = 'java -classpath ../Other/Ilya_Other/sweep/sweepsrc/sweep.jar edu.mit.broad.sweep.Main ExtractAlleleFreqs %(projDir)s/project %(snpInfoFN)s %(pop)s %(chrom)s' % locals(),

                            commandsOld = 'java -classpath ../Other/Ilya_Other/sweep/sweepsrc/sweep/target/sweep-1.0-SNAPSHOT-jar-with-dependencies.jar edu.mit.broad.sweep.Main ExtractAlleleFreqs %(projDir)s/project %(snpInfoFN)s %(pop)s %(chrom)s' % locals(),                            
                            depends_on = ( ancestralImportedFN, genotypesImportedFN, genMapImportedFN ),
                            creates = snpInfoFN )

    chr2dihhFN = {}

    for chrom in chroms:
        with pr.settingAttrs( 'chrom' ):
      
            chrom_s = 'chr' + str( chrom )
            chromDir = os.path.join( cmsDir, chrom_s )

            xpopScoresFN = os.path.join( chromDir, AddFileSfx( 'max_xpop.tsv', chrom_s, selPop, pops ) )

            pr.addInvokeRule( invokeFn = gatherXPOPscores,
                              invokeArgs = Dict( 'pops chrom selPop sweepDir', outFN = xpopScoresFN ),
                              attrs = dict( pop = allPops, stat = 'max_xpop', piperun_short = True ) )

            ihsFN = getFN_ihs_signif( **Dict( 'sweepDir chrom', pop = selPop ) )
                                      
            ihsScoresFN = os.path.join( chromDir, AddFileSfx( 'iHS.tsv', chrom_s, selPop, pops ) )
            dihhScoresFN = os.path.join( chromDir, AddFileSfx( 'dihh.tsv', chrom_s, selPop, pops ) )

            chr2dihhFN[ chrom ] = dihhScoresFN

            pop2ancFreqFN = os.path.join( cmsDir, chrom_s, AddFileSfx( 'pop2ancFreq.tsv', chrom_s, pops ) )
            pop2sampleSizeFN = os.path.join( cmsDir, chrom_s, AddFileSfx( 'pop2sampleSize.tsv', chrom_s, pops ) )

            pop2snpInfoFN = dict([ ( pop, os.path.join( sweepDir, 'analysis', chrom_s,
                                                        'snps_%(pop)s.tsv' % locals() ) )
                                   for pop in pops ])

            pr.addInvokeRule( invokeFn = gather_snp_info,
                              invokeArgs = Dict( 'pops pop2snpInfoFN pop2ancFreqFN pop2sampleSizeFN' ) )

            pr.addInvokeRule( invokeFn = gather_iHS_scores,
                              invokeArgs = Dict( 'chrom selPop ihsFN pop2ancFreqFN',
#                                                 snpInfoFN = pop2snpInfoFN[ selPop ],
                                                 ihsOutFN = ihsScoresFN, dihhOutFN = dihhScoresFN ),
                              attrs = dict( pop = selPop, stat = ( 'iHS', 'StdDiff' ), piperun_short = True ) )


            freqDiffScoresFN = os.path.join( chromDir, AddFileSfx( 'freqDiff.tsv', chrom_s, selPop, pops ) )
            meanFstScoresFN = os.path.join( chromDir, AddFileSfx( 'meanFst.tsv', chrom_s, selPop, pops ) )

            pr.addInvokeRule( invokeFn = computeMeanFstAndFreqDiffScores,
                              invokeArgs = Dict( 'chrom selPop sweepDir pops pop2ancFreqFN pop2sampleSizeFN',
                                                 outMeanFstFN = meanFstScoresFN,
                                                 outFreqDiffFN = freqDiffScoresFN ),
                              attrs = dict( pop = allPops, stat = ( 'freqDiff', 'meanFst' ), piperun_short = True ) )

            StdDiffScoresFN = os.path.join( chromDir, AddFileSfx( 'StdDiff.tsv', chrom_s, selPop, pops ) )

            rawScoresFN[ chrom ] = dict( iHS = ihsScoresFN, StdDiff = StdDiffScoresFN, meanFst = meanFstScoresFN,
                                         freqDiff = freqDiffScoresFN,
                                         max_xpop = xpopScoresFN )

        # end: with pr.settingAttrs( 'chrom' )
    # end: for chrom in chroms

    #    ihhStdFN = os.path.join( cmsDir, 'dihhstd.tsv' )

    dihhGlobalStdFN = os.path.join( cmsDir, AddFileSfx( 'dihh_global_std.tsv', selPop, pops ) )
    dihhBinMeansFN = os.path.join( cmsDir, AddFileSfx( 'dihh_bin_means.tsv', selPop, pops ) )

    pr.addInvokeRule( invokeFn = normalizeByFreq_getMeanStd_tsv,
                      invokeArgs = dict( iHHDiffFNs = [ chr2dihhFN[k] for k in chroms ],
                                         globalStatFN = dihhGlobalStdFN, binsStatFN = dihhBinMeansFN ),
                      name = 'compute_dihh_meanstd' )
    
    # pr.addInvokeRule( invokeFn = computeMeanStd_binned_tsvs,
    #                   invokeArgs = dict( inFNs = chr2dihhFN.values(), valCol = 'iHHDiff',
    #                                      binCol = 'normingFreqs', binMin = 0.05, binMax = 1.05, binStep = .05,
    #                                      outFN = ihhStdFN ),
    #                   name = 'compute_dihh_std' )

    for chrom in chroms:
        with pr.settingAttrs( 'chrom' ):
            chrom_s = 'chr' + str( chrom )
            chromDir = os.path.join( cmsDir, chrom_s )
            
            StdDiffScoresFN = os.path.join( chromDir, AddFileSfx( 'StdDiff.tsv', chrom_s, selPop, pops ) )
            dbg( 'chrom chr2dihhFN[chrom]' )
            pr.addInvokeRule( invokeFn = normalizeByFreq_compute_normed_tsv,
                              invokeArgs = dict( iHHDiffFN = chr2dihhFN[ chrom ],
                                                 globalStatFN = dihhGlobalStdFN,
                                                 binsStatFN = dihhBinMeansFN,
                                                 StdDiffFN = StdDiffScoresFN ) )

    statFNs = {}
    statLikesRatioFNs = {}

    for stat in CMSBins.CMSstats:
        with pr.settingAttrs( stat = stat, pop = ( selPop, ) if stat in ( 'iHS', 'StdDiff' ) else allPops, piperun_short = True ):
            if stat not in CMSBins.nonNormedStats:
                rawFNs = [ rawScoresFN[ chrom ][ stat ] for chrom in chroms ]
                meanStdFN = os.path.join( cmsDir, AddFileSfx( 'meanStd.tsv', stat, selPop, pops ) )

                # DefineRulesTo_computeMeanStd( pr, inFNs = rawFNs, colNum = 1,
                #                               outFN = meanStdFN,
                #                               addRuleArgs = \
                #                               dict( name = 'computeMeanStd_for_stat',
                #                                     attrs = dict( chrom = chroms ) ) )

#                meanStdBzFN = os.path.join( cmsDir, stat + '_meanStdForStat.tsv' )
                pr.addInvokeRule( invokeFn = computeMeanStd,
                                  invokeArgs = dict( inFNs = rawFNs, colName = stat, outFN = meanStdFN ) )
                
            # end: if stat not in CMSBins.nonNormedStats

            for chrom in chroms:
                with pr.settingAttrs( 'chrom' ):
                    statFN = rawScoresFN[ chrom ][ stat ]

                    if stat not in CMSBins.nonNormedStats:
                        normedFN = AddFileSfx( statFN, 'normed' )
                        
                        DefineRulesTo_normalizeOneColumn( pr, inFN = statFN,
                                                          meanStdFN = meanStdFN,
                                                          colName = stat,
                                                          outFN = normedFN,
                                                          addRuleArgs = dict( attrs = Dict( 'chrom' ) ) )
                        statFN = normedFN
                        
                    bins_beg = CMSBins.stat_start[ stat ]
                    bins_end = CMSBins.stat_end[ stat ]
                    bins_n = CMSBins.stat_nbin[ stat ]

                    statFNs[ ( chrom, stat ) ] = statFN

                    statLikesRatioFN = AddFileSfx( rawScoresFN[ chrom ][ stat ], 'likesRatio' )
                    statLikesRatioFNs[ ( chrom, stat ) ] = statLikesRatioFN
                    
                    pr.addInvokeRule( invokeFn = computeLikeRatioForStat,
                                      invokeArgs = dict( stat = stat,
                                                         statValsFN = statFN,
                                                         hitLikesFN = '../Data/Common_Data/sim/likes/hitsLikes_toneutFixed_1.tsv',
                                                         missLikesFN = '../Data/Common_Data/sim/likes/missLikes_toneutFixed_1.tsv',
                                                         stat_start = bins_beg,
                                                         stat_end = bins_end,
                                                         stat_nbin = bins_n,
                                                         statLikesRatioFN = statLikesRatioFN ) )

                                      
                # end: with pr.settingAttrs( 'chrom' )
            # end: for chrom in chroms
        # end: with pr.settingAttrs( stat = stat, piperun_short = True )
    # end: for stat in CMSBins.CMSstats

    for chrom in chroms:
        with pr.settingAttrs( chrom = chrom, stat = CMSBins.CMSstats ):
            chrom_s = 'chr' + str( chrom )
            chromDir = os.path.join( cmsDir, chrom_s )
            
            likesRatioFN = os.path.join( chromDir, AddFileSfx( 'likesRatio.tsv', CMSBins.CMSstats, selPop, pops ) )
            pr.addInvokeRule( invokeFn = addLikesRatios,
                              invokeArgs = dict( inFNs = [ statLikesRatioFNs[ ( chrom, stat ) ] for stat in CMSBins.CMSstats ],
                                                 colNames = [ colName + 'likeRatio' for colName in CMSBins.CMSstats ],
                                                 outFN = likesRatioFN ) )

            joinStatsFN = os.path.join( chromDir, AddFileSfx( 'joinStats.tsv', CMSBins.CMSstats, selPop, pops ) )
            snpInfoFN = os.path.join( sweepDir, 'analysis/chr%(chrom)s/snps_%(selPop)s.tsv' % locals() )
            pr.addInvokeRule( invokeFn = joinStats,
                              invokeArgs = dict( snpInfoFN = snpInfoFN,

                                                 statLikesFNs = [ statLikesRatioFNs[ ( chrom, stat ) ] for stat in CMSBins.CMSstats ],
                                                 likesRatioFN = likesRatioFN,
                                                 outFN = joinStatsFN ),
                              attrs = dict( stat = CMSBins.CMSstats, chrom = chrom ) )
            

# end: def DefineRulesTo_fastCMS( pr, pops, chroms, selPop, popSampleSizesFN, sweepDir, cmsDir )


if __name__ == '__main__':
#    gatherXPOPscores( pops = ( 'ASN', 'BEB', 'YRI' ), selPop = 'CEU', chrom = 22, sweepDir = '/idi/sabeti-scratch/ilya/gsvn/Data/Elinor_Data/BNG_CMS_1KG_P3_gw/sweepanalysis', outFN = 'cmp.tsv' )

    computeMeanFstAndFreqDiffScores( pops = ( 'CEU', 'ASN', 'BEB', 'YRI' ), selPop = 'CEU', chrom = 22, sweepDir = '/idi/sabeti-scratch/ilya/gsvn/Data/Elinor_Data/BNG_CMS_1KG_P3_gw/sweepanalysis', outMeanFstFN = 'cmp.tsv', outFreqDiffFN = 'cmp2.tsv' )

    

    

    

    