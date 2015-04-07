from __future__ import with_statement, division

__all__ = ( 'DefineRulesTo_normalizeColumnsWithinGroups', )

import operator, math, logging, sys, os, itertools, functools
import numpy as np
import matplotlib
import matplotlib.pyplot as pp
import pandas as pd
from Operations.MiscUtil import Dict, MakeSeq, AddFileSfx, AddFileSubdir, Sfx, flatten, dbg, SlurpFile, DumpFile, Str, IsSeq, dbg, \
    Histogrammer, ReplaceFileExt
from Operations.IDotData import IDotData
from Operations.Ilya_Operations.SnpStats import GraphHistograms
from Classes.DotData import DotData
from Operations.Ilya_Operations.PipeRun.python.PipeRun import GetCreates


def computeSumsWithinGroups( inFN, cols, groupCols, groupsAreContiguous = True, outFN = None, getio = None ):
  """For a tsv file, compute sums, sumsquares and counts for each of the given columns within groups
  defined by groupCols.

  >>> z = IDotData( names = ( 'a', 'b' ), Records = ( ( 1, 2 ), ( 1, 3 ), ( 2, 4 ), ( 2, 5 ) ) )
  >>> computeSumsWithinGroups( inFN = z, cols = 'b', groupCols = 'a', outFN = sys.stdout )
  ... # doctest: +NORMALIZE_WHITESPACE
  a	b_count	b_sum	b_sumSq	b_numNaN
  1	2	5.0	13.0	0
  2	2	9.0	41.0	0

  """

  cols = tuple( MakeSeq( cols ) )
  groupCols = tuple( MakeSeq( groupCols ) )
  if outFN is None: outFN = AddFileSubdir( 'stats', AddFileSfx( inFN, 'sums', *( cols + groupCols ) ) )

  def combiner( inFNs, outFN ): IDotData.mergeColumnSummaries( iDotDatas = inFNs, cols = cols, groupCols = groupCols ).save( outFN )
  
  if getio: return dict( depends_on = inFN, creates = outFN,
                         splitByCols = { inFN: dict( keyCols = () ) },
                         combiner = { outFN: combiner }  )

  IDotData( inFN ).summarizeColumnsWithinGroups( **Dict( 'cols groupCols groupsAreContiguous' ) ).save( outFN )

def computeMeanStdWithinGroups( inFN, cols, groupCols, groupsAreContiguous = True, outFN = None, getio = None ):
    """Add columns representing mean and std within each group.
    """

    sumsFN = GetCreates( computeSumsWithinGroups, **Dict( 'inFN cols groupCols groupsAreContiguous' ) )[0]
    if outFN is None: outFN = AddFileSubdir( 'stats', AddFileSfx( inFN, 'meanStd', *( cols + groupCols ) ) )
    if getio: return dict( depends_on = sumsFN, creates = outFN, attrs = dict( piperun_short = True ) )

    return IDotData( sumsFN ).addMeanStdCols( cols = cols ).save( outFN )
  
def normalizeColumnsWithinGroups( inFN, cols, groupCols, outFN, groupsAreContiguous = True, getio = None ):
    """Normalize the specified columns of a table within groups.

    Params:

       inFN - the input table
       cols - the columns to be normalized
       groupCols - the columns that define the groups: rows that have the same combination of values
          in the group columns, are in the same group.
       outFN - the output table
       groupsAreContiguous - if True, rows belonging to the same group must be contiguous in the table;
          if False, no such assumption is made.
    """

    cols = tuple( MakeSeq( cols ) )
    groupCols = tuple( MakeSeq( groupCols ) )

    meansFN = GetCreates( computeMeanStdWithinGroups, **Dict( 'inFN cols groupCols groupsAreContiguous' ) )[0]
    
    if getio: return dict( depends_on = ( inFN, meansFN ), creates = outFN,
                           splitByCols = { inFN: dict( keyCols = () ) } )

    inFile = IDotData( inFN )
    means = IDotData( meansFN )

    inFile.normalizeColumnsWithinGroups_using_means( **Dict( 'cols groupCols groupsAreContiguous means' ) ).save( outFN )
        
#    IDotData( inFN ).mapRecords( lambda r: [ r[ c ] if c not in cols else   for c in r.headings      ]      )


def DefineRulesTo_meanStdWithinGroups( pr, inFN, cols, groupCols, groupsAreContiguous = True, nameSfx = '' ):
  """Adds rules to create a version of a table with given columns normalized within groups."""

  pr.addInvokeRule( invokeFn = computeSumsWithinGroups, invokeArgs = Dict( 'inFN cols groupCols groupsAreContiguous' ),
                    name = 'computeSumsWithinGroups' + Sfx( nameSfx ) )
  pr.addInvokeRule( invokeFn = computeMeanStdWithinGroups, invokeArgs = Dict( 'inFN cols groupCols groupsAreContiguous' ),
                    name = 'computeMeanStdWithinGroups' + Sfx( nameSfx ) )
    
def DefineRulesTo_normalizeColumnsWithinGroups( pr, inFN, cols, groupCols, groupsAreContiguous = True, nameSfx = '', outFN = None ):
  """Adds rules to create a version of a table with given columns normalized within groups."""

  cols = tuple( MakeSeq( cols ) )
  groupCols = tuple( MakeSeq( groupCols ) )
  
  DefineRulesTo_meanStdWithinGroups( **Dict( 'pr inFN cols groupCols groupsAreContiguous nameSfx' ) )
  pr.addInvokeRule( invokeFn = normalizeColumnsWithinGroups,
                    invokeArgs = Dict( 'inFN cols groupCols outFN groupsAreContiguous' ),
                    name = 'normalizeColumnsWithinGroups' + Sfx( nameSfx ) )
                    
def recordCount( inFN, outFN = None, getio = None ):
  """Count the lines in the given file"""

  if outFN is None: outFN = AddFileSfx( inFN, 'recordCount' )

  def SaveCount( count, outFN ):
    IDotData( names = ( 'recordCount', ), Records = ( count, ) ).save( outFN )
  def combiner( inFNs, outFN ):
    SaveCount( np.sum([ next( iter( IDotData( f ) ) ) for f in inFNs ]), outFN )

  if getio: return dict( depends_on = inFN, creates = outFN,
                         splitByCols = { inFN: dict( keyCols = () ) },
                         combiner = { outFN: combiner } )

  contents = SlurpFile( inFN ).strip()
  SaveCount( contents.count( '\n' ), outFN )
  
  
def checkTableKey( inFN, cols, comparison = 'lt', writeCheckedFile = True,
                   tsvOpts = {}, lineFilter = None, lineFilterCols = (), getio = None ):
  """Check that in the given table, record identifiers increase uniformly.

  Params:

     cols - the columns whose tuple should uniformly inrease
     comparison - this comparison must be true between each record and the next.
       the comparison is the name of a routine in the operator module.
  """

  cols = tuple( MakeSeq( cols ) )
  lineFilterCols = tuple( MakeSeq( lineFilterCols ) )
  checkedFN = Str( '$inFN.checked_${comparison}' ) + Sfx( *cols )
  if getio: return dict( depends_on = inFN, creates = checkedFN if writeCheckedFile else (),
                         attrs = dict( piperun_short = True ) )

  comparisonFunc = getattr( operator, comparison )
  prevRec = None
  loadCols = cols + lineFilterCols

  nskipped = 0
  nchecked = 0
  for i, r in enumerate( IDotData( inFN, ToLoad = loadCols, **tsvOpts ) ):
    if lineFilter and not lineFilter( r ):
      nskipped += 1
      continue
    
    thisRec = r[ cols ] if IsSeq( r ) else ( r, )
    if i > 0 and not comparisonFunc( prevRec, thisRec ):
      logging.error( Str( 'at line $i of $inFN, looking at $cols: $prevRec is not $comparison $thisRec' ) )
      assert False
    else: nchecked += 1
    prevRec = thisRec

  dbg( 'nchecked nskipped' )
  DumpFile( checkedFN, 'checked ok.' )
  
def computeHistograms( inFN, cols, binSizes = None, outFNs = None, getio = None ):
  """Compute histograms of the specified columns of the input"""

  cols = tuple( MakeSeq( cols ) )
  binSizesHere = ( .001, ) * len( cols ) if binSizes is None else tuple( MakeSeq( binSizes ) )
  outFNsHere = outFNs
  if outFNsHere is None: outFNsHere = [ AddFileSubdir( 'stats', AddFileSfx( inFN, 'hist', col ) ) for col in cols ]

  assert len( cols ) == len( binSizesHere ) == len( outFNsHere )
  if getio: return dict( depends_on = inFN, creates = outFNsHere ) 
  # add histogram combiner

  hists = [ Histogrammer( binSize = binSize ) for binSize in binSizesHere ]
  z = IDotData( inFN )
  for h, c, outFN in zip( hists, cols, outFNsHere ):
    h.addVals( z[ c ] )
    h.save( outFN )

def plotHistograms( inFN, cols, outFNs = None, getio = None, **kwargs ):

  histFNs = GetCreates( computeHistograms, **Dict( 'inFN cols outFNs' ) )
  outFN = AddFileSfx( ReplaceFileExt( inFN, '.svg' ), 'hist' )

  if getio: return dict( depends_on = histFNs, creates = outFN )

  GraphHistograms( histFiles = histFNs, outFile = outFN, labels = tuple( MakeSeq( cols ) ),
                   **kwargs )
    

def DotData2TSV( inFN, outFN, readOpts = {}, getio = None ):
  """Convert DotData to TSV"""
  if getio: return dict( depends_on = inFN, creates = outFN )
  DotData( Path = inFN, **readOpts ).saveToSV( outFN )


def DotData2TSV_lowmem( inFN, outFN, readOpts = {}, getio = None ):
  """Convert DotData to TSV"""
  if getio: return dict( depends_on = inFN, creates = outFN )
  IDotData( Path = inFN, **readOpts ).saveToSV( outFN )

def TSV2DotData( inFN, outFN, getio = None ):
  """Convert a TSV file to a DotData dir."""
  if getio: return dict( depends_on = inFN, creates = outFN )
  DotData( SVPath = inFN ).save( outFN )

def tsv2npy( inFN, outFN = None, getio = None ):
  """Convert a TSV file to an .npy file."""
  if outFN is None: outFN = ReplaceFileExt( inFN, '.npy' )
  if getio: return dict( depends_on = inFN, creates = outFN )
  z = DotData( SVPath = inFN )
  np.save( outFN, z )
  
def tsv2npz( inFN, outFN = None, arrayName = None, dotDataArgs = {}, getio = None ):
  """Convert a TSV file to an .npy file."""
  if outFN is None: outFN = ReplaceFileExt( inFN, '.npz' )
  if getio: return dict( depends_on = inFN, creates = outFN, attrs = dict( piperun_short = True ) )
  z = DotData( SVPath = inFN, **dotDataArgs )
  if arrayName is None:
    np.savez_compressed( outFN, z )
  else:
    np.savez_compressed( outFN, **{ arrayName : z } )
  

def sortTableOn( inFN, outFN, keyCols, reverse = False, getio = None ):
  """Sort the given table on the given column(s)."""

  if getio: return dict( depends_on = inFN, creates = outFN )

  result = IDotData( inFN ).sortedOn( *MakeSeq( keyCols ) )
  if reverse:
    d = result.toDotData()
    result = d[ range( len( d )-1, -1, -1 ) ]
  result.save( outFN )


def computeMeanStd_binned_old( inDatas, valCol, binCol, binMin, binMax, binCount ):
  """Compute binned stats for a set of tables"""

  sums = np.zeros( binCount )
  sumsSq = np.zeros_like( sums )
  counts = np.zeros_like( sums )
  bins = np.linspace( binMin, binMax, binCount+1 )
  binSize = ( binMax - binMin ) / binCount
  for d_idx, d in enumerate( inDatas ):
    dbg( 'd_idx d binSize' )
    dbg( 'd[binCol]' )

    for i in range( binCount ):
        binBot = bins[i]
        binTop = bins[i+1]
        dbg( 'binBot binTop' )
#        theIdx = ( (binTop - d[ binCol ]) < binSize ) & ( ( binTop - d[ binCol ] ) > 0 )
        theIdx = ( binBot < d[ binCol ].values ) & ( d[ binCol ].values <= binTop )
        dbg( 'binBot binTop' )
        DotData( names = ('rows',), Columns = theIdx.nonzero() ).saveToSV( 'nz%02d.tsv' % i )
        #rowsStr = ','.join(map(str,list(theIdx.nonzero())))
        #print 'binnedRows=', rowsStr
        hereVals = d[ theIdx ][ valCol ]
        DotData( names = ( 'temp', ), Columns = ( hereVals, ) ).saveToSV( 'temp2%2d.tsv' % i )
        
        dbg( '"BEF" theIdx.sum() i bins[i] bins[i+1] len(hereVals)' )
        counts[i] += len( hereVals )
        sums[i] += np.sum( hereVals )
        sumsSq[i] += np.sum( hereVals * hereVals )
        dbg( '"AFT" i bins[i] bins[i+1] len(hereVals)' )

    if False:
        # fast version
        binsHere = np.digitize( d[ binCol ], bins ) - 1
        dbg( 'len(binsHere) binsHere' )
        np.clip( binsHere, 0, binCount-1, out = binsHere );
        dbg( 'binsHere' )

        counts += np.bincount( binsHere, minlength = binCount )
        sums += np.bincount( binsHere, weights = d[ valCol ], minlength = binCount )
        sumsSq += np.bincount( binsHere, weights = d[ valCol ] * d[ valCol ], minlength = binCount )

  countsOrig = counts.astype( int )
  counts[ counts == 0 ] = np.nan
  means = sums / counts
  stds = sumsSq / counts - means * means

  return pd.DataFrame( dict( binBeg = bins[:-1],
                             binEnd = bins[1:],
                             counts = countsOrig, sums = sums, sumsSq = sumsSq,
                             means = means, stds = stds ) )

# end: def computeMeanStd_binned( inDatas, valCol, binCol, binMin, binMax, binCount ):

def computeMeanStd_binned( inDatas, valCol, binCol, binMin, binMax, binStep ):
  """Compute binned stats for a set of tables"""

  binCount = int( ( binMax - binMin ) / binStep )
  dbg( 'binCount' )
  sums = np.zeros( binCount )
  sumsSq = np.zeros_like( sums )
  counts = np.zeros_like( sums )
  bins = np.arange( binMin, binMax, binStep )
  for d_idx, d in enumerate( inDatas ):
    dbg( 'd_idx d binStep' )
    dbg( 'd[binCol]' )

    binColValues = 1.0 - ( 1.0 - d[ binCol ].values )

    for i in range( binCount ):
#        binBot = bins[i]
        binTop = bins[i]
        theIdx = ( (binTop - binColValues) < binStep ) & ( ( binTop - binColValues ) > 0 )
#        theIdx = ( binBot < d[ binCol ].values ) & ( d[ binCol ].values <= binTop )
 #       DotData( names = ('rows',), Columns = theIdx.nonzero() ).saveToSV( 'nz%02d.tsv' % i )
        #rowsStr = ','.join(map(str,list(theIdx.nonzero())))
        #print 'binnedRows=', rowsStr
        hereVals = d[ theIdx ][ valCol ]
#        DotData( names = ( 'temp', ), Columns = ( hereVals, ) ).saveToSV( 'temp2%2d.tsv' % i )
        
        dbg( '"BEF" theIdx.sum() i bins[i] len(hereVals)' )
        counts[i] += len( hereVals )
        sums[i] += np.sum( hereVals )
        sumsSq[i] += np.sum( hereVals * hereVals )
#        dbg( '"AFT" i bins[i] bins[i+1] len(hereVals)' )

    if False:
        # fast version
        binsHere = np.digitize( d[ binCol ], bins ) - 1
        dbg( 'len(binsHere) binsHere' )
        np.clip( binsHere, 0, binCount-1, out = binsHere );
        dbg( 'binsHere' )

        counts += np.bincount( binsHere, minlength = binCount )
        sums += np.bincount( binsHere, weights = d[ valCol ], minlength = binCount )
        sumsSq += np.bincount( binsHere, weights = d[ valCol ] * d[ valCol ], minlength = binCount )

  countsOrig = counts.astype( int )
  counts[ counts == 0 ] = np.nan
  means = sums / counts
  stds = sumsSq / counts - means * means

  return pd.DataFrame( dict( binBeg = bins - binStep,
                             binEnd = bins,
                             counts = countsOrig, sums = sums, sumsSq = sumsSq,
                             means = means, stds = stds ) )

# end: def computeMeanStd_binned( inDatas, valCol, binCol, binMin, binMax, binStep ):
  

def computeMeanStd_binned_tsvs( inFNs, valCol, binCol, binMin, binMax, binStep, outFN, getio = None ):
  """Compute binned stats for a set of tables"""

  if getio: return dict( depends_on = inFNs, creates = outFN,
                         uses = computeMeanStd_binned )

  computeMeanStd_binned( inDatas = itertools.imap( lambda f: pd.read_table( f, usecols = ( valCol, binCol ) ).dropna(),
                                                   MakeSeq( inFNs ) ),
                         **Dict( 'valCol binCol binMin binMax binStep' ) ).to_csv( outFN, sep = '\t',
                                                                                   index_label = 'binId',
                                                                                   na_rep = 'NaN' )

# end: def computeMeanStd_binned_tsvs( inFNs, valCol, binCol, binMin, binMax, binCount, outFN, getio = None )

def normalizeInBins( inData, valCol, binCol, binMin, binMax, binStep, binMeans, commonStd ):
    """Normalize data within bins, using previously computed bin means"""

    binColValues = 1.0 - ( 1.0 - inData[ binCol ].values )
    binCount = int( ( binMax - binMin ) / binStep )
    bins = np.arange( binMin, binMax, binStep )

    means = np.zeros( len( inData ) )

    for i in range( binCount ):
#        binBot = bins[i]
        binTop = bins[i]
        theIdx = ( (binTop - binColValues) < binStep ) & ( ( binTop - binColValues ) >= 0 )
        means[ theIdx ] = binMeans[ i ]

    result = ( inData[ valCol ].values - means ) / commonStd
    
    if False:
        # Fast version
        bins = np.linspace( binMin, binMax, binCount+1 )
        binsHere = np.digitize( inData[ binCol ], bins ) - 1
        np.clip( binsHere, 0, binCount-1, out = binsHere );
        means = np.take( binMeans, binsHere )
        result = ( inData[ valCol ].values - means ) / commonStd
        
    return result

def normalizeInBins_tsv( inDataFN, valCol, binCol, binMin, binMax, binStep, binsFN, outFN,
                         normedCol,
                         getio = None):
    """Normalize data within bins, using previously computed bin means"""

    if getio: return dict( depends_on = ( inDataFN, binsFN ), creates = outFN, uses = normalizeInBins )

    inData = pd.read_table( inDataFN )
    binStats = pd.read_table( binsFN )
    binMeans = binStats.means
    totCount = float( binStats.counts.sum() )
    totMean = binStats.sums.sum() / totCount
    commonStd = np.sqrt( binStats.sumsSq.sum() / totCount - totMean * totMean )
    dbg( '"CCCCCCCC" commonStd binMeans totCount totMean binStats.sums.sum() binStats.sumsSq.sum()' )
    normed = normalizeInBins( **Dict( 'inData valCol binCol binMin binMax binStep binMeans commonStd' ) )
    inData.insert( len( inData.columns ), normedCol, normed )
    inData.to_csv( outFN, sep = '\t', na_rep = 'NaN', index = False )

def DefineRulesTo_computeMeanStd( pr, inFNs, colNum, outFN, addRuleArgs = {} ):
  """Define rules to compute mean and stddev for a given column in the given tsv files"""

  pr.addRule( commands = ' | '.join(( 'tail -q -n +2 ' + ' '.join( MakeSeq( inFNs ) ),
                                      'cut -f %d' % colNum,
                                      'grep -iv nan',
                                      '../Operations/Ilya_Operations/tblstats' )),
              depends_on = inFNs,
              saveOutputTo = outFN,
              **addRuleArgs )

def normalizeOneColumn( inFN, colName, meanStdFN, outFN, getio = None ):
    """Normalize one column"""

    if getio: return dict( depends_on = ( inFN, meanStdFN ), creates = outFN,
                           attrs = dict( piperun_short = True ) )

    import pandas as pd

    f = pd.read_table( inFN )
    meanStd = pd.read_table( meanStdFN, index_col = 0 )
    f[colName] = ( f[colName] - meanStd.loc['mean','val'] ) / meanStd.loc['std','val']
    f.to_csv( outFN, sep = '\t', header = True, index = False, na_rep = 'NaN' )

def DefineRulesTo_normalizeOneColumn( pr, inFN, colName, meanStdFN, outFN, addRuleArgs = {} ):
    """Define rules to normalize one column"""

    pr.addInvokeRule( invokeFn = normalizeOneColumn,
                      invokeArgs = Dict( 'inFN colName meanStdFN outFN' ),
                      **addRuleArgs )

    
#if __name__ == '__main__':
    
  
  
