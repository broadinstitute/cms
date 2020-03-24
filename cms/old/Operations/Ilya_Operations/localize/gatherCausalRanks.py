from Classes.DotData import DotData
from Operations.Shari_Operations.localize.Scenario import GetSelectionScenarios
from Operations.Shari_Operations.localize.PopConsts import AllPops, AllFreqs, AllAges
from Operations.IDotData import IDotData
from Operations.MiscUtil import Dict, AddFileSfx, dbg, MakeSeq
import os, sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as pp

def gatherCausalRanks( Ddata = None, scenario = None,
                       selpos = 500000, thinSfx = '', complikeSfx = None, likesTableSfx = '',
                       nonNanStats = 'ALL',
                       cmsFileFN = None,
                       causalRankFN = None,
                       getio = None ):
    """For each replica in one scenario, get the rank of the causal SNP by CMS score, and save as a replica statistic.
    """

    assert cmsFileFN or ( Ddata and scenario )
    scenDir = scenario.scenDir() if scenario else 'unknown_scenDir'
    if not cmsFileFN: snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenDir )
    if not cmsFileFN: cmsFileFN = os.path.join( snpStatsDir, AddFileSfx( 'complike.data/', scenario.mutPop, complikeSfx, likesTableSfx ) )

    if not causalRankFN: causalRankFN = os.path.join( Ddata, 'replicastats',
                                                      scenDir, AddFileSfx( 'causalRank.tsv', complikeSfx,
                                                                           'nonNan', *MakeSeq( nonNanStats ) ) )

    if getio: return dict( depends_on = cmsFileFN, creates = causalRankFN,
                           mediumRuleNameSfx = ( scenDir, complikeSfx ) )

    cmsScores = IDotData( cmsFileFN )
    if nonNanStats.upper() == 'ALL': nonNanStats = cmsScores.headings

    with IDotData.openForWrite( causalRankFN, headings = 'replicaNum causalRank causalScore' ) as causalRankFile:
        for replicaNum, cmsScores1, cmsScores2 in cmsScores.groupby( 'Chrom', multiPass = 2 ):
            for r1 in cmsScores1:
                if r1.Pos == selpos:
                    causalScore = r1.complike
                    numHigher = 0
                    for r2 in cmsScores2:
                        if r2.complike > causalScore and all([ np.isfinite( r2[ c ] ) for c in nonNanStats ]):
                            numHigher += 1
                    causalRankFile.writeRecord( int( replicaNum ), numHigher, causalScore )

                if r1.Pos >= selpos: break
    
def gatherCausalStats( Ddata, scenario,
                       selpos = 500000, thinSfx = '', complikeSfx = None, likesTableSfx = '',
                       nonNanStats = 'ALL',
                       getio = None ):
    """For each replica in one scenario, gather causal SNP stats and save them as the replica statistic.
    """

    snpStatsDir = os.path.join( Ddata, 'snpStats'+ thinSfx, scenario.scenDir() )
    replicaStatsDir = os.path.join( Ddata, 'replicastats'+ thinSfx, scenario.scenDir() )
    complikeFN = os.path.join( snpStatsDir,
                               AddFileSfx( 'complike.data/', 'normedLocal', scenario.mutPop, complikeSfx,
                                           'nonNan', *MakeSeq( nonNanStats ) ) )
    
    causalStatsFN = os.path.join( replicaStatsDir, AddFileSfx( 'causalStats.tsv', complikeSfx, likesTableSfx,
                                                               'nonNan', *MakeSeq( nonNanStats ) ) )

    if getio: return dict( depends_on = complikeFN, creates = causalStatsFN,
                           mediumRuleNameSfx = ( scenario.scenDir(), complikeSfx, likesTableSfx ) )

                                                  
    complikeFile = IDotData( complikeFN )
    complikeFile[ complikeFile.Pos == selpos ].addComputedCols( newColNames = 'replicaNum',
                                                                newColFn = lambda r: int( r.Chrom ) ).save( causalStatsFN )
    
def DefineRulesTo_gatherCausalRanks( pr, Ddata, mutPops = AllPops, mutFreqs = AllFreqs, mutAges = AllAges,
                                     selpos = 500000, complikeSfx = None, nonNanStats = 'ALL' ):

    for scenario in GetSelectionScenarios( mutAges = mutAges, mutPops = mutPops, mutFreqs = mutFreqs ):
        pr.addInvokeRule( invokeFn = gatherCausalRanks,
                          invokeArgs = Dict ( 'Ddata scenario selpos complikeSfx nonNanStats' ) )
        pr.addInvokeRule( invokeFn = gatherCausalStats,
                          invokeArgs = Dict ( 'Ddata scenario selpos complikeSfx nonNanStats' ) )
        
        

        
    
            
