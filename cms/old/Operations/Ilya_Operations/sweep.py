
from Operations.MiscUtil import ReserveTmpFileName, SystemSucceed, Str, Dict
from Operations.Shari_Operations.localize.Scenario import GetScenarios
from Operations.Shari_Operations.localize.PopConsts import AllAges, AllPops, AllFreqs, popName
import os

def DefineRulesTo_RunSweep( pr, phasedFilesDir, geneticMapsDir, ancestralDir, outDir, pops, chroms,
                            backgroundSweepDir = None, backgroundPopMap = None, backgroundPopOrder = None,
                            tests = None, runImportsLocally = True,
                            add_iHS_opts = '', add_xpop_opts = '',
                            genomeBuild = None, genMapSfx = None, ihsSfx = '',
                            geneticMapsDir2 = None, geneticMapsName2 = None,
                            parallelizeIHS = True, parallelizeXPOP = True,
                            addLegendHeadings = False, lrh_memreq = None ):
    """Define rules to run Sweep analysis on data in HapMap 2 format."""

    # This routine is just a wrapper for running the hapmap_analysis_pipe.pl script, which
    # does the actual work of defining the pipeline.

    popList = ','.join( pops )
    chromList = ','.join( map( str, chroms ) )

    sweepPipelineFile = ReserveTmpFileName( prefix = 'tmp_sweep_', suffix = '.xml' )

    moreOpts = ''
    if backgroundSweepDir: moreOpts += ' --background-sweep-dir ' + backgroundSweepDir
    if backgroundPopMap: moreOpts += ' --background-pop-map ' + ','.join( k + ':' + v for k, v in list(backgroundPopMap.items()) )
    if backgroundPopOrder: moreOpts += ' --background-pop-order ' + ','.join( backgroundPopOrder )

    if tests: moreOpts += ' --tests ' + ','.join( tests )
    if runImportsLocally: moreOpts += ' --run-imports-locally'
    if add_iHS_opts: moreOpts += ' --add-ihs-opts "' + add_iHS_opts + '"'
    if add_xpop_opts: moreOpts += ' --add-xpop-opts "' + add_xpop_opts + '"'
    if genomeBuild: moreOpts += ' --genome-build ' + genomeBuild
    if lrh_memreq is not None: moreOpts += ' --lrh-memreq ' + str( lrh_memreq ) + ' '
    if addLegendHeadings: moreOpts += ' --add-legend-headings '
    if genMapSfx is None:
        if genomeBuild == 'hg18': genMapSfx = '_b36'
        elif genomeBuild == 'hg19': genMapSfx = '_b37'
        else: genMapSfx = ''
    if genMapSfx: moreOpts += ' --genMapSfx ' + genMapSfx
    if ihsSfx: moreOpts += ' --ihsSfx ' + ihsSfx
    if parallelizeIHS: moreOpts += ' --parallelizeIHS'
    if parallelizeXPOP: moreOpts += ' --parallelizeXPOP'
    if pr.pythonCmd != 'python': moreOpts += ' --pythonCmd ' + pr.pythonCmd

    if geneticMapsDir2 and geneticMapsName2: moreOpts += ' --recomb-root2 ' + geneticMapsDir2 + ' --recomb-name2 ' + geneticMapsName2
    
    SystemSucceed( Str( '../Other/Ilya_Other/sweep/scripts/hapmap_analysis_pipe.pl '
                        ' --only-write-pipeline $sweepPipelineFile --phased-root $phasedFilesDir '
                        ' --chimp-root $ancestralDir '
                        ' --recomb-root $geneticMapsDir --out-root $outDir --pops $popList '
                        ' --chroms $chromList $moreOpts' ) )

    pr.addPipelineFromFile( sweepPipelineFile )


def DefineRulesTo_extractGeneticMapsFromSims( pr, Ddata, mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs ):
    """Define rules to extract genetic map from Sweep import of cosi simulations"""

    for scen in GetScenarios( **Dict( 'mutAges mutPops mutFreqs' ) ):


        projFile = os.path.join( Ddata, 'data', scen.scenDir(), 'project' )
        
        thePop = scen.scenName() + ':' + popName[ mutPops[0] ]

        outFile = os.path.join( Ddata, 'snpStats', scen.scenDir(), 'gdMap.tsv' )
        
        pr.addRule( commands = '../Other/Ilya_Other/sweep/scripts/run-sweep ExtractGeneticMap '
                    '$projFile $thePop $outFile',
                    depends_on = projFile,
                    creates = outFile,
                    name = 'ExtractGeneticMap', mediumRuleNameSfx = scen.scenDir(),
                    comment = 'Extract genetic map for all replicas in one simulation scenario' )

            
    

        

        
    
    
