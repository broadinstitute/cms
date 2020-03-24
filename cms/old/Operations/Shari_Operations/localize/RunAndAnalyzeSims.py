
__all__ = ( 'DefineRulesTo_RunSimsOnly', 'DefineRulesTo_RunSweepOnSims', 'DefineRulesTo_RunSimsAndSweep' )

from Operations.Ilya_Operations.sim.sfs.working.pardis2.RunSimulations import DefineRulesTo_CreateSimulationParams, \
    DefineRulesTo_runSims, WriteSimulationInfo
from Operations.MiscUtil import SystemSucceed, Str, Dict, MergeDicts, WaitForFileToAppear, SysTypeString, dbg, AddFileSfx, Sfx
from Operations.Ilya_Operations.PipeRun.python.PipeRun import PipeRun
from Operations.Shari_Operations.localize.Scenario import GetScenarios
from Operations.Shari_Operations.localize.PopConsts import AllAges, AllPops, AllFreqs, pop2name
import os, operator
from functools import reduce


def DefineRulesTo_RunSimsOnly( pr, mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs, nreplicas = 100,
                               Ddata = '../Data/Shari_Data_3/sim', simsOut = 'simsOut5', suffix = '',
                               inputParamsFiles = '../Data/Ilya_Data/sim2/neutralParams.txt',
                               DdataSeeds = '', useGenMap = None, scen2alternateNeutralParams = {},
                               withGeneConvBug = False, withNewCosi = False, withCosi = None, DdataMimic = None ):
    """Define rules for running simulations.

    Parameters:

       mutAges, mutPops, mutFreqs - parameters defining the selection scenario.
    
    """

    DefineRulesTo_CreateSimulationParams( **Dict( 'pr mutAges mutPops mutFreqs Ddata suffix inputParamsFiles '
                                                  'scen2alternateNeutralParams' ) )
    DefineRulesTo_runSims( **Dict( 'pr mutAges mutPops mutFreqs nreplicas Ddata simsOut suffix DdataSeeds useGenMap withGeneConvBug '
                                   'withNewCosi withCosi DdataMimic' ) )

    
def DefineRulesTo_DoThinning( pr, mutAges, mutPops, mutFreqs, nreplicas, Ddata, simsOut, suffix, thinning, thinExt ):
    """Thin the set of SNPs by removing some fraction of SNPs, to simulate the SNP ascertainment done by
    HapMap.

    Parameters:

       thinning -- either blank, '_hm2thin' or '_hm3thin'.

    """

    assert not thinning or thinning in ( '_hm2thin', '_hm3thin' )

    if thinning:

        thinProbFile = '../Data/Ilya_Data/sim/sfs/working/pardis2/thinprob2' + \
            ( '_hapmap3' if thinning == '_hm3thin' else '' ) + '.txt'

        for scen in GetScenarios( mutAges, mutPops, mutFreqs ):
            for replicaNum in range( nreplicas ):

                
                scenName = 'neutral' if scen.is_neutral() else 'sel%d_%d' % ( scen.mutFreq, scen.mutPop )
                scenFileName = '%d_%s' % ( replicaNum, scenName )
                scenSubdir = ('' if scen.is_neutral() else '%dky/' % scen.mutAge) + scenName

                unthinnedSimsOut = Ddata + '/' + simsOut + suffix
                thinnedSimsOut = unthinnedSimsOut + thinning
                hapPosExts = [ '.%s-%d' % ( hapOrPos, pop ) for hapOrPos in ('pos', 'hap') for pop in (1,4,5) ]
                hapPosExtsOrigThin = [ e + thinExt for e in hapPosExts ]
                modelMutExts = [ '.model', '.mut' ]
                exts = modelMutExts + hapPosExts

                thinnedBase   =   thinnedSimsOut + '/' + scenSubdir + '/' + scenFileName
                unthinnedBase = unthinnedSimsOut + '/' + scenSubdir + '/' + scenFileName

                pr.addRule( comment = "Throw out some fraction of SNPs in each frequency bin, to simulate HapMap SNP " \
                                "ascertainment",
                            targets = [    thinnedBase + ext for ext in exts ],
                            sources = [  unthinnedBase + ext for ext in modelMutExts + hapPosExtsOrigThin ],
                            commands = [ ' '.join([ '../Operations/Ilya_Operations/sim/sfs/working/catfoo/utilities/' \
                                                    'thinMarkers_prob_' + SysTypeString(),
                                                    thinProbFile, '0' ] +
                                                  [  unthinnedSimsOut + '/' + scenSubdir + '/' + scenFileName + ext
                                                     for ext in hapPosExtsOrigThin ]) ] +
                                       [ 'mv ' + unthinnedBase + ext + thinExt + '.thin '
                                         + thinnedBase + ext for ext in hapPosExts ] +
                                       [ 'cp ' + unthinnedBase + ext + ' ' + thinnedBase + ext for ext in modelMutExts ],
                            name = 'thin',
                            mediumRuleName = 'thin_%s_%d' % ( scenName, replicaNum ),
                            attrs = dict( piperun_short = True ) )
    

def DefineRulesTo_RunSweepOnSims( pr, 
                                  Ddata, simsOut, thinExt = '', thinning = '', suffix = '',
                                  mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs, nreplicas = 100,
                                  pop2name = pop2name,
                                  tests = ( 'lrh', 'ihs', 'xpop' ),
                                  acceptExistingSimConfigFiles = False,
                                  setOptions = (), appendOptions = (),
                                  inputParamsFiles = [],
                                  runImportsLocally = True,
                                  noImports = False,
                                  powerSfx = '',
                                  doOnlyStages = None ):
    """Define rules for running simulations and doing Sweep analyses of them.

    Parameters:

       mutAges, mutPops, mutFreqs - parameters defining the selection scenario.
    
    """
    
    # Define the rules to do Sweep analyses of the simulations.
    # These rules are created by running a Perl script, sim_analysis_pipe.pl .
    # Each invocation of the script defines rules for one type of test.
    # The rules for each test are saved to an .xml file, and then all these files are
    # merged into pr.

    # Rather than explicitly invoking the script three times we define a small pipeline to do this.
    # The output of this pipeline is a pipeline definition for analyzing simulation results by each of the tests,
    # saved to an .xml file.
    
    simPipeline = PipeRun( name = 'DefineSims', descr = 'Define simulation and analysis pipeline' )

    if not powerSfx: powerSfx = Sfx( *mutAges )
    
    if not acceptExistingSimConfigFiles:
        simPipeline.addInvokeRule( invokeFn = WriteSimulationInfo,
                                   invokeArgs = Dict( 'Ddata mutAges mutPops mutFreqs nreplicas suffix '
                                                      'inputParamsFiles pop2name powerSfx' ) )

    test2pipeline = {}
    for test in tests:

        simSfx = ''
        if suffix: simSfx += suffix
        if thinning: simSfx += thinning

        thinExtHere = '' if thinning else thinExt

        testPipeline = os.path.join( 'Ilya_Temp', AddFileSfx( 'p.xml', pr.name, test ) )
        test2pipeline[ test ] = testPipeline
        
        if thinning:
            simPipeline.addRule( comment = "Copy config file",
                                 targets = "$Ddata/power_$test$simSfx/config$powerSfx.txt",
                                 sources = "$Ddata/power_$test$suffix/config.txt",
                                 commands = "cp $Ddata/power_$test$suffix/config.txt "
                                               "$Ddata/power_$test$simSfx/config$powerSfx.txt",
                                 mediumRuleName = 'copy_config_' + test, name = 'copy_config' );

        simPipeline.addRule( targets = testPipeline,
                             sources = [ '../Other/Ilya_Other/sweep/sims/scripts/sim_analysis_pipe.pl',
                                         "$Ddata/power_$test$simSfx/config$powerSfx.txt",
                                         "$Ddata/config$suffix/sims$powerSfx.txt",
                                         "$Ddata/config$suffix/scenarios$powerSfx.txt",
                                         "$Ddata/config$suffix/pops$powerSfx.txt" ],
                             commands = "../Other/Ilya_Other/sweep/sims/scripts/sim_analysis_pipe.pl "
                             " --only-write-pipeline $testPipeline " +
                             ( '--run-imports-locally ' if runImportsLocally else '' ) +
                             ( '--no-imports ' if noImports else '' ) + 
                             ( ( '--do-only-stages ' + doOnlyStages + ' ' ) if doOnlyStages else '' ) +
                             ( '--sim-suffix ' + simSfx if simSfx else '' ) +
                             ( '--powerSfx ' + powerSfx if powerSfx else '' )
                             + ( ' --thin-ext ' + thinExtHere if thinExtHere else '' )
                             + reduce( operator.concat, [ ' --set-option ' + setOption[0] + ' ' + setOption[1] for setOption in setOptions  ], '' )
                             + reduce( operator.concat, [ ' --append-option ' + appendOption[0] + " '" + appendOption[1] + "'" for appendOption in appendOptions  ], '' )
                             + " --target-test $test $Ddata/$simsOut$simSfx $Ddata/power_$test$simSfx",
                             name = 'simanal',
                             mediumRuleName = 'simanal_$test$simSfx',
                             comment = 'Define rules for analyzing simulations' )


    dbg('"RUNNING_simPipeline"')
    simPipeline.runForced( aftRuleDelay = 5 )
    dbg('"DONE_RUNNING_simPipeline"')

    for test in tests:
        pr.addPipelineFromFile( test2pipeline[ test ] )

    
def DefineRulesTo_RunSimsAndSweep( pr, Ddata, simsOut = 'simsOut',
                                   mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs, nreplicas = 100,
                                   suffix = '', thinning = '',
                                   thinExt = '',
                                   tests = ( 'lrh', 'ihs', 'xpop' ),
                                   doRunSims = True, doRunThinning = True, doRunSweep = True, inputParamsFiles = [],
                                   acceptExistingSimConfigFiles = False, pop2name = pop2name,
                                   runImportsLocally = True,
                                   doOnlyStages = None, DdataSeeds = '',
                                   setOptions = (), appendOptions = (),
                                   scen2alternateNeutralParams = {},
                                   useGenMap = None,
                                   powerSfx = '',
                                   withGeneConvBug = False, withNewCosi = False, withCosi = None, DdataMimic = None ):
    """Define rules for running simulations and doing Sweep analyses of them.

    Parameters:

       mutAges, mutPops, mutFreqs - parameters defining the selection scenario.
    
    """

    dbg('"Running_doRunSims"')
    if doRunSims: DefineRulesTo_RunSimsOnly( **Dict( 'pr mutAges mutPops mutFreqs nreplicas Ddata simsOut suffix '
                                                     'inputParamsFiles DdataSeeds useGenMap scen2alternateNeutralParams '
                                                     'withGeneConvBug withNewCosi withCosi DdataMimic' )
                                               ) 
    dbg('"Running_doThinning"')
    if doRunThinning: DefineRulesTo_DoThinning( **Dict( 'pr mutAges mutPops mutFreqs nreplicas Ddata simsOut thinning '
                                                        'thinExt suffix' ) )
    dbg('"Running_doRunSweep"')
    if doRunSweep: DefineRulesTo_RunSweepOnSims( **Dict( 'pr mutAges mutPops mutFreqs nreplicas Ddata simsOut thinning '
                                                         'suffix tests inputParamsFiles thinExt setOptions appendOptions '
                                                         'acceptExistingSimConfigFiles pop2name runImportsLocally '
                                                         'doOnlyStages powerSfx' ) )
    dbg('"FINISHED_runSimsAndSweep"')


