
from __future__ import division, with_statement

__all__ = ( 'DefineRulesTo_CreateSimulationParams', 'DefineRulesTo_runSims' )

from Operations.IDotData import IDotData
from Operations.Ilya_Operations.PipeRun.python.PipeRun import PipeRun
from Operations.MiscUtil import SystemSucceed, Dict, cross, Str, DumpFile, SlurpFile, MakeSeq, OutputOfCmd, dbg
from Operations.Shari_Operations.localize.Scenario import GetScenarios, SelectionScenario
from Operations.Shari_Operations.localize.PopConsts import pn_European, pn_EastAsian, pn_WestAfrican, pop2name
from operator import concat
from shutil import copyfile
import itertools, os

def CreateSimsParams_neutral( Ddata, suffix, inputParamsFiles, getio = None ):
	"""Write the neutral parameter file.
	"""

	inputParamsFiles = MakeSeq( inputParamsFiles )
	
	neutralParamsFile = Ddata + '/params_neutral' + suffix

	if getio: return dict( depends_on = inputParamsFiles, creates = neutralParamsFile )

	neutralParams = reduce( concat, map( SlurpFile, inputParamsFiles ) )
	
	DumpFile( neutralParamsFile, neutralParams )
	

def CreateSimsParams_selection( mutAges, mutPops, mutFreqs, Ddata, suffix, scen2alternateNeutralParams = {}, getio = None ):
	"""Create simulation parameter files for the cosi simulator.  Takes as input a list of files for the neutral
	scenario, concatenates them, and creates parameter files for the neutral scenario as well as the selection
	scenarios.
	"""

	neutralParamsFile = Ddata + '/params_neutral' + suffix
	selectionParamsFiles = dict( ( ( mutAge, mutPop, mutFreq ), Ddata + '/params' +
				       suffix + '/%dky/params_sel%d_%d' % ( mutAge, mutFreq, mutPop ) )
				     for mutAge in mutAges for mutPop in mutPops for mutFreq in mutFreqs )

	if getio: return dict( depends_on = [ neutralParamsFile ] + scen2alternateNeutralParams.values(),
			       creates = selectionParamsFiles.values() )

	neutralParams = SlurpFile( neutralParamsFile )

	# Check that the list of pops defined in the param file matches the list of pops we are analyzing
	definedPops = set( [ int( s.split()[1] ) for s in neutralParams.split( '\n' ) if s.startswith( 'pop_define' ) ] )
	assert set( mutPops ) <= definedPops

	for mutAge in mutAges:
		for mutPop in mutPops:

			popSize = [ s for s in neutralParams.split( '\n' )
				    if s.startswith( 'pop_size %d ' % mutPop ) ][0].split()[2]

			for mutFreq in mutFreqs:

				neutralParamsHere = neutralParams
				scen = SelectionScenario( mutAge = mutAge, mutPop = mutPop, mutFreq = mutFreq )
				if scen in scen2alternateNeutralParams:
					neutralParamsHere = SlurpFile( scen2alternateNeutralParams[ scen ] )

				yearsPerGen = 20
				selectionStartGen = int( mutAge * 1000 / yearsPerGen )
				mutFreqFraction = mutFreq / 100

				selCoeff = \
				    OutputOfCmd( Str( 'perl ../Operations/Ilya_Operations/sim/sfs/working/pardis2/calcS.pl '
						      '$popSize $mutFreqFraction $selectionStartGen'
						      ) ).split('\n')[0].split()[2].strip()
				
				with open( selectionParamsFiles[ ( mutAge, mutPop, mutFreq ) ], 'w' ) \
					    as selectionParamsFile:
					selectionParamsFile.write( neutralParamsHere )

					mutPos = .5

					selectionParamsFile.write( Str( 'pop_event sweep "selective sweep" '
									'$mutPop 1.0 $selCoeff $mutPos '
									'$mutFreqFraction\n' ) )
					

def DefineRulesTo_CreateSimulationParams( pr, mutAges, mutPops, mutFreqs,
					  Ddata = '../Data/Ilya_Data/sim/sfs/working/pardis2',
					  suffix = '', inputParamsFiles = None,
					  scen2alternateNeutralParams = {} ):

	"""Create simulation parameters for cross( mutAges, mutPops, mutFreqs ).
	"""

	dbg( '"YYYYYYYYYYYYY" inputParamsFiles' )
	if inputParamsFiles == None:
		inputParamsFiles = [ Ddata + '/simParams' + befAft + suffix + '.txt' for befAft in 'Bef', 'Aft' ]
	dbg( '"ZZZZZZZZZZZZZ" inputParamsFiles' )

	pr.addInvokeRule( invokeFn = CreateSimsParams_neutral, invokeArgs = Dict( 'Ddata suffix inputParamsFiles' ) )
	for mutAge in mutAges:
		for mutPop in mutPops:
			for mutFreq in mutFreqs:
				pr.addInvokeRule( invokeFn = CreateSimsParams_selection,
						  invokeArgs = dict( mutAges = ( mutAge, ), mutPops = (mutPop,),
								     mutFreqs = (mutFreq,),
								     **Dict( 'Ddata suffix scen2alternateNeutralParams' ) ),
						  mediumRuleNameSfx = ( mutAge, mutPop, mutFreq ),
                                                  attrs = Dict( 'mutAge mutPop mutFreq' ) )

	
def DefineRulesTo_runSims( pr, mutAges, mutPops, mutFreqs, nreplicas,
			   allPops = None,
			   Ddata = '../Data/Ilya_Data/sim/sfs/working/pardis2', simsOut = 'simsOut',
			   suffix = '', shortSimTime = True, DdataSeeds = '',
			   useGenMap = None, includeNeutral = True, withGeneConvBug = False, withNewCosi = False,
			   withCosi = None, DdataMimic = None ):
	"""Instantiate, for each combination of ( mutAge, mutPop, mutFreq ),   the script that creates simulation parameters
	for simulations with that selected-mutation-age.
	"""

	assert not ( DdataSeeds and DdataMimic )

	mutPops = MakeSeq( mutPops )
	mutAges = MakeSeq( mutAges )
	mutFreqs = MakeSeq( mutFreqs )
	
	if allPops is None: allPops = mutPops

	Dsims = Ddata + '/' + simsOut + suffix

	for scen in GetScenarios( **Dict( 'mutAges mutPops mutFreqs includeNeutral' ) ):
		if DdataSeeds: seeds = IDotData( os.path.join( DdataSeeds, 'replicastats', scen.scenDir(), 'simSeeds.tsv' ) )

		for replicaNum, seedsLine in zip( range( nreplicas ), seeds if DdataSeeds else itertools.repeat( None, nreplicas ) ):

			assert not DdataSeeds or seedsLine.replicaNum == replicaNum
			
			pfx = os.path.join( Dsims, scen.scenDir(), '%d_%s' % ( replicaNum, scen.scenName() ) )
			recombDir = '../Data/Ilya_Data/sim/sfs/working/pardis2'

			attrs = Dict( 'replicaNum', scenDir = scen.scenDir() )
                        if not scen.is_neutral():
                                attrs.update( mutAge = scen.mutAge, mutPop = scen.mutPop, mutFreq = scen.mutFreq )
                        else: attrs.update( mutAge = 0, mutPop = 0, mutFreq = 0 )
			if shortSimTime: attrs[ 'piperun_short' ] = True 

			mutAge = '%dky' % ( 0 if scen.isNeutral() else scen.mutAge )


			useGenMapFile = os.path.join( DdataMimic, 'simsOut', scen.scenDir(), '%d_%s.model' %
						      ( replicaNum, scen.scenName() ) ) if DdataMimic else ''
			useMutRateFile = os.path.join( DdataMimic, 'simsOut', scen.scenDir(), '%d_%s.mut' %
						       ( replicaNum, scen.scenName() ) ) if DdataMimic else ''
#			dbg( '"GGGGGGGGG" mutPops' )
			pr.addRule( targets = [ pfx + ext for ext in ( [ '.model', '.mut', '.cosiParams' ] +
								       ( [ '.recombParams' ] if not useGenMap else [] ) +
								       [ '.%s-%d' % ( hapOrPos, pop )
									 for hapOrPos in ('hap', 'pos') for pop in allPops ]  +
								       ( [] if ( withNewCosi or withCosi ) else [ os.path.join( Dsims, scen.scenDir(), 'treeinfo',
																'%d_%s.%s' % ( replicaNum, scen.scenName(), which ) )
														  for which in ( 'regions.tsv', 'mutlist.tsv', 'nodes.dat' )
														  + ( () if scen.isNeutral() else ( 'sweepinfo.tsv', ) ) ] ) ) ],
				    sources = [ Ddata + '/' +  ( 'params_neutral' + suffix if scen.isNeutral()
								 else 'params%s/%s/params_%s' % \
									 ( suffix, mutAge, scen.scenName() ) ) ] \
					    + ( [ useGenMap ] if useGenMap else [ recombDir + '/recParams_bestfit_generic', \
											  recombDir + '/autosomes_decode.distr' ] ) + \
					    ( [ useGenMapFile, useMutRateFile ] if DdataMimic else [] ),
				    commands = ' '.join(('perl ../Operations/Ilya_Operations/sim/sfs/working/pardis2/' \
								 'runOneSim.pl' + ( ' --coalSeed %ld --recombSeed %ld --useMutRate %s'
										    % ( long( seedsLine.coalescentSeed ),
											long( seedsLine.recombSeed ),
											seedsLine.GetStrItem( 'mutRate' ) )
										    if DdataSeeds else '' )
							 + ( ' --useGenMap ' + useGenMap if useGenMap else '' )
							 + ( ' --withGeneConversionBug' if withGeneConvBug else '' )
							 + ( ' --withNewCosi' if withNewCosi else '' )
							 + ( ( ' --withCosi ' + withCosi ) if withCosi else '' )
							 + ( ( ' --useGenMap ' + useGenMapFile + ' --useMutRateFile ' + useMutRateFile )
							     if DdataMimic else '' ),
							 scen.scenName(), mutAge,
							 str(replicaNum), Ddata, Dsims, suffix )),
				    name = 'RunOneSim',
				    attrs = attrs,
				    comment = 'Adding simulation', mediumRuleNameSfx = ( scen.scenName(), mutAge, replicaNum ) )

		# end: for each replica	

def WriteSimulationInfo( Ddata, mutAges, mutPops, mutFreqs, nreplicas, suffix, inputParamsFiles,
			 pop2name = pop2name, powerSfx = '', getio = None ):
	"""Write out info files describing what is in the simulation.  These files are used by sim_analysis_pipe.pl
	to define the simulation analysis pipeline.
	"""

	stdDir = '../Data/Shari_Data/sim/stdSimAnalConfig'

	configFiles = dict( [ ( stdDir + '/' + test + '_config.txt', Str( "$Ddata/power_$test$suffix/config$powerSfx.txt" ) )
				for test in ( 'lrh', 'ihs', 'xpop' ) ] ) 


	cfgDir = Ddata + '/config' + suffix

	if getio: return dict( depends_on = configFiles.keys() + list( inputParamsFiles ),
			       creates = [ cfgDir + '/' + f + powerSfx + '.txt' for f in 'scenarios', 'sims', 'pops'  ]
			       + configFiles.values() )

	neutralParams = reduce( concat, map( SlurpFile, inputParamsFiles ), '' )

	# Check that the list of pops defined in the param file matches the list of pops we are analyzing
	assert sorted( [ int( s.split()[1] ) for s in neutralParams.split( '\n' ) if s.startswith( 'pop_define' ) ] ) \
	    == sorted( mutPops )

	assert set( map( int, pop2name.keys() ) ) == set( mutPops )

	DumpFile( cfgDir + '/scenarios%s.txt' % powerSfx, '\n'.join( scen.scenDir()
							      for scen in GetScenarios( mutAges, mutPops, mutFreqs ) ) )
	DumpFile( cfgDir + '/sims%s.txt' % powerSfx, '%d\n%d' % (0, nreplicas-1) )
	DumpFile( cfgDir + '/pops%s.txt' % powerSfx , '\n'.join( '%s\t%d' % ( popName, popNum )
						   for popNum, popName in pop2name.items() ) )
	for fromFile, toFile in configFiles.items():
		copyfile( fromFile, toFile )
        
#if __name__ == '__main__':
#    print 'running sims...'
#    RunSims()
#    print 'finished running sims.'
    


                            
    
        
