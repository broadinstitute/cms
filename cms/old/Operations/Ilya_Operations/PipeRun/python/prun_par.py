import logging, os, getpass
from Operations.MiscUtil import dbg, MakeSeq, Dict, AddFileSfx, DumpFile, AddSlash, SlurpFileLines, RandomString


def doSplit( splitFunc, splitFN, outDir, getio = None ):
    """Do the splitting"""

    chunkListFN = os.path.join( outDir, 'chunks.txt' )

    if getio: return dict( depends_on = splitFN, creates = chunkListFN )

    chunkFNs = splitFunc( splitFN, outDir = outDir )
    DumpFile( chunkListFN, '\n'.join( chunkFNs ) )
    

def runCmdParallelized( commands, depends_on, creates, comment, 
                        splitFunc, joinFunc, saveOutputTo = None, splitFN = None, joinFN = None,
                        name = None, mediumRuleName = None, getio = None ):
    """Run the specified command, using parallelization."""
    from Operations.Ilya_Operations.PipeRun.python.PipeRun import PipeRun

    dbg( '"IN_RUNCMDPAR_EEEEEE" depends_on creates saveOutputTo' )

    if creates is None: creates = ()
    
    commands = MakeSeq( commands )
    depends_on = MakeSeq( depends_on )
    creates = MakeSeq( creates )

    gio = Dict( 'depends_on creates comment name mediumRuleName' )
    dbg( 'gio' )
    if getio: return Dict( 'depends_on creates comment name mediumRuleName saveOutputTo',
                           uses = ( splitFunc, joinFunc ) )

    splitFN = splitFN or list( MakeSeq( depends_on ) )[0]
    joinFN = RandomString(12) if saveOutputTo else ( joinFN or list( MakeSeq( creates ) )[0] )

    assert any([ splitFN in command for command in commands ])
    assert saveOutputTo or any([ joinFN in command for command in commands ])


    logging.info( 'calling ' + str( splitFunc ) + ' to split ' + splitFN )
    outDir = os.path.join( '/broad/hptmp', getpass.getuser(), 'par', os.path.abspath( splitFN )[1:] )

    pr = PipeRun( name = 'splitting', descr = 'splitting' )
    r = pr.addInvokeRule( invokeFn = doSplit, invokeArgs = Dict( 'splitFunc splitFN outDir' ) )
    pr.runSubPipeline()
    
    chunkFNs = SlurpFileLines( r.creates[0] )

    logging.info( 'finished running ' + str( splitFunc ) + ' to split ' + splitFN )
    dbg( '"CHUNKS_ARE" chunkFNs' )

    pr = PipeRun( name = 'parallelizing', descr = 'parallelizing' )
    chunkOutFNs = []
    for chunkFN in chunkFNs:
        chunkOutFN = AddFileSfx( chunkFN, 'out' )
        chunkOutFNs.append( chunkOutFN )

        for command in commands:
            dbg( 'splitFN chunkFN chunkOutFN command command.replace(splitFN,chunkFN)' )
        
        pr.addRule( commands = [ command.replace( splitFN, chunkFN ).replace( joinFN, chunkOutFN ) for command in commands ],
                    depends_on = [ f if f != splitFN else chunkFN for f in depends_on ],
                    creates = [ f if f != joinFN else chunkOutFN for f in creates ],
                    saveOutputTo = None if saveOutputTo is None else chunkOutFN )

    pr.runSubPipeline()

    joinFunc( inFNs = chunkOutFNs, outFN = None if saveOutputTo else joinFN )

