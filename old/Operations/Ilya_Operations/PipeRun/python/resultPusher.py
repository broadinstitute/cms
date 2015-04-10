#!/usr/bin/env python

"""
Script: resultPusher

Checks the queue for completed jobs, and notifies <mqsub> instances that the corresponding job has been completed.
Runs on the same machine from which jobs have been submitted.
"""

#
# Fix: push results only for the local machine.  (the mq task directory should record the machine name on which
# the mqsub script ran.) 
#


from Operations.MiscUtil import SlurpFile, DumpFile, GetHostName, coerceVal, FirstLineOfFile, dbg
import os, errno, time, getpass, multiprocessing, optparse, sys, random, logging, contextlib, signal


def GetTaskAttr( fname, attrName, defaultVal = None ):
    """Return the specified attribute of a task, or the specified default value if the task does not have this attribute."""
    for line in SlurpFile( fname ).rstrip('\n').split('\n'):
        arg, val = line.split('\t')
        if arg == attrName: return coerceVal( val )
    return defaultVal

def GetSubmitHost( submithostFN ):
    """Return the name of the host on which the task runs"""
    return FirstLineOfFile( submithostFN )
            
def PushResults( options ):
    """Push finished results"""

    logging.info( 'Starting result pusher as process %d on host %s with options %s' %
                  ( os.getpid(), GetHostName(), options ) )

    # check that all queues exist
    assert all( os.path.exists( queue ) and os.path.isdir( queue )
                for queue in options.queues.split( os.pathsep ) )

    stopSignal = [ False ]

    def SetStopSignal( sigNum, stkFrm ):
        logging.info( 'Setting stop signal to stop pushers' )
        stopSignal[ 0 ] = True
        dbg( '"aftset" stopSignal' )

    signal.signal( signal.SIGUSR2, SetStopSignal )

    while not stopSignal[0]:
        
        for queue in options.queues.split( os.pathsep ):

            logging.info( 'Pushing results in ' + queue + '...' )

            # find an unclaimed task in this queue, and try to claim it
            try: taskDirs = filter( lambda f: f.startswith('mq'), os.listdir( queue ) )
            except EnvironmentError as e:
                logging.info( 'Error getting list of tasks in queue ' + queue + ': ' + str( e ) )
                # sleep a bit -- maybe it's some transient condition that will resolve itself
                time.sleep( 60 + random.normalvariate( 3.0, 1.0 ) )
                continue

            for taskDir in taskDirs:
                fullTaskDir = os.path.join( queue, taskDir )

                try:
                    pushingFN = os.path.join( fullTaskDir, 'pushing.dat' )
                    submithostFN = os.path.join( fullTaskDir, 'submithost.txt' )
                    if os.path.exists( os.path.join( fullTaskDir, 'completed.dat' ) ) \
                            and os.path.exists( submithostFN ) and GetSubmitHost( submithostFN ) == GetHostName() \
                            and not os.path.exists( pushingFN ):
                        try:
                            fd = os.open( pushingFN, os.O_CREAT|os.O_EXCL|os.O_WRONLY )
                        except EnvironmentError as e:
                            if e.errno not in ( errno.EEXIST, errno.EACCES, errno.EAGAIN ): raise
                            # another resultPusher beat us to this task -- go and check other tasks
                            continue

                        
                        os.write( fd, 'result being pushed by process %d on host %s' % ( os.getpid(), GetHostName() ) )
                        os.close( fd )

                        taskDescr = ''

                        try:
                            attrsFN = os.path.join( fullTaskDir, 'attrs.tsv' )
                            if os.path.exists( attrsFN ):
                                taskDescr += ' output in ' + GetTaskAttr( attrsFN, 'piperun_outputSavedTo' )
                        except EnvironmentError as e:
                            logging.info( 'Could not read attrs for task in ' + fullTaskDir + ': ' + str( e ) )
                        
                        try:
                            infoFNs = [ os.path.join( fullTaskDir, f ) for f in 'command.dat', 'attrs.tsv', 'claimed.dat' ]
                            infoContents = '\n'.join([ SlurpFile( f ) if os.path.exists( f ) else 'missing file: ' + f
                                                       for f in infoFNs ])

                            thePipe = os.open( os.path.join( fullTaskDir, 'getresult' ), os.O_WRONLY | os.O_NONBLOCK )
                            exitCodeReadOk = False
                            writeOk = False
                            closeOk = False
                            exitCode = 'UNKNOWN'
                            try:
                                exitCode = SlurpFile( os.path.join( fullTaskDir, 'exitCode.dat' ) ).strip()
                                exitCodeReadOk = True
                                os.write( thePipe, exitCode )
                                writeOk = True
                            finally:
                                os.close( thePipe )
                                closeOk = True

                            logging.info( 'Pushed result in ' + fullTaskDir + ': ' + ( 'nonzero ' if exitCode != '0' else '' ) + 'exit code ' + exitCode + taskDescr )

                            if not writeOk or not closeOk or not exitCodeReadOk: dbg( 'exitCodeReadOk writeOk closeOk' )
                            if exitCodeReadOk and exitCode != '0': logging.info( infoContents )
                            
                        except EnvironmentError as e:
                            logging.info( 'The task at ' + fullTaskDir + ' seems to have been orphaned: ' + e.strerror )

                except EnvironmentError as e:
                    logging.info( 'Error processing queue ' + queue + ' task ' + taskDir + ': ' + str( e ) )
                    # sleep a bit -- maybe it's some transient condition that will resolve itself
                    time.sleep( 60 + random.normalvariate( 3.0, 1.0 ) )

        # if we pushed at least something, go back and try again.  if not, wait.
        time.sleep( options.sleepInterval + random.normalvariate( 3.0, 1.0 ) )

allPushers = []
        
def StartPushers( use_args = None, as_daemons = False ):
    parser = optparse.OptionParser()
    parser.add_option( '-P', '--num-result-pushers', type='int', dest = 'numResultPushers',
                       help='create NUMPROCS parallel result pushers', metavar='NUMPROCS',
                       default=1 )
    parser.add_option( '-Q', '--queues', dest = 'queues',
                       default = os.path.join( '..', 'Other', 'queues', getpass.getuser() ),
                       help='push results for the specified QUEUES', metavar='QUEUES' )
    parser.add_option( '-S', '--sleepInterval', type='int',
                       help = 'between checks, sleep for SEC seconds', metavar = 'SEC', default = 20 )
    dbg( 'use_args' )
    (options, args) = parser.parse_args( args = sys.argv[1:] if use_args is None else list( use_args ) )
    dbg( 'options args' )
    assert not args
    
    for i in range( options.numResultPushers ):
        p = multiprocessing.Process( target = PushResults, args = ( options, ) )
        allPushers.append( p )
        p.daemon = as_daemons
        p.start()
        time.sleep( min( 1.0, random.normalvariate( 4.0, 1.0 ) ) )


def StopPushers():
    """Stop all runners"""

    for pusher in allPushers:
        if pusher.pid is not None:
            os.kill( pusher.pid, signal.SIGUSR2 )
            pusher.join()
        
@contextlib.contextmanager
def RunPushers( use_args = None ):
    """Do things while running result pushers"""

    StartPushers( use_args )
    yield
    StopPushers()

        
if __name__ == '__main__':
    print 'STARTING PUSHERS FROM MAIN'
    StartPushers()
                    
                    
                
            
        
        
