#!/usr/bin/env python

"""Code for taking tasks from a queue on a local or remote machine, and running them.

questions:

    - if multiple local runners need to _read_ the same input file, and they try to rsync it to the
    same place, will this break things?  should each rule have its own private local dir?
    or, use local locks to avoid this conflict?

"""

from __future__ import division, with_statement
from Operations.MiscUtil import GetHostName, SlurpFile, DumpFile, dbg, SystemSucceed, coerceVal, EnsureDirExists
import os, errno, time, getpass, multiprocessing, atexit, optparse, sys, random, stat, types, logging, socket, signal, traceback, re, shutil

from abc import ABCMeta, abstractmethod

try:
    import fcntl
    have_fcntl = True
except ImportError:
    have_fcntl = False


try:
    raise ImportError
    # If the paramiko library is available, we can use SFTP to interact with a remote queue.
    # If not, we can still function, just the SFTP option will be unavailable.
    import paramiko as pm
    from Crypto import Random
    haveParamiko = True
except ImportError:
    haveParamiko = False

class FileSystem(object):
    """A unified interface to local or remote file systems"""
    __metaclass__ = ABCMeta

    def __init__(self):
        self.path = self  # to allow the usage of os.path routines
        self.pathsep = ':'

    @abstractmethod
    def listdir(self, dir):
        """List a directory"""
        pass

    @abstractmethod
    def exists(self, pathname):
        """Check if a path exists"""
        pass

    @abstractmethod
    def isdir(self, pathname):
        """Check if a path is an existing directory"""

    @abstractmethod
    def remove(self, fileName):
        """Erase a file"""
        pass

    @abstractmethod
    def open(self, fileName, mode):
        """Open a file"""
        pass

    @abstractmethod
    def close(self, fileHandle):
        """Close a file"""
        pass

    @abstractmethod
    def write(self, fileHandle, string):
        """Write a string to a file"""
        pass

    @abstractmethod
    def stat(self, fileName):
        """Return various attributes of a file"""
        pass

    @abstractmethod
    def SlurpFile(self, fname):
        """Slurp entire file into memory"""
        pass

    def disconnect(self):
        """Disconnect from the file system"""
        pass

    def join(self,*args):
        """Call os.path.join on the args, for code compatibility"""
        return os.path.join(*args) 

class LocalFileSystem(FileSystem):

    """An implementation of FileSystem that uses the local file system"""

    def __init__(self):
        super(type(self), self).__init__()
    
    def listdir(self, dir):
        """List a directory"""
        return os.listdir( dir )

    def exists(self, pathname):
        """Check if a path exists"""
        return os.path.exists( pathname )

    def isdir(self, pathname):
        """Check if a path is an existing directory"""
        return os.path.isdir( pathname )

    def remove(self, fileName):
        """Erase a file"""
        return os.remove( fileName )

    def open(self, fileName, mode):
        """Open a file"""
        return os.open( fileName, mode )

    def close(self, fileHandle):
        """Close a file"""
        return os.close( fileHandle )

    def write(self, fileHandle, string):
        """Write a string to a file"""
        os.write( fileHandle, string )

    def stat(self, fileName):
        """Return various attributes of a file"""
        return os.stat( fileName )
        
    def SlurpFile(self, fname):
        """Slurp entire file into memory"""
        return SlurpFile( fname )


class RemoteFileSystem(FileSystem):
    """An implementation of FileSystem that uses remote SFTP access."""

    def __init__(self, remote, pw = None, pkey = None ):
        """Initialize the remote file system"""
        super(type(self), self).__init__()
        if not pw: pw = None
        if not pkey: pkey = None
        self.remote = remote
        self.username, hostAndPath= remote.split( '@' )
        self.hostname, self.root = hostAndPath.split( ':' )
        self.pw = pw
        if not pw and not pkey: pkey = '~/.ssh/identity'
        if isinstance( pkey, types.StringTypes ):
            pkey = pm.RSAKey.from_private_key_file( filename = os.path.expanduser( pkey ) )
        self.pkey = pkey
        self.transport = None
        self.client = None

    def connect(self):
        if self.transport and self.transport.is_authenticated(): return
        self.transport = pm.Transport( self.hostname )
        self.transport.connect( username = self.username, password = self.pw, pkey = self.pkey )
        self.client = self.transport.open_sftp_client()
        self.client.chdir( self.root )

    def disconnect(self):
        self.transport.close()
        self.transport = None
        self.client = None
    
    def listdir(self, dir):
        """List a directory"""
        self.connect()
        return self.client.listdir( dir )
        
    def exists(self, pathname):
        """Check if a path exists"""
        self.connect()
        try:
            self.client.stat( pathname )
            return True
        except EnvironmentError:
            return False

    def isdir(self, pathname):
        """Check if a path is a directory"""
        self.connect()
        try: return stat.S_ISDIR( self.client.stat( pathname ).st_mode )
        except EnvironmentError: return False

    def remove(self, fileName):
        """Erase a file"""
        self.connect()
        return self.client.remove( fileName )

    def open(self, fileName, mode):
        """Open a file"""
        self.connect()
        theMode = ''
        if mode == os.O_RDONLY: theMode += 'r'
        if mode & os.O_WRONLY: theMode += 'w'
        if mode & os.O_APPEND: theMode += 'a'
        if mode & os.O_EXCL: theMode += 'x'
        if mode & os.O_RDWR: theMode += '+'
        if mode & os.O_WRONLY and not mode & os.O_CREAT and not self.exists( fileName ):
            raise IOError( errno.ENOENT, 'File does not exist', fileName )
        return self.client.open( fileName, theMode )

    def close(self, fileHandle):
        """Close a file"""
        fileHandle.close()

    def write(self, fileHandle, string):
        """Write a string to a file"""
        fileHandle.write( string )

    def stat(self, fileName):
        """Return various attributes of a file"""
        return self.client.stat( fileName )
        
    def SlurpFile(self, fname):
        """Slurp entire file into memory"""

        try:
            logging.info( 'slurping file ' + fname )
            f = self.open( fname, os.O_RDONLY )
            return f.read()
        finally:
            self.close( f )


def GetMemReq( fs, fname ):
    """Return max memory requested by a task"""
    for line in fs.SlurpFile( fname ).rstrip( '\n' ).split('\n'):
        arg, val = line.split('\t')
        if arg == 'piperun_memreq': return float( val )
    return 0


def GetProcReq( fs, fname ):
    """Return minimum number of processors requested by a task"""
    for line in fs.SlurpFile( fname ).rstrip( '\n' ).split('\n'):
        arg, val = line.split('\t')
        if arg == 'piperun_procreq': return int( val )
    return 1


def GetTaskAttr( fs, fname, attrName, defaultVal = None ):
    """Return the specified attribute of a task, or the specified default value if the task does not have this attribute."""
    if fs.path.exists( fname ):
        for line in fs.SlurpFile( fname ).rstrip( '\n' ).split('\n'):
            arg, val = line.split('\t')
            if arg == attrName: return coerceVal( val )
    return defaultVal

def RunTasks( options ):
    """Take tasks from the specified queue directory, and run them.

    Params:

       options - see command-line parameter definition in main() below.
       
    """

    if haveParamiko: Random.atfork()

    startClock = time.time()

    logging.info( 'Starting runner (process id %d on host %s) with options %s'
                  % ( os.getpid(), GetHostName(), options ) )

    stopSignal = [ False ]

    def SetStopSignal( sigNum, stkFrm ):
        logging.info( 'Setting stop signal to stop runners' )
        stopSignal[ 0 ] = True
        dbg( '"aftset" stopSignal' )

    signal.signal( signal.SIGUSR1, SetStopSignal )

    fs = RemoteFileSystem( remote = options.remote, pw = options.password, pkey = options.pkey ) \
        if options.remote else LocalFileSystem()

    # check that all queues exist
    assert all( fs.exists( queue ) and fs.isdir( queue )
                for queue in options.queues.split( fs.pathsep ) )

    # register a cleanup routine so that, if we claim a task and then
    # crash midway through it, our lock on the task is erased, so that
    # another runner can pick up the task.  Note that if the _task_ fails
    # with an error code, that's fine -- we just report the error code
    # to the mqsub.sh script instance that submitted the task.
    # The cleanup here happens only if the runner crashes before receiving
    # a proper exit code from the task.
    fileToErase = [ None ]
    @atexit.register
    def DoErase( eraseWhat = fileToErase ):
        if eraseWhat[0] and fs.exists(eraseWhat[0]):
            #print 'runner AtExit: removing ' + eraseWhat[0]
            fs.remove( eraseWhat[0] )

    # var: lastFinish - time when a task last finished.
    lastFinish = time.time()

    queues = options.queues.split( fs.pathsep )
    lastQueueModTime = [ None ] * len( queues )

    skipDirs = set([ 'newtask.dat' ] )

    numTasksRun = 0
    numProcsAvail = int( os.getenv( 'LSB_DJOB_NUMPROC', 1 ) )
    dbg( 'numProcsAvail' )

    for queue in queues:
        EnsureDirExists( os.path.join( queue, 'succ' ) )
        EnsureDirExists( os.path.join( queue, 'fail' ) )
    
    while not stopSignal[0]:

        ranCommand = False

        if options.maxRunHours > 0 and ( time.time() - startClock  ) / 3600.0 > options.maxRunHours:
            logging.info( 'Runner exiting after CPU time of %s hours' % ( time.time() - startClock ) / 3600.0 )
            return

        if stopSignal[ 0 ]:
            logging.info( 'Runner stopped by stop signal' )
            return
        else: dbg( '"chkstop" stopSignal' )

        dbg( 'queues' )
        for queueNum, queue in enumerate( queues ):
            dbg( 'queueNum queue' )

            # do a quick check to see if any tasks have been added to the queue since we last checked
            newTaskFN = os.path.join( queue, 'newtask.dat' )
            try:
                curQueueModTime = fs.stat( newTaskFN ).st_mtime
                if curQueueModTime == lastQueueModTime[ queueNum ]: continue
                lastQueueModTime[ queueNum ] = curQueueModTime
            except EnvironmentError as e:
                if os.path.exists( newTaskFN ):
                    logging.warning( 'ERROR CHECKING FOR NEW TASKS in queue %s: %s' % ( queue, e ) )
                pass
            
            # find an unclaimed task in this queue, and try to claim it
            taskDirs = sorted( fs.listdir( queue ) )
            dbg( 'len(taskDirs)' )
            #random.shuffle( taskDirs )
            dbg( 'os.environ.get("MQ_FIRST_DIR")' )
            if 'MQ_FIRST_DIR' in os.environ and os.environ[ 'MQ_FIRST_DIR' ] in taskDirs:
            
                taskDirs = [ os.environ[ 'MQ_FIRST_DIR' ] ] + taskDirs
                logging.info( 'putting specified dir first' )

            for taskDir in taskDirs:

                if taskDir in skipDirs: continue

                if options.maxRunHours > 0 and ( ( time.time() - startClock ) / 3600.0 ) > options.maxRunHours:
                    logging.info( 'Runner exiting after CPU time of %s hours' % ( ( time.time() - startClock ) / 3600.0 ) )
                    return

                if stopSignal[ 0 ]:
                    logging.info( 'Runner stopped by stop signal' )
                    return
                else: dbg( '"chkstop" stopSignal' )

                try:

                    while fs.path.exists( os.path.join( queue, 'noclaim.dat' ) ):
                        time.sleep( 60 + random.normalvariate( 10.0, 5.0 ) ) 
                
                    fullTaskDir = fs.path.join( queue, taskDir )
                    claimedFN = fs.path.join( fullTaskDir, options.claimedFN )

                    attrsFN = fs.path.join( fullTaskDir, 'attrs.tsv' )
                    cwdFN = fs.path.join( fullTaskDir, 'submitdir.txt' )

                    failedCond = []
                    
                    def saveVal( name, val, fc = failedCond ):
                        if not val: fc.append( name )
                        return val

                    if saveVal( 'ready', fs.path.exists( fs.path.join( fullTaskDir, options.readyFN ) ) ) \
                            and saveVal( 'not claimed', not fs.path.exists( claimedFN ) ) \
                            and saveVal( 'relocatable', ( not options.remote or \
                                                              all([ not f.startswith( '/' ) for which in ( 'sources', 'targets' )  \
                                          for f in fs.SlurpFile( fs.path.join( fullTaskDir, which + '.lst' ) ) ]) ) ) \
                            and saveVal( 'memOk', GetMemReq( fs, attrsFN ) <= options.maxMem ) \
                            and saveVal( 'minMemOk', options.minMem == 0 or GetMemReq( fs, attrsFN ) >= options.minMem ) \
                            and saveVal( 'minProc', GetProcReq( fs, attrsFN ) >= options.minProc ) \
                            and saveVal( 'maxProc', GetProcReq( fs, attrsFN ) <= numProcsAvail ) \
                            and saveVal( 'local', ( options.local_tasks or not GetTaskAttr( fs, attrsFN, 'piperun_run_locally', False ) ) ) \
                            and saveVal( 'onlyLocal',
                                         ( not options.only_local_tasks or GetTaskAttr( fs, attrsFN, 'piperun_run_locally', False ) ) ) \
                            and saveVal( 'short', ( not options.runOnlyShort or GetTaskAttr( fs, attrsFN, 'piperun_short', False ) ) )  \
                            and saveVal( 'long', ( not options.runOnlyLong or not GetTaskAttr( fs, attrsFN, 'piperun_short', False ) ) )  \
                            and saveVal( 'notRequeued', ( not options.noRequeuedTasks or not fs.path.exists( fs.path.join( fullTaskDir, 'requeued.dat' ) ) ) ) \
                            and saveVal( 'notFromHost', ( not options.onlyFromHost or socket.getfqdn() == options.onlyFromHost ) ) \
                            and saveVal( 'notFromPipeline', ( not options.onlyFromPipelineId or  \
                                                                  GetTaskAttr( fs, attrsFN, 'piperun_pipelineId' ) == options.onlyFromPipelineId ) ):

                        # try to claim the task
                        try:
                            fd = fs.open(claimedFN, os.O_CREAT|os.O_EXCL|os.O_WRONLY)
                        except EnvironmentError:
                            # another runner beat us to this task -- go and check other tasks
                            logging.info( 'another job beat us to claiming ' + fullTaskDir )
                            continue

                        try:
                            fs.write( fd, 'locked by process %d on host %s\n' % ( os.getpid(), GetHostName() ) )
                            for v in os.environ.keys():
                                fs.write( fd, '%s=%s\n' % ( v, os.environ[v] ) )
                        finally:
                            fs.close( fd )
                        # Tell our cleanup code to release this task if we crash.
                        fileToErase[0] = claimedFN
                        # get the command to run the task
                        theCMD = fs.SlurpFile( os.path.join( fullTaskDir, 'command.dat' ) ).strip()
                        theCmdDir = fs.SlurpFile( os.path.join( fullTaskDir, 'submitdir.txt' ) ).strip()
                        theCmdEnvFN = os.path.join( fullTaskDir, 'submitenv.txt' )

                        if options.remote:
                            assert have_fcntl
                            SystemSucceed( 'mkdir -p ' + os.path.join( options.localDataDir, fs.root[1:] ) )
                            for needDir in 'Operations', 'Classes', 'System', 'Other':
                                needDirFull = os.path.join( options.localDataDir, fs.root[1:], '..', needDir )
                                if not os.path.exists( needDirFull ):
                                    os.symlink( os.path.realpath( os.path.join( '..', needDir ) ), needDirFull )

                            # copy source files

                            # get exclusive locks on the source files
                            srcFiles = sorted( set( fs.SlurpFile( os.path.join( fs.root, fullTaskDir, 'sources.lst' ) ).rstrip( '\n' ).split( '\n' ) ) )
                            srcLockIds = []
                            srcLockFiles = []
                            for srcFile in srcFiles:
                                lockFile = os.path.join( options.localDataDir, 'mqlocks', srcFile[1:] )
                                if lockFile.endswith('/'): lockFile = lockFile[:-1]
                                lockFile += '.lock'
                                SystemSucceed( 'mkdir -p ' + os.path.dirname( lockFile ) )
                                gotLock = False
                                while not gotLock:
                                    try:
                                        openMode = os.O_CREAT|os.O_EXCL|os.O_WRONLY
                                        logging.info( 'opening ' + lockFile + ' with mode ' + str( openMode ) )
                                        lockId = os.open( lockFile, openMode )
                                        gotLock = True
                                    except EnvironmentError:
                                        logging.info( 'Could not create ' + lockFile + ' , waiting...' )
                                        time.sleep( 10 + random.normalvariate( 3.0, 1.0 ) )
                                fcntl.lockf( lockId, fcntl.LOCK_EX )
                                srcLockIds.append( lockId )
                                srcLockFiles.append( lockFile )
                                logging.info( 'Got lock on ' + lockFile )

                            SystemSucceed( 'rsync -zprv --files-from=:' + os.path.join( fs.root, fullTaskDir, 'sources.lst' ) +
                                           ' ' + fs.username + '@' + fs.hostname + ':/ ' + options.localDataDir )

                            for srcLockId, srcLockFile in zip( srcLockIds, srcLockFiles)[::-1]:
                                fcntl.lockf( srcLockId, fcntl.LOCK_UN )
                                os.close( srcLockId )
                                SystemSucceed( 'rm -rf ' + srcLockFile )

                            targets = fs.SlurpFile( os.path.join( fs.root, fullTaskDir, 'targets.lst' ) ).rstrip( '\n' ).split('\n')
                            targetDirs = set( map( os.path.dirname, filter( None, map( str.strip, targets ) ) ) )
                            dbg( '"DDDDDD" targetDirs' )

                            for targetDir in targetDirs:
                                assert targetDir.startswith( '/' )
                                tdir = os.path.join( options.localDataDir, targetDir[1:] )
                                SystemSucceed( 'mkdir -p ' + tdir + ' ' + os.path.join( tdir, 'makeinfo' ) )

                            theCMD = 'cd ' + os.path.join( options.localDataDir, fs.root[1:] ) + ' && ' + theCMD

                        logging.info( 'Under ' + claimedFN + ' RUNNING: ' + theCMD )
                        # Actually run the task; get its exit code
                        save_cwd = os.getcwd()
                        try:
                            os.chdir( theCmdDir )
                            logging.info( 'CWD=' + os.getcwd() )

                            runScriptFN = os.path.join( fullTaskDir, 'run.sh' )

                            with open( runScriptFN, 'w' ) as out:
                                out.write( '#!/usr/bin/env bash\n' )
                                out.write( 'set -e -o pipefail\n' )
                                with open( theCmdEnvFN ) as envFile:
                                    for line in envFile:
                                        if '=' not in line or line.startswith('module='): break
                                        equalIdx = line.index( '=' )
                                        envVarName = line[ :equalIdx+1 ]
                                        if not ( re.search( r'\W', envVarName ) or envVarName.startswith( 'LSB_' ) or \
                                           envVarName.startswith( 'LSF_' ) or \
                                           envVarName.startswith( 'LS_' ) or envVarName.startswith( 'SLURM' ) or \
                                           envVarName in \
                                           ( 'SYS_TYPE', 'MACHTYPE', 'VENDOR', 'OSTYPE',
                                             'DOMAINNAME', 'HOSTTYPE', 'SHORTHOST', 'SSH_TTY',
                                             'HOST', 'HOSTNAME', 'REMOTEHOST', 'STY' ) ):
                                            out.write( 'export ' + envVarName + "'" + line[ equalIdx+1: -1 ] + "'\n" )
                                out.write( theCMD )

                            os.chmod( runScriptFN, stat.S_IXUSR | stat.S_IRWXU )

                            try:
                                exitCode = os.system( runScriptFN )
                            except ( KeyboardInterrupt, SystemExit ):
                                interruptedFN = os.path.join( fullTaskDir, 'interrupted.dat' )
                                DumpFile( interruptedFN, 'interrupted' );
                                raise
                        finally:
                            os.chdir( save_cwd )
                        logging.info( 'Under ' + claimedFN + ' FINISHED RUNNING: ' + theCMD )
                        logging.info( 'Got exit code %d' % exitCode )

                        if options.remote:
                            # copy the target files and the output log back to the correct dirs on the remote system

                            # first, make sure the files all exist, and are no longer being written to.

                            time.sleep( options.aftTaskDelay )

                            os.system( 'rsync -zprv --files-from=:' + os.path.join( fs.root, fullTaskDir, 'targets.lst' ) +
                                       ' ' + options.localDataDir + ' ' + fs.username + '@' + fs.hostname + ':/' )

                        # If we succeeded in running the task (whether the task itself failed or not),
                        # tell the cleanup code to NOT release this task if we crash.
                        fileToErase[0] = None
                        # Tell the task submitter script that we are done, and what the task's
                        # exit code was.

                        if os.path.exists( os.path.join( fullTaskDir, 'nmq.dat' ) ):

                            time.sleep(3)
                            fd = fs.open( os.path.join( fullTaskDir, 'completed.dat' ), os.O_CREAT|os.O_EXCL|os.O_WRONLY )
                            fs.close(fd)
                            
                            try:
                                shutil.move( fullTaskDir, os.path.join( queue, 'succ' if exitCode == 0 else 'fail' ) )
                            except EnvironmentError as e:
                                logging.warning( 'Error moving ' + fullTaskDir + ' to ' + os.path.join( queue, 'succ' if exitCode == 0 else 'fail' ) + ' : ' + e )
                        else:
                            exitCodeFN = os.path.join( fullTaskDir, 'exitCode.dat' )
                            fd = fs.open( exitCodeFN, os.O_CREAT|os.O_EXCL|os.O_WRONLY )
                            bytesWritten = fs.write( fd, str( exitCode ) )
                            fs.close( fd )

                            time.sleep(3)
                            logging.info( 'Wrote exit code %s to file %s (%s bytes)' % ( exitCode, exitCodeFN, bytesWritten ) )

                            fd = fs.open( os.path.join( fullTaskDir, 'completed.dat' ), os.O_CREAT|os.O_EXCL|os.O_WRONLY )
                            fs.close(fd)

                        # Record that we actually ran a task here.
                        ranCommand = True
                        lastFinish = time.time()
                        numTasksRun += 1

                    else:
                        logging.info( 'did not take task ' + taskDir + ' ; reason: ' + str( failedCond ) );

                except:
                    excInfo = sys.exc_info()
                    logging.warning( 'Error trying to grab task from ' + taskDir + ' (%s), skipping...'
                                     % str( excInfo ) )
                    traceback.print_exc()

        dbg( 'ranCommand lastFinish time.time()-lastFinish' )
        if not ranCommand:
            waitTimeHere = time.time() - lastFinish
            if ( numTasksRun > 0 and options.maxWaitTime > 0 and waitTimeHere > options.maxWaitTime ) \
                    or ( numTasksRun == 0 and options.maxFirstWaitTime > 0 and waitTimeHere > options.maxFirstWaitTime ) :
                logging.info( 'Runner exiting after idle time of %s' % waitTimeHere )
                return
            time.sleep( options.taskCheckInterval + random.normalvariate( 3.0, 1.0 ) )

allRunners = []

def StartRunners( use_args = None, as_daemons = False ):
    parser = optparse.OptionParser()
    parser.add_option( '-P', '--num-runners', type='int', dest = 'numRunners',
                       help='create NUMPROCS parallel runners', metavar='NUMPROCS',
                       default=1 )
    parser.add_option( '-Q', '--queues', dest = 'queues',
                       default = '',
                       help='run processes from the specified QUEUES', metavar='QUEUES' )
    parser.add_option( '-C', '--task-check-interval', type='int', dest='taskCheckInterval',
                       default = 10,
                       help='sleep TASK_CHECK_INTERVAL seconds between checking for new tasks',
                       metavar='TASK_CHECK_INTERVAL' )
    parser.add_option( '-M', '--max-wait-time', type='float', dest='maxWaitTime', default=30,
                       help='if no new tasks found for MAX_WAIT_TIME seconds then quit',
                       metavar='MAX_WAIT_TIME' )
    parser.add_option( '-F', '--max-first-wait-time', type='float', dest='maxFirstWaitTime', default=300,
                       help='if no new tasks found for MAX_WAIT_TIME seconds then quit',
                       metavar='MAX_WAIT_TIME' )
    parser.add_option( '-m', '--max-mem', type='float', dest='maxMem',
                       default=2,
                       help='only run tasks that require at most MEM gigabytes',
                       metavar='MEM' )
    parser.add_option( '-N', '--min-mem', type='float', dest='minMem',
                       default=0,
                       help='only run tasks that require at least MEM gigabytes',
                       metavar='MEM' )
    parser.add_option( '-p', '--min-proc', type='int', dest='minProc',
                       default=1,
                       help='only run tasks that require at _least_ NPROC processors',
                       metavar='NPROC' )
    parser.add_option( '-l', '--local-tasks', action='store_true', default=False,
                       help='allow this runner to pick up tasks marked as piperun_run_locally' )
    parser.add_option( '-k', '--only-local-tasks', action='store_true', default=False,
                       help='allow this runner to pick up only tasks marked as piperun_run_locally' )
    parser.add_option( '-L', '--runtime-limit', type='float', dest='maxRunHours',
                       default=0,
                       help='after running for MAX_RUN_HOURS hours, do not start new jobs; quit.',
                       metavar='MAX_RUN_HOURS' )
    parser.add_option( '-S', '--runOnlyShort', default = False, action = 'store_true',
                       help = 'Run only short tasks' )
    parser.add_option( '-G', '--runOnlyLong', default = False, action = 'store_true',
                       help = 'Run only long tasks' )

    parser.add_option( '-n', '--noRequeuedTasks', default = False, action = 'store_true',
                       help = 'Do not take any tasks requeued after their runner has crashed' )

    parser.add_option( '-H', '--onlyFromHost', default = '',
                       help='Run only jobs submitted from HOST', metavar = 'HOST' )
                      
    parser.add_option( '-i', '--onlyFromPipelineId', default = '',
                       help='Run only jobs submitted from PIPELINEID', metavar = 'PIPELINEID' )
                      
    parser.add_option( '--readyFN', default = 'ready.dat',
                       help='use READYFN to indicate jobs ready to run', metavar = 'READYFN' )
    parser.add_option( '--claimedFN', default = 'claimed.dat',
                       help='use CLAIMEDFN to indicate jobs claimed', metavar = 'CLAIMEDFN' )
                      

    if haveParamiko:
        remoteOpts = optparse.OptionGroup( parser,
                                           'Remote options', 'Options for running tasks from a remote queue via SFTP' )
        remoteOpts.add_option( '-R', '--remote', dest = 'remote',
                               help='use a remote queue at REMOTEDIR, in the form user@host:dir', metavar='REMOTEDIR',
                               default='' )
        remoteOpts.add_option( '-W', '--password', dest = 'password',
                               help='Use PASSWORD for logging in to the remote account', metavar='PASSWORD',
                               default = '' )
        remoteOpts.add_option( '-K', '--pkey', dest='pkey', default='',
                               help='Use private RSA key from FILE for logging in to the remote account', metavar='FILE' )
        remoteOpts.add_option( '-D', '--local-data-dir', dest='localDataDir',
                               default='./mqlocalData',
                               help='copy remote data to specified LOCALDIR', metavar='LOCALDIR' )
        
        remoteOpts.add_option( '-Y', '--aft-task-delay', dest='aftTaskDelay', type='int',
                               default=5,
                               help='wait DELAY seconds after a task finishes before copying results back',
                               metavar='DELAY' )

        remoteOpts.add_option( '-U', '--queue-urls', dest='queueURLs',
                               default='',
                               help='use URLS to check queue status via http (gentler on the server than SFTP);'
                               'one URL for each queue, separated by spaces',
                               metavar='URLS' )

        parser.add_option_group( remoteOpts )
    
    (options, args) = parser.parse_args( args = sys.argv[1:] if use_args is None else list( use_args ) )
    assert not args
    options.remote = None
    if not haveParamiko: options.remote = None

    if not options.queues:
        options.queues = os.path.join( '..', 'Other', 'queues',
                                       getpass.getuser() if not options.remote
                                       else options.remote[ :options.remote.index('@') ] )
    
    for i in range( options.numRunners ):
        p = multiprocessing.Process( target = RunTasks, args = ( options, ) )
        allRunners.append( p )
        p.daemon = as_daemons
        p.start()
        if haveParamiko: Random.atfork()
        time.sleep( min( 1.0, random.normalvariate( 4.0, 1.0 ) ) )

def StopRunners():
    """Stop all runners"""

    for runner in allRunners:
        if runner.pid is not None:
            os.kill( runner.pid, signal.SIGUSR1 )
            runner.join()
        
if __name__ == '__main__': StartRunners()
