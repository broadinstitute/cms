
"""
    Module: PipeRun

    Tools for managing a bioinformatics pipeline.

    Usage:

    There is a set of files and a set of rules.
    Each rules specifies how to build a set of target files
    from a set of source files, by running a list of shell commands.
    You specify the rules, and the system takes care of
    keeping the files up-to-date and distributing the computation
    onto available machines, and of visualizing the pipeline.

    The typical usage sequence is:

    pr = PipeRun.PipeRun()
    pr.addRule( ... )
    pr.addRule( ... )
    pr.addRule( ... )
    pr.run()

    Implementation:

    A PipeRun object accumulates the pipeline description as it is built up;
    then, when you call run(), it writes out the pipeline description as an
    .xml file, and invokes the RunPipeRun program to execute the pipeline.
    
"""

from __future__ import with_statement, division
import platform
assert platform.python_version_tuple() >= ( 2, 6 )
from xml.dom.minidom import Document
import os, sys, inspect, binascii, base64, zlib, pickle, types, logging, hashlib, string, tempfile, pprint, numpy, stat, \
    traceback, functools, itertools, contextlib, operator, socket, __main__
from types import StringType
from urlparse import urlparse
from Operations.MiscUtil import SysTypeString, MakeAlphaNum, Str, ReserveTmpFileName, SystemSucceed, WaitForFileToAppear, \
    relpath, dbg, SlurpFile, Dict, MakeSeq, PV, MergeDicts, Sfx, DumpFile, RegisterTmpFile, AttrDict, GetDefaultArgs, \
    IsValidFileName, MAX_FILE_NAME_LEN, AddFileSfx, EnsureDirExists, DictGet
from Operations.Ilya_Operations.PipeRun.python.FunctionChecksums import FunctionChecksum
from string import Template
from Operations.IDotData import IDotData
from Operations.Ilya_Operations.PipeRun.python.resultPusher import StartPushers, StopPushers
from Operations.Ilya_Operations.PipeRun.python.runner import StartRunners, StopRunners
from Operations.Ilya_Operations.PipeRun.python.prun_par import runCmdParallelized


# turn on debug output
logging.basicConfig( level = logging.DEBUG, format='%(asctime)s %(levelname)-8s %(filename)s:%(lineno)s %(message)s' )

def vstackThese( inFNs, outFN, getio = None ):
    """Vstack these files"""
    if getio: return dict( depends_on = inFNs, creates = outFN )
    IDotData.vstack( *inFNs ).save( outFN )

def do_vstack( inFNs, outFN ):
    inFNs = tuple( inFNs )
    """Vstack these files"""
    for i, inFN in enumerate( inFNs ):
        SystemSucceed( ( 'more +2' if i else 'cat' ) + ' ' + inFN + ( ' >> ' if i else ' > ' )  + outFN )
        SystemSucceed( '../Operations/Ilya_Operations/PipeRun/add_newline.sh ' + inFN + ' >> ' + outFN )
    
def DefineRulesTo_vstack( pr, inFNs, outFN, comment = 'vstack some tsvs', name = 'do_vstack', attrs = dict( piperun_short = True ) ):

    inFNs = tuple( inFNs )
    if len( inFNs ) > 1:
        pr.addRule( commands = [ 'cat ' + inFNs[0], '../Operations/Ilya_Operations/PipeRun/add_newline.sh ' + inFNs[0] ] +
                    reduce( operator.concat, [ [ 'more +2 ' + f, '../Operations/Ilya_Operations/PipeRun/add_newline.sh ' + f ] for f in inFNs[1:] ] ),
                    sources = inFNs, saveOutputTo = outFN,
                    **Dict( 'comment name attrs' ) )
                    
def SplitFN( f, i ):
    """Construct a split file name"""
    dirName, fileName = os.path.split( f )
    return AddFileSfx( os.path.join( dirName, 'split', fileName ), i )


def noIgnoreRuleCond( *args, **kwargs ): return False

class PipeRun:

    """Manages a pipeline.

    Implementation details:

    Fields:

       doc - the XML document we're building
       pipeline - the <pipeline> node at the root of the document.
       rules - the <rules> node within the pipeline
       fileNameMap - the <fileNameMap> node within the pipeline
       pipelineFilesToInclude - xmls that should be merged with the rules
         added to this PipeRun object, before executing.
         For example, rules written out by perl scripts that use
         the Perl interface to PipeRun; or by other pipeline descrs
         written out earlier.

    """

    def __init__(self, name = None, descr = None, pythonCmd = os.getenv( 'PIPERUN_PYTHONCMD', 'python' ) ):
        """Create a PipeRun object. """

        if 'PIPERUN_REDIRECT' in os.environ:
            baseFN = os.path.splitext( os.path.basename( sys.argv[0] ) )[0]
            outDir = os.path.join( os.path.dirname( sys.argv[0] ), 'out', baseFN )
            EnsureDirExists( outDir )
            for i in itertools.count( 1 ):
                outFN = os.path.join( outDir, baseFN + '.out%d' % i )
                if not os.path.exists( outFN ):
                    logging.info( 'Redirecting output to\n' + outFN )
                    out = open( outFN, 'w' )
                    sys.stdout = out
                    sys.stderr = out
                    break
        else: logging.info( 'NOT redirecting output.' )
        
        self.doc = Document()
        self.pipeline = self.__addElt("pipeline", self.doc)
        if not name: name = os.path.splitext( os.path.basename( sys.argv[0] ) )[0]
        if not descr: descr = inspect.getdoc( __main__ )
        dbg( '"NameDescr" name descr' )
        self.__addElt("pipelineName", self.pipeline, name)
        self.name = name
        if descr: self.__addElt("pipelineDescr", self.pipeline, descr) 
        self.rules = self.__addElt("rules")
        self.fileDescrsElt = self.__addElt("fileDescrs")
        self.fileTypesElt = self.__addElt("fileTypes")
        self.fileNameMap = self.__addElt("fileNameMap")
        self.fileTypeMap = self.__addElt("fileTypeMap")
        self.pipelineFilesToInclude = []
        self.allTargets = []
        self.arg2attr = { }
        self.lsfNeedDirs = []
        self.parallelAccessDir = os.getenv( 'PIPERUN_PAR_DIR' )
        if self.parallelAccessDir and not os.path.isdir( self.parallelAccessDir ): self.parallelAccessDir = None
        self.copiedToParallel = set()
        dbg( 'self.parallelAccessDir' )
        self.disableSplitting = True
        self.ignoreRuleNames = set()
        self._ignoreRuleCond = noIgnoreRuleCond
        self.pythonCmd = pythonCmd
        self.defaultAttrs = dict()
        self.ruleGrps = ()

    @contextlib.contextmanager
    def ignoringRuleNames( self, *ruleNames ):
        """Prevent the addition of rules with specified names."""
        save_ignoreRuleNames = self.ignoreRuleNames
        try:
            self.ignoreRuleNames = self.ignoreRuleNames | set( ruleNames )
            yield None
        finally:
            self.ignoreRuleNames = save_ignoreRuleNames

    @contextlib.contextmanager
    def ignoringRuleCond( self, ruleCond ):
        """Prevent the addition of rules with names matching a given condition"""
        save_ignoreRuleCond = self._ignoreRuleCond
        try:
            self._ignoreRuleCond = ruleCond
            yield None
        finally:
            self._ignoreRuleCond = save_ignoreRuleCond

    @contextlib.contextmanager
    def settingAttrs( self, *args, **kwargs ):
        """Add the specified attrs to all rules defined within this context"""
        attrs = args[0] if args else kwargs
        if isinstance( attrs, types.StringTypes ):

            caller = inspect.currentframe()
            thisFileName = caller.f_code.co_filename
            while ( caller.f_code.co_filename == thisFileName ):
                caller = caller.f_back

            attrs = Dict( attrs, frame = caller.f_back )
            del caller
            
        save_defaultAttrs = self.defaultAttrs
        self.defaultAttrs = MergeDicts( self.defaultAttrs, attrs )
        yield None
        self.defaultAttrs = save_defaultAttrs

    @contextlib.contextmanager
    def ruleGrp( self, grpName ):
        """Add all rules within this context to the specified rule group"""
        save_ruleGrps = self.ruleGrps
        self.ruleGrps = self.ruleGrps + ( grpName, )
        yield None
        self.ruleGrps = save_ruleGrps

    def makePipelineId( self ): return socket.getfqdn() + ':' + str( os.getpid() )
        
    def addArgAttrs( self, **args ):
        """Specify that in addInvokeRule, argument arg is to be treated as attribute attr.
        If attr not given, it defaults to arg."""
        self.arg2attr.update( args )

    def addLsfNeedDir( self, dir ):
        """Specify a directory that needs to be available for any rules to run"""
        if dir.endswith( '/' ): dir = dir[ :-1 ]
        if dir not in self.lsfNeedDirs:
            self.lsfNeedDirs.append( dir )

    def addPreExistingFiles( self, preExistingFiles ):
        """Designate the given files as pre-existing files that are not created
        by the pipeline, but are assumed to already exist.

        Note that if you pass targetsWithNoRuleOk = True to
        run(), then any files that are not the target of a rule will
        be assumed to be pre-existing.  If targetsWIthNoRuleOk is False
        (the default), any target file that you try to make
        must be either the target of some rule or explicitly declared
        to be pre-existing by calling this method.
        
        """
        self.__filesToXML( "pre_existing_files", preExistingFiles, self.rules )

    def addInvokeRule( self, invokeFn, invokeArgs = {}, targets = None, sources = None, attrs = {}, binaries = None,
                       mediumRuleName = None, mediumRuleNameSfx = '', name = None, nameSfx = '',
                       creates = None, depends_on = None, uses = None, usesOld = None,
                       comment = None, invokeFnOld = None, fileDescrs = {}, fileTypes = {},
                       ignoreArgs = set(),
                       mapFromStdin = None, mapToStdout = None, regionBeg = None, regionEnd = None,
                       chunkMap = {},
                       commandsReadable = None,
                       saveOutputTo = None, staticMethodClass = None ):
        """Adds a rule for invoking a Python function with given args, in a freshly started
        Python interpreter.

        Implementation notes:

        All rule invocations must be reduced to invoking shell commands, so
        we create a shell command whose invocation will run the specified Python function with
        the specified parameters, in a freshly started Python interpreter.
        """

        if hasattr( invokeFn, 'im_class' ):
            staticMethodClass = invokeFn.im_class
            invokeFn = invokeFn.im_func

        # make sure we can invoke invokeFn using the doOp.py script,
        # by starting a fresh Python interpreter.
        assert os.path.exists( invokeFn.func_code.co_filename )

        # for compatibility with Dan and Elain's system, the 'creates' argument may be used
        # in place of the 'targets' argument, and 'depends_on' argument may be used in place
        # of the 'sources' argument.
        if creates and not targets: targets = creates
        if depends_on and not sources: sources = depends_on 

        # some elements of the rule (inputs, outputs, attrs etc) may be specified inside invokeFn
        # using the special 'getio' parameter (invoking invokeFn with getio==True returns a dictionary
        # giving the values of these rule elements).  this is convenient when you want to specify
        # the rule elements inside the function (close to the code, so easier to keep their correspondence
        # to the function's actual I/O correct), and when the exact files read/written by a function
        # depend on the function's arguments.

        # here, if the function has the special 'getio' parameter, use it to get the values of rule
        # attributes, by invoking the function with getio=True, and examining the returned dictionary.

        args, varargs, varkw, defaults = inspect.getargspec( invokeFn )


        # var: paramMap - map from function parameter to its value during the invocation; includes
        #   the default parameters.
        paramMap = MergeDicts( dict( zip( args[ -len( defaults ): ], defaults ) if defaults else () ),
                               invokeArgs )

        extraInfo = {}

        def applyParams( s ):
            """Replace params with values"""
            return Template( str( s ) ).substitute( paramMap )

        def fixParams( v, filterNone = True ):
            if filterNone: v = filter( None, MakeSeq( v ) )
            return applyParams( v ) if not isinstance( v, ( tuple, list ) ) else tuple( map( applyParams, v ) )

        splitInfo = None
        splitInfoNumChunks = None

        ignoreArgsHere = ignoreArgs

        if 'getio' in args:
            # when invoked with getio = True, invokeFn should immediately return
            # a dictionary specifying the various rule elements.  it should not
            # actually execute the rule.
            ruleElems = invokeFn( getio = True, **invokeArgs )

            # rule elements passed as args directly to addInvokeRule() will override
            # any rule elements defined by the function itself; but if a rule element
            # (e.g. targets) is not given as argument to addInvokeRule(), and is
            # returned by invokeFn as part of its getio=True invocation,
            # then that rule element will be used.

            if not targets and 'creates' in ruleElems: targets = fixParams( ruleElems['creates'] )
            if not targets and 'targets' in ruleElems: targets = fixParams( ruleElems['targets'] )
            if not sources and 'sources' in ruleElems: sources = fixParams( ruleElems['sources'] )
            if not sources and 'depends_on' in ruleElems: sources = fixParams( ruleElems['depends_on'] )

            if not saveOutputTo: saveOutputTo = ruleElems.get( 'saveOutputTo', None )

            if not uses and 'uses' in ruleElems: uses = ruleElems['uses']
            if not comment and 'comment' in ruleElems: comment = fixParams( ruleElems['comment'], filterNone = False )
            if not name and 'name' in ruleElems: name = ruleElems['name']
            if not mediumRuleName and 'mediumRuleName' in ruleElems:
                mediumRuleName = fixParams( ruleElems['mediumRuleName'], filterNone = False )
            if not mediumRuleNameSfx and 'mediumRuleNameSfx' in ruleElems:
                mediumRuleNameSfx = fixParams( ruleElems['mediumRuleNameSfx'], filterNone = False )

            if not invokeFnOld and 'invokeFnOld' in ruleElems: invokeFnOld = ruleElems[ 'invokeFnOld' ]
            if not ignoreArgsHere and 'ignoreArgs' in ruleElems: ignoreArgsHere = ruleElems[ 'ignoreArgs' ]

            if 'extraInfo' in ruleElems: extraInfo = ruleElems[ 'extraInfo' ]

            if not self.disableSplitting and splitInfo is None and 'split' in ruleElems: splitInfo = ruleElems[ 'split' ]
            if splitInfoNumChunks is None and 'splitNumChunks' in ruleElems:
                splitInfoNumChunks = ruleElems[ 'splitNumChunks' ]

            # attrs can be given both via getio=True and as argument to addInvokeRule(): if both are given, the dictionaries
            # are merged; if they share any attributes, the attributes given in the addInvokeRule() call will
            # be used.
            if 'attrs' in ruleElems: attrs = MergeDicts( ruleElems[ 'attrs' ], attrs )
            if 'fileDescrs' in ruleElems:
                fileDescrs = MergeDicts( ruleElems[ 'fileDescrs' ], fileDescrs )

                fileDescrs = dict( ( file, applyParams( descr ) if isinstance( descr, str ) else
                                     ( applyParams( descr[0] ),
                                       tuple( ( col, applyParams( colDescr ) ) for col, colDescr in descr[1] ) ) )
                                     for file, descr in fileDescrs.items() )

            if 'fileTypes' in ruleElems:
                fileTypes = MergeDicts( ruleElems[ 'fileTypes' ], fileDescrs )

        # endif: 'getio' in args

        if 'splitByCols' in attrs and 'piperun_short' not in attrs:
            attrs[ 'piperun_short' ] = True

        if invokeFnOld:
            argsOld, varargsOld, varkwOld, defaultsOld = inspect.getargspec( invokeFnOld )
            if 'getio' in argsOld:
                ruleElemsOld = invokeFnOld( getio = True, **invokeArgs )
                if not usesOld and 'uses' in ruleElemsOld: usesOld = ruleElemsOld[ 'uses' ]

        # OPTFIX: old code versions
        #    - support chaining of code
        #    - check that new version creates at least the old creates, and depends_on at most the old depends_on [?]
        #    - allow saving just the old checksum not the old code


        # make a list of Python functions used by invokeFn.   this includes invokeFn itself,
        # as well as any functions explicitly listed in the uses argument.
        # now that we actually collect the functions that are actually called during execution,
        # this code may be extraneous.   but we'll leave it here, so that the user can always manually
        # specify some functions on which invokeFn depends, in case our automatic determination code
        # fails to pick them up.    
        # var: usedFns - list of Python functions invoked when invokeFn is invoked.
        usedFns = [ invokeFn ]
        if uses: usedFns += ( [ uses ] if isinstance( uses, types.FunctionType ) else list( uses ) )

        if invokeFnOld:
            usedFnsOld = [ invokeFnOld ]
            if usesOld: usedFnsOld += ( [ usesOld ] if isinstance( usesOld, types.FunctionType ) else list( usesOld ) )

        # automatically extract some attrs from invokeArgs, if configured.
        # for example, the user can configure to always extract a 'chr' argument passed to invokeFn
        # as a 'chrom' attribute of the created rule.
        if invokeArgs:
            for arg, val in invokeArgs.iteritems():
                if arg in self.arg2attr:
                    if attrs == None: attrs = { }
                    attrs[ self.arg2attr[ arg ] ] = val

        # to invoke invokeFn from a fresh Python interpreter using the doOp.py script,
        # we need to get the full module name in which invokeFn is defined.
        # then we can import invokeFn from that module inside doOp.py, and call it
        # with invokeArgs.            

        relFilePath = relpath( invokeFn.func_code.co_filename )

        #relFilePath = relFilePath.replace( '../../../../selection/sweep2/nsvn/', '../' )

        assert relFilePath.startswith('..' + os.path.sep ) and os.path.splitext( relFilePath )[1] == '.py'
        modulePath = os.path.splitext( relFilePath )[0][3:].replace( os.sep, '.' )

        commentGeneric = comment if comment else \
            ( inspect.getdoc( invokeFn ) if invokeFn.__doc__ else invokeFn.__name__ )

        defaultArgs = GetDefaultArgs( invokeFn )
        if isinstance( ignoreArgsHere, types.StringTypes ): ignoreArgsHere = ignoreArgsHere.strip().split()
        ignoreArgsHere = set()
        invokeArgsTrimmed = dict( ( arg, val ) for arg, val in invokeArgs.items()
                                  if arg not in ignoreArgsHere and ( arg not in defaultArgs or defaultArgs[ arg ] != val ) )

        func_name = ( '' if staticMethodClass is None else ( staticMethodClass.__name__ + '.' ) ) + invokeFn.func_name 

        def GetEncodedArgs( args, sfx ):

            # to pass args to the doOp.py script via the shell command line,
            # we have to pickle them and encode them in a safe way that does not
            # use any special characters.
            pickledArgs = pickle.dumps( args )
            encodedArgsOld = base64.urlsafe_b64encode( pickledArgs )
            compressedArgs = zlib.compress( pickledArgs, 9 ) 
            encodedArgs = ( '=' + base64.urlsafe_b64encode( compressedArgs ) ) \
                if len( compressedArgs ) < len( pickledArgs ) else encodedArgsOld
            if len(encodedArgs) > 4000:
                # first, try to compress the arguments

                argsFN = os.path.join( '../Other/tmp',
                                       hashlib.sha512( encodedArgs + modulePath +
                                                       func_name ).hexdigest() + '.args' + sfx )
                DumpFile( argsFN, encodedArgs )
                encodedArgs = 'file:' + argsFN
                #RegisterTmpFile( argsFN )

            return encodedArgs, encodedArgsOld

        encodedArgs, encodedArgsOld = GetEncodedArgs( invokeArgsTrimmed, 't' )
        encodedArgs2, encodedArgsOld2 = GetEncodedArgs( invokeArgs, '' )

        ruleName = MakeAlphaNum( invokeFn.__name__ + Sfx( nameSfx ) if not name else name )

        if splitInfo is not None:

            if self.parallelAccessDir:
                parSources = [ os.path.join( self.parallelAccessDir, os.path.abspath( os.path.realpath( source ) )[1:] )
                               + ( '/' if source.endswith( '/' ) else '' )
                               for source in MakeSeq( sources ) ]
                dbg( 'sources parSources' )
                for source, parSource in zip( MakeSeq( sources ), parSources ):
                    if source not in self.copiedToParallel:
                        self.addCpRule( source = source, target = parSource,
                                        name = 'cp_par_' + ruleName, attrs = attrs )
                        self.copiedToParallel.add( source )
                sources = parSources

            numPieces = 10
            chunksSfx = ( ' -c %d' % splitInfoNumChunks ) if splitInfoNumChunks is not None else ''
            for i in range( numPieces ):
                self.addRule( 
                    commands = ' '.join(( self.pythonCmd,
                                          '../Operations/Ilya_Operations/PipeRun/python/doOp.py -i %d -n %d%s'
                                          % ( i, numPieces, chunksSfx ),
                                          modulePath, func_name,
                                          encodedArgs ) ),

                    # the doOp command line is cryptic, because the invokeArgs are encoded.
                    # pass also a readable version of the command, showing the name of the function we are invoking,
                    # as well as the arguments we are passing to it, in readable form.
                    commandsReadable = commandsReadable or \
                        func_name + '(' + ', '.join( [ arg + ' = ' + pprint.pformat( val ) \
                                                                    for arg, val in invokeArgs.iteritems() ]  )+ ')',
                    # the rule comment is normally taken from the Python doc string of invokeFn,
                    # but can be overridden with the comment argument to addInvokeRule().
                    # if there is neither a doc string nor a comment argument, then use the function name
                    # as the comment.
                    commentGeneric = commentGeneric,
                    comment = applyParams( commentGeneric ),
                    # make a hash value of the functions invoked by invokeFn, that could be statically determined.
                    usesTests = [ hashlib.sha512( ''.join( map( FunctionString, usedFns ) ) ).hexdigest() ] \
                        + ( [ hashlib.sha512( ''.join( map( FunctionString, usedFnsOld ) ) ).hexdigest() ] if \
                                invokeFnOld else [] ),
                    name = ruleName,
                    targets = map( functools.partial( SplitFN, i = i ), targets ),
                    **Dict( 'sources binaries mediumRuleName mediumRuleNameSfx attrs fileDescrs fileTypes '
                            'paramMap saveOutputTo' ))

            for target in targets:
                self.addInvokeRule( invokeFn = vstackThese,
                                    invokeArgs = dict( inFNs = [ SplitFN( target, i ) for i in range( numPieces ) ],
                                                       outFN = target ),
                                    name = 'vstack_' + ruleName, attrs = attrs )

            return AttrDict( depends_on = PipeRun.FileTuple( MakeSeq( sources ) ),
                             creates = PipeRun.FileTuple( MakeSeq( targets ) ), extraInfo = extraInfo )

        else:

            return self.addRule( 
                          commands = ' '.join(( self.pythonCmd, '../Operations/Ilya_Operations/PipeRun/python/doOp.py' ) +
                                              ( ( '-I', mapFromStdin ) if mapFromStdin else () ) +
                                              ( ( '-O', mapToStdout ) if mapToStdout else () ) +
                                              ( ( '-b', ','.join( map( str, MakeSeq( regionBeg ) ) ) )
                                                if regionBeg else () ) +
                                              ( ( '-e', ','.join( map( str, MakeSeq( regionEnd ) ) ) )
                                                if regionEnd else () ) +
                                              ( ( '-m', base64.urlsafe_b64encode( zlib.compress( pickle.dumps( chunkMap ) ) ) )
                                                 if chunkMap else () ) +
                                              ( modulePath, func_name,
                                                encodedArgs ) ),
                          # for legacy reasons, pass the old way we used to invoke the invokeFn function.
                          # this is to prevent unnecessary rebuilding of target files constructed under the old
                          # invocation syntax.
                     
        

                          # commandsOld = ' '.join(( self.pythonCmd, '../Operations/Ilya_Operations/PipeRun/python/runOp.py',
                          #                          os.path.realpath( invokeFn.func_code.co_filename ), func_name,
                          #                          binascii.hexlify( pickle.dumps( invokeArgs) ) )),

                          # commandsOld2 = ' '.join(( self.pythonCmd, '../Operations/Ilya_Operations/PipeRun/python/doOp.py',
                          #                           modulePath, func_name,
                          #                           encodedArgsOld2)),

                          # commandsOld3 = ' '.join(( self.pythonCmd, '../Operations/Ilya_Operations/PipeRun/python/doOp.py',
                          #                           modulePath, func_name,
                          #                           encodedArgs2 ) ) if encodedArgs2 != encodedArgs else None,

                          # the doOp command line is cryptic, because the invokeArgs are encoded.
                          # pass also a readable version of the command, showing the name of the function we are invoking,
                          # as well as the arguments we are passing to it, in readable form.
                          commandsReadable = \
                              func_name + '(' + ', '.join( [ arg + ' = ' + pprint.pformat( val ) \
                                                                          for arg, val in invokeArgs.iteritems() ]  )+ ')',
                          # the rule comment is normally taken from the Python doc string of invokeFn,
                          # but can be overridden with the comment argument to addInvokeRule().
                          # if there is neither a doc string nor a comment argument, then use the function name
                          # as the comment.
                          commentGeneric = commentGeneric,
                          comment = applyParams( commentGeneric ),
                          # make a hash value of the functions invoked by invokeFn, that could be statically determined.
                          usesTests = [ hashlib.sha512( ''.join( map( FunctionString, usedFns ) ) ).hexdigest() ] \
                              + ( [ hashlib.sha512( ''.join( map( FunctionString, usedFnsOld ) ) ).hexdigest() ] if \
                                      invokeFnOld else [] ),
                          name = MakeAlphaNum( invokeFn.__name__ + Sfx( nameSfx ) if not name else name ),
                          **Dict( 'targets sources binaries mediumRuleName mediumRuleNameSfx attrs fileDescrs fileTypes paramMap '
                                  'extraInfo saveOutputTo' )
                          )

    def addFetchRule( self, URL, target, comment, name = None, mediumRuleName = None,
                      mediumRuleNameSfx = None, gunzip = True, attrs = None, fileDescrs = {},
                      fileTypes = {}, add_depends_on = () ):
        """Add a rule to fetch a web resource"""

        #assert gunzip
        
        zipExt = ''
        unzipCmd = ''
        unzipCmdOld = ''
        if gunzip and ( URL.endswith( '.gz' ) or URL.endswith( '.tgz' ) ):
            zipExt = '.gz' if URL.endswith( '.gz' ) else '.tgz'
            unzipCmd = 'gunzip -f'
            unzipCmdOld = 'gunzip'
            
        if gunzip and URL.endswith( '.bz2' ):
            zipExt = '.bz2'
            unzipCmd = 'bunzip2 -f'
        
        filename = os.path.basename( urlparse( URL )[2] )

        caller = inspect.currentframe()
        thisFileName = caller.f_code.co_filename
        while ( caller.f_code.co_filename == thisFileName ):
            caller = caller.f_back

        target = Str( target, caller )
        del caller

        origTarget = target
        if target.endswith('/'):
            target +=  ( filename if not zipExt else os.path.splitext( filename )[0] )
#            if URL.endswith( '.tgz' ): target += '.tar'

        outputFile = '' if origTarget.endswith( '/' ) or os.path.basename( origTarget ) + zipExt  == filename else \
            target + zipExt


        
        return self.addRule( sources = MakeSeq( add_depends_on ), # for now; should really always run wget with timestamps, or have an option to.
                             targets = target + ( '.tar' if URL.endswith( '.tgz' ) else '' ),
                             commands =
                             ( 'wget ' + ( ( '--timestamping ' +
                                             ( ( '--directory-prefix=' + os.path.dirname( target ) )
                                             if os.path.dirname( target ) else '' ) )
                                           if not outputFile else '--output-document=' + outputFile + ' ' ) + ' ' + URL, ) +
                             ( ( unzipCmd + ' ' + target + zipExt, ) if zipExt else () ),
                             commandsOld2 = ( 'wget --timestamping --output-document=' + target + zipExt + ' ' + URL, ) +
                             ( ( unzipCmd + ' ' + target + zipExt, ) if zipExt else () ),
                             commandsOld = ( ( 'wget --timestamping --output-document=' + target + zipExt + ' ' + URL, ) +
                             ( ( unzipCmdOld + ' ' + target + zipExt, ) if zipExt else () ) ) if unzipCmdOld else None,
                             **Dict( 'comment name mediumRuleName mediumRuleNameSfx attrs fileDescrs fileTypes' ) )


    def addCpRule( self, source, target, comment = None, attrs = {}, name = 'cp', mediumRuleName = None,
                   mediumRuleNameSfx = None, fileDescrs = {}, fileTypes = {}, symlink = False,
                   hardlink = False, diffAfter = False ):
        """Add a rule to copy one file or directory to another."""

        if 'piperun_short' not in attrs: attrs[ 'piperun_short' ] = True

        #assert not source.endswith( '/' )
        dbg( '"UORIG" source target' )
        if target.endswith( '/' ) and not source.endswith( '/' ): target += os.path.basename( source )
        dirCopy = target.endswith('/') and source.endswith('/')

        assert not ( dirCopy and diffAfter )

        link = symlink or hardlink
        cmd = 'cp ' + ( '-r ' if dirCopy else '' ) + source + ' ' + ( target if not dirCopy else target[:-1] ) \
            if not link else 'ln ' + ( '-s ' if symlink else '') + \
            os.path.abspath( source if not dirCopy else source[:-1] ) + ' ' + ( target if not dirCopy else target[:-1] )
        dbg( '"UUUUU" source target dirCopy cmd' )
        return self.addRule( sources = source, targets = target,
                             commands = cmd if not diffAfter else ( cmd, 'diff ' + source + ' ' + target ),
                             **Dict( 'comment attrs name mediumRuleName mediumRuleNameSfx fileDescrs fileTypes' )  )

    class FileTuple(tuple):

        def containing( self, s ):
            r = filter( lambda v: s in v, self )
            if not r: dbg( 'self s' )
            assert r
            return r[0] if len(r)==1 else r


        def not_containing( self, s ):
            r = filter( lambda v: s not in v, self )
            if not r: dbg( 'self s' )
            assert r
            return r[0] if len(r)==1 else r

        
    def addRule( self, commands, targets = None, sources = None,  comment = None, mediumRuleName = None,
                 mediumRuleNameSfx = None, attrs = {}, usesTests = None,
                 binaries = None, name = None, saveOutputTo = None, commandsOld = None,
                 commandsOld2 = None, commandsOld3 = None, commandsReadable = None,
                 creates = None, depends_on = None, fileDescrs = {}, fileTypes = {}, commentGeneric = None,
                 paramMap = {}, extraInfo = {},
                 splitFunc = None, joinFunc = None, splitFN = None, joinFN = None ):
        """Add a rule to the pipeline.

        The rule specifies how to make a set of target files
        from a set of source files using a sequence of commands.

        Required params:

           targets - list of target files to make; can be a list or a single string (if one file)
           sources - list of source files to make; can be a list or a single string (if one file)
           commands - list of shell commands to invoke, in sequence, in order to make
              the targets from the sources.  can be a list or a single string ( if one command ).
              typically there is just one command.

        Optional params:
              
           comment - a description of what the command(s) do
           mediumRuleName - a string representing this rule; will be used in various
              graph displays.

           fileDescr - description of some output and input files.
              dictionary of the form { fname : descr }, where fname can be a regexp.
              descr is either a string or a pair ( descr, dict ) where dict describes
              the columns in the file.

        """

        #dbg( '"AAAAAAAAAAAAAA" commands targets sources depends_on creates saveOutputTo comment attrs mediumRuleName name' )

        # for compatibility with Dan and Elaine's system, the 'creates' argument may be used
        # in place of the 'targets' argument, and 'depends_on' argument may be used in place
        # of the 'sources' argument.
        if creates and not targets: targets = creates
        if depends_on and not sources: sources = depends_on

        assert saveOutputTo is None or isinstance( saveOutputTo, types.StringTypes )

        if not targets and saveOutputTo: targets = saveOutputTo
        assert targets

        targets = tuple( MakeSeq( targets ) )

        if saveOutputTo and saveOutputTo not in targets: targets += ( saveOutputTo, )

        def IsString( s ): return isinstance( s, types.StringTypes )
        def ChkStrings( strs ):
            if not all( map( IsString, MakeSeq( strs ) ) ):
                dbg( 'MakeSeq(strs) map(IsString,MakeSeq(strs))' )
            assert all( map( IsString, MakeSeq( strs ) ) )
        
        ChkStrings( targets )
        ChkStrings( sources )
        ChkStrings( commands )

        if mediumRuleName and not name: name = mediumRuleName

        if self._ignoreRuleCond( name = name, sources = sources, targets = targets ): return
        
        if name and not mediumRuleName:
            mediumRuleName = name
            if mediumRuleNameSfx:

                mediumRuleName += Sfx( *mediumRuleNameSfx ) if isinstance( mediumRuleNameSfx, tuple ) or \
                   isinstance( mediumRuleNameSfx, list ) else Sfx( mediumRuleNameSfx )

        for f in tuple( MakeSeq( targets ) ) + tuple( MakeSeq( sources ) ):
            if not IsValidFileName( f ): logging.error( Str( "Invalid source or target file name $f for rule $mediumRuleName") )

        assert all( IsValidFileName( f ) for f in targets )
        assert all( IsValidFileName( f ) for f in sources )
        
        caller = inspect.currentframe()
        thisFileName = caller.f_code.co_filename
        while ( caller.f_code.co_filename == thisFileName ):
            caller = caller.f_back
            

        caller_filename, caller_lineno, caller_function = caller.f_code.co_filename, caller.f_lineno, caller.f_code.co_name

        applyStr = lambda s: Str( s, caller )

        def ForceError( showVal = None):
            assert False, 'ERROR'
            return None

        targets, sources, commands, commandsOld, commandsOld2, commandsOld3, commandsReadable, comment, mediumRuleName, \
            name = \
            map( lambda s: s if isinstance( s, types.NoneType ) else \
                     applyStr( s ) if isinstance( s, ( StringType, numpy.chararray ) ) else \
                     map( applyStr, s ),
                 ( targets, sources, commands, commandsOld, commandsOld2, commandsOld3, commandsReadable,
                   comment, mediumRuleName, name ) )

        #del applyStr
        #del caller  # to avoid circular references

        if comment == None: comment = mediumRuleName

        if name not in self.ignoreRuleNames and not self._ignoreRuleCond( name = name, sources = sources, targets = targets ):

            if splitFunc and joinFunc:

                return self.addInvokeRule( invokeFn = runCmdParallelized,
                                           invokeArgs = Dict( 'commands depends_on creates comment saveOutputTo '
                                                              'splitFunc joinFunc splitFN joinFN '
                                                              'name mediumRuleName') )
        
            rule = self.__addElt( "rule", self.rules )

            if saveOutputTo: rule.setAttribute( 'saveOutputTo', saveOutputTo )

            self.__filesToXML( "targets", targets, rule )
            self.allTargets += targets
            self.__filesToXML( "sources", sources, rule )
            self.__stringsToXML( "commands", "command", commands, rule )
            if commandsOld: self.__stringsToXML( "commandsOld", "command", commandsOld, rule )
            if commandsOld2: self.__stringsToXML( "commandsOld2", "command", commandsOld2, rule )
            if commandsOld3: self.__stringsToXML( "commandsOld3", "command", commandsOld3, rule )
            if commandsReadable: self.__stringsToXML( "commandsReadable", "command", commandsReadable, rule )
            self.__addElt( "comment", rule, comment )
            #if commentGeneric: self.__addElt( "commentGeneric", rule, commentGeneric )
            if mediumRuleName:
                self.__addElt( "mediumRuleName", rule, MakeAlphaNum( mediumRuleName ) )
                #                               MakeAlphaNum( mediumRuleName[:( (MAX_FILE_NAME_LEN-5) - max( map( len, targets ) ) )] ) )

            if binaries: self.__filesToXML( "binaries", binaries, rule )

            if name:
                if not attrs: attrs = {}
                attrs[ 'name' ] = name

            if 'piperun_pipelineId' not in attrs: attrs[ 'piperun_pipelineId' ] = self.makePipelineId()

            attrs = MergeDicts( self.defaultAttrs, attrs )
            if self.ruleGrps: attrs = MergeDicts( attrs, dict( ruleGrps = self.ruleGrps ) )

            if attrs:
                attrsElt = self.__addElt( "attrs", rule )
                for a, v in attrs.iteritems():
                    vals = MakeSeq( v )
                    for val in vals:
                        self.__addElt( a, attrsElt, val )

            definedAt = self.__addElt( "definedAt", rule )


            self.__addElt( "file", definedAt, relpath( caller_filename ) )
            self.__addElt( "line", definedAt, caller_lineno )
            self.__addElt( "function", definedAt, caller_function )

            if usesTests: self.__stringsToXML( "usesTests", "usesTest",
                                               MakeSeq( usesTests ), rule )

            if fileDescrs:
                for fileName, fileDescr in fileDescrs.items():
                    fileDescrElt = self.__addElt( "fileDescr", self.fileDescrsElt );
                    self.__addElt( "file", fileDescrElt, MakeSeq( targets )[ fileName ] if isinstance( fileName, int ) else fileName )

                    descrText = fileDescr if isinstance( fileDescr, str ) else fileDescr[0]

                    thisDescrElt = self.__addElt( "descr", fileDescrElt, descrText )
                    if isinstance( fileDescr, tuple ):
                        for col, colDescr in fileDescr[1]:
                            thisColDescrElt = self.__addElt( "columnDescr", thisDescrElt )
                            self.__addElt( "column", thisColDescrElt, col )
                            if isinstance( colDescr, ( StringType, numpy.chararray ) ):
                                self.__addElt( "descr", thisColDescrElt, colDescr )
                            else:
                                self.__addElt( 'sameAsFile', thisColDescrElt, colDescr[0] )
                                self.__addElt( 'sameAsCol', thisColDescrElt, colDescr[1] )

            if fileTypes:
                for fname, fileType in fileTypes.items():
                    for fileName in MakeSeq( fname ):

                        fileTypeElt = self.__addElt( "fileType", self.fileTypesElt );
                        self.__addElt( "file", fileTypeElt, MakeSeq( targets )[ fileName ] if isinstance( fileName, int ) else fileName )
                        self.__addElt( "type", fileTypeElt, fileType )
                            
        return AttrDict( depends_on = PipeRun.FileTuple( MakeSeq( sources ) ),
                         creates = PipeRun.FileTuple( MakeSeq( targets ) ),
                         extraInfo = extraInfo,
                         commands = commands )
        
        # end of: addRule()

    def addPipelineFromFile( self, pipelineFile ):
        """Add a pipeline from a disk file; its rules will be merged with the rules
        added to this PipeRun object. and the merged pipeline will be run
        when run() is called.
        """
        # Sometimes the pipeline file is written out by a call to an external script.
        # Play it safe and make sure it is completely finished writing.
        WaitForFileToAppear( pipelineFile )
        self.pipelineFilesToInclude.append( pipelineFile )

    def programBinary( self, execPath ):
        """Take the base name of a platform-specific binary, and
        ensure that the right binary is chosen based on the platform on which we run."""
        return execPath + self.__programBinarySuffix
        
    def addFileNameMapping( self, fileMapping ):
        """Add substitution rules to be applied to filenames in the pipeline.

        This is needed for web display purposes.  For instance, if you are reading
        a file named ../Data/myfile.tsv, showing the same path over the web
        may not show the file in the pipeline/ directory.  By creating the symlink
        pipeline/Data pointing to ../Data, and adding the file name mapping
        { '../Data' : 'Data' }, references from HTML files to '../Data/myfile.tsv'
        will find the right file.
        
        """
        for fromStr, toStr in fileMapping.iteritems():
            thisMapping = self.__addElt( "fileMapping", self.fileNameMap )
            self.__addElt( "from", thisMapping, fromStr )
            self.__addElt( "to", thisMapping, toStr )

    def addFileTypeMapping( self, fileTypeMapping ):
        """Add mapping of file patterns to file types.  Sometimes e.g. a text file does not
        end in .txt; this lets us add a pattern to say that files matching this pattern
        should be viewed as .txt files.
        """
        for fromStr, toStr in fileTypeMapping:
            thisMapping = self.__addElt( "fileTypeMapping", self.fileTypeMap )
            self.__addElt( "from", thisMapping, fromStr )
            self.__addElt( "to", thisMapping, toStr )
        

    def write( self, fname ):
        """Write out the pipeline description to the given XML file."""
        f=open(fname, "w")
        f.write( self.doc.toprettyxml(indent="  ") )
        f.close()
        WaitForFileToAppear( fname )

    def runForced( self, *args, **kwargs ):
        return self.run( *args, forceRun = True, lsfAbovePar = None, viewOnlyFromArgs = False, maxpar = 4, **kwargs )

    def forceMakeTarget( self, fname ):
        return self.runForced( ruleFilter = Str( 'Bwd( RulesWhere( "$fname" in target ) )' ) )

    def runSubPipeline( self, *args, **kwargs ):
        """Run the pipeline, inheriting values from a parent pipeline"""

        self.run( viewOnly = False, viewOnlyFromArgs = False, maxparFromArgs = False, lsfAbovePar = None,
                  useLsf = os.environ.get( 'USE_LSF', False ), useMQ = os.environ.get( 'USE_MQ', False ),
                  mqQueue = os.environ.get( 'MQ_QUEUE', None ),
                  maxpar = int( os.environ.get( 'MAXPAR', 4 ) ),
                  topScript = os.environ.get( 'PIPERUN_TOP_SCRIPT', None ),
                  shell = os.environ.get( 'PIPERUN_SHELL', None ), isSubPipeline = True, *args, **kwargs )

    def run( self, viewOnly = None, showRuleStatus = True, showLSFRuleStatus = False,
             outputDir = None, showDescrsInNodes = True,
             keepGoing = False, maxpar = None, maxparLocal = 4,
             viewOnlyFromArgs = True, maxparFromArgs = True, forceRun = False,
             forceLocal = False,
             lsfAbovePar = 10,
             lsfMaxPar = 300,
             mqMaxPar = 800,
             useLsf = False, targetsWithNoRuleOk = True, lsfDefaultQueue = 'week',
             lsfShortQueue = None,
             lsfJobGroup = None,
             lsfProject = 'sabeti_localize',
             lsfNeedDir = None,
             lsfResources = None,
             lsfUserJobPriority = None,
             useMQ = False,
             mqQueue = None,
             dryRun = False,
             ruleFilter = None, updateFetched = None, acceptExistingOutputs = None,
             forceRemake = None,
             ignoreCodeChangesFor = None,
             ignoreCmdChangesFor = None,
             ignoreNewDependenciesFor = None,
             ignoreAttrDiffs = None,
             noViewFor = None,
             disableHtmlViewers = True, checkModTimes = False,
             aftRuleDelay = None,
             sanityCheck = True, fnMap = None,
             detailCmdChanges = False,
             shell = None,
             passParentParams = None,
             topScript = None,
             isSubPipeline = False
              ):
        """Run the pipeline, constructing using earlier calls to
        addRule() and other methods.

        Parameters:

           viewOnly - if True, writes an HTML rendering of the pipeline to the directory
              specified.  The rendering shows how the pipeline steps connect together,
              and can also show the state of the pipeline (which steps have completed and
              which have not); it can also be used to browse the output of the steps.
              If viewOnly is False, run() actually runs the pipeline (but does not update
              its HTML rendering -- so you have to call run() a second time, this time with
              viewOnly as True, if you want to update the HTML rendering to reflect the new
              status of the pipeline.)

           showRuleStatus - if True, the HTML pipeline rendering produced by viewOnly=True will
              show the current status of each rule in the pipeline.  if False, only the pipeline
              structure will be shown.

           showLSFRuleStatus - if True, the HTML pipeline rendering produced by viewOnly=True and
              showRuleStatus=True will
              try to show the current LSF run status of rules (not just which rules have finished
              but which are running/pending/etc).

           outputDir - where to put the HTML pipeline rendering produced by viewOnly=True

           showDescrsInNodes - if True, pipeline step descriptions will be shown in nodes of the
              summary graph; if False, only step names will be shown.   (you can shift-doubleclick
              on any node in the summary graph to see the description in the middle frame).

           keepGoing - if True, if an error occurs in some pipeline step the steps that do not depend
              on its result will still be run.  if False, an error in one step will cause the entire
              pipeline to be stopped (currently running steps will be allowed to finish but new steps
              will not be started).

           maxpar - allow parallel processing, with up to this many processes at once.
              if useLsf is False, this specifies how many processes to run locally at once
              (useful if you have a multicore machine); if useLsf is True, this specifies how many
              jobs to submit to LSF at once.  Note that dependencies may prevent you from actually
              running this many jobs at once.

           lsfAbovePar - if an integer, useLsf is automatically set to True if maxpar is at least this value.
              used to automatically require the use of LSF if too many parallel processes are requested
              (since running too many processes locally can overload and freeze the local machine).
              
           useLsf - if True, jobs will be submitted to LSF.  if False, jobs will be run locally.
           
           targetsWithNoRuleOk - if False, any file assumed by the pipeline to already exist must
              be explicitly added by a call to addPreExistingFiles().  if True, any file that is not
              made by a rule is assumed to already exist.

           lsfDefaultQueue - if set, the LSF queue to which LSF jobs will be submitted by default,
              except for rules whose LSF queue designation is individually overridden.

           lsfShortQueue - if set, the LSF queue to which LSF short-running jobs will be submitted.
              
           lsfJobGroup - if set, the LSF job group to which LSF jobs will be added.  lets you use
              LSF commands to operate on all jobs in the named group.

           lsfProject - if set, the LSF project to which the LSF jobs will be assigned.
              this is used to determine which project is "charged" for this computation,
              in terms of allocating dynamic priority to jobs.

           lsfNeedDir - if set, LSF jobs will check that this directory is accessible before
              scheduling the jobs to run.   Needed because sometimes the filesystem with all your
              data becomes intermittently unavailable and if a job starts running at that time,
              it will fail with an error.  This prevents the job from running until the filesystem
              is checked to be available.

           lsfResources - if set, LSF jobs will request this LSF resource string.

           dryRun - if True, will do a "dry run" where the makefile is generated and
              everything is done except actually running the commands.

           ruleFilter - restrict the pipeline to rules matching this condition.
              The condition is a Python expression returning a set of rules.  In this expression,
              the function RulesWhere( cond ) denotes the set of rules matching the condition cond.
              cond is a Python expression in which rule attributes are Python identifiers.
              The most commonly used attribute is name, which denotes rule name.
              Every rule has it.  User-defined attributes can also be used.
              If a user-defined attribute is not defined for a given rule (but defined for at least one other rule),
              it evaluates to None for this rule.
              The expression defined('x') tests if a rule has the attribute 'x'.

              Predefined operations on sets of rules include Fwd(), Bwd(), FwdOnly() and BwdOnly().
              And of course you can use all Python set operators here.

            updateFetched - currently unused; in the future, will tell run() to re-fetch any files it fetches from
               the web, if the web versions have been updated.

            acceptExistingOutputs - for    
              
        Implementation:
        This method writes out the pipeline description as XML, then calls
        the RunPipeRun program to execute the pipeline.

        """

        dbg( '"QQQQQQQQ" useLsf' )

        genScript = os.environ.get( 'PIPERUN_GEN_SCRIPT' )

        if not isSubPipeline and 'PIPERUN_RULE_FILTER' in os.environ:
            ruleFilter = os.environ[ 'PIPERUN_RULE_FILTER' ]

        mqQueue = os.environ.get( 'MQ_QUEUE' )

        if viewOnly is None: viewOnly = False
        else: viewOnlyFromArgs = False

        if maxpar is None: maxpar = 1
        else: maxparFromArgs = False

        if viewOnlyFromArgs: viewOnly = ( len( sys.argv ) == 1 )
        if forceRun: viewOnly = False
        
        if maxparFromArgs:
            if ( len( sys.argv ) > 1 and 'dbg' in sys.argv[1] ): keepGoing = False
            if ( len( sys.argv ) > 1 and 'mq' in sys.argv[1] ): useMQ = True
            
            maxpar = 1 if ( len( sys.argv ) > 1 and 'dbg' in sys.argv[1] ) else \
                ( maxparLocal if ( len( sys.argv ) > 1 \
                             and ( forceLocal or 'loc' in sys.argv[1]  ) ) else ( mqMaxPar if useMQ else lsfMaxPar ) )

        if lsfAbovePar != None:
            if isinstance( maxparLocal, int ): lsfAbovePar = max( lsfAbovePar, maxparLocal+1 )
            useLsf = ( maxpar >= lsfAbovePar )

        #assert not sanityCheck or not ( maxpar > 10 and not ( useLsf or useMQ ) )

        if self.lsfNeedDirs:
            assert len( self.lsfNeedDirs ) == 1
            assert not lsfNeedDir or lsfNeedDir == self.lsfNeedDirs[0]
            lsfNeedDir = self.lsfNeedDirs[0]

        if lsfJobGroup is None:
            lsfJobGroup = '/ilya/' + os.path.splitext( os.path.basename( sys.argv[0] ) )[0]

        dbg( '"LLLLLL" lsfJobGroup' )

        if fnMap: self.addFileNameMapping( fnMap )

        pipelineFileName = ReserveTmpFileName( prefix = 'PipeRun', suffix = 'pipeline.xml',
                                               tmpDir = genScript or '../Other/tmp' )
        self.write( pipelineFileName )
        dbg( 'self.fileTypeMap pipelineFileName' )
        os.system( 'cp ' + pipelineFileName + ' pipeline.xml' )  # for reading
        cmd = ( ( os.getenv( 'PIPERUN_DIR' ) if os.getenv( 'PIPERUN_DIR' ) else ( ( os.getenv( 'PIPERUN_ROOT' ) if os.getenv( 'PIPERUN_ROOT' ) else '..' ) + '/Operations/Ilya_Operations/PipeRun/binaries/' + SysTypeString() ) ) + '/' + os.environ.get( 'PIPERUN_BINARY', 'RunPipeRun' ) + ' PIPELINE='
                + ( '"{' + ','.join( [ relpath( pipelineFileName ) ] + self.pipelineFilesToInclude ) + '}"' ) )
        logging.info( 'GOT COMMAND %s' % cmd )
        dbg( 'cmd' )
        if viewOnly: cmd += " VIEW_PIPELINE_AND_QUIT=True"
        if not outputDir: outputDir = 'pipeline_' + self.name
        if outputDir: cmd += " OUTPUT_DIR=" + outputDir
        if keepGoing: cmd += " KEEP_GOING=True"
        if maxpar > 1: cmd += " MAXPAR=%d" % maxpar
        if useLsf: cmd += " USE_LSF=True"
        if useMQ: cmd += " USE_MQ=True"
        if mqQueue: cmd += " MQ_QUEUE=" + mqQueue
        if lsfDefaultQueue != None: cmd += " LSF_DEFAULT_QUEUE=" + lsfDefaultQueue
        if lsfShortQueue != None: cmd += " LSF_SHORT_QUEUE=" + lsfShortQueue
        if lsfJobGroup != None: cmd += " LSF_JOB_GROUP=" + lsfJobGroup
        if lsfProject != None: cmd += " LSF_PROJECT=" + lsfProject
        if lsfNeedDir != None: cmd += " LSF_NEED_DIR=" + lsfNeedDir
        if lsfResources != None: cmd += " LSF_RESOURCES=" + lsfResources
        if lsfUserJobPriority is not None: cmd += " LSF_USER_JOB_PRIORITY=%d" % lsfUserJobPriority
        if targetsWithNoRuleOk: cmd += " TARGETS_WITH_NO_RULE_OK=True"
        if not showDescrsInNodes: cmd += " SHOW_DESCR_IN_NODES=False"
        if dryRun: cmd += " DRY_RUN=True"
        if disableHtmlViewers: cmd += " DISABLE_HTML_VIEWERS=True"
        if showRuleStatus: cmd += " SHOW_RULE_STATUS=True"
        if not showLSFRuleStatus: cmd += " SHOW_LSF_RULE_STATUS=False"
        if checkModTimes: cmd += " CHECK_MOD_TIMES=True"
        if aftRuleDelay != None: cmd += " AFT_RULE_DELAY=%d" % aftRuleDelay
        if detailCmdChanges: cmd += " DETAIL_CMD_CHANGES=True"
        if passParentParams: cmd += " PASS_PARENT_PARAMS=True";
        if shell: cmd += " SHELL=" + shell
        if self.pythonCmd != 'python': cmd += " PYTHON_CMD='%s'" % self.pythonCmd
        if ignoreAttrDiffs: cmd += " IGNORE_ATTR_DIFFS='%s'" % ','.join( MakeSeq( ignoreAttrDiffs ) )
        
        cmd += " TOP_SCRIPT='%s'" % ( topScript or sys.argv[0] )

        useNQ = bool( int( os.environ.get( 'PIPERUN_USENQ', '0' ) ) )

        if useMQ and not genScript and not useNQ:
            queueSpec = () if not mqQueue else ( '-Q', mqQueue )
            if not os.getenv( 'PIPERUN_NO_OWN_PUSHER' ): StartPushers( use_args = ( '-S', '40' ) + queueSpec )
            memSpec = ()
            if os.getenv( 'PIPERUN_MEM' ): memSpec = ( '--max-mem', os.getenv( 'PIPERUN_MEM' ) )
            if not os.getenv( 'PIPERUN_NO_OWN_RUNNER' ):
                extraArgs = ()
                if not os.getenv( 'PIPERUN_RUN_ANY' ): extraArgs += ( '-i', self.makePipelineId() )
                if not os.getenv( 'PIPERUN_NOTLOCAL' ): extraArgs += ( '-l', )
                StartRunners( use_args = ( '-M', 20000, '-F', 20000 ) + memSpec + queueSpec + extraArgs )
        
        tmpFiles = []
        def addOpt( optName, optVal, tmpFiles ):
            """Add an option the value of which may need to be passed through a tempfile"""

            if not optVal: return ''
            if optVal.find('"') == -1  and  optVal.find("'") == -1:
                opt = '"' + optVal + '"'
            else:
                tmp_fd, tmp_name = tempfile.mkstemp( dir = genScript or '../Other/tmp',
                                                     prefix = 'tmp_prun_', suffix='.txt', text = True)
                f = os.fdopen(tmp_fd, 'w')
                f.write( optVal )
                f.close()
                tmpFiles.append( tmp_name )
                opt = 'file:' + tmp_name
            return ' ' + optName + '=' + opt

        for optName, optVal in ( ( 'RULE_FILTER', ruleFilter ), ( 'FORCE_REMAKE', forceRemake ),
                                 ( 'ACCEPT_EXISTING_OUTPUTS', acceptExistingOutputs ),
                                 ( 'IGNORE_CODE_CHANGES_FOR', ignoreCodeChangesFor ),
                                 ( 'IGNORE_CMD_CHANGES_FOR', ignoreCmdChangesFor ),
                                 ( 'IGNORE_NEW_DEPENDENCIES_FOR', ignoreNewDependenciesFor ),
                                 ( 'NO_VIEW_FOR', noViewFor ), ( 'UPDATE_FETCHED', updateFetched ) ):
            cmd += addOpt( optName, optVal, tmpFiles )

        cmd += ' 1>&2'
        logging.debug( 'PR command: ' + cmd )

        if genScript:
            runScriptFN = os.path.join( genScript, 'run.sh' )
            with open( runScriptFN, 'w' ) as out:
                out.write( '#!/usr/bin/env bash\n' )
                out.write( 'set -e\n' )
                out.write( 'set -o pipefail\n' )
                out.write( cmd );
            os.chmod( runScriptFN, stat.S_IXUSR | stat.S_IRWXU )

        else:
        
            try: SystemSucceed( cmd )
            finally:
                logging.debug( 'Finished PR command: ' + cmd )

                if 'DEBUG_PIPERUN' not in os.environ:
                    for f in tmpFiles: os.remove( f )
                    os.remove( pipelineFileName )

                if useMQ and not genScript:
                    if not os.getenv( 'PIPERUN_NO_OWN_RUNNER' ): StopRunners()
                    if not os.getenv( 'PIPERUN_NO_OWN_PUSHER' ): StopPushers()
        
    ####################################
    # private methods of class PipeRun #       
    ####################################

    def __addElt(self, eltName, addTo = None, addTxt = None ):
        """Add a blank XML element of given type to the
        pipeline descr, and return the element.

        Parameters:

           eltName - XML element type
           addTo - the parent element (created earlier by self.doc.createElement() )
           addTxt - the text to put into the body of the element (i.e. between the
             <eltName> and <eltName/> tags).

        """

        if addTo == None: addTo = self.pipeline

        assert isinstance( eltName, StringType )
        elt = self.doc.createElement(eltName)
        if addTxt != None:
            elt.appendChild( self.doc.createTextNode( str( addTxt ) ) )
        return addTo.appendChild(elt)

    def __stringsToXML(self, eltName, eltChildName, strings, addTo):
        """Convert a list of strings to an XML element containing these
        strings.

        Parameters:

           eltName - name of the outer XML element containing the strings.
           eltChildName - name of each child of <eltName>, each containing
              one string in its body
           strings - the list of strings, or one string.
           addTo - the XML node to which to add the resulting <eltName>.

        """
        assert isinstance( eltName, StringType )
        stringsElt = self.doc.createElement( eltName )
        if isinstance( strings, StringType ): strings = [ strings ]
        for s in strings:
            assert isinstance( s, StringType )
            self.__addElt( eltChildName, stringsElt, s )

        return addTo.appendChild( stringsElt )               

    def __filesToXML( self, eltName, files, addTo ):
        """Convert a list of files to an XML representation.

        Parameters:

           eltName - the name of the outer XML element containing
              <file> elements representing the individual files
              in the list
           files - the list of files in the list, as strings.
           addTo - the XML node to which the file list will be
              added as a child.

        """
        return self.__stringsToXML( eltName, "file", files, addTo )

    __programBinarySuffix = 'zZz_PLATFORM'

    def getAllTargets(self):
        """Return all targets of all rules, as an iterable"""

        return self.allTargets



########################
# end of class PipeRun #
########################

def FunctionString( fn ):
    """Return a string representing the code of the function
    and any functions that it calls.  A checksum of the string
    can tell us whether the function has materially changed.
    """

    # allow dependencies on a class; if anything in 
    if isinstance( fn, (types.ClassType, types.ModuleType) ): return inspect.getsource( fn )
    if isinstance( fn, types.MethodType ): fn = fn.im_func

    assert isinstance( fn, types.FunctionType )
    code = fn.func_code
    result = ''.join( FunctionChecksum( code, fn ) )

    # Because identifiers in Python are resolved dynamically at runtime,
    # it's hard to automatically identify from a function's source code
    # all the functions that it calls.  We therefore rely on the user
    # to manually tell us the callees of a given function, by
    # giving a special 'calls' argument and listing the callees in its
    # default value.
    
    argspec = inspect.getargspec( fn )
    args, defaults = argspec[0], argspec[3]
    if defaults and 'calls' in args[-len(defaults):] :
        for called in defaults[ args[-len(defaults):].index('calls') ]:
            result += FunctionString( called )

    return result        
    

def GetCreates( fn, *args, **kwargs ):
    """Return the creates filename(s)"""
    return PipeRun.FileTuple( MakeSeq( fn( getio = True, *args, **kwargs )[ 'creates' ] ) )

def GetDependsOn( fn, *args, **kwargs ):
    """Return the depends_on filename(s)"""
    return PipeRun.FileTuple( MakeSeq( DictGet( fn( getio = True, *args, **kwargs ), 'depends_on', () ) ) )

def GetAttrs( fn, *args, **kwargs ):
    """Return the depends_on filename(s)"""
    return tuple( DictGet( fn( getio = True, *args, **kwargs ), 'attrs', {} ))

def GetIO( fn, *args, **kwargs ):
    """Return the depends_on filename(s)"""
    ioInfo = fn( getio = True, *args, **kwargs )
    if 'attrs' not in ioInfo: ioInfo[ 'attrs' ] = {}
    if 'depends_on' not in ioInfo: ioInfo[ 'depends_on' ] = ()
    ioInfo[ 'depends_on' ] = MakeSeq( ioInfo[ 'depends_on' ] )
    ioInfo[ 'creates' ] = MakeSeq( ioInfo[ 'creates' ] )
    
    return ioInfo

#pr = PipeRun()

#pr.addPreExistingFiles( [ 's1', 's2' ] ) 
#pr.addRule( [ 't1', 't2' ], [ 's1', 's2' ], [ 'cmd' ], 'cmnt' )

#pr.addRule( [ "ls.out" ], [], [ "touch ls.out" ], "create a file" )

#pr.run( maxpar = 3 )


