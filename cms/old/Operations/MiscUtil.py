"""Miscellaneous general-purpose utilities, that are not specific to our particular datasets.

For general utilities related to our datasets, see Operations.Ilya_Operations.DataUtil.py .
"""


import platform
from functools import reduce
assert platform.python_version_tuple() >= ( 2, 6 )

import sys, os, subprocess, operator, itertools, functools, inspect, logging, traceback, shutil, math, copy, types, \
    atexit, time, glob, numbers, pickle, base64, collections, re, tempfile, string, gzip, bz2, socket, platform, random, \
    bisect, zipfile, tarfile, errno, zlib, hashlib

#if 'MISCUTIL_NO_MATPLOTLIB' not in os.environ:
#  try:
#    import matplotlib
#    matplotlib.use('agg')
#  except ImportError: pass

try:
    import numpy as np
    haveNumpy = True
    from Classes.DotData import DotData
except ImportError:
    haveNumpy = False

try:
    import fcntl
    have_fcntl = True
except ImportError:
    have_fcntl = False

try:
  from Bio import SeqIO
  haveBio = True
except ImportError:
  haveBio = False

from subprocess import Popen, PIPE
from contextlib import contextmanager, nested

def dump_args(func, fname = None):
    "This decorator dumps out the arguments passed to a function before calling it"
    argnames = func.__code__.co_varnames[:func.__code__.co_argcount]
    if not fname: fname = func.__name__
    def echo_func(*args,**kwargs):
        callStr = fname + ":" + ', '.join(
            '%s=%r' % entry
            for entry in list(zip(argnames,args)) + list(kwargs.items()))
        logging.info( 'calling ' + str(callStr) )
        retVal = func(*args, **kwargs)
        logging.info( 'from ' + str(callStr) + ' returned ' + str(retVal) )
        return retVal
    return echo_func

#######################################################################

ioMap = {}
ioMapActive = False
file2chunk = {}

def AddPathMapping( fromFile, toFile ):
    """Record a mapping"""
    dirName, fileName = os.path.split( fromFile )
    if dirName and not fileName: fromFile = dirName
    ioMap[ fromFile ] = toFile

def AddFileView( fileName, beg, end, **tableReadOpts ):
  """Restrict the specified file to a given view."""
  file2chunk[ fileName ] = ( beg, end, tableReadOpts )

@dump_args    
def MapPath( fn ):
    """Map a filename as specified by the current remapping."""
    logging.info( 'ioMapActive=' + str(ioMapActive) + ' ioMap=' + str(ioMap) )
    if not ioMapActive: return fn
    if fn in ioMap: return ioMap[ fn ]
    dirName, fileName = os.path.split( fn )
    if dirName in ioMap: return os.path.join( ioMap[ dirName ], fileName )
    return fn

# var: orig_fns - map from ( module, function_name ) to the original definition
#    of that function in that module, before we replaced it.
orig_fns = {}

def GetOrigFn( module, fnName ):
    """Return the original function."""
    return DictGet( orig_fns, ( module, fnName ), getattr( module, fnName ) )

class NonClosingStream(object):
    """Wraps a standard stream with a decorator that ignores attempts to
    close the stream.
    """

    def __init__(self, stdstream):
        self.stdstream = stdstream

    def close(self):
        """Prevent the closing of the stream"""
        logging.info( 'NOT closing' )
        pass
    
    def __getattr__(self, attr):
        """Delegate everything but close to the stream"""

        logging.info( 'delegating %s' % attr )
        return getattr(self.stdstream, attr)

class FileChunk(object):
    """A decorator on a file-like object, that filters out all but a specified window of the file.
    Preserves the comments, headers and line structure."""

    def __init__(self, fileObj, beg, end, fileName = None, **tableReadOpts ):

        useHeadings = tableReadOpts.get( 'useHeadings', True )
        commentPrefix = tableReadOpts.get( 'commentPrefix', '#' )

        if beg is None or beg < 0: beg = 0
        if end is None: end = -1
        
        self.fileObj = fileObj
        self.fileName = fileName

        fileObj.seek( 0 )
        self.headerLines = []

        dataStart = 0
        if commentPrefix or useHeadings:
          while True:
              line = fileObj.readline()
              isComment = commentPrefix and line.startswith( commentPrefix )
              if isComment or useHeadings:
                self.headerLines.append( line )
                dataStart = fileObj.tell()
              if not isComment: break

        self.beg = max( beg, dataStart )
        fileObj.seek( self.beg )
        self.end = end

        self.did_readline = False
        self.did_read = False

    def readline( self ):
        assert not self.did_read
        self.did_readline = True
        if self.headerLines:
          result = self.headerLines.pop( 0 )
          logging.info( "rrreadline_returning_header %s" % result )
          return result
        if self.end >= 0 and self.fileObj.tell() >= self.end:
          logging.info( "rrreadline_returning_EOF %s" % self.end )
          return ''
        result = self.fileObj.readline()
        return result

    def read( self ):

        assert not ( self.did_read or self.did_readline )
        self.did_read = True

        if self.end < 0: mainChunk = self.fileObj.read()
        else:
          readSize = self.end - self.beg
          logging.info( 'read file chunk: %s:%d-%d; pos=%d, reading=%d' % ( self.fileName, self.beg, self.end,
                                                                            self.fileObj.tell(), readSize) )
          mainChunk = self.fileObj.read( readSize ) if readSize > 0 else ''

          logging.info( 'read file chunk: %s:%d-%d' % ( self.fileName, self.beg, self.end ) )

        return ''.join( self.headerLines ) + mainChunk

    def close( self ): self.fileObj.close()

    def __iter__( self ):
        while True:
            line = self.readline()
            if not line: break
            yield line

    def __enter__( self ): return self
    def __exit__( self, *excinfo): self.close()

def remap( orig_func, n_args = 1, module = None ):
    """Decorate the specified function so that it remaps the input files as specified by ioMap."""


    @functools.wraps( orig_func )
    def newFunc( *args, **kwargs ):
        newPaths = list(map( MapPath, args[ :n_args ] ))
        if str( newPaths ) != str( args[ :n_args ] ):
          logging.info( 'MiscUtil.remap: applying ' + orig_func.__name__ + ' to ' + ','.join( newPaths ) + ' instead of ' + ','.join( args[ :n_args ] ) )
        specialPaths = dict( stdin = sys.stdin, stdout = sys.stdout )

        if n_args == 1 and newPaths[0] in specialPaths:
          assert orig_func == orig_open
          result = NonClosingStream( specialPaths[ newPath ] )
        else:
          result = orig_func( *( tuple ( newPaths ) + tuple( args[n_args:] ) ), **kwargs )

        if orig_func == orig_open and newPaths[0] in file2chunk:
          beg, end, tableReadOpts = file2chunk[ newPaths[0] ]
          if newPaths[0].endswith( '.vcf' ) or newPaths[0].endswith( '.vcf.gz' ):
            tableReadOpts[ 'commentPrefix' ] = '##'
          result = FileChunk( fileObj = result, beg = beg, end = end, fileName = newPaths[0], **tableReadOpts )

        return result

    newFunc.__doc__ += '\n\nMapped by ioremap module.'

    orig_module = module or inspect.getmodule( orig_func )

    orig_fns[ ( orig_module, orig_func.__name__ ) ] = newFunc
    setattr( orig_module, orig_func.__name__, newFunc )

def InitPathMapping():
    """Give effect to path mappings specified by AddPathMapping()."""

    if ioMap or file2chunk:
  
      try:
          import posix

          dummy = list(map( functools.partial( remap, module = posix ),
               ( posix.access, posix.chdir, posix.chmod, posix.chown, posix.listdir, posix.link, posix.lstat,
                 posix.stat, posix.mknod, posix.mkdir, posix.pathconf, posix.remove, posix.readlink,
                 posix.rmdir, posix.stat, posix.statvfs, posix.unlink, posix.utime ) ))
          
          remap( posix.rename, module = posix, n_args = 2 )
          logging.info( 'REMAPPED os.rename' )
          
      except ImportError: pass

      try:
          import posixpath

          dummy = list(map( remap, ( posixpath.walk, ) ))
      except ImportError: pass

      global orig_open
      orig_open = open

      remap( open )
      
      dummy = list(map( functools.partial( remap, module = os ),
           ( os.access, os.chdir, os.chmod, os.chown, os.listdir, os.link, os.lstat,
             os.stat, os.mknod, os.mkdir, os.makedirs, os.open, os.pathconf, os.remove, os.readlink, os.removedirs,
             os.rmdir, os.stat, os.statvfs, os.unlink, os.utime, os.walk, os.path.walk,
             os.path.exists ) ))

      remap( os.rename, module = os, n_args = 2 )
      logging.info( 'REMAPPED os.rename' )

      global ioMapActive
      ioMapActive = True


    
#DoRemap()
    
#######################################################################


def SysTypeString( ):
    """Returns a string identifying the processor type and OS type.
    Useful for choosing among platform-specific binaries to execute.
    """

    if sys.platform == 'darwin':
        machType, osType = 'i386', 'darwin'
    else:    
        machType = platform.uname()[4]
        if machType.startswith('x86_64'): machType = 'x86_64'
        osType = platform.system().lower()
        if osType.startswith( 'darwin' ): osType = 'darwin'

    # linuxType = ''
    # if platform.system() == 'Linux':
    #   linuxType += '-' + MakeAlphaNum( platform.linux_distribution()[0].strip() )
        
    return machType + '-' + osType


def delete(ToDelete):
	'''
	unified "strong" version of delete that uses os.remove for a file 
	and shutil.rmtree for a directory tree
	'''
	if os.path.isfile(ToDelete):
		os.remove(ToDelete)
	elif os.path.isdir(ToDelete):
		shutil.rmtree(ToDelete)		

def MakeDir(DirName,creates = ()):
	'''
	is a "strong" directory maker -- if DirName already exists, this deletes it first 
	'''
        DirName = MapPath( DirName )
	if os.path.exists(DirName): delete(DirName)
	os.mkdir(DirName)


def SlurpFile( fname ):
    """Read entire file into one string"""
    if not os.path.isfile( fname ): raise IOError( 'File not found: %s' % fname )
    with open( fname ) as f:
        return f.read()

def SlurpFileLines( fname ):
    """Read entire file and return a newly allocated array of its lines"""
    return SlurpFile( fname ).strip( '\n' ).split( '\n' )


def NthLineOfFile( fname, n = 0 ):
    """Return the nth line of file (zero-indexed)"""
    with open( fname ) as f:
      while n > 0:
        f.readline()
        n -= 1
      return f.readline().strip()
    
FirstLineOfFile = NthLineOfFile
SecondLineOfFile = functools.partial( NthLineOfFile, n = 1 )

def LastLineOfFile( fname ):
  """Returns the last line of the file"""
  with open( fname ) as f:
    for line in f:
      strippedLine = line.strip()
      if strippedLine:
        lastLine = strippedLine
  return lastLine
  
    
def NumLinesInFile( fname ):
    cnt = 0
    with open( fname ) as f:
        for line in f:
            if line.strip(): cnt += 1
    return cnt

def OutputOfCmd( cmd, exitCodesOk = ( 0, ) ):
    """Run a shell command and return its output as a string"""
    #return Popen( cmd.split(), stdout = PIPE ).communicate()[0]

    tmpFN = ReserveTmpFileName()
    SystemSucceed( cmd + ' >& ' + tmpFN, exitCodesOk = exitCodesOk )
    time.sleep(3)
    WaitForFileToAppear( tmpFN )
    return SlurpFile( tmpFN )

def DumpFile( fname, contents, dbg = False, getio = None ):
    """Write out a file whose entire contents is the given string"""
    if getio: return dict( depends_on = (), creates = fname, attrs = dict( piperun_short = True ) )
    logging.info( 'Writing ' + fname  )
    with open( fname, 'wt' ) as f:
        f.write( str( contents ) )
    logging.info( 'Wrote ' + fname  )

def DumpFileIfChanged( fname, contents ):
    """Write out a file whose entire contents is the given string; if
    file already contains this exact string, its modtime is left unchanged. 
    """
    try:
        if os.path.isfile( fname ) and SlurpFile( fname ) == contents: return
    except:
        dbg( 'DumpFileIfChanged: error loading existing file' )
        
    DumpFile( fname, contents )


def CmdOutput( cmd ):
    """Run a shell command, capture its entire output to a string, and return the string.
    If the command does not complete successfully, raises an exception.
    """
    p = subprocess.Popen( cmd, shell = True, stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT, universal_newlines = True )
    output = p.communicate()[0]
    if p.returncode != 0: raise subprocess.CalledProcessError(p.returncode,cmd)
    return output


def SystemSucceed( cmd, dbg = False, exitCodesOk = ( 0, ) ):
    """Run a shell command, and raise a fatal error if the command fails."""
    logging.info( 'Running command ' + cmd + ' ; called from ' + sys._getframe(1).f_code.co_filename + ':' +
                  str( sys._getframe(1).f_lineno ) )
    if int( os.environ.get( 'DISABLE_SYSTEM_SUCCEED', '0' ) ):
        logging.info( 'SKIPPING COMMAND ' + cmd )
        exitCode = 0
    else:
        exitCode = os.system( cmd )
    logging.info( 'Finished command ' + cmd + ' with exit code ' + str( exitCode ) )
    if exitCodesOk != 'any' and exitCode not in exitCodesOk:
        raise IOError( "Command %s failed with exit code %d" % ( cmd, exitCode ) )
    return exitCode


def joinstr( sep, *args ):
    """Join the arguments after converting them to a string"""
    return sep.join( map( str, args ) )

def tabjoin( *args, **kwargs ):
    """Join args by tab"""
    sep = DictGet( kwargs, 'sep', '\t' )
    return sep.join( map( str, args ) )

def tabwriten( f, *args, **kwargs ):
    """Write a tab-separated line to a file"""
    f.write( tabjoin( *args, **kwargs ) )

def tabwrite( f, *args, **kwargs ):
    """Write a tab-separated line to a file, followed by a newline"""
    tabwriten( f, *args, **kwargs )
    f.write( '\n' )


def simple_decorator(decorator):
    """This decorator can be used to turn simple functions
    into well-behaved decorators, so long as the decorators
    are fairly simple. If a decorator expects a function and
    returns a function (no descriptors), and if it doesn't
    modify function attributes or docstring, then it is
    eligible to use this. Simply apply @simple_decorator to
    your decorator and it will automatically preserve the
    docstring and function attributes of functions to which
    it is applied.

    Taken from http://wiki.python.org/moin/PythonDecoratorLibrary

    """
    def new_decorator(f):
        g = decorator(f)
        g.__name__ = f.__name__
        g.__doc__ = f.__doc__
        g.__dict__.update(f.__dict__)
        return g
    # Now a few lines needed to make simple_decorator itself
    # be a well-behaved decorator.
    new_decorator.__name__ = decorator.__name__
    new_decorator.__doc__ = decorator.__doc__
    new_decorator.__dict__.update(decorator.__dict__)
    return new_decorator


def ApplyToResult( func ):
    """Creates a decorator that applies func to the result of the decorated function.
    """

    @simple_decorator
    def wrap( f ):
        def new_function(*args, **kw):
            return func( f( *args, **kw ) )
        return new_function
        
    return wrap

def func_sig(func):

    """func_sig(func)

       Returns the signature of a Python function/method as string.
       Keyword initializers are also shown using
       repr(). Representations longer than 100 bytes are truncated.

       XXX Anonymous arguments ((a,b,c)=(1,2,3)) are not supported and
           probably never will be since they require disassembling the
           byte code which is bound to fail once byte code optimizers find
           their way into every Pythoneers home...


       Adapted from http://www.koders.com/python/fid3235F0CCCAF64A3B9CEFB8E89BAF352823A9B3FC.aspx?s=cdef%3Atimer#L375 and
       http://www.faqts.com/knowledge_base/view.phtml/aid/5666 .

    """
    if hasattr(func,'im_func'):
        # func is a method
        func = func.__func__
    code = func.__code__
    fname = code.co_name
    callargs = code.co_argcount
    # XXX Uses hard coded values taken from Include/compile.h
    args = list(code.co_varnames[:callargs])
    if func.__defaults__:
        i = len(args) - len(func.__defaults__)
        for default in func.__defaults__:
            try:
                r = repr(default)
            except:
                r = '<repr-error>'
            if len(r) > 100:
                r = r[:100] + '...'
            arg = args[i]
            if arg[0] == '.':
                # anonymous arguments
                arg = '(...)'
            args[i] = '%s=%s' % (arg,r)
            i = i + 1
    if code.co_flags & 0x0004: # CO_VARARGS
        args.append('*'+code.co_varnames[callargs])
        callargs = callargs + 1
    if code.co_flags & 0x0008: # CO_VARKEYWORDS
        args.append('**'+code.co_varnames[callargs])
        callargs = callargs + 1
    return '%s(%s)' % (fname,string.join(args,', '))

def flatten(*args):
    """flatten(sequence) -> tuple

    Returns a single, flat tuple which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, tuple((8,9,10))])
    (1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10)

    Taken from http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks

    """

    x = args[0] if len( args ) == 1 else args

    result = []
    for el in x:
        #if isinstance(el, (list, tuple)):
        if IsSeq( el ): result.extend(flatten(el))
        else: result.append(el)
    return tuple( result )

def obj2path( obj, braceDepth = 0, parenDepth = 0 ):
    """Encode an object as a pathname."""

    if isinstance( obj, collections.Mapping ):
        return os.path.join( '_BRACE%d' % braceDepth,
                             os.path.join(*[  os.path.join( MakeAlphaNum( str( k ) ),
                                                            obj2path( obj = v, braceDepth = braceDepth+1,
                                                                      parenDepth = parenDepth ) )
                                              for k, v in sorted( list( obj.items() ) ) if v ]),
                             'BRACE%d_' % braceDepth )
    elif IsSeq( obj ):
        return os.path.join( '_PAREN%d' % parenDepth,
                             os.path.join(*[ obj2path( obj = v, braceDepth = braceDepth,
                                                       parenDepth = parenDepth+1 )
                                             for v in sorted( obj ) ]),
                             'PAREN%d_' % parenDepth )
    else: return MakeAlphaNum( str( obj ) )

def Dict( vars, frame = None, **kwargs ):
    """Construct a dictionary mapping each var to its value.  Var names are
    given in one string, separated by whitespace.  Useful for passing specified arguments through to a function you
    call.

    Has the same effect as invoking RestrictDict( MergeDicts( globals(), locals() ), vars.split() ),
    where globals() and locals() refer to names visible at the point where Dict() is invoked.
    """
    if frame is None: frame = inspect.currentframe().f_back
    a_map = frame.f_globals
    a_map.update(frame.f_locals)
    return MergeDicts( dict([ (v, a_map[v]) for v in vars.split() ]), kwargs )

tmpFiles = []

@atexit.register
def RemoveTmpFiles():
    """Remove any tempfiles"""
    
    if 'DEBUG_PIPERUN' not in os.environ:
        for f in tmpFiles:
            if os.path.exists( f ):
                logging.info( 'removing tmpfile ' + f )
                os.remove( f )
    else:
        if tmpFiles: logging.info( 'preserving tmpfiles:\n' + '\n'.join( tmpFiles ) )

def RegisterTmpFile( fname ):
    """Remember a file to be deleted at exit"""
    tmpFiles.append( fname )

def ReserveTmpFileName( prefix = 'mytmp', suffix = '.tmp', text = True, tmpDir = '../Other/tmp' ):
    """Create a unique temp file, close it, and return its name.
    The caller then would typically overwrite this file,
    being assured that no other ReserveTmpFileName() call can return
    the same filename to another process or thread.
    The file is automatically deleted at program exit.
    """
    EnsureDirExists( tmpDir )
    fileDescr, fileName = tempfile.mkstemp( prefix = prefix, suffix = suffix, dir = tmpDir, text = text )
    os.close( fileDescr )
    WaitForFileToAppear( fileName )
    tmpFiles.append( fileName )
    return fileName


def EnsureDirExists( d ):
    """Create the given dir (and any parent dirs) if it does not already exist"""
    if not os.path.isdir( d ): os.makedirs( d )


def WaitForFileToAppear( fname ):
    """Wait until the given file appears, then return."""

    caller = inspect.currentframe()
    thisFileName = caller.f_code.co_filename
    while ( caller.f_code.co_filename == thisFileName ):
        caller = caller.f_back
    
    SystemSucceed( "perl ../Operations/Ilya_Operations/PipeRun/wait_for_file_to_appear.pl '%s' 1>&2" % MapPath( Str( fname, caller ) ) )
    del caller

class Struct(dict):

    def __init__(self, *args, **kwargs):
        assert not args
        for k, v in list(kwargs.items()):
            setattr(self, k, v)

    def __getattr__(self,name):

        try:
            val=self[name]
        except KeyError:
            val=getattr(super(Struct,self), name)
            
        return val

    def __setattr__(self,name,val):
        self[name]=val

from operator import itemgetter as _itemgetter
from keyword import iskeyword as _iskeyword
import sys as _sys
        
def namedtuple(typename, field_names, verbose=False):
    """Returns a new subclass of tuple with named fields.

    >>> Point = namedtuple('Point', 'x y')
    >>> Point.__doc__                   # docstring for the new class
    'Point(x, y)'
    >>> p = Point(11, y=22)             # instantiate with positional args or keywords
    >>> p[0] + p[1]                     # indexable like a plain tuple
    33
    >>> x, y = p                        # unpack like a regular tuple
    >>> x, y
    (11, 22)
    >>> p.x + p.y                       # fields also accessable by name
    33
    >>> d = p._asdict()                 # convert to a dictionary
    >>> d['x']
    11
    >>> Point(**d)                      # convert from a dictionary
    Point(x=11, y=22)
    >>> p._replace(x=100)               # _replace() is like str.replace() but targets named fields
    Point(x=100, y=22)

    """

    # Parse and validate the field names.  Validation serves two purposes,
    # generating informative error messages and preventing template injection attacks.
    if isinstance(field_names, str):
        field_names = field_names.replace(',', ' ').split() # names separated by whitespace and/or commas
    field_names = tuple(field_names)
    for name in (typename,) + field_names:
        if not min(c.isalnum() or c=='_' for c in name):
            raise ValueError('Type names and field names can only contain alphanumeric characters and underscores: %r' % name)
        if _iskeyword(name):
            raise ValueError('Type names and field names cannot be a keyword: %r' % name)
        if name[0].isdigit():
            raise ValueError('Type names and field names cannot start with a number: %r' % name)
    seen_names = set()
    for name in field_names:
        if name.startswith('_'):
            raise ValueError('Field names cannot start with an underscore: %r' % name)
        if name in seen_names:
            raise ValueError('Encountered duplicate field name: %r' % name)
        seen_names.add(name)

    # Create and fill-in the class template
    numfields = len(field_names)
    argtxt = repr(field_names).replace("'", "")[1:-1]   # tuple repr without parens or quotes
    reprtxt = ', '.join('%s=%%r' % name for name in field_names)
    dicttxt = ', '.join('%r: t[%d]' % (name, pos) for pos, name in enumerate(field_names))
    template = '''class %(typename)s(tuple):
        '%(typename)s(%(argtxt)s)' \n
        __slots__ = () \n
        _fields = %(field_names)r \n
        def __new__(cls, %(argtxt)s):
            return tuple.__new__(cls, (%(argtxt)s)) \n
        @classmethod
        def _make(cls, iterable, new=tuple.__new__, len=len):
            'Make a new %(typename)s object from a sequence or iterable'
            result = new(cls, iterable)
            if len(result) != %(numfields)d:
                raise TypeError('Expected %(numfields)d arguments, got %%d' %% len(result))
            return result \n
        def __repr__(self):
            return '%(typename)s(%(reprtxt)s)' %% self \n
        def _asdict(t):
            'Return a new dict which maps field names to their values'
            return {%(dicttxt)s} \n
        def _replace(self, **kwds):
            'Return a new %(typename)s object replacing specified fields with new values'
            result = self._make(map(kwds.pop, %(field_names)r, self))
            if kwds:
                raise ValueError('Got unexpected field names: %%r' %% kwds.keys())
            return result \n\n''' % locals()
    for i, name in enumerate(field_names):
        template += '        %s = property(itemgetter(%d))\n' % (name, i)
    if verbose:
        logging.info( template )

    # Execute the template string in a temporary namespace
    namespace = dict(itemgetter=_itemgetter)
    try:
        exec(template, namespace)
    except SyntaxError as e:
        raise SyntaxError(e.message + ':\n' + template)
    result = namespace[typename]

    # For pickling to work, the __module__ variable needs to be set to the frame
    # where the named tuple is created.  Bypass this step in enviroments where
    # sys._getframe is not defined (Jython for example).
    if hasattr(_sys, '_getframe'):
        result.__module__ = _sys._getframe(1).f_globals['__name__']

    return result


        
def MakeAlphaNum( str ):
    """Return a version of the argument string, in which all non-alphanumeric chars have been replaced
    by underscores.
    """
    return re.sub( '\W+', '_', str )

def IsAlphaNum( str ):
    """Check whether the string consists of only alphanumeric chars or underscores"""
    return str.replace( '_', 'A' ).isalnum()

class AtomicForIsSeq(object):
    """Instances of classes derived from this class will be called _not_ sequences
    by IsSeq(), even if they otherwise look like a sequence."""
    pass

def IsSeq( val ):
    """Test if the value is a sequence but not a string or a dict.
    Derive your class from AtomicForIsSeq to force its instances to be called
    _not_ sequences by this function, regardless of anything else.
    """
    return ( isinstance( val, ( collections.Sequence, types.GeneratorType ) ) or
             ( hasattr( val, '__getitem__' ) and hasattr( val, '__len__' ) ) )  \
             and not isinstance( val, ( (str,), collections.Mapping, AtomicForIsSeq ) )

def MakeSeq( val ):
    """If val is a sequence (tuple or list) then return val, else return a new singleton tuple
    containing val."""
    return val if IsSeq( val ) else ( val, )

def IsScalar( val ):
    """Check if the value is a scalar value"""
    return isinstance( val, ( numbers.Number, (str,) ) )

def ToString( s ):
    """Convert s to string; handle numpy.charray correctly."""
    
    if haveNumpy and isinstance( s, np.chararray ):
        # a numpy.chararray, when converted to a string, looks like
        # ['string_contents']
        # we want just the actual string contents, so that we get
        # the original simple string to work with.
        return str( s )[2:-2]
    else:
        return str( s )

def Str( s, frame = None, **substs ):
    """Perl-style variable interpolation"""

    s = ToString( s )
    if not frame: frame = inspect.currentframe().f_back
    a_map = frame.f_globals
    a_map.update(frame.f_locals)
    del frame
    return string.Template( s ).substitute(a_map, **substs)

def MergeDicts( *dicts ):
    """Construct a merged dictionary from the given dicts.
    If two dicts define the same key, the key from the dict later in the list is chosen."""
    return dict( reduce( operator.add, list(map( dict.items, dicts )), [] ) )

def RestrictDict( aDict, restrictSet ):
    """Return a dict which has the mappings from the original dict only for keys in the given set"""
    restrictSet = frozenset( restrictSet )
    return dict( item for item in list(aDict.items()) if item[0] in restrictSet )


def DictExcept( aDict, exceptSet ):
    """Return a dict which has the mappings from the original dict except for keys in the given set"""
    if not isinstance( exceptSet, collections.Set ):
        exceptSet = frozenset( MakeSeq( exceptSet ) )
    return dict( item for item in list(aDict.items()) if item[0] not in exceptSet )


def RestrictDictValues( aDict, restrictSet ):
    """Return a dict which has the mappings from the original dict only for values in the given set"""
    return dict( item for item in list(aDict.items()) if item[1] in restrictSet )


def InvertDict( aDict ):
    """Return an inverse mapping for the given dictionary"""
    assert len( set( aDict.values() ) )  == len( aDict )
    return dict( ( v, k ) for k, v in list(aDict.items()) )

def DictGet( d, k, dflt ):
    """Get a value from the dictionary; if not such value, return the specified default"""
    return d.get( k, dflt )

def DictGetNotNone( d, k, dflt ):
    """Get a value from the dictionary; if not such value, return the specified default"""

    return d[k] if k in d and d[k] is not None else dflt
    
class AttrDict(dict):

    """A dictionary that lets you access its keys using Python
    attribute notation.

    From http://parand.com/say/index.php/2008/10/13/access-python-dictionary-keys-as-properties/
    """
    def __getattr__(self, name):
        try:
            return self.__getitem__(name)
        except KeyError:
            return super(type(self),self).__getattribute__(name)
        
    
def cross(*sequences):
    """Compute the cross product of a sequence of sequences.  Returns an iterator yielding the cross-product
    elements.
    
    >>> tuple(cross(('a','b'),(1,2)))
    (('a', 1), ('a', 2), ('b', 1), ('b', 2))
    
    Can be used for compactly writing nested loops, as well as loops where
    the number of nested loops varies dynamically.

    Taken from the "Python recipes" site.
    """
    # visualize an odometer, with "wheels" displaying "digits"...:
    wheels = list(map(iter, sequences))
    digits = [next(it) for it in wheels]
    while True:
        yield tuple(digits)
        for i in range(len(digits)-1, -1, -1):
            try:
                digits[i] = next(wheels[i])
                break
            except StopIteration:
                wheels[i] = iter(sequences[i])
                digits[i] = next(wheels[i])
        else:
            break

logging.basicConfig( level = logging.DEBUG, format='%(process)d %(asctime)s %(levelname)-8s %(filename)s:%(lineno)s %(message)s' )

def dbgnt( *exprs ):
    """Print the values of the given items.  Each item is either a string of Python expressions
    separated by whitespace (then the expressions themselves may not contain whitespace),
    or a string starting with a hash (#) and followed by a single Python expression which may
    contain whitespace.   Free variables in the expressions are interpreted as they would be
    in the caller.
    """
    caller = inspect.currentframe()
    thisFileName = caller.f_code.co_filename
    while ( caller.f_code.co_filename == thisFileName ):
        caller = caller.f_back
            

    caller_filename, caller_lineno, caller_function = caller.f_code.co_filename, caller.f_lineno, caller.f_code.co_name

    callerLocals = caller.f_locals
    callerGlobals = caller.f_globals
    del caller

    logging.debug( 'at ' + caller_filename + ':' + str(caller_lineno) + ' in ' + caller_function + ', ' +
                   ', '.join([ e + ' = ' + str(eval( e, callerGlobals, callerLocals )) for ex in exprs for e in ( ex.strip().split(' ') if not ex.startswith('#') else (ex[1:],) ) ] ) )

    del callerGlobals, callerLocals


def dbg( *exprs, **kwargs ):
    """Print the values and types of the given items.  Each item is either a string of Python expressions
    separated by whitespace (then the expressions themselves may not contain whitespace),
    or a string starting with a hash (#) and followed by a single Python expression which may
    contain whitespace.   Free variables in the expressions are interpreted as they would be
    in the caller.
    """
    caller = inspect.currentframe()
    thisFileName = caller.f_code.co_filename
    while ( caller.f_code.co_filename == thisFileName ):
        caller = caller.f_back
            

    caller_filename, caller_lineno, caller_function = caller.f_code.co_filename, caller.f_lineno, caller.f_code.co_name

    callerLocals = caller.f_locals
    callerGlobals = caller.f_globals
    del caller

    
    strs = []
    for ex in exprs:
        for e in ( ex.strip().split(' ') if not ex.startswith('#') else (ex[1:],) ):
            evalResult = eval( e, callerGlobals, callerLocals )
            strs.append( e + ' (' + str(type(evalResult)) + ') = ' + str(evalResult) )

    logging.debug( 'at ' + caller_filename + ':' + str(caller_lineno) + ' in ' + caller_function + ', ' +
                   ', '.join(strs) )

    del callerGlobals, callerLocals

    if DictGet( kwargs, 'tb', False ): traceback.print_stack()


def assertEq( x, y ):
    """Check that the two values are equal, if not then print them and abort"""

    if x != y:
        print('x (', type(x), ') = ', x)
        print('y (', type(y), ') = ', y)
        assert x == y
    
def PV( valName, val, printVal = None ):
    """Print value and return it"""
    if printVal != None:
        print(valName, ' = ', printVal, ' type=', type( printVal ))
    else:
        print(valName, ' = ', val, ' type=', type( val ))
    return val
    
def ListIdx( lst ):
    """Build a map from list value to its index in the list.
    If the value appears more than once, the index of the last occurrence is used."""
    return dict( list(map( reversed, enumerate( lst ) )) )

def iter_ith( it, item ):
  """Return i'ith item from iterator"""
  for i, v in enumerate( it ):
    if i == item: return v
  raise IndexError( "iter_ith: iterator does not have item number " + str( item ) )
    
def relpath(target, base=os.curdir):
    """
    Return a relative path to the target from either the current dir or an optional base dir.
    Base can be a directory specified either as absolute or relative to current dir.

    Taken from http://code.activestate.com/recipes/302594/
    """

    if not os.path.exists(target):
        raise OSError('Target does not exist: '+target)

    if not os.path.isdir(base):
        raise OSError('Base is not a directory or does not exist: '+base)

    base_list = (os.path.abspath(base)).split(os.sep)
    target_list = (os.path.abspath(target)).split(os.sep)

    # On the windows platform the target may be on a completely different drive from the base.
    if os.name in ['nt','dos','os2'] and base_list[0] != target_list[0]:
        raise OSError('Target is on a different drive to base. Target: '+target_list[0].upper()+', base: '+base_list[0].upper())

    # Starting from the filepath root, work out how much of the filepath is
    # shared by base and target.
    for i in range(min(len(base_list), len(target_list))):
        if base_list[i] != target_list[i]: break
    else:
        # If we broke out of the loop, i is pointing to the first differing path elements.
        # If we didn't break out of the loop, i is pointing to identical path elements.
        # Increment i so that in all cases it points to the first differing path elements.
        i+=1

    rel_list = [os.pardir] * (len(base_list)-i) + target_list[i:]
    return os.path.join(*rel_list)

scriptsCopied = set()

def HTMLScript( dir, scriptName ):
    """Install a given script in the given dir, and return the code to include the script in an HTML file."""
    if ( dir, scriptName ) not in scriptsCopied:
        if not os.path.isdir( dir ): os.makedirs( dir )
        shutil.copy( '../Operations/Ilya_Operations/PipeRun/auxfiles/' + scriptName + '.js', dir )
        WaitForFileToAppear( dir + '/' + scriptName + '.js' )
        os.system( 'chmod 777 ' + dir + '/' + scriptName + '.js' )
        scriptsCopied.add( ( dir, scriptName ) )
    return "<SCRIPT type=\"text/javascript\" src=\"" + scriptName + ".js\"> </SCRIPT>\n"

def isnumeric(x):
    """Test if the value is convertible to a float"""
    try: 
        float(x)
        return True
    except ValueError:
        return False

def progress( name, i, total, freq = 10 ):
    """Print progress reports, every freq percent of the way."""
    if total == 0:
        logging.info( 'Error: total == 0; %s: %d of %d ...' % ( name, i, total ) )
        return

    chunk = int( total / ( 100 / freq ) )
    if ( chunk == 0  or  i % chunk ) == 0:
        logging.info( '%s: %d of %d ( %d %% ) ...' % ( name, i, total, int( round( i / total * 100 ) ) ) )

        
def Enum(*names):
    """Create an enumerated type.
    
    Usage examples:
    
    print '\n*** Enum Demo ***'
    print '--- Days of week ---'
    Days = Enum('Mo', 'Tu', 'We', 'Th', 'Fr', 'Sa', 'Su')
    print Days
    print Days.Mo
    print Days.Fr
    print Days.Mo < Days.Fr
    print list(Days)
    for each in Days:
    print 'Day:', each
    print '--- Yes/No ---'
    Confirmation = Enum('No', 'Yes')
    answer = Confirmation.No
    print 'Your answer is not', ~answer    
    
    
    Taken from Python recipe 413486-1 .
    
    """
    assert names, "Empty enums are not supported" # <- Don't like empty enums? Uncomment!

    class EnumClass(object):
       __slots__ = names
       def __iter__(self):        return iter(constants)
       def __len__(self):         return len(constants)
       def __getitem__(self, i):  return constants[i]
       def __repr__(self):        return 'Enum' + str(names)
       def __str__(self):         return 'enum ' + str(constants)

    class EnumValue(object):
       __slots__ = ('__value')
       def __init__(self, value): self.__value = value
       Value = property(lambda self: self.__value)
       EnumType = property(lambda self: EnumType)
       def __hash__(self):        return hash(self.__value)
       def __cmp__(self, other):
          # C fans might want to remove the following assertion
          # to make all enums comparable by ordinal value {;))
          assert self.EnumType is other.EnumType, "Only values from the same enum are comparable"
          return cmp(self.__value, other.__value)
       def __invert__(self):      return constants[maximum - self.__value]
       def __bool__(self):     return bool(self.__value)
       def __repr__(self):        return str(names[self.__value])

    maximum = len(names) - 1
    constants = [None] * len(names)
    for i, each in enumerate(names):
       val = EnumValue(i)
       setattr(EnumClass, each, val)
       constants[i] = val
    constants = tuple(constants)
    EnumType = EnumClass()
    return EnumType

def TabSep( *items ):
    """Return a tab-separated, newline-terminated string containing string representations of the arguments"""

    return '\t'.join( map( str, items ) ) + '\n'


class EarlyExit(Exception):
    """An exception that can be raised inside a loop when we want to exit the loop's body early.
    Using an exception instead of a continue statement lets us define some code that gets run regardless
    of whether we exited the loop early."""

    def __init__(self, value = None):
        self.value = value
    def __str__(self):
        return repr(self.value)

def AddFileSfx( fullFileName, *sfx ):
    """Add the specified suffix to a filename, inserting it before the file extension.
    So, if the file was named ../Data/myfile.tsv and the suffix is 'verA',
    this returns '../Data/myfile_verA.tsv' .
    """

    fileName, fileExt = os.path.splitext( fullFileName if not fullFileName.endswith( '/' ) else fullFileName[:-1] )
    return ( ( fileName + Sfx( sfx ) ) if fileExt not in ( '.gz', 'bz2' ) else \
               AddFileSfx( fileName, *sfx ) ) + fileExt + ( '' if not fullFileName.endswith( '/' ) else '/' )

def AddExts( baseFN, *exts ):
  """Return tuple of filenames obtained by adding specified extensions to the base name"""
  return tuple( [ baseFN + ( '.' if ext and not ( baseFN.endswith('.') or ext.startswith( '.' ) ) else '' ) + ext
                                           for extList in exts for ext in ( extList.strip().split() if extList.strip() else ('',) ) ] )

def AddFileSubdir( subdir, fname ):
  """Add a given subdir before the given filename.

  >>> AddFileSubdir( 'stats', '/home/mydir/myfile.tsv' )
  '/home/mydir/stats/myfile.tsv'

  """

  return os.path.join( os.path.dirname( fname ), subdir, os.path.basename( fname ) )

def Sfx( *vals ):
    """Return a suffix string based on the value: if the value is non-empty, return '_' + val; but if it is
    already empty, then just return the empty string.
    """

    def underscorify( val ):
        """prepend undescrore if needed"""
        if not val and val != 0: return ''
        noPrefix = isinstance( val, (str,)) and val.startswith('#')
        alphaVal = MakeAlphaNum( str( val ) )
        return alphaVal[1:] if noPrefix else ( alphaVal if alphaVal.startswith( '_' ) else '_' + alphaVal )
    
    return ''.join( map( underscorify, flatten( vals ) ) )

matplotlib_linestyles = ( '-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', 'H',
                          '+', 'x', 'D', 'd', '|', '_' )

def MakeExt( ext ):
  """Returns a file extension, starting with a dot if the given one does not."""
  return ( '' if not ext or ext.startswith( '.' ) else '.' ) + ext

def ReplaceFileExt( fullFileName, newExt ):
    """Replace file extension with a new one"""
    fileName, fileExt = os.path.splitext( fullFileName )
    return fileName + MakeExt( newExt )

if haveBio:
    
  def LoadFastaAsStr( fname ):
      """Loading a single-sequence FASTA file as a string.
      Note that position 0 in the string will correspond to position 1
      in the fasta file, where positions are 1-offset!"""
      logging.info( 'Loading fasta from ' + fname )
      with open( fname ) as f:
          sequences = list( SeqIO.parse( f, 'fasta' ) )
      assert len( sequences ) == 1
      seqRecord = sequences[0]

      logging.info( 'Loaded fasta from ' + fname )
      return seqRecord.seq.tostring()

#    def WordStartingAt( str, x ):
#        return str[ x : str.index( ' ', x+1 ) ]


def compile_expr( expr ):
    """Return a code object for the given expression.  Shorthand for compile( expr, 'this_filename', 'eval' )."""
    return compile( expr, sys._getframe(1).f_code.co_filename, 'eval' )

class Histogrammer(object):
    """Incrementally builds a histogram, from values given one at a time.  You only need to specify bin size,
    but not min and max bin bounds in advance: Histogrammer divides the entire real line into
    bins of binSize, and keeps track of only the non-empty bins.   You can later coarsen the bins, or clamp them
    to a given range, as needed.
    """

    #
    # Invariant: 
    #

    # to-do:

    #   - set how many bins you want, and/or min bin size, and have bins automatically merge to achieve this.

    #   - remove dependence on numpy where we need a simple thing

    #   - allow non-uniformly spaced bins (but support specification of uniform bins also)

    #   - operation to merge given bins, e.g. merge all bins above 100 into "and higher" bin.
    
    def __init__(self, binSize, binIntType = int, binShift = 0.0, clampLo = float('nan'),
                 clampHi = float('nan') ):
        """Create a histogrammer with bins of the form:

        [binShift-binSize, binShift), [binShift, binShift+binSize), [binShift+binSize, binShift+2*binSize), ...
        
        """
        assert binSize > 0.0
        self.binSize = binSize
        self.binShift = binShift % binSize
        self.bin2count = collections.defaultdict( binIntType )
        self.valSum = 0.0
        self.numVals = 0
        self.numNaNs = 0
        self.numInfs = 0
        self.clampLo = clampLo
        self.clampHi = clampHi

    def binVal( self, val ):
        """Return bin number for given value"""
        return int( math.floor( ( val - self.binShift ) / self.binSize ) )

    def addVal( self, x ):
        """Add a value.  The value itself is not kept. NaN and inf values are ignored (but counted)."""

        if math.isnan(x):
            self.numNaNs += 1
            return
        if math.isinf(x):
            self.numInfs += 1
            return
        bin = self.binVal( x )
        self.bin2count[ bin ] += 1
        self.valSum += float(x)
        self.numVals += 1

    def addVals( self, vals ):
        """Add multiple values to the histogram"""
        for v in vals: self.addVal( v )

    def getNumVals( self ):
        """Return the total number of values added to the histogram"""
        result = sum( self.bin2count.values() ) #self.numVals
        assert result == self.numVals
        return result

    def getNumNaNs( self ):
        """Return the number of NaN values added to the histogram"""
        return self.numNaNs
        
    def getNumInfs( self ):
        """Return the number of inf values added to the histogram"""
        return self.numInfs

    def getValSum( self ):
        """Return the sum of all values in the histogram"""
        return self.valSum

    def getValAvg( self ):
        """Return the average of all values added to the histogram"""
        nvals = self.getNumVals()
        return self.getValSum() / float( nvals ) if nvals > 0 else np.nan

    def getNumBins( self ):
        """Return the number of non-empty bins"""
        return len( bin2count )

    def getMin( self ):
        """Return the left side of the smallest non-empty bin"""
        return min( self.getBinLefts() )
    
    def getMax( self ):
        """Return the right side of the largest non-empty bin"""
        return max( self.getBinRights() )

    def getBinLefts( self ):
        """Return the left boundaries of all non-empty bins"""
        return  [ self.binShift  +  bin * self.binSize  for bin, count in sorted( self.bin2count.items() ) ]

    def getBinRights( self ):
        """Return the right boundaries of all non-empty bins"""
        return [ leftBound + self.binSize for leftBound in self.getBinLefts() ]


    def getAllBinIds( self, minBinId = None, maxBinId = None ):
        """Return a sequence of all bin IDs, including empty bins."""
        keys = list(self.bin2count.keys())
        if not keys: keys = ( 0, )
        return list(range( min( keys ) if minBinId is None else minBinId ,
                      max( keys )+1 if maxBinId is None else maxBinId, 1))
    
    def getAllBinLefts( self, minBinId = None, maxBinId = None ):
        """Return the left boundaries of all bins"""
        return  [ self.binShift  +  bin * self.binSize  for bin in self.getAllBinIds( minBinId = minBinId, maxBinId = maxBinId ) ]
    
    def getBinBounds( self ):
        """Return the left boundaries of all non-empty bins"""
        return  [ ( self.binShift  +  bin * self.binSize, self.binShift  +  ( bin + 1 ) * self.binSize )
                    for bin, count in sorted( self.bin2count.items() ) ]

    def getAllBinCounts( self, normed = False, cumulative = False,
                         minBinId = None, maxBinId = None ):
        """Return the counts of all non-empty bins"""
        return self.__fixBinCounts( binCounts = [ DictGet( self.bin2count, binId, 0 )
                                                  for binId in self.getAllBinIds( minBinId = minBinId, maxBinId = maxBinId) ],
                                    **Dict( 'cumulative normed' ) )

    def getCumulativeBinFor( self, cumulativeUpTo ):
        """Return the bin id such that the bins to the left of this bin, inclusive, contain fraction cumulativeUpTo of the values"""
        bc = self.getBinCounts( normed = True, cumulative = True )
        if not bc: return 0
        for bin, fraction in zip( sorted( self.bin2count.keys() ), bc ):
            if fraction >= cumulativeUpTo: return bin
        return max( self.bin2count.keys() )

    def __fixBinCounts( self, binCounts, cumulative, normed ):
        """Convert bin counts to cumulative or normed as needed."""
        if cumulative:
            sumSoFar = 0
            for i in range( len( binCounts ) ):
                thisBinCount = binCounts[ i ]
                binCounts[ i ] += sumSoFar
                sumSoFar += thisBinCount

        if normed:
            binSum = self.getNumVals()
            binCounts = [ binCount / binSum for binCount in binCounts ]

        return binCounts

    def getBinCounts( self, normed = False, cumulative = False ):
        """Return the counts of all non-empty bins"""
        return self.__fixBinCounts( binCounts = list(map( operator.itemgetter( 1 ), sorted( self.bin2count.items() ) )),
                                    **Dict( 'cumulative normed' ) )

    def coarsenBy( self, factor ):
        """Coarsen the histogram by merging adjacent bins together.
        Each group of 'factor' adjacent bins will be merged together.
        This method creates a new histogram and returns it.
        """

        newHist = Histogrammer( binSize = self.binSize * factor, binShift = self.binShift )
        for bin, count in list(self.bin2count.items()):
            newBin = int( math.floor( bin / factor ) )
            newHist.bin2count[ newBin ] += count

        newHist.valSum = self.valSum
        newHist.numVals = self.numVals
        newHist.numNaNs = self.numNaNs
        newHist.numInfs = self.numInfs
        
        return newHist

    def __str__( self ):
        """Make a printable representation of the histogram"""
        return str( self.__dict__ )

    def __iadd__( self, v ):
        """Add counts from the other hist into this.
        The bin sizes must be the same."""

        if isinstance( v, numbers.Number ): self.addVal( v )
        else:
            assert abs( self.binSize - v.binSize ) < 1e-12
            assert abs( self.binShift - v.binShift ) < 1e-12

            for bin, count in list(v.bin2count.items()):
                self.bin2count[ bin ] += count

            self.valSum += v.valSum
            self.numVals += v.numVals
            self.numNaNs += v.numNaNs
            self.numInfs += v.numInfs

        return self

    def __add__( self, hist ):
        newHist = Histogrammer( binSize = self.binSize, binShift = self.binShift )
        newHist += self
        newHist += hist
        return newHist

    def saveToStream( self, out, html = False ):
        """Save the histogram to an output stream"""

        headings = ( 'bin', 'binStart', 'binEnd', 'binGap', 'count', 'countNormed', 'countCumul', 'countNormedCumul',
                     'hist', 'histCumul' )
                
        binRecords = []
        if not sum( self.bin2count.values() ):
            # if there are no data points, at least record
            # the bin size.
            binRecords = [ [ 0, self.binShift, self.binShift + self.binSize, 0, 0, 0.0, 0, 1.0, '*', '*' ] ]
            assert len( binRecords[0] ) == len( headings )
        else:  # if we have at least one non-empty bin
            numVals =self.getNumVals()
            countCumulative = 0
            prevBinStart = None
            maxBin = max( self.bin2count.values() )
            for bin, count in sorted( self.bin2count.items() ):
                if not count: continue
                binStart = self.binShift + bin * self.binSize
                countCumulative += count
                numStars = int( 100.0 * ( count / maxBin ) )
                binRecord = [ bin, binStart, binStart + self.binSize,
                                     ( binStart - prevBinStart ) if prevBinStart != None else 0,
                                     count, '%.2f' % ( count / numVals ), countCumulative,
                                     '%.2f' % ( countCumulative / numVals ), '*' * numStars,
                                   '*' * int( 100.0 * ( countCumulative / numVals ) ) ]
                assert len( binRecord ) == len( headings )
                binRecords.append( binRecord )
                prevBinStart = binStart
                                   

        if html:
            with htag('table', border = 1):
                with htag('thead'):
                    with htag('tr'):
                        for heading in headings:
                            with htag('th'):
                                print(heading)
                with htag('tbody'):
                    for r in binRecords:
                        with htag('tr'):
                            for val in map( str, r ):
                                with htag('td'):
                                    print(val)
        else:    

            tabwrite( out, *headings )
            isFirst = True
            for r in binRecords:
                if not isFirst: out.write( '\n' )
                tabwriten( out, *r )
                isFirst = False

    def save( self, fname ):
        """Save the histogram to a .tsv file"""
        with open( fname, 'w' ) as f:
            self.saveToStream( f )

        self.saveStats( AddFileSfx( fname, 'stats' ) )

    def saveStats( self, fname ):
        """Save various histogram statistics to the given file"""

        with open( fname, 'w' ) as f:
            tabwrite( f, 'stat', 'val', 'encodedVal' )
            isFirst = True
            for stat, val in sorted( self.__dict__.items() ):
                if isinstance( val, ( numbers.Number, (str,) ) ):
                    if not isFirst: f.write( '\n' )
                    tabwriten( f, stat, val, base64.urlsafe_b64encode( pickle.dumps( val ) ) )
                    isFirst = False
                
    @staticmethod
    def load( fname ):
        """Load histogram from .tsv file"""

        logging.info( 'loading histogram from ' + fname )

        attrs = {}
        with open( AddFileSfx( fname, 'stats' ) ) as f:
            f.readline()  # skip header
            for line in f:
                stat, val, encodedVal = line.strip().split()
                attrs[ stat ] = pickle.loads( base64.urlsafe_b64decode( encodedVal ) ) 

        result = Histogrammer( binSize = attrs[ 'binSize' ] )
        for attr, val in list(attrs.items()):
            setattr( result, attr, val )

        with TableIter( fname ) as f:
            for r in f:
                assert abs( r.binEnd - r.binStart - result.binSize ) < 1e-12
                if r.count > 0:
                    result.bin2count[ r.bin ] = r.count

        assert result.numVals == sum( result.bin2count.values() )
                    
        return result
        
    def printVals( self ):
        """Print the histogram"""

        print('binSize=', self.binSize, ' bin2count=', self.bin2count, ' avgVal=', self.getValAvg())

    __hash__ = None

    def printBins(self, minBin, maxBin):
        """Print the bins in a given range"""
        pass

    def __str__(self):
        return '[Histogrammer: binSize=%f binShift=%f bin2count=%s' % ( self.binSize, self.binShift,
                                                                        sorted( self.bin2count.items() ) )

class Histog(object):
  """Incrementally builds a histogram from a stream of values, automatically choosing bins so that they contain
  roughly similar numbers of values."""

  pass
        

@contextmanager
def htag(x, **kwargs):
    """Write out an html tag"""
    print('<%s' % x)
    for k, v in list(kwargs.items()):
        print(' %s="%s"' % ( k, v ))
    print('>')
    yield
    if x != 'col': print('</%s>' % x)

def make_htag( out ):
    """Create an htag function for writing out HTML tags to a given outfile"""

    @contextmanager
    def htag(x, **kwargs):
        """Write out an html tag"""
        out.write( '<%s' % x )
        for k, v in list(kwargs.items()):
            out.write( ' %s="%s"' % ( k, v ) )
        out.write( '>' )
        yield
        if x != 'col': out.write( '</%s>' % x )

    return htag

@contextmanager
def hparen(x = '()'):
    """Print open and close parentheses"""
    print(x[0])
    yield
    print(x[1])
    
    
def chomp( s ):
    """Return string without terminal newline if present"""
    return s[:-1] if s.endswith('\n') else s


def coerceVal( s ):
    """Coerce a string s representing a value, to the most specific value type if possible; else leave as string.
    Non-string values are returned as-is.
    The strings True and False are converted to the corresponding bool values; ints are converted to int;
    floats are converted to float.  Strings not interpretable as a Boolean, int or float are returned as-is,
    and no exception is generated.
    """
    if not isinstance( s, (str,) ): return s
    if s == 'True': return True
    if s == 'False': return False
    try: return int(s)
    except ValueError:
        try: return float(s)
        except ValueError: return s

class TableReaderBase(object):

    """Abstract base class for reading tables"""

    class TSVRecord(object):

        def __init__(self, lineVals, colName2num, fname, headings):

            self.lineVals = tuple( lineVals )
            self.colName2num = colName2num
            self.fname = fname
            self.headings = headings
        
        def __getattr__(self, name):
            """Get value of column by name"""
            if name.startswith( '__' ): return object.__getattr__( self, name )
            try: return coerceVal( self.lineVals[ self.colName2num[ name ] ] )
            except KeyError:
                print('lineVals=', self.lineVals, ' fname=', self.fname, ' headings=', self.headings)
                raise
                

        def __getitem__(self, key):
            """Get value of column by name or number"""

            if isinstance(key,(tuple,list)): return tuple( self.__getitem__(k) for k in key )

            v = self.lineVals[ self.colName2num[ key ] if isinstance( key, str ) else key ]

            retVal = tuple( map( coerceVal, v ) ) if isinstance(v, (tuple, list)) else coerceVal(v)
            
            return retVal

        def __iter__(self):
            return map( coerceVal, self.lineVals )

        def __str__(self):
            return str( list(zip( self.headings, self.lineVals )) )

        def __repr__(self):
            return self.__str__()

        def __len__(self):
            return len(self.headings)

    def __init__(self, headings, fname = None, allFiles = (), closeWhenDone = False,
                 valueFixer = None):
        """Create a TableReaderBase"""

        self.headings = tuple( headings )

        if not len( set( headings ) ) == len( headings ):
            print('headings=', headings)
            seen = set()
            for h in headings:
                if h in seen:
                    print('duplicate heading: ', h)
                seen.add( h )
            
        assert len( set( headings ) ) == len( headings )
        
        self.fname = fname
        self.allFiles = allFiles
        self.colName2num = dict( ( colName, colNum ) for colNum, colName in enumerate( self.headings ) )
        self.closeWhenDone = closeWhenDone
        self.valueFixer = valueFixer
        self.lineNum = 0
        
    def __iter__(self):
        return self

    def __getitem__(self,item):
        """Returns an TblIter that gets only specified columns"""

        return IterReader( headings = tuple( self.headings[h] if isinstance(h,int) else h for h in MakeSeq(item) ),
                           recordsIter = map( operator.itemgetter(*MakeSeq(item)), self ) )

    def __getattr__(self,name):
        """Return a given column"""
        if name not in self.headings: raise AttributeError
        return self[name]

    def numCols(self):
        """Return number of columns"""
        return len(self.headings)

    def takewhile(self, pred):
        return IterReader( headings = self.headings,
                           recordsIter = itertools.takewhile( pred, self ) )

    def renameCols(self, renamings):
        """Get a version of this TblIter with some columns renamed"""
        return IterReader( headings = [ h if h not in renamings else renamings[h] for h in self.headings ],
                           recordsIter = self )
    
    def __next__(self):
        try:
            lineVals = self._nextLine()
            #print '***lineVals: ', lineVals

            self.lineNum += 1

            if self.lineNum % 20000  ==  0:
                logging.info( joinstr(' ', self.lineNum, ' of ', \
                                  self.fname, ': ', list(zip( self.headings, lineVals )) ) )

            if len( lineVals ) != len( self.headings ):
                print('lineVals=', lineVals, ' self.headings=', self.headings, ' \nzip=', list(zip( self.headings, lineVals )))
                
            assert len( lineVals ) == len( self.headings ), 'lineVals=%d self.headings=%d: %s' % ( len( lineVals ), len( self.headings ), self.headings )

            if self.valueFixer: lineVals = list(map( self.valueFixer, lineVals ))

            return TableReaderBase.TSVRecord( lineVals, self.colName2num, self.fname, self.headings )
        except StopIteration:
            if self.closeWhenDone: self.close()
            raise


class TSVReader(TableReaderBase):

    """Class for reading .tsv files"""

    def __init__(self, f, fname, sep = '\t', closeWhenDone = False, headings = None,
                 valueFixer = None, skipFirstLines = 0 ):

        self.headerLines = []
        if headings == None:
            headings = f.readline()
            while headings.startswith( '#' ) and not ( fname.endswith( '.vcf' ) and headings.startswith( '#CHROM' ) ):
                self.headerLines.append( headings.strip() )
                headings = f.readline()
            headings = chomp( headings )
            headings = headings.split( sep )

        TableReaderBase.__init__(self, headings = headings, fname=fname, allFiles=(fname,),
                                 closeWhenDone = closeWhenDone, valueFixer = valueFixer )
        
        self.sep = sep
        self.freader = f.__iter__()
        self.f = f
        self.skipFirstLines = skipFirstLines


    def _nextLine(self):
        """Get next line"""
        rawLine = chomp( next(self.freader) )
        if not rawLine: raise StopIteration

        # skip any comment lines.  normally these only appear at the beginning.
        while rawLine.startswith( '#' ):
            rawLine = chomp( next(self.freader) )
            if not rawLine: raise StopIteration

        while self.skipFirstLines > 0:
            rawLine = chomp( next(self.freader) )
            if not rawLine: raise StopIteration
            self.skipFirstLines -= 1

        line = rawLine.split( self.sep )
        if not line: raise StopIteration
        return line

    def close(self):
        """Close the file"""
        self.f.close()
        
        
class DotDataReader(TableReaderBase):

    """Class for reading .tsv files"""

    def __init__(self, fname, fs, cols, allFiles, closeWhenDone = False,
                 valueFixer = None):

        TableReaderBase.__init__(self, headings = cols, fname=fname, allFiles=allFiles,
                                 closeWhenDone = closeWhenDone, valueFixer = valueFixer )
        
        fs = tuple( fs )
        self.fs = fs
        self.freaders = [ f.__iter__() for f in fs ]
        self.headerLines = []


    def _nextLine(self):
        """Get next line"""
        nextLine =  tuple( chomp( next(freader) ) for freader in self.freaders )
        if not nextLine: raise StopIteration
        return nextLine

    def close(self):
        """Close all open files we had"""
        for f in self.fs: f.close()
        

class IterReader(TableReaderBase):

    """Class for reading records from an iterator"""

    def __init__(self, headings, recordsIter):
        TableReaderBase.__init__(self, headings = headings)

        self.recordsIter = filter( lambda x: x != None, iter( recordsIter ) )
        if len(headings)==1:
            self.recordsIter = map( lambda v: (v,) if isinstance( v, (int, float, str)) else v,
                                               self.recordsIter )


    def _nextLine(self): return next(self.recordsIter)
        
@contextmanager
def TableIter( fname, sep = '\t', fileType = None ):
    """Returns an iterator over SV files.  The returned object f has the field f.headings which gives the column names.
    Iterating over this object yields line objects.   Each line object is like a named tuple.
    """

    logging.info( 'TabelIter: creating iter over ' + fname )
    
    if fileType == '.ssv': sep = ' '
    if sep == ' ': fileType = '.tsv'
    
    if fname.endswith( '.tsv' ) or fname.endswith( '.tsv.gz' ) or fileType == '.tsv':
        with OpenForRead( fname ) as f:
            yield TSVReader( f = f, fname = fname, sep = sep )
    else:
        if not fname.endswith('/'): fname += '/'
        assert fname.endswith('.data/')

        headerFN = fname + os.path.basename(fname[:-len('.data/')]) + '.header.txt'
        if not os.path.exists( headerFN ):
            headerFNs = glob.glob( fname + '*.header.txt' )
            assert len( headerFNs ) == 1
            headerFN = headerFNs[0]

        cols = tuple( SlurpFileLines( headerFN ) )

        tsvFiles = [ f for f in reduce( operator.concat, list(map( operator.itemgetter(2), os.walk( fname ) )) ) if f.endswith('.csv') ]

        col2file = dict( ( os.path.splitext( os.path.splitext( tsvF )[0] )[0], tsvF ) for tsvF in tsvFiles )
        allFiles = [ fname + col2file[ col ] for col in cols ]
        
        with nested( *list(map( open, allFiles )) ) as fs:
            yield DotDataReader( fname = fname, fs = fs, cols = cols, allFiles = allFiles )
            

def TblIter( fname, sep = '\t', fileType = None, headings = None, valueFixer = None, skipFirstLines = 0 ):
    """Returns an iterator over SV files.  The returned object f has the field f.headings which gives the column names.
    Iterating over this object yields line objects.   Each line object is like a named tuple.
    """

    logging.info( 'TblIter: creating iter over ' + fname )

    if fileType == None  and  not fname.endswith('.data') and not fname.endswith('.data/'): fileType = '.tsv'
    if sep == ' ': fileType = '.tsv'
    
    if fname.endswith( '.tsv' ) or fname.endswith( '.tsv.gz' ) or fileType == '.tsv':
        return TSVReader( f = OpenForRead( fname ), fname = fname, sep = sep, closeWhenDone = True, headings = headings,
                          valueFixer = valueFixer, skipFirstLines = skipFirstLines )
    else:
        if not fname.endswith('/'): fname += '/'
        assert fname.endswith('.data/')

        headerFN = fname + os.path.basename(fname[:-len('.data/')]) + '.header.txt'
        if not os.path.exists( headerFN ):
            headerFNs = glob.glob( fname + '*.header.txt' )
            assert len( headerFNs ) == 1
            headerFN = headerFNs[0]

        cols = tuple( SlurpFileLines( headerFN ) )

        tsvFiles = [ f for f in reduce( operator.concat, list(map( operator.itemgetter(2), os.walk( fname ) )) ) if f.endswith('.csv') ]

        col2file = dict( ( tsvF.split('.')[0], tsvF ) for tsvF in tsvFiles )
        allFiles = [ fname + col2file[ col ] for col in cols ]
        
        return DotDataReader( fname = fname, fs = list(map( open, allFiles )), cols = cols, allFiles = allFiles, closeWhenDone = True,
                              valueFixer = valueFixer )

            
def SaveTableIterToSV( tableIter, fname ):
    """Save the output of TableIter to a .tsv file"""

    with OpenForWrite( fname ) as f:
        f.write( '\t'.join( tableIter.headings ) )
        for line in tableIter:
            f.write( '\n' + '\t'.join( map( str, line ) ) )


def SaveTableIterToDotData( tableIter, fname ):
    """Save the output of TableIter to a .tsv file"""

    if not fname.endswith('/'): fname += '/'
    assert fname.endswith('.data/')

    MakeDir( fname )
    
    DumpFile( fname + os.path.basename( fname[:-len('.data/')] ) + '.header.txt', '\n'.join( tableIter.headings ) )

    recordsIter = iter( tableIter )
    
    firstRec = next(recordsIter)
    colTypes = [ type(coerceVal(val)) for val in firstRec ]

    type2name = { float: 'float', int: 'int', str: 'str' }

    typeOrder = ( int, float, str )

    newColTypes = copy.copy( colTypes )
    
    with nested( *[ open(fname + heading + '.' + type2name[colType] + '.csv', 'w')
                    for heading, colType in zip( tableIter.headings, colTypes ) ] ) as fs:
        isFirstLine = True
        for r in itertools.chain( (firstRec,), recordsIter ):
            for colNum, ( v, f ) in enumerate( zip( r, fs) ):
                if not isFirstLine: f.write( '\n' )
                f.write( str(v) )

                newType = type( coerceVal( v ) )
                oldType = newColTypes[ colNum ]
                if typeOrder.index( newType ) > typeOrder.index( oldType ): newColTypes[ colNum ] = newType

            isFirstLine = False
                

    for heading, colType, newColType in zip( tableIter.headings, colTypes, newColTypes ):
        if newColType != colType:
            # we have to rename the file to a new correct name.
            # make sure it has been fully written and closed.
            time.sleep(20)

            oldName = fname + heading + '.' + type2name[colType] + '.csv'
            WaitForFileToAppear( oldName )
            newName = fname + heading + '.' + type2name[newColType] + '.csv'
            os.rename( oldName, newName )
            time.sleep(10)
            WaitForFileToAppear( newName )
            
                
def PrintTableIter( tableIter ):
    """Save the output of TableIter to a .tsv file"""

    sys.stdout.write( '\t'.join( tableIter.headings ) )
    for line in tableIter:
        sys.stdout.write( '\n' + '\t'.join( map( str, line ) ) )

    sys.stdout.write( '\n' )


            
def makeIdxPrepender( idx ):
    """Creates a function that prepends a given index to its argument, returning a tuple."""
    return lambda v: ( idx, v )

class attrDbg(object):

    def __getattribute__(self, attr):
        print('getting attr ', attr)
        return attr
    
def makeKeyGetter( k ):
    """Creates a function that gets a key"""
    def myFunc( v ):
        return k( v[1] )
    print('making key getter for k=', k)
    return myFunc

def makeKeyApplier( k ):
    """Make a function that applies a key"""
    def myApplier( v ):
        r = ( k(v), v )
        #print 'applying ', k, ' to ', v, ' got ', r
        return r
        
    return myApplier


def itermerge( iters, keys = None, ids = None, includeKeys = False ):
    '''Merge multiple sorted inputs into a single sorted output.

    Equivalent to:  sorted(itertools.chain(*iterables))

    Params:

      keys - if given, then for each iter gives a function to extract the ordering key
        from elements yielded by that iter

      ids - if True, the yielded sequence is a sequence of pairs with the first element of the tuple being
        the index of the iterator from which the original element comes, and the second being the original element.
        if a sequence, instead of the numeric index of the iterator, a corresponding value from this sequence is used.

    >>> list(itermerge(iters=([1,3,5,7], [0,2,4,8], [5,10,15,20], [], [25])))
    [0, 1, 2, 3, 4, 5, 5, 7, 8, 10, 15, 20, 25]

    Adapted from http://code.activestate.com/recipes/535160/

    '''

    logging.debug( 'in itermerge: keys=' + str( keys ) )

    iters = tuple( iters )

    if keys:

        def fixKey( k ):
            if isinstance( k, int ): return operator.itemgetter( k )
            if isinstance( k, str ): return operator.attrgetter( k )
            if isinstance( k, ( tuple, list ) ):
                return functools.partial( lambda keys, vals: tuple( a_key( a_val ) for a_key, a_val in zip( keys, vals ) ),
                                          list(map( fixKey, k )) )

                return fixTuple( list(map( fixKey, k )) )

            return k
        
        keys = [ ( operator.itemgetter( k ) if isinstance( k, (int,tuple,list) ) else
                   ( operator.attrgetter( k ) if isinstance( k, str ) else k ) ) for k in keys ]
    
    if ids:
        return itermerge( iters = [ map( makeIdxPrepender( idx ), i ) for idx, i in
                                    ( enumerate( iters ) if ids == True else list(zip( ids, iters )) ) ],
                          keys = ( operator.itemgetter( 1 ), ) * len( iters ) if keys == None else
                          list(map( makeKeyGetter, keys )),
                          includeKeys = includeKeys )
    
    if keys:
        print('keys are: ', keys)
        result = itermerge( iters = [ map( makeKeyApplier( k ), i )
                                      for k, i in zip( keys, iters ) ] )
        if not includeKeys: result = map( operator.itemgetter( 1 ), result )
        return result
    
    def merge( i1, i2 ):
        next1 = iter( i1 ).__next__
        next2 = iter( i2 ).__next__
        try:        
            v1 = next1()
        except StopIteration:
            while True:
                yield next2()
        try:
            v2 = next2()
        except StopIteration:
            while True:
                yield next1()
        while True:
            if v1 < v2:
                yield v1
                try:        
                    v1 = next1()
                except StopIteration:
                    yield v2
                    while True:
                        yield next2()
            else:
                yield v2
                try:
                    v2 = next2()
                except StopIteration:
                    yield v1
                    while True:
                        yield next1()
    iters_cnt = len( iters )
    if iters_cnt == 1:
        return iter( iters[0] )
    if iters_cnt == 2:
        return merge( iters[0], iters[1] )
    bisect = int( iters_cnt / 2 )
    return merge( itermerge( iters = iters[:bisect] ), itermerge( iters = iters[bisect:] ) )

class GeneratorWrap(object):

    """Wrap a generator to create an object that can have attributes
    assigned to it.
    """

    def __init__(self, gen):
        self.gen = gen

    def __iter__(self): return self
    def __next__(self):
        v =  next(self.gen)
        #print 'yiiiield ', v
        return v 
    
def TableMaker( fn ):
    """Create a proper table maker"""

    def tmaker( *args, **kwargs):

        gen = fn( *args, **kwargs )
        headings = next(gen)

        logging.info( 'creating IterReader with headings %s' % headings )
        return IterReader( headings = headings,
                           recordsIter = gen )
                           
    return tmaker

    
def TableIterInnerJoin( tableIters, cols, suffixes = None, blanks = None, concat = True ):
    """Create an iterator that yields an inner join of the specified tables"""

    @TableMaker
    def TableIterInnerJoinAux( tableIters, cols, headings, blanks ):


        yield headings if concat else [ 'h_%d' % i for i in range( len( tableIters ) ) ]

        headingLens = [ len( ti.headings ) for ti in tableIters ]
        assert sum( headingLens ) == len( headings )
        blanks = [ blank if blank is None or hasattr(blank,'__len__') else (blank,)*headingLen
                   for blank, headingLen in zip(blanks, headingLens) ]
        
        prevKey = None
        for k, g in itertools.groupby( itermerge( iters = tableIters, keys = cols, ids = True,
                                                  includeKeys = True ),
                                       key = operator.itemgetter( 0 ) ):

            # check that the keys are sorted in strictly increasing order -- important for correct operation of join 
            if not( prevKey==None or k > prevKey ):
                logging.info( 'prevKey=' + str(prevKey) + ' key=' + str(k) )
            assert prevKey==None or k > prevKey
            prevKey = k

            records = tuple( g )

            origins = [ r[1][0] for r in records ]
            if not is_sorted( origins, strict = True ):
                print('records are ', records)
                print('origins are ', origins)
            assert is_sorted( origins, strict = True )

            recordsList = [ None ] * len( tableIters )
            positionFilled = [ False ] * len( tableIters )
            for r in records:
                recordsList[ r[1][0] ] = r[1][1]
                positionFilled[ r[1][0] ] = True
                
            for i in range( len( tableIters ) ):
                if not positionFilled[ i ] and blanks[ i ] != None:
                    recordsList[ i ] = blanks[ i ]
                    positionFilled[ i ] = True

            if all( positionFilled ):
                rec = list(map( tuple, recordsList ))
                assert list(map( len, rec )) == headingLens
                if concat: rec = reduce( operator.concat, rec )
                yield rec

    tableIters = tuple( tableIters )
    if blanks == None: blanks = (None,) * len( tableIters )
    if suffixes == None: suffixes = [ '_%d' % i for i in range( len( tableIters ) ) ]

    assert all([ not hasattr(blank,'__len__') or  len( blank ) == len( tableIter.headings )
                 for tableIter, blank in zip( tableIters, blanks ) ])

    assert len( suffixes ) == len( tableIters )
    logging.info( 'suffixes are ' + str(suffixes) )

    sharedColNames = set( [] )
    allColNames = set( [] )
    for tableIter in tableIters:
            for heading in tableIter.headings:
                    ( sharedColNames if heading in allColNames else allColNames ).add( heading )

    logging.info( 'sharedColNames=%s' % sharedColNames )
                    
    newHeadings = [ n if n not in sharedColNames else n + sfx for tableIter, sfx in zip( tableIters, suffixes )
                    for n in tableIter.headings ]

    logging.info( 'newHeadings are: %s' % newHeadings )

    return TableIterInnerJoinAux( tableIters = tableIters, cols = cols, headings = newHeadings, blanks = blanks )


def mergeRecs( *vals ):
    return reduce( operator.concat, list(map( tuple, vals )) )


def ihstack( *tableIters ):
    """Create a new table iterator that yields these tables stacked horizontally"""

    return IterReader( headings = reduce( operator.concat, list(map( operator.attrgetter( 'headings' ), tableIters )) ),
                       recordsIter = map( mergeRecs, *tableIters ) )


def ivstack( *tableIters ):
    """Create a new table iterator that yields these tables stacked vertically"""

    headings = tableIters[0].headings
    assert all([ t.headings == headings for t in tableIters ])
    
    return IterReader( headings = headings,
                       recordsIter = itertools.chain( *tableIters ) )

@TableMaker
def ivstack_iterable( tableIters ):
    """Create a new table iterator that yields these tables stacked vertically"""

    for i, it in enumerate(tableIters):
        if i == 0:
            headings = it.headings
            yield headings
        else:
            assert it.headings == headings
        for r in it: yield r
    
#    firstIter = tableIters.__iter__().next()
    
#    headings = firstIter.headings
    
#    return IterReader( headings = headings,
#                       recordsIter = itertools.chain.from_iterable( itertools.chain( (firstIter,), tableIters ) ) )


def TblIterFromDotData( dotData ):
    """Create a TblIter from DotData"""

    return IterReader( headings = dotData.dtype.names,
                       recordsIter = dotData if dotData.numCols() > 1 else map( (lambda x: (x,)), dotData ) )


def TblIterFilter( pred, tableIter ):
    """Filter a TableIter with a predicate"""

    return IterReader( headings = tableIter.headings,
                       recordsIter = filter( pred, tableIter ) )
    

def is_sorted( s, strict = True ):
    """Test if a sequence is sorted"""

    prev_elem = None
    for x in s:
        if prev_elem != None and x < prev_elem or ( x == prev_elem and strict ):
            return False
        prev_elem = x
    return True

def is_iterable( x ):
    """Check if x is an iterable object"""
    return hasattr( x, '__iter__' )


IUPAC_ambig_codes = { 'A': 'A',
                      'C': 'C',
                      'G': 'G',
                      'T': 'T',
                      'AC': 'M',
                      'AG': 'R',
                      'AT': 'W',
                      'CG': 'S',
                      'CT': 'Y',
                      'GT': 'K',
                      'ACG': 'V',
                      'ACT': 'H',
                      'AGT': 'D',
                      'CGT': 'B',
                      'ACGT': 'N' }

IUPAC_ambig_codes_set = dict( ( frozenset( k ), v ) for k, v in list(IUPAC_ambig_codes.items()) )

def ambiguousIUPACdnaCode( codes ):
    """Given a set of DNA codes, return the one-letter IUPAC code representing them;
    from http://www.ncbi.nlm.nih.gov/SNP/iupac.html"""

    return IUPAC_ambig_codes_set[ frozenset( v.upper() for v in codes ) ]

def BreakString( s, maxLen ):
    """Break string into lines of at most given length"""

    quot, rem = divmod( len( s ), maxLen )
    r = [ s[i*maxLen:(i+1)*maxLen] for i in range( quot ) ]
    if len( s ) % maxLen:
        r.append( s[-rem:] )

    return r

class namespace: pass

def GetDefaultArgs( fn ):
    """Get the default args of the function and their values, as a dictionary"""

    argInfo = inspect.getargspec( fn )
    if not argInfo.defaults: return {}
    return dict( list(zip( argInfo.args[ -len( argInfo.defaults ): ], argInfo.defaults )) )

def IsValidFileName( fname ):
    """Return true if fname is a valid file name"""
    if not hasattr( os, 'pathconf' ): return True
    maxNameLen = os.pathconf( '.', 'PC_NAME_MAX' ) if not hasattr( os, 'pathconf_names' ) or 'PC_NAME_MAX' in os.pathconf_names else 255
    maxPathLen = os.pathconf( '.', 'PC_PATH_MAX' ) if not hasattr( os, 'pathconf_names' ) or 'PC_PATH_MAX' in os.pathconf_names else 4096
    if len( fname ) > maxPathLen: return False
    while True:
        d, n = os.path.split( fname )
        if len( n ) > maxPathLen: return False
        if len( d ) < 2 or d.endswith( '/' ): return True
        fname = d

MAX_FILE_NAME_LEN = 120

class DbgIter(collections.Iterator):

    """Takes a given iterator, and yields its values, also keeping track
    of them for debug purposes."""

    def __init__(self, it):
        self.it = it
        self.valNum = 0
        self.lastVal = None

    def __next__(self):
        self.valNum += 1
        self.lastVal = next( self.it )
        return self.lastVal


def IsConvertibleTo( x, aType ):
    """Test if a value is convertible to a given type"""
    try:
        dummy = aType(x)
        return True
    except ValueError: return False

def IsPosInt( x ):
    """Test if a number is an integer"""
    try:
        val = int(x)
        return val > 0
    except ValueError: return False

class DupChecker(object):

  """Check for duplicate objects, assert if one found"""

  def __init__(self):
    self.seen = set()

  def add( self, obj ):
    """Add and check for duplicates"""
    assert obj not in self.seen, 'duplicate object %s' % str( obj )
    self.seen.add( obj )

  def __len__( self ): return len( self.seen )
  def __str__( self ): return str( self.seen )

# end class DupChecker    
    
class NonClosingContext( object ):

  def __init__(self, f):
    self.f = f

  def __enter__(self): return self.f
  def __exit__(self, *exc_info): pass

  def __iter__(self): return iter( self.f )

  def close(self): pass

  def __getattr__(self, attr):
    """Delegate everything but enter and exit to the stream"""
    return getattr(self.f, attr)

class ClosingContext( NonClosingContext ):

  def __init__(self, f): super(type(self),self).__init__( f )
  def __exit__(self, *exc_info): self.f.close()
  
def OpenForRead( fname ):
    """Open a normal or compressed file for reading"""
    if fname == sys.stdin: return NonClosingContext( fname )
    return ClosingContext( gzip.open( fname, 'rb' ) if fname.endswith( '.gz' ) else
                           ( bz2.BZ2File( fname ) if fname.endswith( '.bz2' ) else open( fname ) ) )

def OpenForWrite( fname ):
    """Open a normal or compressed file for writing"""

    if fname in ( sys.stdout, sys.stderr ):
      return NonClosingContext( fname )

    return ClosingContext( gzip.open( fname, 'wb' ) ) if fname.endswith( '.gz' ) else open( fname, 'w' )

def IsFileType( f, *ftypes, **kwargs ):
  """Returns True if the given file is one of the given types or its compressed variants.

  >>> IsFileType( 'f.tsv', 'tsv' )
  True
  >>> IsFileType( 'f.tsv.gz', 'tsv' )
  True
  >>> IsFileType( 'f.tsv.', 'tsv' )
  False
  
  """

  compressedSfxs = tuple( MakeSeq( kwargs.get( 'compressedSfxs', ( '.gz', ) ) ) )
  return hasattr( f, 'endswith' ) and any( f.endswith( MakeExt( t ) + MakeExt( c ) ) for t in ftypes
                                           for c in ('',) + compressedSfxs )


class SumKeeper(object):
    """A class for keeping an accurate running sum for a stream of numbers, without loss of precision, by keeping
    intermediate sums.

    Adapted from: http://code.activestate.com/recipes/393090/

    """

    def __init__(self):
        self.clear()

    def add( self, x ):
        """Add a value to the SumKeeper"""

        if math.isnan(x): self.numNaNs += 1
        elif math.isinf(x): self.numInfs += 1
        else:
            self.numVals += 1

            i = 0
            for y in self.partials:
                if abs(x) < abs(y):
                    x, y = y, x
                hi = x + y
                lo = y - (hi - x)
                if lo:
                    self.partials[i] = lo
                    i += 1
                x = hi
            self.partials[i:] = [x]
            
        return self

    def __iadd__( self, x ):
        """Add a value to the SumKeeper"""
        return self.add( x )

    def addVals( self, vals ):
        """Add multiple values to the SumKeeper"""
        for x in vals: self.add( x )

    def clear( self ):
        """Reset this SumKeeper to zero."""
        self.partials = []
        self.numNaNs = 0
        self.numInfs = 0
        self.numVals = 0

    def getSum( self, dtype = None ):
        """Return the sum of values added so far"""
        return sum( self.partials, 0.0 ) if dtype is None else np.sum( self.partials, dtype = dtype )

class StatKeeper(object):
    """A class for keeping accurate stats for a stream of numbers."""

    def __init__(self):
      self.clear()

    def add(self, x ):
      """Add a value to the StatKeeper"""
      if not np.isnan( x ):
          self.sum.add( x )
          self.sumSq.add( x * x )
          self.numVals += 1
      else: self.numNaNs += 1
          
      return self

    def __iadd__( self, x ):
      """Add a value to the StatKeeper"""
      return self.add( x )

    def addVals( self, vals ):
      """Add multiple values to the SumKeeper"""
      for x in vals: self.add( x )

    def clear( self ):
        """Reset this StatKeeper to zero."""
        self.sum = SumKeeper()
        self.sumSq = SumKeeper()
        self.numVals = 0
        self.numNaNs = 0

    def getCount( self ):
        """Return the number of values added so far, not counting NaNs."""
        return self.numVals

    def getNumNaNs( self ):
        """Return the number of NaNs seen so far"""
        return self.numNaNs
        
    def getSum( self, dtype = None ):
        """Return the sum of values added so far"""
        return self.sum.getSum( dtype )

    def getMean( self, dtype = None ):
        """Return the mean of values added so far"""
        return self.getSum( dtype ) / self.getCount() if self.numVals else np.nan

    def getStd( self, dtype = None ):
        """Return the stddev of values added so far"""
        if not self.numVals: return np.nan
        meanSoFar = self.getMean( dtype )
        return math.sqrt( self.sumSq.getSum( dtype ) / self.numVals - ( meanSoFar * meanSoFar ) )

def DirUnder( parentDir, subDir ):
    """If subdir is either an absolute path or a path relative to . or ..,
    return subDir and ignore parentDir; otherwise, return os.pathjoin( parentDir, subDir )."""

    return subDir if subDir.startswith( '/' ) or subDir.startswith( '.' ) else os.path.join( parentDir, subDir )

def tmap( *args ):
    """Make a map and then wrap it in a tuple"""
    return tuple( map( *args ) )

def clamp( v, v_min, v_max ):
    """Return v clamped to the given interval"""
    if v < v_min: v = v_min
    if v > v_max: v = v_max
    return v

def PrintDictDiff( d1, d2 ):
    """Print the differences between dictionaries"""

    d1only = [ k1 for k1 in list(d1.keys()) if k1 not in d2 ]
    d2only = [ k2 for k2 in list(d2.keys()) if k2 not in d1 ]
    diffVals = [ k for k in list(d1.keys()) if k in list(d2.keys()) and d1[k] != d2[k] ]

    print('d1only=', d1only, ' d2only=', d2only, ' diffVals=', diffVals)
    
        
def GetHostName():
    """Return the full hosname of the machine on which we are running"""
    hn = 'unknown'
    try:
      hn = socket.gethostname()
      hn = socket.gethostbyaddr(hn)[0]
    except socket.herror as err:
      logging.warning( 'error getting hostname: %s; hostname=%s' % ( err, hn ) )
    return hn
        
def WordStartingAt( str, x ):
  """Return the word in the given string starting at the given position"""
  return str[ x : str.index( ' ', x+1 ) ]

def AddSlash( dirName ):
  """Add a slash to the dir name if not there already"""
  return dirName if dirName.endswith('/') else dirName + '/'

def RemoveSlash( dirName ):
  """Remove a slash from the end of dir name if present"""
  return dirName[:-1] if dirName.endswith('/') else dirName 
    
def parseNumList( s ):
    """Parse list of numbers, possibly with ranges"""

    return reduce( operator.concat,
                   [list(range(*list(map( int, term.split( '-' ) )))) if '-' in term else [ int( term ) ] for term in s.split( ',' )],
                   [] )

    
def NoDups( s ):
  """Return true if s has no duplicates"""
  return len( set( s ) ) == len( s )
    
def mapValues( d, f ):
  """Return a new dict, with each value mapped by the given function"""
  return dict([ ( k, f( v ) ) for k, v in list(d.items()) ]) 
        

class EventCounter(object):

  """Helper class for counting how many events of various types happened during a processing step"""

  def __init__(self):
    self.eventCounts = collections.defaultdict( int )
    self.key2descr = {}

  def __call__(self, key, descr = None ):
    """Count the given event"""
    if isinstance( key, Exception ):
      exc = key
      key = exc.args[0]
      if len( exc.args ) > 1: descr = exc.args[1]
    self.eventCounts[ key ] += 1
    self.key2descr[ key ] = descr

  def setStat(self, key, val, descr = None ):
    """Explicitly set the value of a statistic"""
    self.eventCounts[ key ] = val
    self.key2descr[ key ] = descr


  def save( self, outFN, header = None ):
    """Save the counts to a file"""
    with OpenForWrite( outFN ) as f:
      if header: f.write( '# ' + header + '\n' )
      tabwrite( f, 'event', 'count', 'descr' )
      for key, count in sorted( self.eventCounts.items() ):
        tabwrite( f, key, count, DictGetNotNone( self.key2descr, key, key ) )

  @staticmethod
  def combine( inFNs, outFN, getio = None ):
    """Adds up the counters in the input files and writes them to output file""" 
    if getio: return dict( depends_on = inFNs, creates = outFN )
    inFNs = MakeSeq( inFNs )
    event2count = collections.defaultdict( int )
    event2descr = {}
    headerLines = []
    for i, inFN in enumerate( inFNs ):
        headerSeen = False
        for line in SlurpFileLines( inFN ):
          if line.startswith( '#' ):
            if i == 0: headerLines.append( line )
            continue
          if not headerSeen:
            headerSeen = True
            continue

          key, count, descr = line.split( '\t' )
          event2count[ key ] += int( count )
          event2descr[ key ] = descr
          
    with OpenForWrite( outFN ) as outFile:
      if headerLines: outFile.write( '\n'.join( headerLines ) + '\n' )
      tabwrite( outFile, 'event', 'count', 'descr' )
      for key, count in sorted( event2count.items() ):
        tabwrite( outFile, key, count, DictGetNotNone( event2descr, key, key ) )
        
def xor( a, b ):
  """Compute the boolean xor of two values"""
  return bool( a ) != bool( b )
    
    
    
class TabularFileIdx(object):
  """An index of a tabular file that lets one quickly get to a given line or value"""

  def __init__(self, fileName, keyCols, commentPrefix = '#'):
    self.commentPrefix = commentPrefix
    self.keyCols = keyCols

  
def vstack_tsvs( inFNs, outFN, commentPrefix = '#', blockSize = 67000, useHeader = True,
                 headings = None, getio = None ):
  """Vertically stack tsv files"""

  if getio: return dict( depends_on = inFNs, creates = outFN, attrs = dict( piperun_short = True )  )

  with OpenForWrite( outFN ) as outF:

    class FileTailKeeper(object):

      """Keeps track of the tail of data written to a file.
      This lets us know whether a newline at the end is needed."""

      def __init__(self, f ):
        self.f = f
        self.tail = ''

      def write( self, data ):
        if data:
          self.f.write( data )
          self.tail = data[-1]

    ftk = FileTailKeeper( outF )

    if headings: ftk.write( '\t'.join( MakeSeq( headings ) ) + '\n' )
    
    for i, inFN in enumerate( inFNs ):

      logging.info( 'appending %s to %s' % ( inFN, outFN ) )

      with OpenForRead( inFN ) as f:

        if useHeader:
          while True:
            line = f.readline()
            if i == 0: ftk.write( line )
            if commentPrefix and line.startswith( commentPrefix ):
              continue
            break

        while True:
          chunk = f.read( blockSize )
          if not chunk: break
          ftk.write( chunk )

      if ftk.tail and ftk.tail != '\n': ftk.write( '\n' )


def vstack_tsvs_lenient( inFNs, outFN, commentPrefix = '#', blockSize = 67000, useHeader = True,
                         headings = None, getio = None ):
  """Vertically stack tsv files"""

  inFNs = list(filter( os.path.exists, inFNs ))

  if getio: return dict( depends_on = inFNs, creates = outFN, attrs = dict( piperun_short = True )  )

  with OpenForWrite( outFN ) as outF:

    class FileTailKeeper(object):

      """Keeps track of the tail of data written to a file.
      This lets us know whether a newline at the end is needed."""

      def __init__(self, f ):
        self.f = f
        self.tail = ''

      def write( self, data ):
        if data:
          self.f.write( data )
          self.tail = data[-1]

    ftk = FileTailKeeper( outF )

    if headings: ftk.write( '\t'.join( MakeSeq( headings ) ) + '\n' )
    
    for i, inFN in enumerate( inFNs ):

      logging.info( 'appending %s to %s' % ( inFN, outFN ) )

      with OpenForRead( inFN ) as f:

        if useHeader:
          while True:
            line = f.readline()
            if i == 0: ftk.write( line )
            if commentPrefix and line.startswith( commentPrefix ):
              continue
            break

        while True:
          chunk = f.read( blockSize )
          if not chunk: break
          ftk.write( chunk )

      if ftk.tail and ftk.tail != '\n': ftk.write( '\n' )

      
def FirstVal( *vals ):
    """Returns the first non-None value"""
    for v in vals:
        if v is not None: return v
    return None
      
def SplitStr( s, sep ):
  """Split a string into tokens.

  Params:
     s - the string to split
     sep - the separator; if equal to 'whitespace', split on whitespace.
  """
  return s.split( sep if sep != 'whitespace' else None )
  

def groupIntoTuples(lst, n):
    """group([0,3,4,10,2,3], 2) => iterator

    Group an iterable into an n-tuples iterable. Incomplete tuples
    are discarded e.g.

    >>> list(groupIntoTuples(range(10), 3))
    [(0, 1, 2), (3, 4, 5), (6, 7, 8)]

    From http://code.activestate.com/recipes/303060-group-a-list-into-sequential-n-tuples/
    """
    return zip(*[itertools.islice(lst, i, None, n) for i in range(n)])

def IsEven( x ):
  """Tests if x is even"""
  return ( x % 2 ) == 0

def ExtractOpts( opts, *optNamesAndDefaults ):
  """Returns a statement that when exec'ed extracts the given
  options from the given dictionary, with specified defaults
  when the options are not in the dictionary."""
  assert IsEven( len( optNamesAndDefaults ) )
  return '\n'.join( optName + ' = ' + repr( opts.get( optName, dfltVal ) )
                    for optName, dfltVal in groupIntoTuples( optNamesAndDefaults, 2 ) )


def RandomString( N ):
  """Return a random string of given length"""
  return ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(N))


#
# Bisect-related routines from http://docs.python.org/2/library/bisect.html
#

def bisect_index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def bisect_find_lt(a, x):
    'Find rightmost value less than x'
    i = bisect_left(a, x)
    if i:
        return a[i-1]
    raise ValueError

def bisect_find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError

def bisect_find_gt(a, x):
    'Find leftmost value greater than x'
    i = bisect_right(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def bisect_find_ge(a, x):
    'Find leftmost item greater than or equal to x'
    i = bisect_left(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def SymLinkRel( fromFN, toFN, getio = None ):
  """Create a symlink from fromFN to toFN using relative paths"""

  if getio: return dict( depends_on = toFN, creates = fromFN,
                         attrs = dict( piperun_short = True ) )
  
  os.symlink( os.path.relpath( toFN, os.path.dirname( fromFN ) ),
              fromFN )

def concatFiles( inFNs, outFN, getio = None ):
    """Concatenate specified files, save result to outFN"""

    inFNs = MakeSeq( inFNs )
    if getio: return dict( depends_on = inFNs, creates = outFN,
                           attrs = dict( piperun_short = True ) )
    SystemSucceed( 'rm -f ' + outFN )
    SystemSucceed( 'touch ' + outFN )
    for f in inFNs:
        SystemSucceed( 'cat ' + f + ( ' | bzip2 -c ' if outFN.endswith( 'bz2' ) else '' ) + ' >> ' + outFN )

# end: def concatFiles        

def zipFiles( inFNs, outFN, pathRelTo = None, getio = None ):
    """Zip specified files to given output archive"""

    inFNs = MakeSeq( inFNs )
    if getio: return dict( depends_on = inFNs, creates = outFN,
                           attrs = dict( piperun_short = True ) )
    SystemSucceed( 'rm -f ' + outFN )
    with zipfile.ZipFile( outFN, 'w' ) as out:
        for f in inFNs:
            arcName = f if pathRelTo is None else os.path.relpath( f, pathRelTo )
            out.write( f, arcName )

# end: def zipFiles

def tarFiles( inFNs, outFN, pathRelTo = None, getio = None ):
    """Tar specified files to given tar archive."""

    inFNs = MakeSeq( inFNs )
    if getio: return dict( depends_on = inFNs, creates = outFN,
                           attrs = dict( piperun_short = True ) )
    
    SystemSucceed( 'rm -f ' + outFN )
    with tarfile.open( outFN, 'w:bz2' ) as out:
        for f in inFNs:
            arcName = f if pathRelTo is None else os.path.relpath( f, pathRelTo )
            out.add( f, arcName )

# end: def tarFiles    

def DefineRulesTo_concatFiles( pr, inFNs, outFN, outBaseFN = None,
                               fromIdx = None, toIdx = None, ruleName = 'concatFiles' ):
    """Define rules to concatenate files"""

    if fromIdx is None: fromIdx = 0
    if toIdx is None: toIdx = len( inFNs )
    if outBaseFN is None: outBaseFN = outFN

    assert toIdx >= fromIdx

    if toIdx == fromIdx:
        pr.addRule( comment = 'concat', name = ruleName + Sfx( fromIdx, toIdx ),
                    depends_on = (), creates = outFN,
                    commands = 'touch ' + outFN )
    else:
        midIdx = int( ( fromIdx + toIdx ) / 2 )
        deps = []
        if midIdx > fromIdx:
            if midIdx == fromIdx+1:
                half1FN = inFNs[ fromIdx ]
            else:
                half1FN = AddFileSfx( outBaseFN, fromIdx, midIdx )
                DefineRulesTo_concatFiles( **Dict( 'pr inFNs fromIdx outBaseFN',
                                                   toIdx = midIdx, outFN = half1FN ) )
            deps.append( half1FN )

        if toIdx > midIdx:
            if toIdx == midIdx+1:
                half2FN = inFNs[ midIdx ]
            else:
                half2FN = AddFileSfx( outBaseFN, midIdx, toIdx )
                DefineRulesTo_concatFiles( **Dict( 'pr inFNs toIdx outBaseFN',
                                                   fromIdx = midIdx, outFN = half2FN ) )
            deps.append( half2FN )
            
        pr.addRule( comment = 'concat',
                    name = ruleName + Sfx( fromIdx, toIdx ),
                    depends_on = deps,
                    saveOutputTo = outFN,
                    commands = ' '.join( [ 'cat' ] + deps ) )
    
# end: def DefineRulesTo_concatFiles

if have_fcntl:
    class SimpleFlock(object):
       """Provides the simplest possible interface to flock-based file locking. Intended for use with the `with` syntax. It will create/truncate/delete the lock file as necessary.

       Adapted from https://github.com/derpston/python-simpleflock.git
       """

       def __init__(self, path, timeout = None, minCheckInterval = 0.1, maxCheckInterval = 10 ):
          self._path = path
          self._timeout = timeout
          self._fd = None
          self._minCheckInterval = minCheckInterval
          self._maxCheckInterval = maxCheckInterval
          self._rand = random.Random()
          self._rand.seed()

       def __enter__(self):
          self._fd = os.open(self._path, os.O_CREAT)
          start_lock_search = time.time()
          checkInterval = self._minCheckInterval
          while True:
             try:
                fcntl.flock(self._fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                # Lock acquired!
                return
             except IOError as ex:
                if ex.errno != errno.EAGAIN: # Resource temporarily unavailable
                   raise
                elif self._timeout is not None and time.time() > (start_lock_search + self._timeout):
                   # Exceeded the user-specified timeout.
                   raise

             # TODO It would be nice to avoid an arbitrary sleep here, but spinning
             # without a delay is also undesirable.
             if checkInterval < self._maxCheckInterval: checkInterval *= 2
             time.sleep(checkInterval + self._rand.random() )

       # end: def __enter__() 

       def __exit__(self, *args):
          fcntl.flock(self._fd, fcntl.LOCK_UN)
          os.close(self._fd)
          self._fd = None

          # Try to remove the lock file, but don't try too hard because it is
          # unnecessary. This is mostly to help the user see whether a lock
          # exists by examining the filesystem.
          # try:
          #    os.unlink(self._path)
          # except:
          #    pass

    # end: class SimpleFlock
    
def EncodeObjAsText( obj ):
    """Encode object as a text string"""

    pickledObj = pickle.dumps( obj )
    encodedObjOld = base64.urlsafe_b64encode( pickledObj )
    compressedObj = zlib.compress( pickledObj, 9 ) 
    encodedObj = ( '=' + base64.urlsafe_b64encode( compressedObj ) ) \
                  if len( compressedObj ) < len( pickledObj ) else encodedObjOld
    return encodedObj

def DecodeObjFromText( obj ):
    """Decode object from text string"""
    
    compressed = ( obj.startswith( '=' ) )
    if compressed: obj = obj[1:]
    unhexed = base64.urlsafe_b64decode( obj )
    if compressed: unhexed = zlib.decompress( unhexed )
    unpickled = pickle.loads( unhexed )
    return unpickled

    
def TupleIf( v1, v2 ):
    return () if not v1 else MakeSeq( v2 )

def tup( *args ): return tuple( args )    

def ShortenWithHash( s, maxLen ):
    """Shorten a long string by replacing its suffix beyond maxLen with a hash of the whole string.
    Useful e.g. when constructing a directory name, to avoid making it too long."""

    digest = '_' + hashlib.sha224( s ).hexdigest()
    if len( s ) + len( digest ) > maxLen:
        logging.debug( 'shortening: s=%s len(s)=%d maxLen=%d' % ( s, len(s), maxLen ) )
        s = s[ :max( 0, maxLen - len( digest ) ) ] + digest
        logging.debug( 'shortened: s=%s len(s)=%d maxLen=%d' % ( s, len(s), maxLen ) )
    return s
    
    
def ShortenPath( s ):
    """Shorten a path to comply with filesystem limits"""

    maxNameLen = os.pathconf( '.', 'PC_NAME_MAX' ) if 'PC_NAME_MAX' in os.pathconf_names else 255
    maxPathLen = os.pathconf( '.', 'PC_PATH_MAX' ) if 'PC_PATH_MAX' in os.pathconf_names else 4096

    components = []

    while True:
        d, n = os.path.split( s )
        n = ShortenWithHash( n, maxNameLen )
        components.append( n )
        s = d
        if len( s < 2 ) or s.endswith( '/' ):
            break

    return os.path.join( s, *reversed( components ) )
