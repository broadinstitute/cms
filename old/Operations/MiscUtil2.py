
import inspect, logging, sys, os

def SlurpFile( fname ):
    """Read entire file into one string"""
    if not os.path.isfile( fname ): raise IOError( 'File not found: %s' % fname )
    with open( fname ) as f:
        return f.read()

def SlurpFileLines( fname ):
    """Read entire file and return a newly allocated array of its lines"""
    return SlurpFile( fname ).strip( '\n' ).split( '\n' )

def chomp( s ):
    """Return string without terminal newline if present"""
    return s[:-1] if s.endswith('\n') else s


def dbg( *exprs ):
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

    logging.debug( 'at ' + caller_filename + ':' + str(caller_lineno) + ' in ' + caller_function + ', ' +
                   ', '.join([ e + ' (' + str(type(eval( e, callerGlobals, callerLocals ))) + ') = ' + str(eval( e, callerGlobals, callerLocals )) for ex in exprs for e in ( ex.strip().split(' ') if not ex.startswith('#') else (ex[1:],) ) ] ) )

    del callerGlobals, callerLocals

class GeneratorWrap(object):

    """Wrap a generator to create an object that can have attributes
    assigned to it.
    """

    def __init__(self, gen):
        self.gen = gen

    def __iter__(self): return self
    def next(self):
        return self.gen.next()



def SystemSucceed( cmd, dbg = False ):
    """Run a shell command, and raise a fatal error if the command fails."""
    logging.info( 'Running command ' + cmd + ' ; called from' + sys._getframe(1).f_code.co_filename + ':' +
                  str( sys._getframe(1).f_lineno ) )
    exitCode = os.system( cmd )
    logging.info( 'Finished command ' + cmd + ' with exit code ' + str( exitCode ) )
    if exitCode != 0:
        raise IOError( "Command %s failed with exit code %d" % ( cmd, exitCode ) )

def AddSlash( dirName ):
  """Add a slash to the dir name if not there already"""
  return dirName if dirName.endswith('/') else dirName + '/'
