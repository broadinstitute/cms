import hashlib, base64, inspect

from Operations.MiscUtil import dbg


def CodeContents( code ):
    """Return the contents of a code object as a string.  The contents includes the bytecode and all the constants used in it,
    but deliberately does not include any whitespace or comments.   So you can use this to see if the code of an operation has
    changed in some substantial way.
    """

    # One complication: the tuple of constants used in the code, code.co_consts, may include lambda objects.
    # These are converted to string as <code object <lambda> at 0x2a96c8af30 ... > rather than as their actual code.
    # So we must recursively convert any such "constants" to strings, when constructing a string representation
    # of this code object.
    
    constsStr = str( [x if not inspect.iscode( x ) else CodeContents( x ) for x in code.co_consts] )
    return code.co_code + constsStr + str( code.co_argcount ) + str( code.co_flags ) + str( code.co_nlocals ) + str( code.co_names ) + str( code.co_varnames ) 
    

def FunctionChecksum( code, func_obj ):
    """Take a code object and (optionally) its function object, and return two checksums.
    The first represents the code, and the second represents the default args.
    """

    assert not func_obj or func_obj.__code__ is code

    funcContents = CodeContents( code )

    funcDefaults = ( str( func_obj.__defaults__ ) + str( func_obj.__dict__ ) ) if func_obj else ''
    funcContentsChecksum = base64.urlsafe_b64encode( hashlib.sha512( funcContents ).digest() )
    funcDefaultsChecksum = base64.urlsafe_b64encode( hashlib.sha512( funcDefaults ).digest() )

    return funcContentsChecksum, funcDefaultsChecksum 

