from __future__ import with_statement

from contextlib import contextmanager, nested
import os, operator, itertools, glob
from Operations.MiscUtil2 import SlurpFileLines, chomp, dbg, SystemSucceed, AddSlash
from Operations.IDotData import IDotData


"""Code for manipulating large DotData files without loading them fully into memory.
"""

def VStackDotDataFiles( inFiles, outFile, getio = None ):
    """Vertically stack the specified DotData or tsv files"""

    if getio: return dict( depends_on = inFiles, creates = outFile )

    IDotData.vstackFromIterable( itertools.imap( IDotData, inFiles ) ).save( outFile )

def VStackTableFiles( inFiles, outFile, fileIdColumn = None, fileIdVals = None, getio = None ):
    """Vertically stack the specified DotData or tsv files"""

    if getio: return dict( depends_on = inFiles, creates = outFile,
                           attrs = dict( piperun_short = True ) )

    assert not fileIdColumn or len( fileIdVals ) == len( inFiles )

    with open( outFile, 'w' ) as outF:

        needNewline = False
        theHeader = None
        
        for inFileNum, inFile in enumerate( inFiles ):
            
            assert inFile.endswith( '.tsv' ) # for now
            
            with open( inFile ) as inF:

                dbg( 'inFile' )
                header = inF.readline()
                if not header.endswith( '\n' ): header += '\n'
                if inFileNum == 0:
                    if fileIdColumn: outF.write( fileIdColumn + '\t' )
                    outF.write( header )
                    theHeader = header
                else:
                    if header != theHeader:
                        print 'theHeader=', theHeader
                        print 'header=', header
                    assert header == theHeader
                    
                if needNewline: outF.write( '\n' )
                needNewline = False
                for inLine in inF:
                    if fileIdColumn: outF.write( fileIdVals[ inFileNum ] + '\t' )
                    outF.write( inLine )
                    needNewline = not inLine.endswith( '\n' )

                    
def linkDotDatas( source, target, getio = None ):
    """Create a hard link to an existing DotData under a new name"""

    source = os.path.abspath( source )
    target = os.path.abspath( target )

    if getio: return dict( depends_on = AddSlash( source ),
                           creates = AddSlash( target ) )

    if source.endswith('/'): source = source[:-1]
    if target.endswith('/'): target = target[:-1]

    print 'linking ' + source + ' to ' + target
    SystemSucceed( 'mkdir -p ' + target )
    for f in os.listdir( source ):
        os.link( os.path.join( source, f ),
                 os.path.join( target, os.path.basename( f if not f.endswith( '.header.txt' ) else os.path.splitext( os.path.basename( target ) )[0] + '.header.txt' ) ) )
