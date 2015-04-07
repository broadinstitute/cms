
"""Code for low-memory processing of tabular data.
This code works with IDotDatas, which are immutable streams of records.   An IDotData has a fixed number of named columns,
and lets you take out iterators to iterate sequentially over the records (table rows).  Unlike DotData, it does not provide
efficient random access to the rows, and does not allow changing the table.  In return, IDotData operations require only
a small constant amount of memory, compared to DotData's linear-memory requirements.

An IDotData represents an immutable stream of records from a table with a fixed number of named columns.
It must have the following attributes:

   headings - a tuple of names of the table columns.
   __iter__() - returns an iterator over the records. 

      It must yield immutable records, which are named tuples of field values.

notes:

  comparisons between IDotDatas: equality comparisons produce an IDotData for a pairwise comparisons of rows.
  Boolean testing tests if there are any rows.
  so, in particular, bool( a == b ) is true as long as they're both nonempty (even if different),
  and is false even if they're both equally empty (because the answer is an empty sequence of booleans).
  the headings also do not play a role in comparisons.

  given an iter you can make a routine that returns this iter once and fails on any other attempt,
  which turns an iterator into a container.

  - so, imap does not work here.

  we need to create a new IDotDataRoot() object,
  with the specified headings,

  which, whenever asked, for an iterator, will ask _us_ for an iterator,
  and return an object which takes items from our iterator and selects from them.

  spokojno, ehto nuzhno sdelat' akkuratno.


  - so, forget abstraction, implement what you want, what each method should do, even if with code/object duplication;
  then, when everything is correct, re-abstract -- it will be more clear how to abstract.
  but you already know what each method should _do_.

  >>> z = IDotData( names = ( 'a', 'b' ), Records = ( ( 1, 2 ), ( 1, 3 ), ( 2, 4 ), ( 2, 5 ) ) )
  ... # doctest: +NORMALIZE_WHITESPACE
  >>>

  """

from __future__ import with_statement, division

__all__ = ( 'IDotData', )

import sys, os, logging, itertools, operator, copy, contextlib, time, numbers, types, glob, math, collections, \
    heapq, inspect, abc, __builtin__
import traceback as tb
from abc import ABCMeta, abstractmethod
from Operations.MiscUtil import chomp, dbg, coerceVal, is_sorted, IsSeq, MakeSeq, MakeDir, DumpFile, SlurpFileLines, \
    WaitForFileToAppear, joinstr, Sfx, IsScalar, DictGet, flatten, SystemSucceed, ReplaceFileExt, BreakString, joinstr, \
    RestrictDict, MergeDicts, tabwriten, AtomicForIsSeq, DbgIter, Dict, OpenForRead, OpenForWrite, SumKeeper, tmap, \
    MapPath, DictGetNotNone, tabwrite, AddFileSfx, SplitStr, FirstVal, ExtractOpts, IsFileType, iter_ith

try:
    import numpy as np
    haveNumpy = True
except ImportError:
    haveNumpy = False

try:
    from Classes.DotData import DotData
    haveDotData = True
except ImportError:
    haveDotData = False

class IDotDataRecord(collections.Sequence):
    
    """Represents an immutable named tuple (i.e. where you can access tuple components by name as well
    as by index).  Very similar to collections.namedtuple, but column names can be arbitrary strings and
    need not be valid Python identifiers.
    IDotDataRecord represents one record (table row) of an IDotData, returned by an iterator taken out
    over IDotData.  
    
    Fields:

    colName2num - map from column name to integer column number.  Used for accessing tuple
        elements by column name.  This map is created once for each .tsv or .data/ file,
        and is shared by all IDotDataRecords read from that file.

    lineVals - string values of the fields of this record.  We do not convert them
      to non-string datatypes until asked -- because often, they're just passed around as-is,
      so not converting them from string to int or float and then back saves time.

    fname - name of the file we are reading; needed for error reporting only.
    
    headings - names of the columns.  needed for error reporting.
    
    """
    
    def __init__(self, lineVals, colName2num, headings, fname = None ):
        
        self.lineVals = tuple( lineVals )
        self.colName2num = colName2num
        self.headings = headings
        assert len( self.lineVals ) == len( self.headings ) == len( self.colName2num )
        
    def __getattr__(self, name):
        """Get value of one element of this record, by its field name"""

        if name.startswith( '__' ): return object.__getattr__( self, name )
        try: return coerceVal( self.lineVals[ self.colName2num[ name ] ] )
        except KeyError:
            raise AttributeError(name)
                
    def __getitem__(self, key):
        """Get value of one field of this named tuple, either by field name or by position (column number).
        The key can be a sequence of keys, in which a tuple of field values is returned.
        """

        if IsSeq( key ): return tuple( self[k] for k in key )
        return  coerceVal( self.lineVals[ self.colName2num[ key ] if isinstance( key, types.StringTypes ) else key ] )

    def GetStrItem(self, key):
        """Get value of one field of this named tuple, either by field name or by position (column number).
        The key can be a sequence of keys, in which a tuple of field values is returned.
        The field is left as a string and is not coerced.
        """

        if IsSeq( key ): return tuple( self[k] for k in key )
        return  self.lineVals[ self.colName2num[ key ] if isinstance( key, types.StringTypes ) else key ]

    
    def __iter__(self): return itertools.imap( coerceVal, self.lineVals )
    def __repr__(self): return str( zip( self.headings, self.lineVals ) )
    def __len__(self): return len(self.headings)

    def __add__(self,other):
        """Return a concatenation of the two records"""
        return IDotDataRecord( lineVals = self.lineVals + other.lineVals,
                               colName2num = MergeDicts( self.colName2num, other.colName2num ),
                               headings = self.headings + other.headings )
    
    def __eq__(self,other):
        other = MakeSeq( other )
        return len(self) == len(other) and all([ a == b for a,b in itertools.izip( self, other ) ])

    def __ne__(self,other): return not self == other
    def __lt__(self,other):
        other = MakeSeq( other )
        return any( a < b for a,b in itertools.izip( self, other ) ) \
            or len(self) < len(other) and all([ a == b for a,b in itertools.izip( self, other ) ])

    def __le__(self,other): return self == other or self < other
    def __gt__(self,other): return not self <= other
    def __ge__(self,other): return not self < other

    def __hash__(self): return hash(tuple(lineVals))

    def asDict(self): return dict( ( n, v ) for n, v in zip( self.headings, self ) )  


class my_dtype(object):

    """For compatibility with numpy arrays which have a dtype field, we define a class that mimics the dtype interface
    and add a field of this type to IDotData."""
    
    def __init__(self, names ):
        self.names = names

    
class IDotDataRoot(collections.Iterable, AtomicForIsSeq):

    __metaclass__ = ABCMeta

    """Abstract base class for various IDotData implementations.
    The implementation (a sublcass of IDotDataRoot) provides the tuple of headings when it calls
    IDotDataRoot's __init__() method, and provides the _iter() method which returns an iterator yielding
    the records as tuples.  IDotDataRoot then adds many useful operations for manipulating IDotDatas.

    Fields:

      headings - names of columns of this IDotData.  the names can be arbitrary strings and need not be valid Python
        identifiers.

      colName2num - map from column name to its position in headings.  used for fast lookup of fields by their name.

      parents - IDotDatas from which this IDotData's data is derived, if any.

    notes:

      - should we _require_ that our inputs are repeatable collections and not just iters?
      that's not the common case though...

      - have the method, IDotData_mergeOnKeyCols() and manage other static methods that way.

      alternately, assign this as an attribute to the IDotData function object.

      
    
    """

    # Mark IDotData as not hashable
    __hash__ = None

    isIDotData = True

    @abstractmethod
    def _iter(self):
        """Return the next record, as either a tuple or a single value or an IDotDataRecord """
        pass

    def __init__(self, headings, parents = () ):

        headings = tuple( MakeSeq( headings ) )
        #assert headings
        if len( set( tuple( headings ) ) ) != len( headings ):
            dbg( 'headings' )
        assert len( set( tuple( headings ) ) ) == len( headings ), ( 'headings=%s' % str( headings ) )
        self.headings = headings
        self.colName2num = dict( ( colName, colNum ) for colNum, colName in enumerate( self.headings ) )
        self.parents = parents
        self.isInfinite = any( p.isInfinite for p in parents )
        self.len = None

        self.dtype = my_dtype( names = headings )

        self.creationTraceback = tb.extract_stack()

    def __getitem__(self,item):
        """Returns a sub-rectangle of this IDotData.   What is returned depends on the type of item:

           - if item is a string, it must be a column name, and we return a one-column IDotData consisting of just that
           column from this IDotData

           - if the item is a sequence of strings, each string must be a column name, and we return an IDotData
           consisting of these columns of this IDotData

           - if the item is a sequence of Booleans, we return an IDotData consisting of those rows of this IDotData
           for which item is True.  
        """

        if all([ isinstance( it, types.StringTypes ) for it in MakeSeq( item ) ]):
            for it in MakeSeq( item ):
                if it not in self.headings: raise AttributeError( item )
            return IDotDataGetColumns( parent = self, item = item )
        if isinstance( item, ( types.IntType, types.LongType ) ): return iter_ith( self, item )
        if isIDotData( item ): return IDotDataFilterBool( parent = self, filter = item )
        if hasattr( type( item ), 'isDotData' ): return self[ IDotData.fromDotData( item ) ] 
        asDotData = self.toDotData()[ item ]
        if hasattr( type( asDotData ), 'isDotData' ): return IDotData.fromDotData( asDotData )
        return IDotData( headings = self.headings, Columns = ( asDotData, ) )
                
    def __getattr__(self,name):
        """Return a given column"""
        try: return collections.Iterable.__getattr__(self,name)
        except AttributeError: return self[name] if name != '__setstate__' else None

    def __iter__(self):

        """Returns an iterator over the records of this IDotData.  If this IDotData consists of exactly one column,
        yields values from that column; otherwise, yields named tuples of values.
        """

        singleColumn = ( len( self.headings ) == 1 )
        
        for r in self._iter():
            r = MakeSeq( r )
            if len( r ) != len( self.headings ):
                dbg( 'len(r) len(self.headings) r tuple(r) self.headings zip(self.headings,r)' )
            assert len( r ) == len( self.headings )
            if singleColumn: yield coerceVal( r[0] )
            else:
                yield r if isinstance( r, IDotDataRecord ) else \
                    IDotDataRecord( lineVals = r, colName2num = self.colName2num, headings = self.headings )

    def recordsIter(self):

        """Returns an iterator over the records of this IDotData.  This iterator always yields IDotDataRecords,
        even when the IDotData consists of exactly one column.
        """

        return iter( self ) if len( self.headings ) > 1 else self.wrapIter()

    def wrapIter(self):
        assert len( self.headings ) == 1
        for r in self:
            yield IDotDataRecord( lineVals =  ( r, ), colName2num = self.colName2num, headings = self.headings )

    def flatIter(self):
        """Returns an iterator that iterates over all values in all records of this IDotData, returning a flat
        stream of just the values."""
        for r in self:
            for v in MakeSeq( r ):
                yield v

#    def argsort(self, *cols):
#        """Returns the tuple that would sort this IDotData by the given columns.
#        Note that this currently requires reading the whole thing into memory."""

#        assert len( cols ) == 1
#        return np.array( tuple( self[ cols[0] ] ) ).argsort()

    def argsort(self, *args, **kwargs): return np.argsort( np.asarray( self ), *args, **kwargs )

    def getColNum( self, col ):
        """Returns the column number from column name"""
        return self.colName2num[ col ]
    
    def sortedOn(self, *cols):
        """Return a version of this IDotData that is sorted on the specified columns
        (the first column being the primary key, the second column the secondary key, etc).
        Note that at the moment, the results in instantiating the full IDotData in memory.
        Later we'll add options for external-memory sorting.
        """
        return self[ np.lexsort( [ self[c] for c in reversed( cols ) ] ) ]

    sortedBy = sortedOn

    def numCols(self):
        """Return number of columns in this IDotData."""
        return len(self.headings)

    def takewhile(self, pred):
        """Returns an IDotData consisting of those first rows of this IDotData for which the callable pred returns True.
        """
        return IDotDataTakeWhile( parent = self, pred = pred )

    def oneSlice( self, colNames, partNum, totalParts, numChunks = None):
        """Return one slice of an IDotData split by specified column(s)"""
        return IDotDataOneSlice( parent = self, **Dict( 'colNames partNum totalParts numChunks' ) )

    def oneRegion( self, keyCols, beg, end ):
        """Return one region of an IDotData split by specified column(s)"""
        return IDotDataOneRegion( parent = self, **Dict( 'keyCols beg end' ) )

    
    # versions for all other itertools methods, esp. groupby.

    def hstack(self, *iDotDatas):
        """Return the horizontal (side-by-side) stacking of this IDotData and specified IDotDatas"""
        return IDotDataHStack( self, *iDotDatas )

    def addCols(self, names, cols):
        """Return this IDotData with horizontally stacked columns added."""
        return self.hstack( IDotData( names = names, Columns = cols ) )

    def addCol(self, name, col): 
        """Return this IDotData with one horizontally stacked column added."""
        return self.addCols( names = (name,), cols = (col,) )

    def vstack(self, *iDotDatas, **kwargs):
        """Return the vertical stacking of this IDotData and specified IDotDatas"""
        return IDotDataVStack( self, *iDotDatas, **kwargs )

    def vstackFromIterable(self, iDotDatas, headings = None):
        """Return the vertical stacking of this IDotData and specified IDotDatas"""
        return IDotDataVStackFromIterable( self, iDotDatas = iDotDatas, headings = headings )
    
    def filter(self, pred):
        """Return an IDotData consisting of those rows of this IDotData for which the specified predicate returns True"""
        return IDotDataFilter( self, pred )

    def removeDups(self, keyCols = None ):
        """Return an IDotData consisting of non-duplicate rows of this IDotData.  Rows are considered
        duplicate if they're adjacent and equal at column(s) keyCols.  Of each group of duplicate rows, the first
        is selected to be part of the new IDotData."""

        return IDotDataRemoveDups( self, keyCols )
        

    def firstWhere(self, pred):
        """Return the first row, that matches the given condition.
        If no such row, raises IndexError."""
        for r in self:
            if pred( r ): return r
        raise IndexError( 'No row matches the specified condition!' )

    def mapRecords(self, func, headings = None):
        """Return an IDotData whose records (rows) are obtained by applying a transformer function to the records (rows)
        of this IDotData."""
        if not headings: headings = self.headings
        return IDotData.mapRecords( func = func, iDotDatas = (self,), headings = headings )

    def mapVals(self, f):
        """Return an IDotData with all values mapped by the given function"""
        return self.mapRecords( ( lambda x: f(x) ) if self.numCols() == 1 else ( lambda r: map( f, r ) ) )

    def addComputedCols( self, newColNames, newColFn):
        """Add a column computed from the records according to specified function"""

        if isinstance( newColNames, types.StringTypes ) and ( ' ' in newColNames or '\t' in newColNames ):
            newColNames = tuple( newColNames.split( '\t' if '\t' in newColNames else None ) )
        
        return self.mapRecords( func = lambda r: tuple( r ) + tuple( MakeSeq( newColFn( r ) ) ),
                                headings = self.headings + tuple( MakeSeq( newColNames ) ) )

    def replaceCol( self, col, newColFn ):
        """Create a new IDotData where values in given col are replaced according to a given function on records"""
        colNum = self.colName2num[ col ]
        def mapRec( r ):
            t = tuple( r )
            return tuple( t[ :colNum ] ) + tuple( MakeSeq( newColFn( r ) ), )  + tuple( t[ colNum+1: ] )
        return self.mapRecords( func = mapRec )

    def starmap(self, func, headings = None):
        if not headings: headings = self.headings
        return IDotData.starmap( func = func, iDotData = self, headings = headings )
    
    def renameCols(self, renamings):
        """Get a version of this IDotData with some columns renamed.
        
        Parameters:

           renamings - map from old column name to new column name.  Columns not renamed
              keep their old name.
        """
        return IDotDataRenameCols( parent = self, renamings = renamings )

    def columns(self, cols = None):
            """Return the columns, as a list.  Convenient for mapping a function over the columns."""
            if cols == None: cols = self.headings
            return [ self[ colName ] for colName in cols ]
    
    def __doOp(self, other, func, flip = False, headings = ( 'V', ) ):
        """Return an IDotData whose rows are obtained by applying the specified binary function to the corresponding
        rows of this IDotData and another IDotData.

        If 'other' is a scalar value, first promote it to a one-column IDotData consisting of that value repeated.
        
        """

        if IsScalar( other ): return self.__doOp( IDotData.repeat( heading = 'N', value = other ), func, flip, headings )
        if not isIDotData( other ): return NotImplemented
        return IDotData.mapRecords( func = func, iDotDatas = ( self, other ) if not flip else ( other, self ),
                                    headings = headings )

    def __doROp(self, other, func, headings = ( 'V', ) ):
        return self.__doOp( other, func, flip = True, headings = headings )
    
    def __ge__(self, other): return self.__doOp( other, operator.ge )
    def __gt__(self, other): return self.__doOp( other, operator.gt )
    def __le__(self, other): return self.__doOp( other, operator.le )
    def __lt__(self, other): return self.__doOp( other, operator.lt )
    def __eq__(self, other): return self.__doOp( other, operator.eq )
    def __ne__(self, other): return self.__doOp( other, operator.ne )
    
    def __add__(self, other): return self.__doOp( other, operator.add )
    def __sub__(self, other): return self.__doOp( other, operator.sub )
    def __mul__(self, other): return self.__doOp( other, operator.mul )
    def __div__(self, other): return self.__doOp( other, operator.div )
    def __truediv__(self, other): return self.__doOp( other, operator.truediv )
    def __floordiv__(self, other): return self.__doOp( other, operator.floordiv )
    def __and__(self, other): return self.__doOp( other, operator.and_ )
    def __or__(self, other): return self.__doOp( other, operator.or_ )
    def __xor__(self, other): return self.__doOp( other, operator.xor )
    def __lshift__(self, other): return self.__doOp( other, operator.lshift )
    def __rshift__(self, other): return self.__doOp( other, operator.rshift )
    def __mod__(self, other): return self.__doOp( other, operator.mod )
    def __divmod__(self, other): return self.__doOp( other, divmod, headings = ( 'div', 'mod' ) )
    def __pow__(self, other): return self.__doOp( other, operator.pow )
    
    def __radd__(self, other): return self.__doROp( other, operator.add )
    def __rsub__(self, other): return self.__doROp( other, operator.sub )
    def __rmul__(self, other): return self.__doROp( other, operator.mul )
    def __rdiv__(self, other): return self.__doROp( other, operator.div )
    def __rtruediv__(self, other): return self.__doROp( other, operator.truediv )
    def __rfloordiv__(self, other): return self.__doROp( other, operator.floordiv )
    def __rmod__(self, other): return self.__doROp( other, operator.mod )
    def __rdivmod__(self, other): return self.__doROp( other, divmod, headings = ( 'div', 'mod' ) )
    def __rpow__(self, other): return self.__doROp( other, operator.pow )
    def __rand__(self, other): return self.__doROp( other, operator.and_ )
    def __ror__(self, other): return self.__doROp( other, operator.or_ )
    def __rxor__(self, other): return self.__doROp( other, operator.xor )
    def __rlshift__(self, other): return self.__doROp( other, operator.lshift )
    def __rrshift__(self, other): return self.__doROp( other, operator.rshift )

    def __neg__(self): return self.mapRecords( operator.neg )
    def __pos__(self): return self.mapRecords( operator.pos )
    def __abs__(self): return self.mapRecords( operator.abs )
    def __invert__(self): return self.mapRecords( operator.invert )

    def __len__(self):
        """Return the number of data rows in this IDotData.   Note that the first time this is called,
        it may take linear time to run."""
        if self.len is None: self.len = sum( itertools.imap( lambda x: 1, self ) )
        return self.len

    def __nonzero__(self):
        """Return True if self has at least one row, False otherwise"""

        try:
            next( iter( self ) )
            return True
        except StopIteration: return False
    
    #def __contains__(self, val):
    #    return 1

    #def __neg__(self): return self.

    def dbg(self, msg = None):
        print '=============================='
        if msg is not None:
            print '\n' + msg + '\n'
            print '=============================='
        tabwrite( sys.stdout, *self.headings )
        nlines = 0
        for line in self:
            tabwrite( sys.stdout, *line )
#            sys.stdout.write( '\n' + '\t'.join( map( str, line )
#                                                if ( hasattr(line,'__iter__') and not isinstance(line,types.StringTypes) )
#                                                else (str(line),) ) )
            nlines += 1
        sys.stdout.write( '\n%d rows\n' % nlines )
        print '=============================='

    def dbgGraph(self, fname):
        """Write a graph representation of the current object"""

        with open( fname, 'w' ) as f:
            f.write( 'digraph A {\n' )

            nodes = [ self ]

            while nodes:
                n = nodes.pop()

                f.write( 'n%d' % id( n ) );
                f.write( '[ shape = "box", label = "%s" ];\n' % '\\n'.join( BreakString( str( n ), 20 ) ) )
                for p in n.parents:
                    nodes.append( p )
                    f.write( 'n%d -> n%d\n' % ( id( p ), id( n ) ) )
                
            f.write( '}\n' )

        SystemSucceed( 'dot -Tsvg -o' + ReplaceFileExt( fname, '.svg' ) + ' ' + fname )

    def save(self, fname, useHeadings = True, fileType = None, sep = '\t', comments = () ):
        """Save this IDotData to either a .tsv file or to a .data/ directory, depending on the extension"""
        with IDotData.openForWrite( fname, self.headings if useHeadings else None, **Dict( 'sep fileType comments' ) ) as f:
            for line in self: f.writeRecord( line )
        return self

    saveToSV = save
    saveToDotData = save


    def toDotData(self):
        """Convert to DotData"""
        assert haveDotData
        return DotData.fromIDotData( self )
    
    def __array__(self):
        """Convert this to a numpy.ndarray.   This lets many numpy methods take
        IDotDatas as arguments."""
        return np.array( tuple( self ) ) if self.numCols() == 1 else self.toDotData()

    def nanmax(self):
        isFirst = True
        for v in self.flatIter():
            if isFirst or v > result: result = v
            isFirst = False
        return result if not isFirst else float('nan')

    def groupby(self, *cols, **kwargs):
        """Yield groups of rows of this IDotData, grouped by values in the given columns, as individual IDotDatas.

        Params:

          cols - column name or sequence of column names of this IDotData.
             Adjacent rows for which these columns have a given
             combination of values will be grouped together in the output.

          multiPass - if True ( default ), the IDotDatas yielded by the returned iterator (see 'Returns' below)
             may be iterated over multiple times; if False, they can be iterated over only once (but this takes less
             memory).  If an integer k, we yield tuples of the form (key, idd_1, idd_2, ..., idd_k) where each of
             idd_i is a single-pass IDotData over the grouped rows.

        Returns:

           an iterator yielding pairs of values; the first value is the value or tuple of values of the key columns;
           the second is an IDotData representing the next group of adjacent rows with a particular combination of
           values in the key columns.
         
        """

        multiPass = DictGet( kwargs, 'multiPass', True )
        if not isinstance( multiPass, bool ):
            for v in itertools.izip( *[ self.groupby( *cols, multiPass = False ) for i in range( multiPass ) ] ):
                yield ( v[0][0], ) + tmap( operator.itemgetter( 1 ), v )
            return

        for k, g in itertools.groupby( self, operator.itemgetter( cols ) ):
            yield ( k[0] if isinstance( k, tuple ) and len( k ) == 1 else k ), \
                IDotData.fromIterable( headings = self.headings,
                                       iterable = tuple( g ) if multiPass else g )


    def selectRowsUsingOneCriterion(self, expr):
        """Return selected rows.

        If expr is a string or a code object, returns rows for which expr is True.  Expr is
        either a string or a code object; if a string, it is compiled into a code object.  The
        expression is evaluated for each row, in an environment in which column names are
        bound to column values in that row. A DotData containing copies of rows for which
        the expression evaluates to True is returned.

        If expr is a map of column names to values, returns rows for which the specified
        columns have the specified values. Column names mapped to None do not
        participate in the filtering.

        You can get the same effect by using self[expr], but the [] operator is getting very
        overloaded...


        See also: selectRows(), which can select rows according to a list of criteria.

        """

        if isinstance(expr, types.CodeType) or isinstance(expr, types.StringTypes):
            exprCode = expr if isinstance(expr, types.CodeType) \
                else compile( expr, sys._getframe(0).f_code.co_filename, 'eval' )

            # the expression may reference only the column names (as python variables),
            # or it may reference the caller's local and global vars.
            # check which is the case here, and if needed get the caller's
            # local and global variable mappings.
            callerGlobals = { }
            callerLocals = None
            exprVars = frozenset( exprCode.co_names )
            if not exprVars <= frozenset( self.dtype.names ):
                callerFrame = inspect.currentframe()
                thisFileName = callerFrame.f_code.co_filename
                while ( callerFrame.f_code.co_filename == thisFileName ):
                        callerFrame = callerFrame.f_back
                callerLocals = RestrictDict( callerFrame.f_locals, exprVars )
                callerGlobals = RestrictDict( callerFrame.f_globals, exprVars )
                del callerFrame

            def evalRow( rec, cGlobals = callerGlobals, cLocals = callerLocals ):
                r = dict( ( n, v ) for n, v in zip( self.headings, rec ) )
                return eval( exprCode, cGlobals,
                             r if not cLocals else MergeDicts( cLocals, r ) )

            result = self.filter( evalRow )

            del callerGlobals, callerLocals
            return result

        elif callable(expr):
            return self.filter( expr )
        elif isinstance(expr, dict):
            exprItems = frozenset( [ (columnName, val) for columnName, val in expr.items() if val != None ] )
            return self.filter( lambda r: exprItems <= frozenset( r.asDict().items() ) )
        else:
            raise TypeError('selectRows() requires either an expression (string or code object) or a dictionary')

    def selectRows(self, expr, *exprs):
        """Returns a copy of rows matching all criteria in a list.  See selectRowsAux()
        for documentation of the criteria.

        This method takes a variable number of arguments; each argument is a
        separate criterion, and rows are returned which match ALL the criteria.

        See also: __getitem__().
        """
        result = self;
        for e in (expr,) + exprs:
            result = result.selectRowsUsingOneCriterion(e)
        return result	


    def getColumnStats(self, *cols):
        """For each of the specified columns, compute some statistics.
        Return the results as a new IDotData with one row per column in 'cols',
        and the headings 'col', 'mean', 'std', 'n', 'ntot'.  The meanings of the
        columns are:

           n - the number of non-NaN values in the column; these were the values
              used for the mean and stddev computation.
              
           ntot - the total number of values seen, including NaNs.  same as
              self.numRows().
           """

        if len( cols ) == 1 and IsSeq( cols[0] ): cols = cols[0] 

        meansStds = imeanstd_plusStats( self[ cols ] )
        if len( cols ) == 1: meansStds = [ meansStds ]

        return IDotData( names = 'col mean std n ntot',
                         Records = [ ( col, mean, std, n, ntot )
                                     for col, ( mean, std, n, ntot ) in zip( cols, meansStds ) ] )
            
    def normalizeColumns(self, *cols, **kwargs):
        """Return a version of this IDotData with the specified columns normalized and the remaining columns unchanged"""

        cols = flatten( cols )

        if 'meansStds' in kwargs: meansStds = kwargs[ 'meansStds' ]
        else:
            meansStds = imeanstd( self[ cols ] )
            if len( cols ) == 1: meansStds = [ meansStds ]
        
        col2meanStd = dict( zip( cols, meansStds ) )

        colMeanStds = [ DictGet( col2meanStd, c, ( None, None ) ) for c in self.headings ]

        return self.mapRecords( lambda r: [ v if mean is None else ( ( v - mean ) / std if std > 0 else np.nan )
                                            for v, ( mean, std ) in zip( r, colMeanStds ) ] )

    def normalizeColumnsWithinGroups(self, cols, groupCols, normedSfx = None, col2normed = {}, meanSfx = None,
                                     stdSfx = None, countSfx = None ):
        """Return a version of this IDotData with the specified columns normalized within row groups defined by groupCols
        (each group of adjacent rows that have an identical combination of values in columns 'groupCols' is one row group),
        and the remaining columns unchanged.  If col2normed is given, it specifies which columns' normalized versions
        should be added as new columns rather than replacing the old columns."""

        cols = flatten( MakeSeq( cols ) )
        groupCols = flatten( MakeSeq( groupCols ) )

        colNums = [ self.colName2num[ c ] for c in cols ]

        # Set up a place to keep the running sum and sum-of-squares for each
        # column to be normalized.  These are allocated once and then reused
        # within each row block.
        sums = [ None ] * len( cols )
        sumsSq = [ None ] * len( cols )
        for i in range( len( cols ) ):
            sums[ i ] = SumKeeper()
            sumsSq[ i ] = SumKeeper()
        counts = np.zeros( len( cols ), dtype = int )

        colsSet = frozenset( cols )

        # var: isNormeds - for each of the original columns, whether
        #   a normalized version of this column will be computed.
        isNormeds = [ h in colsSet for h in self.headings ]

        colMeans = np.zeros( len( self.headings ) )
        colStds = np.zeros( len( self.headings ) )

        extraCols = []
        extraColsFrom = []
        for colNum, c in zip( colNums, cols ):
            customName = ( c in col2normed and col2normed[ c ] != c )
            if normedSfx or customName:
                extraCols.append( col2normed[ c ] if customName else c + Sfx( normedSfx ) )
                extraColsFrom.append( colNum )
         
        def GetNormedRows():

            yield self.headings + tuple( extraCols ) \
                + ( tuple( [ c + Sfx( meanSfx ) for c in cols ] ) if meanSfx else () ) \
                + ( tuple( [ c + Sfx( stdSfx ) for c in cols ] ) if stdSfx else () ) \
                + ( tuple( [ c + Sfx( countSfx ) for c in cols ] ) if countSfx else () )
            
            for groupId, rows1, rows2 in self.groupby( *groupCols, multiPass = 2 ):

                #
                # Determine the mean and stddev for each column on which we're normalizing
                #

                for i in range( len( cols ) ):
                    sums[ i ].clear()
                    sumsSq[ i ].clear()
                counts.fill( 0 )

                for r in rows1:
                    for i, colNum in enumerate( colNums ):
                        x = r[ colNum ]
                        if not ( math.isnan( x ) or math.isinf( x ) ):
                            sums[ i ] += x
                            sumsSq[ i ] += x*x
                            counts[ i ] += 1

                for i, colNum in enumerate( colNums ):
                    if counts[ i ]:
                        colMeans[ colNum ] = sums[ i ].getSum() / counts[ i ]
                        colStds[ colNum ] = \
                            math.sqrt( sumsSq[ i ].getSum() / counts[ i ] - colMeans[ colNum ]*colMeans[ colNum ] )
                    else:
                        colMeans[ colNum ] = np.nan
                        colStds[ colNum ] = np.nan
                        
                # Now normalize the rows in this group using the means and stds determined above.

                for r in rows2:
                    r_new = [ x if not isNormed
                              else ( np.nan if math.isnan( colMean ) or colStd == 0.0
                                     else ( x - colMean ) / colStd ) for x, isNormed, colMean, colStd in
                              zip( r, isNormeds, colMeans, colStds ) ]
                    newVals = []
                    for ecFrom in extraColsFrom:
                        newVals.append( r_new[ ecFrom ] )
                        r_new[ ecFrom ] = r[ ecFrom ]
                    yield r_new + newVals + ( [ colMeans[i] for i in colNums ] if meanSfx else [] ) + \
                        ( [ colStds[j] for j in colNums ] if stdSfx else [] ) + \
                        ( [ counts[k] for k in range( len( colNums ) ) ] if countSfx else [] )
                                    
        return IDotData.fromFn( GetNormedRows )


    def summarizeColumnsWithinGroups(self, cols, groupCols, groupsAreContiguous = True ):
        """Return an IDotData with the specified columns summarized within row groups defined by groupCols
        (each group of adjacent rows that have an identical combination of values in columns 'groupCols' is one row group).

        >>> z = IDotData( names = ( 'a', 'b' ), Records = ( ( 1, 2 ), ( 1, 3 ), ( 2, 4 ), ( 2, 5 ) ) )
        >>> z.summarizeColumnsWithinGroups( cols = 'b', groupCols = 'a' ).dbg()
        ... # doctest: +NORMALIZE_WHITESPACE
        ==============================
        a	b_count	b_sum	b_sumSq	b_numNaN
        1	2	5.0	13.0	0
        2	2	9.0	41.0	0
        <BLANKLINE>
        2 rows
        ==============================
        >>> z.summarizeColumnsWithinGroups( cols = 'b', groupCols = () ).dbg()
        ... # doctest: +NORMALIZE_WHITESPACE
        ==============================
        b_count	b_sum	b_sumSq	b_numNaN
        4	14.0	54.0	0
        <BLANKLINE>
        1 rows
        ==============================
        >>> z = IDotData( names = ( 'a', 'b' ), Records = ( ( 1, 2 ), ( 2, 4 ), ( 2, 5 ), ( 1, 3 ) ) )
        >>> z.summarizeColumnsWithinGroups( cols = 'b', groupCols = 'a', groupsAreContiguous = False ).dbg()
        ... # doctest: +NORMALIZE_WHITESPACE
        ==============================
        a	b_count	b_sum	b_sumSq	b_numNaN
        1	2	5.0	13.0	0
        2	2	9.0	41.0	0
        <BLANKLINE>
        2 rows
        ==============================

        """

        groupCols = flatten( MakeSeq( groupCols ) )
        if not groupsAreContiguous and groupCols:
            return self.sortedOn( *groupCols ).summarizeColumnsWithinGroups( cols = cols, groupCols = groupCols,
                                                                             groupsAreContiguous = True )

        cols = flatten( MakeSeq( cols ) )

        colNums = [ self.colName2num[ c ] for c in cols ]

        # Set up a place to keep the running sum and sum-of-squares for each
        # column to be normalized.  These are allocated once and then reused
        # within each row block.
        sums = [ None ] * len( cols )
        sumsSq = [ None ] * len( cols )
        for i in range( len( cols ) ):
            sums[ i ] = SumKeeper()
            sumsSq[ i ] = SumKeeper()
        counts = np.zeros( len( cols ), dtype = int )
        nanCounts = np.zeros( len( cols ), dtype = int )

        def GetSummarizedRows():

            #summarizedRowsHeadings = tuple( groupCols ) + tuple([ col + '_' + sfx for col in cols for sfx in ( 'count', 'sum', 'sumSq', 'numNaN' ) ])
            #dbg( 'summarizedRowsHeadings' )
            
            yield tuple( groupCols ) + tuple([ col + '_' + sfx for col in cols for sfx in ( 'count', 'sum', 'sumSq', 'numNaN' ) ])

            prevGroupId = None
            for groupId, rows1 in self.groupby( *groupCols ):

                thisGroupId = tuple( MakeSeq( groupId ) )
                assert prevGroupId is None or prevGroupId < thisGroupId

                for i in range( len( cols ) ):
                    sums[ i ].clear()
                    sumsSq[ i ].clear()
                counts.fill( 0 )

                for r in rows1:
                    for i, colNum in enumerate( colNums ):
                        x = r[ colNum ]
                        if math.isnan( x ) or math.isinf( x ):
                            nanCounts[ i ] += 1
                        else:
                            sums[ i ] += x
                            sumsSq[ i ] += x*x
                            counts[ i ] += 1

                yield thisGroupId + tuple( reduce( operator.concat, [ ( counts[ i ], sums[ i ].getSum(), sumsSq[ i ].getSum(), nanCounts[ i ] ) for i in range( len( cols ) ) ] ) ) 

        return IDotData.fromFn( GetSummarizedRows )
    
    @staticmethod
    def mergeColumnSummaries( iDotDatas, cols, groupCols ):
        """Merge column summaries produced by summarizeColumnsWithinGroups"""

        dbg( '"mergeColumnSummaries" cols groupCols' )
        
        def GetMergedRows():

            ourHeadings = tuple( groupCols ) + tuple([ col + '_' + sfx for col in cols for sfx in ( 'count', 'sum', 'sumSq', 'numNaN' ) ])
            dbg( 'ourHeadings' )
            yield tuple( groupCols ) + tuple([ col + '_' + sfx for col in cols for sfx in ( 'count', 'sum', 'sumSq', 'numNaN' ) ])

            sumCols = ( 'count', 'sum', 'sumSq', 'numNaN' )

            keyHeadings = tuple( 'grp_' + g for g in groupCols )
            dbg( 'keyHeadings' )
            for r in IDotData.mergeOnKeyCols( iDotDatas = iDotDatas, cols = ( groupCols, ) * len( iDotDatas ),
                                              blanks = ( ( np.nan, ) * len( groupCols ) +
                                                         ( 0, ) * len( cols ) * len( sumCols ), ) * len( iDotDatas ),
                                              keyHeadings = keyHeadings ):

                dbg( '"rrrrrr" r.headings' )
                
                yield tuple( list( r[ keyHeadings ] ) + [ sum([ r[ col + '_' + sfx + ( ( '_%d' % i ) if len( iDotDatas ) > 1 else '' ) ]
                                                                for i in range( len( iDotDatas ) ) ])
                                                          for col in cols for sfx in sumCols ])

        return IDotData.fromFn( GetMergedRows )


    def addMeanStdCols( self, cols ):
        """Add columns representing mean and std within each group.

        >>> z = IDotData( names = ( 'a', 'b' ), Records = ( ( 1, 2 ), ( 1, 3 ), ( 2, 4 ), ( 2, 5 ) ) )
        >>> z.summarizeColumnsWithinGroups( cols = 'b', groupCols = 'a' ).addMeanStdCols( cols = 'b' ).dbg()
        ... # doctest: +NORMALIZE_WHITESPACE
        ==============================
        a	b_count	b_sum	b_sumSq	b_numNaN	b_mean	b_std
        1	2	5.0	13.0	0	2.5	0.5
        2	2	9.0	41.0	0	4.5	0.5
        <BLANKLINE>
        2 rows
        ==============================
        
        """

        cols = MakeSeq( cols )

        def getMeanStd( a_sum, sumSq, count ):
            if count < 1: return ( np.nan, np.nan )
            avg = a_sum / count
            return ( avg, math.sqrt( np.max( ( 0, sumSq / count - avg * avg ) ) ) )

        return self.addComputedCols( newColNames = [ c + sfx for c in cols for sfx in ( '_mean', '_std' ) ],
                                     newColFn = lambda r: \
                                         reduce( operator.concat,
                                                 [ getMeanStd( r[ c + '_sum'], r[ c + '_sumSq' ], r[ c + '_count' ] ) for c in cols  ], () ) )


    def normalizeColumnsWithinGroups_using_means( self, cols, groupCols, means, groupsAreContiguous = True ):
        """Normalize the specified columns of a table within groups.

        Params:

           cols - the columns to be normalized
           groupCols - the columns that define the groups: rows that have the same combination of values
              in the group columns, are in the same group.
           means - the means and stds computed e.g. with addMeanStdCols   
           groupsAreContiguous - if True, rows belonging to the same group must be contiguous in the table;
              if False, no such assumption is made.

        >>> z = IDotData( names = ( 'a', 'b' ), Records = ( ( 1, 2 ), ( 1, 3 ), ( 2, 4 ), ( 2, 5 ) ) )
        >>> means = z.summarizeColumnsWithinGroups( cols = 'b', groupCols = 'a' ).addMeanStdCols( cols = 'b' )
        >>> z.normalizeColumnsWithinGroups_using_means( cols = 'b', groupCols = 'a', means = means ).dbg()
        ... # doctest: +NORMALIZE_WHITESPACE
        ==============================
        a	b
        1	-1.0
        1	1.0
        2	-1.0
        2	1.0
        <BLANKLINE>
        4 rows
        ==============================
        """

        if not groupsAreContiguous: grp2means = dict([ ( r[ groupCols ], r ) for r in means ])

        cols = flatten( MakeSeq( cols ) )
        groupCols = flatten( MakeSeq( groupCols ) )

        isNormeds = [ h in cols for h in self.headings ]

        meanCols = [ means.colName2num[ h + '_mean' ] if isNormed else -1  for h, isNormed in zip( self.headings, isNormeds ) ] 
        stdCols = [ means.colName2num[ h + '_std' ] if isNormed else -1  for h, isNormed in zip( self.headings, isNormeds ) ]

        def GetRows():
          yield self.headings

          if groupsAreContiguous: meansIter = iter( means )

          curGroup = None
          for r in self:
            if not groupsAreContiguous:
              newGroup = True
              curGroup = grp2means[ r[ groupCols ] ]
            else:
              newGroup = False
              while curGroup is None or curGroup[ groupCols ] < r[ groupCols ]:
                curGroup = next( meansIter )
                newGroup = True

            if newGroup:
                curMeans = [ curGroup[ meanCol ] if meanCol >= 0 else np.nan for meanCol in meanCols  ]
                curStds = [ curGroup[ stdCol ] if stdCol >= 0 else np.nan for stdCol in stdCols  ]

            yield tuple([ ( ( val - mean_val ) / std_val if std_val > 0 else 0.0 ) if isNormed else val
                          for val, isNormed, mean_val, std_val in zip( r, isNormeds, curMeans, curStds ) ])

        return IDotData.fromFn( GetRows )
            
    
    def nlargest(self, n, *cols ):
        """Return a new IDotData consisting of n largest rows of this IDotData, *sorted from largest to smallest*
        (i.e. NOT in original order), as given by comparison on columns cols.
        Note that this forces the immediate reading of all values of this IDotData as soon as the first row is requested;
        however, only O(n) memory will be used.
        """

        if not cols: cols = self.headings
        
        def GetResult():
            yield self.headings
            for r in heapq.nlargest( n, self.recordsIter(), key = operator.itemgetter( cols ) ): yield r
            
        return IDotData.fromFn( GetResult )

    def nsmallest(self, n, *cols ):
        """Return a new IDotData consisting of n smallest rows of this IDotData, *sorted from smallest to largest*
        (i.e. NOT in original order), as given by comparison on columns cols.
        Note that this forces the immediate reading of all values of this IDotData as soon as the first row is requested;
        however, only O(n) memory will be used.
        """

        if not cols: cols = self.headings
        
        def GetResult():
            yield self.headings
            for r in heapq.nsmallest( n, self.recordsIter(), key = operator.itemgetter( cols ) ): yield r
            
        return IDotData.fromFn( GetResult )

    def isPrimaryKey(self, columnName):
        """Test whether the given column is a primary key, i.e. that each row has a
        unique value in this column, different from the value in any other row."""
        seen = set()
        for x in self[ columnName ]:
            if x in seen: return False
            seen.add( x )
        return True

    def numRows(self):
        i = 0
        for r in self: i += 1
        return i

    def __str__(self): return str( type( self ) )

# end: class IDotDataRoot    

    
class TableReaderBase(collections.Iterator):
    
    """Abstract base class for reading tables, whether in .tsv or .data/ format."""

    def __init__(self, headings, fname = None, allFiles = (), 
                 valueFixer = None):
        """Create a TableReaderBase"""

        self.headings = tuple( headings )

        if not len( set( headings ) ) == len( headings ):
            print 'headings=', headings
            seen = set()
            for h in headings:
                if h in seen:
                    print 'duplicate heading: ', h
                seen.add( h )
            
        assert len( set( headings ) ) == len( headings )
        
        self.fname = fname
        self.allFiles = allFiles
        self.colName2num = dict( ( colName, colNum ) for colNum, colName in enumerate( self.headings ) )
        self.valueFixer = valueFixer
        self.lineNum = 0
        
    def next(self):
        try:
            lineVals = self._nextLine()

            self.lineNum += 1

            if self.lineNum % 200000  ==  0:
                logging.info( joinstr(' ', self.lineNum, ' of ', \
                                  self.fname, ': ', zip( self.headings, lineVals ) ) )

            if len( lineVals ) != len( self.headings ):
                logging.error( 'invalid table format: header has %d columns, line %d has %d items'
                               % ( len( self.headings ), self.lineNum, len( lineVals ) ) ) 
                print 'lineVals=', lineVals, ' self.headings=', self.headings, ' \nzip=', zip( self.headings, lineVals )
                
            assert len( lineVals ) == len( self.headings )

            if self.valueFixer: lineVals = map( self.valueFixer, lineVals )

            return IDotDataRecord( lineVals = lineVals,
                                   colName2num = self.colName2num,
                                   headings = self.headings, fname = self.fname )
        except StopIteration:
            self.close()
            raise


class IDotDataGetColumns(IDotDataRoot):

    """Given an IDotData, constructs a new IDotdata consisting of some subset of its columns, possibly reordered.
    """

    def __init__(self, parent, item):
        item = MakeSeq( item )
        badColumns = [ h for h in item if h not in parent.headings ]
        if badColumns: raise AttributeError( ','.join( badColumns ) ) 
        super(type(self), self).__init__( headings = item, parents = (parent,) )
        self.parent = parent
        itemCols = [ parent.colName2num[ c ] for c in item ]
        igetter = operator.itemgetter( itemCols )

        def itemGetter(rec):
            return IDotDataRecord( lineVals = igetter(rec),
                                   colName2num = self.colName2num,
                                   headings = self.headings,
                                   fname = None )
        self.getter = itemGetter
        self.item = item

    def _iter(self):
        return itertools.imap( self.getter, self.parent.recordsIter() )

    def __str__(self):
        return 'IDotDataGetColumns(' + ','.join( self.item ) + ')'

class IDotDataFilterBool(IDotDataRoot):

    """An IDotData that consists of selected rows of a parent IDotData.  The rows are specified by an index
    IDotData which is a one-column stream of Booleans.  IDotDataFilterBool represents an IDotData
    consisting of those rows of the parent IDotData for which the corresponding item in filter IDotData is True."""

    def __init__(self, parent, filter):
        assert len( filter.headings ) == 1
        super(type(self), self).__init__( headings = parent.headings, parents = (parent,) )
        self.parent = parent
        self.filter = filter

    def _iter(self):
        for p, f in itertools.izip( self.parent, self.filter.flatIter() ):
            assert isinstance( f, bool )
            if f: yield p


class IDotDataRemoveDups(IDotDataRoot):

    """Return an IDotData consisting of non-duplicate rows of this IDotData.  Rows are considered
    duplicate if they're adjacent and equal at column(s) keyCols.  Of each group of duplicate rows, the first
    is selected to be part of the new IDotData."""

    def __init__(self, parent, keyCols = None):
        super(type(self), self).__init__( headings = parent.headings, parents = (parent,) )
        self.parent = parent
        self.keyCols = keyCols if keyCols is not None else parent.headings

    def _iter(self):
        lastRec = None
        for r in self.parent.recordsIter():
            if lastRec is None or r[ self.keyCols ] != lastRec[ self.keyCols ]:
                yield r
                lastRec = r
            
class IDotDataTakeWhile(IDotDataRoot):

    def __init__(self, parent, pred):
        super(type(self), self).__init__( headings = parent.headings, parents = (parent,) )
        self.parent = parent
        self.pred = pred

    def _iter(self):
        return itertools.takewhile( self.pred, iter( self.parent ) )

    def __str__(self):
        try: s = inspect.getsource( self.pred )
        except IOError:
            s = joinstr( ',', inspect.getsourcemodule( self.pred ), inspect.getsourcefile( self.pred ), inspect.getsourcelines( self.pred ) )

        s = s.replace( '\n', '\\n' )
        return 'IDotDataTakeWhile(%s)' % s

class IDotDataRenameCols(IDotDataRoot):

    def __init__(self, parent, renamings):
        renamings = dict( renamings )
        assert set( renamings.keys() ) <= set( parent.headings ), \
            'keys are %s headings are %s' % ( renamings.keys(), parent.headings )
        
        super(type(self), self).__init__( headings =
                                          [ DictGet( renamings, h, h) for h in parent.headings ],
                                          parents = (parent,) )
        self.parent = parent

    def _iter(self): return self.parent.recordsIter()

class IDotDataOneSlice(IDotDataRoot):
    """One slice of an IDotData split by specified column(s)"""

    def __init__(self, parent, colNames, partNum, totalParts, numChunks = None ):

        super(type(self),self).__init__( headings = parent.headings, parents = ( parent, ) )

        if not colNames: colNames = ()
        self.parent = parent
        self.colNames = colNames


        if numChunks is None:
            numChunks = 0
            lastId = None
            for i, r in enumerate( parent ):
                thisId = r[ colNames ]
                if thisId != lastId:
                    numChunks += 1
                lastId = thisId

        self.fromChunk = int( ( partNum / totalParts ) * numChunks )
        self.toChunk = int( ( ( partNum+1 ) / totalParts ) * numChunks )

        dbg( 'numChunks self.fromChunk self.toChunk' )

    def _iter(self):
        lastId = None
        chunkNum = 0
        for r in self.parent:
            thisId = r[ self.colNames ]
            if lastId is not None and thisId != lastId: chunkNum += 1
            lastId = thisId
            if chunkNum >= self.toChunk: break
            if chunkNum >= self.fromChunk: yield r


class IDotDataOneRegion(IDotDataRoot):
    """One region of an IDotData split by specified column(s)"""

    def __init__(self, parent, keyCols, beg, end ):

        super(type(self),self).__init__( headings = parent.headings, parents = ( parent, ) )

        self.parent = parent
        self.keyCols = MakeSeq( keyCols or () )
        self.beg = None if beg is None else tuple( map( coerceVal, MakeSeq( beg ) ) )
        self.end = None if end is None else tuple( map( coerceVal, MakeSeq( end ) ) )

        dbg( 'self.beg self.end' )

    def _iter(self):
        it = iter( self.parent )
        if self.beg is not None: it = itertools.dropwhile( lambda r: tuple( r[ self.keyCols ] ) < self.beg, it )
        if self.end is not None: it = itertools.takewhile( lambda r: tuple( r[ self.keyCols ] ) < self.end, it )
        return it

            
class IDotData_tsv(IDotDataRoot):
    """An IDotData that reads a .tsv file"""

    def __init__(self, fname, **tableReadOpts ):
        self.fname = fname
        self.tableReadOpts = tableReadOpts
        exec ExtractOpts( tableReadOpts,
                          'headings', None,
                          'headingSep', None,
                          'commentPrefix', '##' )
        self.comments = []

        if headings is None:
            # read headings from file.
            with OpenForRead( fname ) as f:
                headings = chomp( f.readline() )
                while headings.startswith( commentPrefix ):
                    self.comments.append( headings )
                    headings = chomp( f.readline() )
                headings = SplitStr( headings, FirstVal( headingSep, sep ) )

        dbg( '"IDotData_tsv" headings' )
                
        super(type(self),self).__init__( headings = headings, parents = () )
                
    def _iter(self):
        return TSVReader( fname = self.fname, **self.tableReadOpts )

    def __str__(self):
        return "IDotData_tsv('%s')" % self.fname


class IDotData_dotData(IDotDataRoot):
    """An IDotData that reads a .data/ directory"""

    def __init__(self, fname, valueFixer = None):

        if not fname.endswith('/'): fname += '/'
        assert fname.endswith('.data/')
        if not os.path.exists( fname ): raise IOError( 'Directory not found: ' + fname )
        self.fname = fname
        self.valueFixer = valueFixer

        headerFN = fname + os.path.basename(fname[:-len('.data/')]) + '.header.txt'
        if not os.path.exists( headerFN ):
            headerFNs = glob.glob( fname + '*.header.txt' )
            assert len( headerFNs ) == 1
            headerFN = headerFNs[0]

        super(type(self),self).__init__( headings = tuple( SlurpFileLines( headerFN ) ) )

        tsvFiles = [ f for f in reduce( operator.concat,
                                        map( operator.itemgetter(2), os.walk( fname ) ) ) if f.endswith('.csv') ]
        colNames = [ os.path.splitext( os.path.splitext( tsvF )[0] )[0] for tsvF in tsvFiles ]

        assert len( set( colNames ) ) == len( colNames )
        assert sorted( colNames ) == sorted( self.headings )
        
        col2file = dict( zip( colNames, tsvFiles ) )
        self.allFiles = [ fname + col2file[ col ] for col in self.headings ]

    def _iter(self):
        return DotDataReader( fname = self.fname, headings = self.headings, valueFixer = self.valueFixer,
                              allFiles = self.allFiles )
    

    def __str__(self):
        return "IDotData_dotData('%s')" % self.fname
    
class TSVReader(TableReaderBase):

    """An iterator that yields records from a TSV file as named tuples."""

    def __init__(self, fname, sep = '\t', headingSep = None, headings = None,
                 valueFixer = None, skipFirstLines = 0, lineFixer = None,
                 commentPrefix = '#' ):

        f = OpenForRead( fname )
        
        if headings == None:
            headings = f.readline()
            while headings.startswith( commentPrefix ): headings = f.readline()
            headings = chomp( headings )
            headings = SplitStr( headings, sep if headingSep is None else headingSep )

        super(type(self),self).__init__(headings = headings, fname=fname, allFiles=(fname,), valueFixer = valueFixer )
        
        self.sep = sep if sep != 'whitespace' else None
        self.freader = iter( f )
        self.f = f
        self.skipFirstLines = skipFirstLines
        self.lineFixer = lineFixer

    def _nextLine(self):
        """Get next line"""
        rawLine = chomp( self.freader.next() )
        if not rawLine:

            # possible bug: if there is only one heading and the line is
            # a blank string, we'll stop here.

            # we could have a "blank value allowed, and here is what it should be"
            # option for each column.
            
            raise StopIteration

        # skip any comment lines.  normally these only appear at the beginning.
        while rawLine.startswith( commentPrefix ):
            rawLine = chomp( self.freader.next() )
            if not rawLine: raise StopIteration

        while self.skipFirstLines > 0:
            dbg( '"skipping" rawLine' )
            rawLine = chomp( self.freader.next() )
            if not rawLine: raise StopIteration
            self.skipFirstLines -= 1

        line = rawLine.split( self.sep )
        if not line: raise StopIteration
        if self.lineFixer: line = self.lineFixer( line )
        return line

    def close(self):
        """Close the file"""
        if self.f is not None:
            self.f.close()
            self.f = None
            self.freader = None
        
class DotDataReader(TableReaderBase):

    """Class for reading .tsv files"""

    def __init__(self, fname, headings, allFiles, valueFixer = None):

        super(type(self),self).__init__( headings = headings, fname=fname, allFiles=allFiles,
                                         valueFixer = valueFixer )
        
        self.fs = tuple( map( open, allFiles ) )
        self.freaders = map( iter, self.fs )

    def _nextLine(self):
        """Get next line"""
        result =  tuple( chomp( freader.next() ) for freader in self.freaders )
        if not result: raise StopIteration
        if len(result) != len( self.fs ):
            print 'result=', result
        assert len(result) == len( self.fs )
        return result

    def close(self):
        """Close all open files we had"""
        if self.fs is not None:
            for f in self.fs: f.close()

# var: iddSlice - map from IDotData file name to slicing params
iddSliceInfo = {}
iddRegionInfo = {}

def CanonFN( fname ): return fname[:-1] if fname.endswith( '.data/' ) else fname

def SetSliceInfo( fname, *args ):
    """Specify slice info for the given file"""
    iddSliceInfo[ CanonFN( fname ) ] = args

def SetRegionInfo( fname, keyCols, beg, end ):
    """Specify slice info for the given file"""
    iddRegionInfo[ CanonFN( fname ) ] = Dict( 'keyCols beg end' )
    
def IDotData( fname = None, sep = '\t', headingSep = None, fileType = None, headings = None, valueFixer = None, lineFixer = None,
              skipFirstLines = 0, names = None, commentPrefix = '#',
              Path = None, SVPath = None, SVValueFixer = None, SVLineFixer = None, SVDelimiter = None,
              SVSkipFirstLines = None, Header = True, Records = None, Columns = None, multiPass = True,
              ToLoad = None, lookupSlice = True, lookupRegion = True, useHeadings = False ):
    """Creates an IDotData from one of several possible sources: a .tsv file, a .data/ directory, a list of
    records, or a list of columns.
    """

    if IDotData.isA( fname ): return fname
    if hasattr( type( fname ), 'isDotData' ): return IDotData.fromDotData( fname )

    assert sum( map( bool, ( fname, Path, SVPath, Records is not None, Columns is not None ) ) ) == 1
    
    # handle legacy arguments for backwards compatibility with DotData
    if not fname and Path: fname = Path
    if not fname and SVPath: fname = SVPath

    if fname and lookupSlice and CanonFN( fname ) in iddSliceInfo:
        return IDotData( lookupSlice = False,
                         **Dict( 'fname sep fileType headings valueFixer lineFixer skipFirstLines names commentPrefix '
                                 'Path SVPath SVValueFixer SVLineFixer SVDelimiter SVSkipFirstLines Header Records '
                                 'Columns multiPass ToLoad' ) ).oneSlice( *( iddSliceInfo[ CanonFN( fname ) ] ) )

    if fname and lookupRegion and CanonFN( fname ) in iddRegionInfo:
        return IDotData( lookupRegion = False,
                         **Dict( 'fname sep fileType headings valueFixer lineFixer skipFirstLines names commentPrefix '
                                 'Path SVPath SVValueFixer SVLineFixer SVDelimiter SVSkipFirstLines Header Records '
                                 'Columns multiPass ToLoad' ) ).oneRegion( **( iddRegionInfo[ CanonFN( fname ) ] ) )
    
    if not headings and names: headings = names
    if isinstance( headings, types.StringTypes ) and ( ' ' in headings or '\t' in headings ):
        headings = tuple( headings.split( '\t' if '\t' in headings else None ) )
    if not skipFirstLines and SVSkipFirstLines: skipFirstLines = SVSkipFirstLines
    if SVDelimiter: sep = SVDelimiter
    if SVValueFixer: valueFixer = SVValueFixer
    if SVLineFixer: lineFixer = SVLineFixer

    assert ( fname and ( Header or headings ) ) or ( ( Records is not None or Columns ) and headings )

    assert not Columns or len( Columns ) == len( headings )

    if Records is not None: return IDotData.fromIterable( headings = headings, iterable = Records )

    if Columns:
        # make a one-column IDotData from each column, then hstack them.
        return IDotDataHStack( *[ IDotData.fromIterable( headings = (h,), iterable = c, multiPass = multiPass )
                                  for h, c in zip( headings, Columns ) ])

    logging.info( 'IDotData: creating iter over ' + fname )

    if fileType is None and IsFileType( fname, '.data', '.data/' ): fileType = '.tsv'
    if sep == ' ': fileType = '.tsv'

    if IsFileType( fname, '.vcf' ): commentPrefix = '##'
    
    if fname.endswith( '.data' ) or fname.endswith( '.data/' ):
        result = IDotData_dotData( fname = fname, valueFixer = valueFixer )
    else:
        result =  IDotData_tsv( **Dict( 'fname sep headingSep headings valueFixer skipFirstLines lineFixer commentPrefix' ) )

    if ToLoad:
        dbg( '"BBBBEF" result.headings ToLoad' )
        result = result[ ToLoad ]
        dbg( '"AAAAFT" result.headings' )
    return result

def makeIdxPrepender( idx ):
    """Creates a function that prepends a given index to its argument, returning a tuple."""
    return lambda( v ): ( idx, v )

class attrDbg(object):

    def __getattribute__(self, attr):
        print 'getting attr ', attr
        return attr
    
def makeKeyGetter( k ):
    """Creates a function that gets a key"""
    def myFunc( v ):
        return k( v[1] )
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

    Equivalent to:  sorted(itertools.chain(*iters))
    except that iterables need not be finite.

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

    iters = tuple( iters )

    if keys is not None:
        keys = [ ( operator.itemgetter( k ) if isinstance( k, (int,tuple,list) ) else
                   ( operator.attrgetter( k ) if isinstance( k, types.StringTypes ) else k ) ) for k in keys ]
    
    if ids:
        return itermerge( iters = [ itertools.imap( makeIdxPrepender( idx ), i ) for idx, i in
                                    ( enumerate( iters ) if ids == True else zip( ids, iters ) ) ],
                          keys = ( operator.itemgetter( 1 ), ) * len( iters ) if keys is None else
                          map( makeKeyGetter, keys ),
                          includeKeys = includeKeys )
    
    if keys is not None:
        result = itermerge( iters = [ itertools.imap( makeKeyApplier( k ), i )
                                      for k, i in zip( keys, iters ) ] )
        if not includeKeys: result = itertools.imap( operator.itemgetter( 1 ), result )
        return result
    
    def merge( i1, i2 ):
        next1 = iter( i1 ).next
        next2 = iter( i2 ).next
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
    if iters_cnt == 0:
        return iter( () )
    if iters_cnt == 1:
        return iter( iters[0] )
    if iters_cnt == 2:
        return merge( iters[0], iters[1] )
    bisect = int( iters_cnt / 2 )
    return merge( itermerge( iters = iters[:bisect] ), itermerge( iters = iters[bisect:] ) )

def TableIterInnerJoinAux( tableIters, cols, headings, blanks, headingLens, colName2num, keyHeadings = () ):

    cols = MakeSeq( cols )
    assert sum( headingLens ) == len( headings )
    blanks = [ blank if blank is None or IsSeq( blank ) else (blank,)*headingLen
               for blank, headingLen in zip(blanks, headingLens) ]

    prevKey = None
    for k, g in itertools.groupby( itermerge( iters = tableIters, keys = cols, ids = True,
                                              includeKeys = True ),
                                   key = operator.itemgetter( 0 ) ):


        # check that the keys are sorted in strictly increasing order -- important for correct operation of join 
        if not( prevKey==None or k > prevKey ):
            dbg( 'prevKey k tuple(g)' )
        assert prevKey==None or k > prevKey
        prevKey = k

        records = tuple( g )

        origins = [ r[1][0] for r in records ]
        if not is_sorted( origins, strict = True ):
            dbg( 'records origins' )
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
            rec = map( tuple, recordsList )
            assert map( len, rec ) == headingLens
            rec = reduce( operator.concat, rec )
            yield IDotDataRecord( lineVals = rec + ( k if keyHeadings else () ), colName2num = colName2num,
                                  fname = None, headings = headings + keyHeadings )


def mergeHeadings( iDotDatas, suffixes = None ):

    if suffixes == None: suffixes = [ '_%d' % i for i in range( len( iDotDatas ) ) ]
    assert len( suffixes ) == len( iDotDatas )

    sharedColNames = set( [] )
    allColNames = set( [] )
    for iDotData in iDotDatas:
            for heading in iDotData.headings:
                    ( sharedColNames if heading in allColNames else allColNames ).add( heading )

    return tuple([ n if n not in sharedColNames else n + sfx for iDotData, sfx in zip( iDotDatas, suffixes )
                   for n in iDotData.headings ])
    

def GetIDotDatas(iDotDatas):
    return tuple( IDotData(iDotData) if isinstance(iDotData, types.StringTypes) else
                  ( IDotData.fromDotData( iDotData ) if hasattr( type( iDotData ), 'isDotData' ) else iDotData )
                  for iDotData in iDotDatas  )

class IDotDataJoin(IDotDataRoot):

    def __init__(self, iDotDatas, cols, suffixes = None, blanks = None, keyHeadings = () ):


        iDotDatas = GetIDotDatas( iDotDatas )
        
        super(type(self),self).__init__( headings = mergeHeadings( iDotDatas, suffixes ) + keyHeadings, parents = iDotDatas )
        
        if blanks == None: blanks = (None,) * len( iDotDatas )

        assert all([ not hasattr(blank,'__len__') or  len( blank ) == len( iDotData.headings )
                     for iDotData, blank in zip( iDotDatas, blanks ) ])

        self.iDotDatas = iDotDatas
        self.cols = cols
        self.blanks = blanks
        self.headingLens = [ len( iDotData.headings ) for iDotData in iDotDatas ]
        self.keyHeadings = keyHeadings

    def _iter(self):
        return TableIterInnerJoinAux( tableIters = [ idd.recordsIter() for idd in self.iDotDatas ],
                                      cols = self.cols, headings = self.headings if not self.keyHeadings else self.headings[ :-len( self.keyHeadings ) ], blanks = self.blanks,
                                      headingLens = self.headingLens,
                                      colName2num = self.colName2num,
                                      keyHeadings = self.keyHeadings
                                      )


def IDotDataAsTuples( iDotDatas, cols, blanks = None ):

    return TableIterInnerJoinAuxAsTuples( tableIters = map( iter, iDotDatas ), cols = cols,
                                          blanks = blanks if blanks is not None else [ (None,)*len(idd.headings)
                                                                                       for idd in iDotDatas ],
                                          headingLens = [ len( idd.headings ) for idd in iDotDatas ] )
    
    
def TableIterInnerJoinAuxAsTuples( tableIters, cols, blanks, headingLens ):
    """Merge the outputs of several sorted iterators into one unified sorted iterator.
    The unified iterator yields _tuples_ with each tuple position corresponding to one
    input iterator.   Where some input iterator does not have a value for some key,
    that position is set to None.
    """

    blanks = [ blank if blank is None or IsSeq( blank ) else (blank,)*headingLen
               for blank, headingLen in zip(blanks, headingLens) ]

    numJoinsSkipped, numJoinsAllowed = 0, 0
    prevKey = None
    for k, g in itertools.groupby( itermerge( iters = tableIters, keys = cols, ids = True,
                                              includeKeys = True ),
                                   key = operator.itemgetter( 0 ) ):

        # check that the keys are sorted in strictly increasing order -- important for correct operation of join 
        if not( prevKey==None or k > prevKey ):
            print 'prevKey=', prevKey, ' key=', k, ' g is ', tuple( g )
        assert prevKey==None or k > prevKey
        prevKey = k

        records = tuple( g )

        origins = [ r[1][0] for r in records ]
        if not is_sorted( origins, strict = True ):
            print 'records are ', records
            print 'origins are ', origins
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
            yield recordsList
            numJoinsAllowed += 1
        else: numJoinsSkipped += 1

    dbg( 'cols numJoinsAllowed numJoinsSkipped' )


class IDotDataVStack(IDotDataRoot):

    def __init__(self, *iDotDatas, **kwargs):

        print 'befFlatten: ', iDotDatas
        iDotDatas = flatten( iDotDatas )

        print 'aftFlatten: ', iDotDatas
        assert len( iDotDatas ) > 0

        iDotDatas = map( IDotData, iDotDatas )
        if 'sourceCol' in kwargs:
            ids = DictGetNotNone( kwargs, 'sourceIds', map( operator.attrgetter( 'fname' ), iDotDatas ) )
            dbg( 'kwargs["sourceCol"] ids' )
            iDotDatas = [ iDotData.hstack( IDotData.repeat( kwargs[ 'sourceCol' ], idd_id ) )
                          for idd_id, iDotData in zip( ids, iDotDatas ) ]


        if 'sourceLabels' in kwargs:
            # hstack to each IDotData some identifying columns

            sourceLabels = IDotData( kwargs[ 'sourceLabels' ] )

            iDotDatas = [ iDotData.hstack( IDotData.repeat( headings = sourceLabels.headings,
                                                            values = sourceLabel ) )
                          for sourceLabel, iDotData in zip( sourceLabels, iDotDatas ) ]
            
        headings = DictGet( kwargs, 'headings', iDotDatas[ 0 ].headings )
            
        super(type(self), self).__init__( headings = headings )
        self.iDotDatas = iDotDatas
        self.comments = reduce( operator.concat, [ idd.comments if hasattr( idd, 'comments' ) else []
                                                   for idd in iDotDatas ] )

    def _iter(self):
        return itertools.chain( *self.iDotDatas )

class IDotDataVStackFromIterable(IDotDataRoot):

    def __init__(self, iDotDatas, headings = None):

        isInfinite = False
        if not headings:
            iDotDatas = iter( iDotDatas )
            firstIDotData = next( iDotDatas )
            headings = firstIDotData.headings
            isInfinite = firstIDotData.isInfinite
            iDotDatas = itertools.chain.from_iterable( ( (firstIDotData,), iDotDatas ) )
        
        super(type(self), self).__init__( headings = headings )
        self.iDotDatas = iDotDatas
        self.iDotDatasIterUsedUp = False
        self.isInfinite = isInfinite

    def _iter(self):
        assert not self.iDotDatasIterUsedUp
        if isinstance( self.iDotDatas, collections.Iterator ): self.iDotDatasIterUsedUp = True
        return itertools.chain.from_iterable( self.iDotDatas )

    
class IDotDataHStack(IDotDataRoot):

    def __init__(self, *iDotDatas ):

        iDotDatas = GetIDotDatas( iDotDatas )
        self.iDotDatas = tuple( iDotDatas )
        super(type(self), self).__init__( headings = mergeHeadings( iDotDatas ), parents = iDotDatas )

    def _iter(self):

        lineValsGetter = operator.attrgetter( 'lineVals' )
        def mergeRecs( *vals ):
            return  IDotDataRecord( lineVals = reduce( operator.concat, map( lineValsGetter, vals ) ),
                                    colName2num = self.colName2num,
                                    headings = self.headings )

        return itertools.imap( mergeRecs, *[ idd.recordsIter() for idd in self.iDotDatas ] )


class IDotDataFilter(IDotDataRoot):

    def __init__(self, parent, pred):
        super(type(self),self).__init__( headings = parent.headings, parents = (parent,) )
        self.parent = parent
        self.pred = pred

    def _iter(self):
        for r in self.parent:
            if self.pred( r ):
                yield r
        #return itertools.ifilter( self.pred, iter( self.parent ) )

    
class IDotDataMapRecords(IDotDataRoot):

    """Create an IDotData whose rows (records) are obtained by applying a specified function to values from
    specified IDotDatas."""

    def __init__(self, func, iDotDatas, headings):
        super(type(self),self).__init__( headings = headings, parents = iDotDatas )
        self.iDotDatas = iDotDatas
        self.func = func

    def _iter(self): return itertools.imap( self.func, *self.iDotDatas )

class IDotDataStarMap(IDotDataRoot):

    def __init__(self, func, iDotData, headings):
        super(type(self),self).__init__( headings = headings, parents = (iDotData,) )
        self.iDotData = iDotData
        self.func = func

    def _iter(self):
        return itertools.starmap( self.func, self.iDotData )
    

IDotData.mergeOnKeyCols = IDotDataJoin
IDotData.merge = IDotDataJoin
IDotData.mergeAsTuples = IDotDataAsTuples
IDotData.mapRecords = IDotDataMapRecords
IDotData.starmap = IDotDataStarMap
IDotData.mergeColumnSummaries = IDotDataRoot.mergeColumnSummaries

IDotData.TableIterInnerJoinAuxAsTuples = TableIterInnerJoinAuxAsTuples

IDotData.setSliceInfo = SetSliceInfo
IDotData.setRegionInfo = SetRegionInfo

class IDotDataRepeat(IDotDataRoot):

    """A one-column IDotData consisting of a specified number of repeats
    of a given constant value"""

    def __init__(self, heading = None, value = None, headings = None, values = None, times = None):
        assert heading is not None and value is not None and headings is None and values is None \
            or heading is None and value is None and headings is not None and values is not None
        super(type(self),self).__init__( headings = ( heading, ) if headings is None else headings )
        self.values = value if values is None else values
        self.times = times
        if times is None:
            self.repeater = itertools.repeat( self.values )
            self.isInfinite = True

    def _iter(self):
        return self.repeater if self.times is None else itertools.repeat( self.values, self.times )

IDotData.repeat = IDotDataRepeat
IDotData.ones = lambda n = None: IDotData.repeat( heading = 'v', value = 1, times = n )
IDotData.zeros = lambda n = None: IDotData.repeat( heading = 'v', value = 0, times = n )

IDotData.vstack = IDotDataVStack
IDotData.vstackFromIterable = IDotDataVStackFromIterable
IDotData.hstack = IDotDataHStack

class IDotDataWhere(IDotDataRoot):

    """Creates an IDotData which takes each row from one of two given IDotDatas, depending on the value of a dispatcher
    IDotData.

    Fields:

       which - a one-column IDotData whose values specify which values to take.
    
    """

    def __init__(self, which, ifTrue, ifFalse):

        assert len( which.headings ) == 1
        assert ifTrue.headings == ifFalse.headings
        
        super(type(self),self).__init__( headings = ifTrue.headings, parents = ( which, ifTrue, ifFalse ) )
        
        self.which = which
        self.ifTrue = ifTrue
        self.ifFalse = ifFalse

    def _iter(self):
        for whichVal, ifTrueVal, ifFalseVal in itertools.izip( self.which, self.ifTrue, self.ifFalse ):
            yield ifTrueVal if whichVal else ifFalseVal

IDotData.where = IDotDataWhere


class IDotDataChoose(IDotDataRoot):

    """Creates an IDotData which takes each row from one of the given IDotDatas, depending on the value of a dispatcher
    IDotData.

    Fields:

       which - a one-column IDotData whose values specify which values to take.
    
    """

    def __init__(self, which, choices):

        assert len( which.headings ) == 1
        assert all([ c.headings == choices[0].headings for c in choices ])
        
        super(type(self),self).__init__( headings = choices[0].headings, parents = ( which, ) + tuple( choices ) )
        
        self.which = which
        self.choices = choices

    def _iter(self):
        for r in itertools.izip( self.which, *self.choices ):
            yield r[ r[0] + 1 ]

IDotData.choose = IDotDataChoose

class IDotDataFromFn(IDotDataRoot):

    def __init__(self, fn, *args, **kwargs):
        super(type(self),self).__init__( headings = next( fn(*args, **kwargs) ) )
        self.fn = fn
        self.args = args
        self.kwargs = kwargs

    def _iter(self):
        it = self.fn( *self.args, **self.kwargs)
        next(it)  # skip the headings
        return it

IDotData.fromFn = IDotDataFromFn

class IDotDataFromIterable(IDotDataRoot):

    """An IDotData that takes its records from a specified iterable.   The iterable may be a one-pass iterator,
    in which case this IDotData may be iterated over at most once.
    """

    def __init__(self, headings, iterable, multiPass = False):
        if isinstance( headings, types.StringTypes ) and ( ' ' in headings or '\t' in headings ):
            headings = tuple( headings.split( '\t' if '\t' in headings else None ) )
        super(type(self),self).__init__( headings = headings )
        if multiPass and isinstance( iterable, collections.Iterator ): iterable = tuple( iterable )
        self.iterable = iterable
        self.iteratorUsedUp = False

    def _iter(self):
        assert not self.iteratorUsedUp
        iterHere = iter(self.iterable)
        if iterHere is self.iterable or isinstance( self.iterable, collections.Iterator ): self.iteratorUsedUp = True
        return iterHere

IDotData.fromIterable = IDotDataFromIterable

def IDotDataFromDotData( d ):
    """Create an IDotData from a DotData"""

    if not hasattr( type( d ), 'isDotData' ):
        if hasattr( d, 'shape' ) and len( d.shape ) == 1:
            return IDotData( names = ( 'val', ), Columns = ( tuple( d ), ) )
        assert False
    
    result = IDotData.fromIterable( headings = d.dtype.names,
                                    iterable = d )
    result.len = len( d )
    return result

IDotData.fromDotData = IDotDataFromDotData

class IDotDataWriterRoot(object):
    __metaclass__ = ABCMeta

    """Abstract base class for classes that help you incrementally write out an IDotData file record-by-record."""

    def __init__(self, headings):
        #assert headings

        if isinstance( headings, types.StringTypes ) and ( ' ' in headings or '\t' in headings ):
            headings = tuple( headings.split( '\t' if '\t' in headings else None ) )
        self.headings = headings

    @abstractmethod
    def writeRecord(self, *line): pass

    def writeRecords(self, lines):
        """Write all lines from an iterable"""
        for line in lines: self.writeRecord( line )

    @abstractmethod
    def close(self): pass

    def __enter__(self): return self
    def __exit__(self, exc_type, exc_value, traceback):

        if exc_type is not None:
            dbg( '"Error_while_writing_IDotData:" exc_type exc_value traceback' )
            tb.print_exception( exc_type, exc_value, traceback )
        
        self.close()
        return False   # don't suppress any exceptions
    

class IDotDataSVWriter(IDotDataWriterRoot):
    """Writes records to SV file"""

    def __init__(self, fname, headings, sep = '\t', comments = () ):

        logging.info( 'opening tsv for write: ' + str( fname ) )
        self.f = OpenForWrite( fname )
        super(type(self), self).__init__( headings = headings )
        for c in comments: self.f.write( c + '\n' )
        if self.headings:
            tabwrite( self.f, *self.headings, sep = sep )
            logging.info( 'headings are: ' + str( self.headings ) )
        self.sep = sep

    def writeRecord(self, *line):
        """Write one record to the file"""

        if isinstance( line, collections.Sized ) and len( line ) == 1 and IsSeq( line[0] ): line = line[0]
        if isinstance( line, types.GeneratorType ): line = tuple( line )
        line = MakeSeq( line )
        assert not self.headings or len( line ) == len( self.headings )
        
        tabwrite( self.f, *MakeSeq( line ), sep = self.sep )
        self.isFirstLine = False

    def close(self): self.f.close()

class IDotDataDotDataWriter(IDotDataWriterRoot):
    """Class for incrementally writing to a .data/ directory."""
    
    def __init__(self, fname, headings, comments = ()):

        if not fname.endswith('/'): fname += '/'
        assert fname.endswith('.data/')

        self.fname = fname
        super(type(self), self).__init__( headings = headings )

        MakeDir( fname )

        DumpFile( fname + os.path.basename( MapPath( fname )[:-len('.data/')] )
                  + '.header.txt', '\n'.join( self.headings ) )
        self.isFirstRecord = True

    type2name = { float: 'float', int: 'int', str: 'str', bool: 'bool', long: 'long' }
    if haveNumpy: type2name.update( ( ( np.int64, 'long' ), (np.bool_, 'bool' ), ( np.float64, 'float' ) ) )
    typeOrder = ( bool, ) + ( ( np.bool_, ) if haveNumpy else () ) + \
        ( int, long ) + ( ( np.int64, np.float64 ) if haveNumpy else () ) +  ( float, str )

    def __initColTypes(self, colTypes):

        self.colTypes = tuple( colTypes )
        self.newColTypes = list( self.colTypes )

        self.fs = [ open(self.fname + heading + '.' + type(self).type2name[colType] + '.csv', 'w')
                    for heading, colType in zip( self.headings, self.colTypes ) ]

    
    def writeRecord( self, *line ):
        """Append a record to a .data/"""

        if isinstance( line, collections.Sized ) and len( line ) == 1 and IsSeq( line[0] ): line = line[0]
        if isinstance( line, types.GeneratorType ): line = tuple( line )
        
        line = MakeSeq( line )

        assert len( line ) == len( self.headings )

        if self.isFirstRecord: self.__initColTypes( type(coerceVal(val)) for val in line )

        typeOrder = type(self).typeOrder
            
        for colNum, ( v, f ) in enumerate( zip( line, self.fs) ):
            if not self.isFirstRecord: f.write( '\n' )
            f.write( str(v) )

            newType = type( coerceVal( v ) )
            oldType = self.newColTypes[ colNum ]
            if newType not in typeOrder: dbg( 'newType' )
            if typeOrder.index( newType ) > typeOrder.index( oldType ): self.newColTypes[ colNum ] = newType

        self.isFirstRecord = False
                    
    def close(self):

        if self.isFirstRecord: self.__initColTypes( ( str, ) * len( self.headings ) )
        
        for f in self.fs: f.close()
        
        type2name = type(self).type2name
        for heading, colType, newColType in zip( self.headings, self.colTypes, self.newColTypes ):
            if type2name[ newColType ] != type2name[ colType ]:
                # we have to rename the file to a new correct name.
                # make sure it has been fully written and closed.
                time.sleep(20)

                oldName = self.fname + heading + '.' + type2name[colType] + '.csv'
                WaitForFileToAppear( oldName )
                newName = self.fname + heading + '.' + type2name[newColType] + '.csv'
                os.rename( oldName, newName )
                time.sleep(10)
                WaitForFileToAppear( newName )

        logging.info( 'saved IDotData to ' + self.fname )


@contextlib.contextmanager
def IDotDataOpenForWrite( fname, headings, sep = '\t', fileType = None, comments = () ):
    """Context manager for writing .data/ directories"""

    openFn = IDotDataDotDataWriter if IsFileType( fname, '.data', '.data/' ) else IDotDataSVWriter
    with openFn( fname = fname, headings = headings, comments = comments ) as f: yield f

IDotData.openForWrite = IDotDataOpenForWrite

def isIDotData(x):
    """Test if the given object is an IDotData.
    We could test if it inherits from IDotDataRoot, but we also want to
    allow unrelated objects to implement the iDotData interface.
    If their class has an isIDotData attribute, we assume the object
    is an IDotData.  If a class inherits from IDotDataRoot then it automatically
    inherits this attribute."""
    return hasattr( type( x ), 'isIDotData' )

IDotData.isA = isIDotData

IDotData.rootClass = IDotDataRoot

def imean(iterable):
    """Return the mean value of an iterable, or nan if
    there are no values."""
    sum = SumKeeper()
    count = 0
    for x in iterable:
        sum += x
        count += 1
    return sum.getSum() / count if count > 0 else np.nan

def imeanstd_old( iterable ):
    """Return the mean and stddev of an iterable,
    or (nan,nan) if there are no values.  """

    sum = SumKeeper()
    sumSq = SumKeeper()
    n = 0
    for x in iterable:
        if not ( math.isnan( x ) or math.isinf( x ) ):
            sum += x
            sumSq += x*x
            n += 1

    if n == 0: return np.nan, np.nan

    seqMean = sum.getSum() / n
    seqStd = math.sqrt( sumSq.getSum() / n - seqMean * seqMean )

    return seqMean, seqStd

def imeanstd( iterable ):
    """Return the mean and stddev of an iterable,
    or (nan,nan) if there are no values.  If iterable yields
    sequence values, then return (mean,std) for each column
    of the sequence."""

    sums = None
    sumSqs = None
    numCols = None
    ns = None
    for x in iterable:

        xSeq = MakeSeq( x )
        if numCols is None:
            numCols = len( xSeq )
            sums = [ SumKeeper() for i in range( numCols ) ]
            sumSqs = [ SumKeeper() for i in range( numCols ) ]
            ns = [ 0 ] * len( xSeq )
        else: assert numCols == len( xSeq )

        for dim in range( numCols ):

            xVal = xSeq[ dim ]
        
            if not ( math.isnan( xVal ) or math.isinf( xVal ) ):
                sums[ dim ] += xVal
                sumSqs[ dim ] += xVal*xVal
                ns[ dim ] += 1

    if numCols is None: return np.nan, np.nan

    seqMeans = [ sum.getSum() / n for sum, n in zip( sums, ns ) ]
    seqStds = [ math.sqrt( sumSq.getSum() / n -
                           seqMean * seqMean )
                for sumSq, seqMean, n in zip( sumSqs, seqMeans, ns ) ]
    seqMeansStds = zip( seqMeans, seqStds )

    return seqMeansStds if numCols > 1 else seqMeansStds[ 0 ]

def imeanstd_plusStats( iterable ):
    """Return the mean and stddev of an iterable, as well as
    counts of values (total and the non-nan values on which
    this is based),
    or (nan,nan) if there are no values.  If iterable yields
    sequence values, then return (mean,std) for each column
    of the sequence."""

    sums = None
    sumSqs = None
    numCols = None
    ns = None
    ntots = None
    for x in iterable:

        xSeq = MakeSeq( x )
        if numCols is None:
            numCols = len( xSeq )
            sums = [ SumKeeper() for i in range( numCols ) ]
            sumSqs = [ SumKeeper() for i in range( numCols ) ]
            ns = [ 0 ] * numCols
            ntots = [ 0 ] * numCols
        else: assert len( xSeq ) == numCols

        for dim in range( numCols ):

            xVal = xSeq[ dim ]

            ntots[ dim ] += 1
            if not ( math.isnan( xVal ) or math.isinf( xVal ) ):
                sums[ dim ] += xVal
                sumSqs[ dim ] += xVal*xVal
                ns[ dim ] += 1

    if numCols is None: return np.nan, np.nan, 0, 0

    seqMeans = [ sum.getSum() / n for sum, n in zip( sums, ns ) ]
    seqStds = [ math.sqrt( sumSq.getSum() / n -
                           seqMean * seqMean )
                for sumSq, seqMean, n in zip( sumSqs, seqMeans, ns ) ]
    seqMeansStds = zip( seqMeans, seqStds, ns, ntots )

    return seqMeansStds if numCols > 1 else seqMeansStds[ 0 ]

if haveNumpy:

    #
    #  Extend numpy routines to handle IDotData arguments.
    #

    def extendNumpy(f):

        """A decorator that replaces a function f with a new function which calls f
        if any arguments are IDotDatas, otherwise calls a numpy function of the same name.
        """
        
        orig_fn = eval( 'np.' + f.func_name )

        def new_f( *args, **kwargs ):
            return f( *args, **kwargs ) if any( IDotData.isA( a ) for a in args ) \
                else orig_fn( *args, **kwargs )

        new_f.__doc__ = orig_fn.__doc__

        setattr( np, f.func_name, new_f )
        
        return f
        
    numpy_isnan = np.isnan
    @extendNumpy    
    def isnan(arg): return arg.mapVals( numpy_isnan )

    numpy_isfinite = np.isfinite
    @extendNumpy    
    def isfinite(arg): return arg.mapVals( numpy_isfinite )

    numpy_isinf = np.isinf
    @extendNumpy    
    def isinf(arg): return arg.mapVals( numpy_isinf )

    numpy_isneginf = np.isneginf
    @extendNumpy    
    def isneginf(arg): return arg.mapVals( numpy_isneginf )

    numpy_isposinf = np.isposinf
    @extendNumpy    
    def isposinf(arg): return arg.mapVals( numpy_isposinf )

    numpy_isscalar = np.isscalar
    @extendNumpy    
    def isscalar(arg): return arg.mapVals( numpy_isscalar )

    numpy_invert = np.invert
    @extendNumpy
    def invert(arg): return arg.mapVals( numpy_invert )

    numpy_exp = np.exp
    @extendNumpy
    def exp(arg): return arg.mapVals( numpy_exp )

    numpy_log = np.log
    @extendNumpy
    def log(arg): return arg.mapVals( numpy_log )
    
    @extendNumpy
    def nanmax(a): return a.nanmax()

    @extendNumpy
    def mean(a): return imean( a )

    @extendNumpy
    def min(a): return __builtin__.min( a )

    @extendNumpy
    def max(a): return __builtin__.max( a )

    numpy_abs = np.abs
    @extendNumpy
    def abs(arg): return arg.mapVals( numpy_abs )

#    @extendNumpy
#    def asarray( a, dtype = None):
#        assert a.numCols() == 1
#        return np.fromiter( a, dtype = float if dtype is None else dtype )

 #   @extendNumpy
#    def asanyarray( a, dtype = None):
#        assert a.numCols() == 1
#        return np.fromiter( a, dtype = float if dtype is None else dtype )
    
    numpy_atleast_1d = np.atleast_1d

    def atleast_1d(*arys):
        if len( arys ) == 1:
            a = arys[0]

            if IDotData.isA( a ):
                return np.asanyarray( a )
            else:
                return numpy_atleast_1d( *arys  )
        else:
            return map( atleast_1d, arys )

    np.atleast_1d = atleast_1d

    @extendNumpy
    def where( cond, ifTrue, ifFalse ): return IDotData.where( cond, ifTrue, ifFalse )

    @extendNumpy
    def choose( a, choices ): return IDotData.choose( a, choices )

