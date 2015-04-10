'''
Classes and functions pertaining to the DotData class, a structured tabular data object. 
'''

from __future__ import with_statement
import numpy, types, copy, sys, itertools, operator, System.Utils
from Operations.DotDataFunctions.SaveDotDataAsSV import *
from Operations.DotDataFunctions.SaveDotData import *
from Operations.DotDataFunctions.SaveColumnsInDotData import *
from Operations.DotDataFunctions.DotDataListFromDirectory import *
from Operations.DotDataFunctions.DotDataListFromSV import *
import Operations.DotDataFunctions.datavstack as dv
from Operations.DotDataFunctions.datahstack import *
from Operations.DotDataFunctions.DictionaryOps import MergeDicts, RestrictDict
from System.Colors import GrayScale
from Operations.MiscUtil import dbg, TableIterInnerJoin, TblIterFromDotData, tabwriten

class DotData(numpy.core.records.recarray):
	"""A numpy recarray (a table with named columns where each column is of a uniform Python type),
	with added functionality and ability to define named groups of columns.

	Invariants:

	The names of all columns are distinct (unique) within one DotData.

	"""

	isDotData = True
		
	def __new__(subtype,Array = None, Records = None, Columns = None, shape = None, dtype = None, formats = None,
		    names = None, ToLoad = None,coloring = None, Path = '', PathList = '' , SVPath = '', 
		    SVDelimiter = None, SVDelimiterRegExp = None, SVLineBreak = None, SVHeader = True, SVHash = None,
		    SVLineFixer = None, SVValueFixer = None, Wrap = None, SVSkipFirstLines = 0,rowdata = None):
		"""Construct a new DotData from a copy of the data that's passed in.   
		The new DotData is a separate object from all parts of the original data.

		There are several ways of creating a DotData.  You need to specify
		the table data, and the column headings.

		Specifying the data:

		   - from a numpy.ndarray (Array arg)
		   - from list of records (Records arg), where each record is represented
		   	as a sequence of elements
		   - from list of columns (Columns arg), by passing each column as 
		   	list of uniform Python values.
		   - from a .data directory (Path arg);
		   - from an SV (separated-values) file (SVPath, SVDelimiter,
		   	SVDelimiterRegExp, SVLineBreak, SVHeader, SVHash, 
		   		SVSkipFirstLines, SVLineFixer, SVValueFixer args)

		Specifying the column names:

		  Column names can be inferred from the input data:
		     - if you're constructing from a numpy.ndarray, the column names 
		     	are taken from its dtype attribute
		     - if you're constructing from .data directory, the column names
		     	are taken from the names of .csv
		     files in that directory.
		     - if you're constructing from an SV file, the column names 
		     	can be taken from the headers in the file
		     (by default they're assumed to be there).

  		  Or you can specify the column names as arguments:   

		     - you can specify column names as the 'names' argument 
		     	(list of strings), or as the 'dtype' argument (numpy.dtype object).
		   
		  If you specify your own column names, and you're constructing from a numpy.recarray,
		  they will override any column names from the numpy.recarray.
		   
		You can specify a _subset_ of the columns to use, by specifying 
		the list of column names as the ToLoad argument.

		Specifying the colors (named column groups):

	   Colorings can be inferred from the input data:

		  If constructing from a .data directory, colorings will be automaticallly
		  inferred from the directory tree.

	   Colorings can be passed as argument:

		  In the 'coloring' argument, pass a dictionary each color names
		  (a string) to list of column names in that color.

	   If colorings are passed as argument, they override any colorings
	   inferred from the input data.
	   
	   The Wrap argument adds a color with name Wrap (thus Wrap is a string) 
	   listing all column names.
	   (When this DotData is saved to a .data directory, all columns will be
	   nested in an additional directory, Wrap.data.)

		You can also specify some rowdata:
		
		The rowdata argument must either be None (no rowdata) or a record
		array the same length as the data set.     This rowdata is meant 
		to represent "row-by-row" metadata.   These are columns that
		you don't really want to appear in displays of the data set, but which 
		are meant to "travel with" no matter what subsetting you do on the 
		"real data" columns -- so that the rowdata columns will still be there 
		without your explicitly having to remember them.   
	
		One use of .rowdata columns is for formatting and communicating "side" 
		information to other DotData methods. For instance, various specially 
		designated columns in the .rowdata information can be used to tell other 
		applications that use DotDatas how to interpret the rows in a way that
		would be tedious for the user to have to remember to supply. Two instances of this are:
			
		-- A '__color__' column present in a DotData's rowdata is specially 
		interpreted by the browser's DotData- Html representation. The '__color__' column
		is expected in each row to contain a web-safe hex triplet color specification, 
		e.g a string of the form
			'#XXXXXX' (see http://en.wikipedia.org/wiki/Web_colors).  
		-- The 'Aggregates' columns is used to disambiguate rows that 
		are aggregates of data in other sets of rows for the .aggregate_in method (see comments on that below). 
		
		Rowdata information can also be used to specify arbitrary higher-level 
		groups of rows, in analogy to how the "coloring" attribute 
		specifies groupings of columns.   This would work either by: 
			-- storing in a .rowdata column whose name specifies group name,
				a boolean in each row as to whether the row belongs to that group, or
			-- for a "type" of grouping consisting of several nonintersecting row 
				groups, a single column specifying by some string or integer code 
				which group the row belongs to.  (An example of this is the "Aggregates"
				column used by the .aggregate_in method, see below for info about this.) 

		Parameters:

		  subtype - the class object for the actual type of the newly created 
		  DotData object.  (will be either type(DotData) or the type of a subclass).


		   Arguments for recarray: shape, formats (see numpy docs)
		   
		"""

		assert not ( Path and SVPath )

		if isinstance( Array, types.StringTypes ):
			if Array.endswith( '.data' ) or Array.endswith( '.data/' ):
				Path = Array
			else:
				SVPath = Array
			Array = None

		if not Path and SVPath and ( SVPath.endswith('.data') or SVPath.endswith('.data/') ):
			Path, SVPath = SVPath, None
		if not SVPath and Path and Path.endswith('.tsv'):
			SVPath, Path = Path, None
			

		# Ensure that any iterable can be used as arguments: e.g. tuples, not just lists
		Records, Columns, names, ToLoad = map( lambda x: x if x == None or isinstance(x, types.ListType)
						       else list(x),
						       ( Records, Columns, names, ToLoad ) )

		###
		
		if not Array is None:
			if len(Array) > 0:
				DataObj = numpy.rec.fromrecords(Array, dtype = dtype, formats = formats, names = names)
			else:
				DataObj = numpy.rec.fromarrays([[]]*len(Array.dtype), dtype = dtype, formats = formats, names = names)
		elif not Records is None:
			assert not names or all( len( r ) == len( names ) for r in Records )
			DataObj = numpy.rec.fromrecords(Records, dtype = dtype, formats = formats, names = names)
		elif Columns != None:
			DataObj = numpy.rec.fromarrays(Columns, dtype = dtype, formats = formats, names = names)
		elif Path != '':
			System.Utils.chkExists(Path)
			[Columns,givennames,coloring,rowdata] = DotDataListFromDirectory(Path,ToLoad=ToLoad)
			if names == None:
				names = givennames
			DataObj = numpy.rec.fromarrays(Columns,names = names)
		elif PathList != '':
			[System.Utils.chkExists(x) for x in PathList]
			result = DotDataFromPathList(PathList)
			result.headings = result.dtype.names
			return result
		elif SVPath != '':
			System.Utils.chkExists( SVPath )
			[Columns,givennames] = DotDataListFromSV(SVPath, SVDelimiter, SVDelimiterRegExp, SVLineBreak, SVHeader, SVHash, ToLoad = ToLoad, SkipFirstLines = SVSkipFirstLines, LineFixer = SVLineFixer, ValueFixer = SVValueFixer, SVNames = names )
			if SVHeader and names == None:
				names = givennames
			DataObj = numpy.rec.fromarrays(Columns,names = names)
			DataObj.loadedFrom = SVPath
		else:
			DataObj = numpy.core.records.recarray.__new__(numpy.core.records.recarray,shape,dtype = dtype,formats = formats,names=names)

		self = numpy.core.records.recarray.__new__(subtype,DataObj.shape,dtype = DataObj.dtype,names = names)
		if coloring != None:
			coloringsInNames = list(set(coloring.keys()).intersection(set(self.dtype.names)))
			if len(coloringsInNames) == 0:
				self.coloring = coloring				
			else:
				print "Warning: the following coloring keys,", coloringsInNames, ", are also attribute (column) names in the DotData. This is not allowed, and so these coloring keys will be deleted. The corresponding columns of data will not be lost and will retain the same names."
				for c in coloringsInNames:
					coloring.pop(c)
				self.coloring = coloring
		else:
			self.coloring = {}		
			
		self.rowdata = rowdata

		if Wrap != None:
			self.coloring[Wrap] = self.dtype.names
	
		if self.dtype.names != None:
			for attrib in self.dtype.names:
				self[attrib] = DataObj[attrib]
			self.headings = self.dtype.names
			return self
		else:
			return DataObj
		
	def extract(self):
		"""Creates a copy of this DotData in the form of a numpy.ndarray.

		Useful if you want to do math on array elements, e.g. if you have a subset of the columns
		that are all numerical, you can construct a numerical matrix and do matrix operations.
		
		"""
		return numpy.vstack([self[r] for r in self.dtype.names]).T.squeeze()
	
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

		if isinstance(expr, types.CodeType) or isinstance(expr, types.StringType):
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
			result = self[ numpy.array( [ eval( exprCode, callerGlobals,
							    r if not callerLocals else MergeDicts( callerLocals, r ) )
						      for r in self.recordsAsDicts() ] ) ]
			del callerGlobals, callerLocals
			return result
		
		elif callable(expr):
			return self[ numpy.array( [ expr(r) for r in self.recordsAsDicts() ]  ) ]
		elif isinstance(expr, dict):
			exprItems = frozenset( [ (columnName, val) for columnName, val in expr.items() if val != None ] )
			return self[ numpy.array( [ exprItems <= frozenset( r.items() ) for r in self.recordsAsDicts() ]  ) ]
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

	def __getitem__(self,attribute):
		"""Returns a subrectangle of the table.

		The representation of the subrectangle depends on type(attribute). Also,
		whether the returned object represents a new independent copy of the
		subrectangle, or a "view" into this self object, depends on type(attribute).

		- if you pass the name of an existing coloring, you get a DotData consisting of
		copies of columns in that coloring

		- if you pass a list of existing coloring names and/or column names, you get a
		DotData consisting of copies of columns in the list (name of coloring is
		equivalent to list of names of columns in that coloring; duplicate columns are
		deleted).

		- if you pass a Python expression as a code object (see compile_expr()), where
		columns names are used in the expression as variables, you get subset of rows
		for which the expression is true; more specifically, you get a DotData consisting
		of copies of rows for which the expression evaluates to True, in an environment
		where column names are bound to the corresponding values for each row.

		- if you pass a dictionary mapping column names to values, you get back rows
		where the specified columns have the specified values (more specifically, you
		get a new DotData containing copies of these rows). Dictionary entries mapping
		to None are ignored.

		- if you pass an ndarray, you get a DotData consisting a subrectangle of the
		table, as handled by recarray's: * if you pass a 1D ndarray of booleans of
		len(DotData), the rectangle contains copies of the rows for which the
		corresponding entry is True. * if you pass a list of row numbers, you get a
		DotData containing copies of these rows.


		See also: selectRows(), which lets you specify several row-selection criteria at
		once, and lets you pass an expression for selecting the rows as a string rather
		than only as a Python code object.

		"""

		if 'coloring' in dir(self) and attribute in self.coloring.keys():
			restrictedcoloring = dict([(a,self.coloring[a]) for a in self.coloring.keys() if set(self.coloring[a]) < set(self.coloring[attribute])])
			return DotData(Columns = [numpy.core.records.recarray.__getitem__(self.view(),attrib) for attrib in self.coloring[attribute]],dtype = numpy.dtype([(a,self.dtype[a].str) for a in self.coloring[attribute]]),coloring = restrictedcoloring,rowdata = self.rowdata)
		elif isinstance(attribute,(list,tuple)) and all( isinstance( x, types.StringType ) for x in attribute ):
			if not ( set(attribute) <= set(self.coloring.keys()).union(self.dtype.names) ):
				print "Unknown column name(s):", set(attribute) - set(self.coloring.keys()).union(self.dtype.names)
			assert set(attribute) <= set(self.coloring.keys()).union(self.dtype.names)
			attribset = []
			for att in attribute:
				if att in self.dtype.names and att not in attribset:
					attribset += [att]
				elif att in self.coloring.keys():
					attribset += [j for j in self.coloring[att] if j not in attribset]
			#return DotData(Columns = [numpy.core.records.recarray.__getitem__(self.view(),a) for a in attribset], names = attribset, coloring = dict([(a,self.coloring[a]) for a in attribute if a in self.coloring.keys()]))
			return DotData(Columns = [numpy.core.records.recarray.__getitem__(self.view(),a) for a in attribset], names = attribset, coloring = dict([(a,list(set(self.coloring[a]).intersection(set(attribute)))) for a in self.coloring.keys() if len(set(self.coloring[a]).intersection(set(attribute))) > 0]),rowdata=self.rowdata)
		elif isinstance(attribute,numpy.ndarray) or (isinstance(attribute,list)):
			if self.rowdata != None:
				rowdata = self.rowdata[attribute]
			else:
				rowdata = None
			if isinstance(self[0],numpy.record):
				return DotData(Array = numpy.core.records.recarray.__getitem__(self,attribute),dtype = self.dtype, coloring = self.coloring,rowdata=rowdata)
			else:
				D = numpy.core.records.recarray.__getitem__(self,attribute)
				D.coloring = self.coloring
				D.rowdata = rowdata
				return D
		elif isinstance(attribute,types.CodeType) or isinstance(attribute,types.FunctionType) or isinstance(attribute,dict):
			return self.selectRows( attribute )
		else:
			obj =  numpy.core.records.recarray.__getitem__(self.view(),attribute)
			if isinstance(obj,numpy.record) and len(obj) == 1:
				return obj[0]
			return obj
			
	def __getslice__(self,attribute1,attribute2):
		"""Returns a slice into the array: a _contiguous_ range of its rows.  
		This is not a copy but a mutable reference into the array.
		"""
		obj = numpy.core.records.recarray.__getslice__(self,attribute1,attribute2)
		try:	
			obj.coloring = self.coloring
			if self.rowdata != None:
				obj.rowdata = self.rowdata[attribute1:attribute2]
			else:
				obj.rowdata = None
		except:
			pass
		return obj
	

	def renameColumns(self, newnames):
		""" Create a new DotData with column names given by new names.
		"""
		return DotData(Columns = [self[n] for n in self.dtype.names], names = newnames, coloring = self.coloring, rowdata = self.rowdata)
	
	
	def vstack(self,new):
		""" Create a new DotData that has rows of self followed by rows of new,
		with the headers, colorings and rowdata correctly merged.
		"""
		return dv.datavstack([self,new])
	
	def hstack(self,new):
		""" Create a new DotData that has columns of self followed by columns of
		new, with the headers and colorings correctly merged.  rowdata is merged by
		simple SafeColumnStack (if there's repeated column names from different
		arrays, no exception is raised and the data from the first array that has that
		rowdata is taken)
		"""
		return datahstack([self,new])
		
	def addrecords(self,new):
		"""Append one or more records to the end of the array. Can take a single
		record or tuple, or a list of records/tuples.
		"""
		if isinstance(new,numpy.record) or isinstance(new,tuple):
			new = [new]
		
		if self.rowdata != None:
			newrowdata = numpy.rec.fromarrays([[nullvalue(self.rowdata[l][0])]*len(new) for l in self.rowdata.dtype.names], names = self.rowdata.dtype.names)
			rowdata = SimpleStack1([self.rowdata,newrowdata])
		else:
			rowdata = None
		
		return DotData(Array = numpy.append(self,numpy.rec.fromrecords(new,dtype = self.dtype),axis = 0),dtype = self.dtype,coloring = self.coloring,rowdata =rowdata)

	def copy(self):
		x = numpy.core.records.recarray.copy(self)
		x.coloring = self.coloring
		x.rowdata = self.rowdata
		x.headings = self.headings
		return x
	
	def recordsAsDicts(self):
		"""Return an interator over the records, which yields each record as a 
		dictionary mapping column name to the value of that column in that record.

		"""
		for r in self:
			yield dict( [ (n,v) for n, v in zip( self.dtype.names, r ) ] )

	def view(self,obj=numpy.core.records.recarray):
		return numpy.core.records.recarray.view(self,obj)
				
	def save(self,TargetFolderName,HeaderOn = 1):
		"""Save the DotData to a folder, preserving the colorings.  Can later 
		be loaded back by passing the folder as the Path argument to DotData 
		constructor.
		"""
		if TargetFolderName.endswith( '.tsv' ): self.saveToSV( TargetFolderName )
		else: SaveDotData(self,TargetFolderName,HeaderOn)

	def ColIdx(self, colName):
		"""Return the index of the column with the given name"""
		return list( self.dtype.names ).index( colName )


	@staticmethod
	def mergeOnKeyCols( dotDatas, primaryKeyCols, blanks = None, suffixes = None, innerJoin = False, verbose = False  ):
		"""Do a relational join"""

		if blanks == None and not innerJoin: blanks = [ (None,) * d.numCols() for d in dotDatas ]
		if suffixes == None: suffixes = [ '_%d' % i for i in range( len( dotDatas ) ) ]
		if innerJoin: blanks = (None,) * len(dotDatas)
            
		assert len( dotDatas ) == len( primaryKeyCols ) == len( blanks ) == len( suffixes )
		assert len( suffixes ) == len( set( suffixes ) )   # ensure that all suffixes differ.
		#assert all( len( blank ) == dotData.numCols() for dotData, blank in zip( dotDatas, blanks ) )


		return DotDataFromTblIter( TableIterInnerJoin( tableIters = map( TblIterFromDotData, dotDatas ),
							       cols = primaryKeyCols,
							       suffixes = suffixes,
							       blanks = blanks ) )



	@staticmethod
	def mergeOnKeyColsOld( dotDatas, primaryKeyCols, blanks = None, suffixes = None, innerJoin = False, verbose = False  ):
		"""Merge several dotDatas based on primary key columns.

		Parameters:

		dotDatas - sequence of DotDatas to join primaryKeyCols - sequence of column
		names, one for each DotData in dotDatas. these columns should all have the
		same type of data, typically an identifier telling us for what entity (e.g. for
		what SNP) that row gives information.

		If each item in the sequence is a tuple of column names rather than a single
		column name, then the primary key from the corresponding DotData is the
		tuple of values from these columns, rather than a single value from one
		column.

		blanks - when a DotData does not have a value for a given primary key, this
		tuple is filled in. suffixes - sequence of suffixes, one for each DotData, to
		append to names of columns from that DotData to make them unique, if
		needed.

		For each primary key value existing in at least one DotData, we create one
		output record in the resulting merged DotData; these records are in order of
		primary key.   The columns in the output are obtained by concatenating
		column name lists of the input DotDatas; any column names shared between
		two or more dotDatas are made unique by appending the corresponding suffix
		from 'suffixes'.   The record corresponding to a given primary key is obtained
		by concatenating the records of the input dotDatas; when some input dotData
		does not have a row with the given primary key, the corresponding tuple from
		'blanks' is used instead.

		"""

		if blanks == None: blanks = [ (None,) * d.numCols() for d in dotDatas ]
		if suffixes == None: suffixes = [ '_%d' % i for i in range( len( dotDatas ) ) ]
		
		assert len( dotDatas ) == len( primaryKeyCols ) == len( blanks ) == len( suffixes )
		assert len( suffixes ) == len( set( suffixes ) )   # ensure that all suffixes differ.
		assert all( len( blank ) == dotData.numCols() for dotData, blank in zip( dotDatas, blanks ) )

		

		# construct the merged list of column names, by concatenating the column name lists of input dotDatas.
		# if any column name appears in two or more dotDatas, make these columns unique in the output by appending
		# the corresponding suffix.
		sharedColNames = set( [] )
		allColNames = set( [] )
		for d in dotDatas:
			for n in d.dtype.names:
				( sharedColNames if n in allColNames else allColNames ).add( n )

		newNames = [ n if n not in sharedColNames else n + sfx for d, sfx in zip( dotDatas, suffixes ) \
				     for n in d.dtype.names ]

		#
		# Construct the records (rows) of the merged DotData: one record for each distinct primary key value
		# appearing in at least one of the input dotDatas.
		#

		newRecords = [ ]

		# Create a flat list that has a triple (primaryKeyVal, dotDataId, row) for each occurrence of a primary key value in a dotData row.
		vals = sorted( [ (primaryKeyVal if isinstance( primaryKeyCol, types.StringType ) or \
					  len( primaryKeyCol ) == 1 else tuple( primaryKeyVal ), dotDataId, row ) \
				  for dotDataId, ( dotData, primaryKeyCol ) in enumerate( zip( dotDatas, primaryKeyCols ) ) \
					 for row, primaryKeyVal in enumerate( dotData[ primaryKeyCol ] ) ] )
		iteration = 0
		
		# Group the triples by primary key value.
		recordType = None
		for primaryKeyVal, primaryKeyOccs in itertools.groupby( vals, operator.itemgetter(0) ):

			if verbose:
				if ( iteration % 100 ) == 0:
					iteration = 0
					print "merging: primary key ", primaryKeyVal
				iteration += 1
					
			# For the primary key value 'primaryKeyVal', make a map from dotDataId to the row in that dotData where this primary key value occurs.
			# For those dotDatas where primaryKeyVal does not occur in any row, there is no entry in the map.
			dotDataId2row = dict( [ ( dotDataId, row ) for pkval, dotDataId, row in primaryKeyOccs ] )
			# Map each dotDataId to either the record from the corresponding input dotData (if that dotData has a row with primaryKeyVal)
			# or to the corresponding blank record (if that dotData does _not_ have a row with primaryKeyVal).
			# Concatenate these records to obtain the output record for primaryKeyVal.
			if innerJoin and len( dotDataId2row ) < len( dotDatas ): continue
			newRecord = reduce( operator.concat, \
						   [ tuple( dotDatas[ dotDataId ][ dotDataId2row[ dotDataId ] ] \
								    if dotDataId2row.has_key( dotDataId ) \
								    else blanks[ dotDataId ] ) \
						     for dotDataId in range( len( dotDatas ) ) ] )
			newRecordType = map( type, newRecord )
			if recordType == None: recordType = newRecordType
			else: assert newRecordType == recordType
			newRecords.append( newRecord )

		return DotData( Records = newRecords, names = newNames )


	@staticmethod
	def fromIDotData( idd ):
		"""Create a DotData from an IDotData"""
		return DotData( names = idd.headings, Records = idd.recordsIter() )
	
	def numCols(self):
		"""Return the number of columns in this DotData"""
		return len( self.dtype.names )

	def numRows(self):
		"""Return the number of rows in this DotData"""
		return len( self )

	def columns(self, cols = None):
		"""Return the columns, as a list.  Convenient for mapping a function over the columns."""
		if cols == None: cols = self.dtype.names
		return [ self[ colName ] for colName in cols ]
	
	def saveToSV(self,TargetFileName,sep = None,linesep='\n', UseHeader = True):
		"""Save the DotData to a single flat file.  Column headers are kept, but 
		colorings are lost.
		"""
		SaveDotDataAsSV(self,TargetFileName,sep=sep,linesep=linesep,UseHeader=UseHeader)
		
	def saveColumns(self,TargetFolderName):
		"""Save the DotData to a set of flat .csv files in .data format (i.e. .int.csv, 
		.float.csv, .str.csv). Colorings are lost.
		"""
		SaveColumnsInDotData(self, TargetFolderName)

	def withoutNans( self, cols = None ):
		"""Return a new DotData obtained from this one by removing rows that have nan values
		in any of the given columns"""

		if cols == None: cols = self.dtype.names

		dbg( '#zip( cols, map( type, [ numpy.isnan( self[ c ] ) for c in cols ] ) )' )
		
		return self[ numpy.invert( reduce( operator.or_, [ numpy.isnan( self[ c ] ) for c in cols ] ) ) ]

	def replace(self,old,new,strict=True,cols=None,rows = None):
		"""
			Replace value "old" with "new" everywhere it appears in self. 
			If strict = True, replace only exact occurences of "old"; if strict = False, assume 'old' and 'new' are strings and replace all occurences of substrings (e.g. like str.replace)
			cols = columns to make replacements in, if none, make replacements everywhere
			rows = rows to make replacements in, if none, make replacements everywhere
			
			NOTE: this function does in-place replacements.  thus there are issues handling datatypes here when replacement dtype is larger than original dtype.  Can be resolved later by making new array when necessary ... 
		"""
		
		if cols == None:
			cols = self.dtype.names
			
		if rows == None:
			rows = numpy.ones((len(self),), bool)
			
		if strict:
			new = numpy.array(new)
			for a in cols:
				if self.dtype[a] < new.dtype:
					print 'WARNING: dtype of column', a, 'is inferior to dtype of ', new, 'which may cause problems.'
				self[a][(self[a] == old)[rows]] = new
		else:
			for a in cols:
				QuickRep = True
				if rows != None:
					colstr = ''.join(self[a][rows])
				else:
					colstr = ''.join(self[a])		
				avoid = [ord(o) for o in uniqify(old + new + colstr)]
				ok = set(range(256)).difference(avoid)
				if len(ok) > 0:
					sep = chr(list(ok)[0])
				else:
					ok = set(range(65536)).difference(avoid)
					if len(ok) > 0:
						sep = unichr(list(ok)[0])
					else:
						print 'All unicode characters represented in column', a ,', can\t replace quickly.'
						QuickRep = False
						
				if QuickRep:
					newrows = numpy.array(sep.join(self[a][rows]).replace(old,new).split(sep))			
				else:
					newrows = numpy.array([aa.replace(old,new) for aa in self[a][rows]])
				self[a][rows] = numpy.cast[self.dtype[a]](newrows)
					
				if newrows.dtype > self.dtype[a]:
						print 'WARNING: dtype of column', a, 'is inferior to dtype of its replacement which may cause problems.'

	def isPrimaryKey(self, columnName):
		"""Test whether the given column is a primary key, i.e. that each row has a
		unique value in this column, different from the value in any other row."""
		return len( frozenset( self[ columnName ] ) ) == self.numRows()

	def dbg(self, msg = None):
		print '=============================='
		if msg is not None:
		    print '\n' + msg + '\n'
		    print '=============================='
		tabwriten( sys.stdout, self.headings )
		nlines = 0
		for line in self:
		    sys.stdout.write( '\n' + '\t'.join( map( str, line )
							if ( hasattr(line,'__iter__') and not isinstance(line,types.StringTypes) )
							else (str(line),) ) )
		    nlines += 1
		sys.stdout.write( '\n%d rows\n' % nlines )
		print '=============================='

	def aggregate(self,On=None,AggFuncDict = None,AggFunc=None):
		"""
		Aggregate a dataset, on specifed columns, using specified aggregation functions.  See commends for DotDataAggregate function. 
		"""
		
	
		return DotDataAggregate(self,On,AggFuncDict=AggFuncDict,AggFunc=AggFunc)
		
	def pivot(self,a,b,Keep=None):
		"""
			Pivot a DotData table  on columns a,b.  
			See comments for Pivot function, and http://en.wikipedia.org/wiki/Pivot_table. 
		"""
		
		return Pivot(self,a,b,Keep)

	def aggregate_in(self,On=None,AggFuncDict = None,AggFunc = None,interspersed = True):
		"""
		Take aggregate of data set on specified columns, then add the 
		resulting rows back into data set to make a composite object 
		containing both original non-aggregate data rows as well as the 
		aggregate rows. 
		
		First read commends for aggregate method.  Now:
		ARGUMENTS:
		--On, AggFuncDict -- same as arguments for aggregate method. 
		-- interspersed : boolean, if true aggregate rows are interleaved 
		with the data of which they are aggregates, if false, all aggregate 
		rows placed at the end of the array. 
		
		RETURNS:
		A DotData, with number of rows equaling:
				len(self) + len(A)
		where A is the the result of self.aggregate(On,AggFuncDict).
		A represents the aggregate rows and the other 
		rows were the original data rows. 
				
		This function supports _multiple_ aggregation, meaning that one
		can first aggregate on one set of factors, then repeat
		aggregation on the result for another set of factors, without the
		results of the first aggregation interfering the second.  To achieve 
		this, the method adds (or augments, if already present), some
		.rowdata information.  (See comments on __new__ of DotData 
		for more about rowdata information in general).     
		The specific rowdata information added by the aggregate_in 
		method is a column called "Aggregates" specifying on which 
		factors the rows that are aggregate rows were aggregated.  
		Rows added by aggregating on factor 'A' (a column in the original 
		data set) will have 'A' in the 'Aggregates' column of the self.rowdata array.  
		When multiple factors 'A1', 'A2' , ... are aggregated on, the notation
		is a comma-separated list 'A1,A2,...'.    This way, when you call
		.aggregate_in again, the function only aggregates on the columns 
		that have the empty char '' in their "Aggregates" rowdata.   
		
		The function also adds (or adds rows to) a '__color__' column 
		to self.rowdata, specifying Gray-Scale colors for aggregated
		rows that will be used by the Data Environment system's 
		browser for colorizing the  data.   When there are multiple
		levels of aggregation, the coarser aggregate groups 
		(e.g. on fewer factors) get darker gray color then those on
		finer aggregate groups (e.g. more factors).  
		
		"""
		#see if there's an Aggregates column in rowdata.  If so, strip off all those that atre non trivial 
		
		if self.rowdata != None and 'Aggregates' in self.rowdata.dtype.names:
				X = self[self.rowdata['Aggregates'] == ''][:]
				OldAggregates = self[self.rowdata['Aggregates'] != ''][:]
				AggVars = uniqify(ListUnion([x.split(',') for x in OldAggregates.rowdata['Aggregates']]))
		else:
			X = self
			OldAggregates = self[0:0]
			AggVars = []

		if On == None:
			On = []
			
		NewAggregates = DotDataAggregate(X,On,AggFuncDict=AggFuncDict,AggFunc=AggFunc)
		on = ','.join(On)
		NewAggregates.rowdata = numpy.rec.fromarrays([[on]*len(NewAggregates)],names = ['Aggregates'])
		AggVars = uniqify(AggVars + On)
		
		Aggregates = dv.datavstack([OldAggregates,NewAggregates])
		
		ANLen = numpy.array([len(x.split(',')) for x in Aggregates.rowdata['Aggregates']])
		U = numpy.array(uniqify(ANLen)); U.sort()
		[A,B] = System.Utils.fastequalspairs(ANLen,U)
		Grays = numpy.array(GraySpec(len(U)))
		AggColor = numpy.rec.fromarrays([Grays[A]],names = ['__color__'])
		
		Aggregates.rowdata = AddOrReplaceColumns(Aggregates.rowdata,AggColor)
		
		#s = len(ANLen) - ANLen.argsort() - 1
		s = ANLen.argsort()
		Aggregates = Aggregates[s[range(len(Aggregates)-1,-1,-1)]]
		
		if not interspersed or len(AggVars) == 0:
			return dv.datavstack([X,Aggregates])
		else:
			X.sort(order = AggVars)
			Diffs = numpy.append(numpy.append([0],1 + (X[AggVars][1:] != X[AggVars][:-1]).nonzero()[0]),[len(X)])		
			DiffAtts = ([[t for t in AggVars if X[t][Diffs[i]] != X[t][Diffs[i+1]]] for i in range(len(Diffs) - 2)] if len(Diffs) > 2 else []) + [AggVars]
			
			HH = {}
			for l in uniqify(Aggregates.rowdata['Aggregates']):
				Avars = l.split(',')
				print Avars
				HH[l] = System.Utils.FastRecarrayEqualsPairs(X[Avars][Diffs[:-1]],Aggregates[Avars])
			
			Order = []
			for i in range(len(Diffs)-1):
				t = time.time()
				Order.extend(range(Diffs[i],Diffs[i+1]))
				
				Get = []
				for l in HH.keys():
					Get += [len(X) + j  for j in HH[l][2][range(HH[l][0][i],HH[l][1][i])] if len(set(DiffAtts[i]).intersection(Aggregates.rowdata['Aggregates'][j].split(','))) > 0 and set(Aggregates.rowdata['Aggregates'][j].split(',')) == set(l.split(','))]
	
				Order.extend(Get)
			
			return dv.datavstack([X,Aggregates])[Order]
				
		
def AddOrReplaceColumns(X,Cols):
	"""
	technical dependency of .aggregate_in method
	"""
	return numpy.rec.fromarrays([X[a] for a in X.dtype.names if a not in Cols.dtype.names] + [Cols[a] for a in Cols.dtype.names], names = [a for a in X.dtype.names if a not in Cols.dtype.names] + list(Cols.dtype.names))


		
def GraySpec(k):
	"""
	For integer argument k, returns list of k gray-scale colors, increasingly light, 
	linearly in the HSV color space, as web hex triplets, .   
	Technical dependency of .aggregate_in method. 
	"""
	ll = .5
	ul = .8
	delta = (ul - ll)/k
	return [GrayScale(t) for t in numpy.arange(ll,ul,delta)]
	
		
def DotDataAggregate(X,On=None,AggFuncDict = None,AggFunc=None,returnsort = False):
	"""
	used to aggregate a DotData on a set of specified factors, using
	specified aggregation functions. 
	
	ARGUMENTS:
	-- X, DotData record array
	-- On, list of column names in X
	-- AggFuncDict -- dictionary where:
		-- keys are some (or all) column names of X that are NOT in 'On'  
		-- values are functions that can be applied to lists or numpy arrays
	-- returnsort, boolean:  because the aggregation function sorts 
	the original data, there's an argsort that can be returned.  If 
	this boolean is true the return value is	[A,s] where A is the 
	aggregated dotdata and s is the sorting permutation on the
	original data set.  If false (the default) the return value is just A
	
	Intuitively, this function will aggregate the dataset X on the 
	columns listed in 'On', so that the resulting aggregate data set 
	has one record for each unique tuples of values in those columns.   
	The more factors listed in On argument, the "finer" is the 
	aggregation, the fewer factors, the "coarser" the aggregation.    
	For example, if On = ['A','B'], the resulting data set will have 
	one record for each unique value of pairs (a,b) in X[['A','B']].   
	
	The AggFuncDict argument specifies how to aggregate the factors
	_not_ listed in 'On' -- the so-called "Off" columns. For example, if
	On = ['A','B'] and 'C' is some other column, then AggFuncDict['C'] 
	is the function that will be used to reduce to a single value
	the (potentially multiple) values in the 'C' column 
	corresponding to unique values in the 'A', 'B' columns.   
	For instance, if
		AggFuncDict['C'] = numpy.mean
	then the result will be that the values in the 'C' column 
	corresponding to a single 'A','B' value will be averaged.  
	
	If an "Off" column is _not_ provided as a key in AggFuncDict, 
	a default aggregator function will be used:  the sum 
	function for numerical columns, concatenation for string columns.    
		
	"""
	
	names = X.dtype.names
	
	if On == None:
		On = []
	elif isinstance(On,str):
		On = On.split(',')
		
	assert all([o in names for o in On]), "Axes " + str([o for o in On if o not in names]) + " can't be found."
	Off = set(names).difference(On)
	
	if AggFuncDict == None:
		AggFuncDict = {}
		
	if AggFunc != None:
		AggFuncDict.update( dict([(o,AggFunc) for o in Off if o not in AggFuncDict.keys()]) )
		
	NotProvided = Off.difference(AggFuncDict.keys()) if AggFuncDict else Off
	if len(NotProvided) > 0:
		print 'No aggregation function provided for axes: ' , NotProvided, 'so assuming "sum".'
		AggFuncDict.update(dict([(o,sum) for o in NotProvided]))

	if len(On) > 0:	
		if len(On) == 1:
			[D,s] = System.Utils.FastRecarrayUniqify(X[On[0]])
		else:
			[D,s] = System.Utils.FastRecarrayUniqify(X[On])
		X = X[s]
		Diffs = numpy.append(numpy.append([-1],D[1:].nonzero()[0]),[len(D)])
	else:
		Diffs = numpy.array([-1,len(X)])

	ColDict = dict([(o,X[o][Diffs[:-1]+1]) for o in On])
	ColDict.update(dict([(o,[AggFuncDict[o](X[o][Diffs[i]+1:Diffs[i+1]+1]) for i in range(len(Diffs) - 1)]) for o in Off]))

	Columns = [ColDict[n] for n in names]
	
	if returnsort:
		return [DotData(Columns = Columns,names=names,coloring=X.coloring),s]
	else:
		return DotData(Columns = Columns,names=names,coloring=X.coloring)


def DotDataFromPathList(PathList):
	"""
	Opens dotdatas from a list of dot-data paths, assuming they have disjoint 
	columns and identical numbers of rows;  
	then stacks them horizontally, e.g. adding columns side-by-side, 
	aligning the rows. 
	
	"""
	Data = DotData(Path = PathList[0])
	
	for path in PathList[1:]:
		Data = Data.hstack(DotData(Path = path))
	
	return Data

	
def Pivot(X,a,b,Keep=None):
	'''
	Implements pivoting on dotdatas.  
	See http://en.wikipedia.org/wiki/Pivot_table for information about pivot tables.
	
	ARGUMENTS:
	X -- dotdata 
	a,b -- columns in X
	Keep -- list of columns in X
	
	RETURNS:
	--X pivoted on (a,b) with a as the row axis and b values 
	as the column axis. 
	
	So-called "nontrivial columns relative to b" in X are added as
	color-grouped sets of columns, and "trivial columns relative to b"  
	are also retained as cross-grouped sets of columns if they're 
	listed in 'Keep' argument.   
	
	(A column 'c' in X is "trivial relative to b" if for all rows i, X[c][i]  
	can be determined from X[b][i], e.g the elements in X[c] are in 
	many-to-any correspondence with the values in X[b].)
		
	The function will raise an exception if the list of pairs of values
	in X[[a,b]] is not the product of the individual columns values, e.g. 
		X[[a,b]]   ==   set(X[a])   x   set(X[b])   ,
	in some ordering. 
	'''

	for c in [a,b]:
		assert c in X.dtype.names, 'Column ' + c + ' not found'
		
	assert len(X) == len(set(X[[a,b]].tolist())) , 'Pairs of values in columns ' + a + ' and ' + b +  ' must be unique'
	Da = len(set(X[a]))
	Db = len(set(X[b]))
	assert len(X) == Da * Db, 'The set of pairs of values in columns ' + a + ' and ' + b +  ' must be the product of the two sets of values.'	
	X.sort(order = [a,b])
	Bvals = X[b][:Db]
	bnames = [str(bv).replace(' ','') for bv in Bvals]

	othernames = [o for o in X.dtype.names if o not in [a,b]]
	
	assert len(set(othernames).intersection(bnames)) == 0 and a not in bnames, 'Processed values of column ' + b + ' musn\'t intersect with other column names.'

	acol = X[a][::Db]

	Cols = [acol]
	names = [a]
	Trivials = []
	NonTrivials = []
	for c in othernames:
		Z = X[c].reshape((Da,Db)) 
		
		if all([len(set(Z[:,i])) == 1 for i in range(Z.shape[1])]):
			Trivials.append(c)
		else:
			NonTrivials.append(c)
			Cols += [Z[:,i] for i in range(Z.shape[1])]
			names += [bn + '_' + c for bn in bnames]
	D = DotData(Columns = Cols,names = names)
	
	if Keep != None:
		Trivials = set(Trivials).intersection(Keep)
		for c in Trivials:
			X.sort(order=[c])
			cvals = numpy.array(uniqify(X[c]))
			[AA,BB] = System.Utils.fastequalspairs(cvals,X[c])
			
			for (i,cc) in enumerate(cvals):
				print time.time() - t, i, cc
				blist = [str(bv).replace(' ','') for bv in Bvals if bv in X[b][AA[i]:BB[i]]]
				D.coloring[str(cc)] = [a] + [bn + '_' + d for bn in blist for d in NonTrivials]
				for d in NonTrivials:
					D.coloring[str(cc) + '_'  + d]  = [a] + blist
				
	for c in NonTrivials:
		D.coloring[c] = [a] + [bn + '_' + c for bn in bnames]
	for bn in bnames:
		D.coloring['_' + bn] = [a] + [bn + '_' + c for c in NonTrivials]
		
	return D
	
	
def nullvalue(test):
	"""
		Returns a null value for each of various kinds of test values. 
	"""
	return False if isinstance(test,bool) else 0 if isinstance(test,int) else 0.0 if isinstance(test,float) else ''


def SVIter( SVPath, SVDelimiter = '\t', ToLoad = None ):
	"""Iterate over tuples of values from a .tsv file.
	"""

	with open( SVPath ) as f:
		headings = f.readline()
		headings = headings.strip().split( SVDelimiter )

		if ToLoad:
			useColumn = numpy.array( [ h in ToLoad for h in headings ] )

		for line in f:
			vals = line.strip().split( SVDelimiter )

			if not vals: break
			if ToLoad: vals = list( numpy.array( vals )[ useColumn ] )
			yield vals

			
def DotDataFromTblIter( tableIter ):
	"""Construct a DotData from a TableIter"""
	return DotData( names = tableIter.headings, Records = map( tuple, tableIter ) )

def SaveToSV( dotData, outFN, getio = None ):
	"""Save given DotData to a .tsv file"""

	if getio: return dict( depends_on = (), creates = outFN, attrs = dict( piperun_short = True ) )
	dotData.saveToSV( outFN )
	


			
				

			

