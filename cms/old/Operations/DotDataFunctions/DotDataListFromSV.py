import numpy, types, logging, csv
from System.Utils import *
from Operations.TypeInfer import TypeInfer

def DotDataListFromSV(Path, Delimiter = None, DelimiterRegExp = None, LineBreak = None, Header = True, SVHash = None, ToLoad = None, SkipFirstLines = 0, ValueFixer = None, LineFixer = None,
		      returnasrecords = False, SVNames = None ):
	"""Takes a tabular file with the specified delimeter and end-of-line character,
	and returns a list of columns (i.e. each list item is a list of data in one column, 
	converted to native Python types) and a list of column names (column headers).

	These two are returned as [attributes_list, attribute_names].

	Parameters:

	    Path - filename
	    Delimiter - if given, the column separator in the file (if not given, automatically inferred from file extension)
	    LineBreak - if given, the line separator in the file (defaults to newline)
	    Header - whether the first line of the file is a header line
	    SVHash - whether we expect one or more hash-lines (comment lines starting
	    with a hash) at top of file.
	       If both Header and SVHash are True, then we assume the header line also 
	       begins with a hash,
	       and is the last line beginning with a hash.
	    ToLoad - list of attribute names (column names) to load.   each name in this 
	    list must match the name
	       of a column in the file (and you must also have Header == True )
	    SkipFirstLines - this many first lines of the file will be ignored   

	"""

	logging.info( 'Loading DotData from ' + Path + '...' )
	logging.info( 'ToLoad = ' + str( ToLoad ) )

	if Path[-4:] == '.csv':
		if Delimiter == None:
			print("Warning: no Delimiter argument specified, assuming ',' because of file extension.")
			Delimiter = ','
	elif Path[-4:] == '.tsv':
		if Delimiter == None:
			print("Warning: no Delimiter argument specified, assuming '\\t' because of file extension.")
			Delimiter = '\t'

	if Delimiter == None:
		print("Warning: no Delimiter argument specified, defaulting to '\\t'.")
		Delimiter = '\t'
	
	if LineBreak == None or LineBreak == '\n':
		if LineBreak == None:
			print("Warning: no LineBreak argument specified, using '\\n'.")
			LineBreak = '\n'
		F = open(Path,'rU').read().strip(LineBreak).split(LineBreak)
		
	else:		
		F = open(Path,'r').read().strip(LineBreak).split(LineBreak)

	logging.info( 'DotData file contents read in from ' + Path )

	F = F[SkipFirstLines:]

	if LineFixer: F = list(map( LineFixer, F ))
	if DelimiterRegExp:
		if isinstance( DelimiterRegExp, bytes ):
			import re
			DelimiterRegExp = re.compile( DelimiterRegExp )
		F = [DelimiterRegExp.sub( Delimiter, line ) for line in F]
	
	if SVHash != None:
		headerlines = 0
		for line in F:
			if line[:len(SVHash)] == SVHash:
				headerlines += 1
			else:
				break
		if Header:
			print("Assuming that the last line in ", Path, " beginning with '", SVHash, "' has attribute names.")
			F = F[headerlines-1:]
			F[0] = F[0][1:]
		else:
			F = F[headerlines:]				

	# attribute_names - list of column names, or None
	# records_list - list of lines in the file as strings, sans any comment or header lines.
	if Header:
		#Extract attribute names (single header line)
		attribute_names =  list(csv.reader(F[0:1],delimiter=Delimiter))[0]
		#Extract records to a list of strings (skip header)
		records_list = F[1:]
	else:
		#No header in SV
		attribute_names = SVNames
		#Extract records to a list of strings
		records_list = F

	print('attribute_names=', attribute_names)
		
	# Parse records into a list of lists
	# records_parsed - list of lists of strings, where each inner list of strings is the parsing of one file line
	#    into individual values.
	#records_parsed =  [record.strip(LineBreak).split(Delimiter) for record in records_list]
	records_parsed = list(csv.reader(records_list,delimiter=Delimiter))

	if returnasrecords:
		return [records_parsed,attribute_names]

	if attribute_names:
		for r in records_parsed:
			if not len( r ) >= len( attribute_names ):
				print('r=', r)
			assert len( r ) >= len( attribute_names )
	
	# Type the columns	and save to a list of (attribute) columns
	# attributes_list - list of columns, where each column is a list of values of one Python type.
	attributes_list = [TypeInfer([ record[column] if not ValueFixer else ValueFixer( record[column] ) for record in records_parsed]) for column in range(len(attribute_names  if attribute_names else records_parsed[0]))]
	
	if Header and ToLoad != None:
		# Restrict the columns to the subset requested by the user
		for TL in ToLoad:
			if TL not in attribute_names:
				print('column ', TL, ' not in ', attribute_names)
			
		load_ind = [attribute_names.index(TL) for TL in ToLoad]
		attribute_names = [ attribute_names[ i ] for i in load_ind ]
		attributes_list  = [ attributes_list[ i ] for i in load_ind ]

	logging.info( 'Loaded DotData from ' + Path + ': ' + str(len( attributes_list[0] )) + ' rows, ' + str( len( attributes_list ) ) + ' columns.' )
	
	return [attributes_list, attribute_names]
