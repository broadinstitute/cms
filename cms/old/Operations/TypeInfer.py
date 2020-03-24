import numpy

TrueStrings = ( 'True', 'true' )
FalseStrings = ( 'False', 'false' )
BoolStrings = TrueStrings + FalseStrings

def TypeInfer(column):
	"""Take a list of strings, and attempts to infer a numeric datatype that fits them all.
	If the strings are all integers, returns a list of corresponding Python integers.
	If the strings are all floats, returns a list of corresponding Python floats.
	Otherwise, returns the original list of strings.

	Used to determine the datatype of a column read from a tab- or comma-separated text file.
	"""
	try:
		return [int(x) if x != '' else numpy.nan for x in column]
	except:
		try:
			return [float(x) if x != '' else numpy.nan for x in column]
		except:
			if len(column) > 0 and column[0] in BoolStrings and all([ x in BoolStrings for x in column]):
				return [ x in TrueStrings for x in column ]
			else:
				return column
