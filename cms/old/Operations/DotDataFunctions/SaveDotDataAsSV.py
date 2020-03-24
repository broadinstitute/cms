from System.Utils import *
import csv, logging, os, time, itertools

def SaveDotDataAsSV(Data, TargetFileName,sep = None,linesep=None, UseHeader = True):

	if os.path.exists( TargetFileName ):
		os.remove( TargetFileName )
		time.sleep( 5 )
	
	if TargetFileName[-4:] == '.csv':
		if sep == None:
			sep = ','
		elif sep != ',':
			print("WARNING: You've specified a .csv ending for the file", TargetFileName , "but have given something other than a comma as the delimiter.")
	elif TargetFileName[-4:] == '.tsv':
		if sep == None:
			sep = '\t'
		elif sep != '\t':
			print("WARNING: You've specified a .tsv ending for the file", TargetFileName , "but have given something other than a tab as the delimiter.")

	if sep == None:
		#print 'NOTICE:  No delimiter specified, using tab.' 
		sep = '\t'

	if linesep == None:
		print('NOTICE:  No lineseparator specified, using \\n.') 
		linesep = '\n'

	Header = sep.join(Data.dtype.names)
	typevec = []
	ColStr = []
	UseComplex = True

	if False:
		for attribute_name in Data.dtype.names:
			typevec.append(Data.dtype[attribute_name].name.strip('0123456789').rstrip('ing'))
			D = Data[attribute_name]
			if D.ndim > 1:
				D = D.flatten()	
			if typevec[-1] == 'str':
				if sum([sep in d for d in D]) > 0:
					print("\nWARNING: An entry in the '" + attribute_name +"' column contains at least one instance of the delimiter '" + sep + "', and therefore will use Python csv module quoting convention (see online documentation for Python's csv module).  You may want to choose another delimiter not appearing in records, for performance reasons.\n")
					UseComplex = True
					break
				else:
					ColStr.append(D)
			else:
				ColStr.append(str(D.tolist()).strip('[]').split(', '))


	with open(TargetFileName,'wb') as F:
		if UseHeader: F.write(Header + linesep)
		if UseComplex:
			csv.writer(F,delimiter = sep, lineterminator = linesep).writerows(Data if len( Data.dtype ) > 1 else \
												  map( lambda x: (x,), Data ))
		else:
			if len( D ) > 0:
				F.write(linesep.join([sep.join([col[i] for col in ColStr]) for i in range(len(ColStr[0]))]))

	logging.info( 'Saved DotData with %d rows to %s' % ( len( Data ), TargetFileName ) )
