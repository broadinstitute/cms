from System.Utils import *
from numpy import int64
import logging
import cPickle
from Operations.MiscUtil import dbg, chomp

def DotDataListFromDirectory(path, data_list = None, attribute_names = None, rootpath = '', rootheader = None, coloring = None, ToLoad = None, Nrecords = None):
	"""Load DotData from a .data/ directory.  
	Each column is stored in a separate file 
		(for column named myColumn the file is myColumn.datatype.csv
		where datatype is the data type of the column data ( int, int64, float, str, ... ), 
		and the ordered list of columns is in a separate header file.)

	See also: SaveDotData()
	"""

	logging.info( 'Loading DotData from ' + path + '...' )

	if IsDir(path) and path[-1] != '/':
		path = path + '/'
	if data_list == None:
		data_list = []
	if attribute_names == None:
		attribute_names = []
	if rootpath == '':
		rootpath = path
	if rootheader == None:
		#Will use any .header.txt files to order attributes (this is not required) 
		rootheader = []
		if IsDir(path):
			H = [h for h in listdir(path) if h[-11:] == '.header.txt']	
			expectedheader = rootpath.strip('/').split('/')[-1][:-5] + '.header.txt'
			if len(H) == 1:
				if H[0] != expectedheader:
					print "Warning: the file, ", rootpath + H[0], " is being used to determine the order of the attribute names, even though ", rootpath + expectedheader, " was expected."
				rootheader =  open_for_read(path + H[0])[0].read().strip('\n').split('\n')
			elif len(H) > 1:
				if expectedheader in H:
					rootheader = open_for_read(path + H[0])[0].read().strip('\n').split('\n')
					print "Warning: the file, ", rootpath + expectedheader, " is being used to determine the order of the attribute names. Multiple .header.txt files were provided."
				else:
					print "Warning: there are multiple .header.txt files in the directory ", rootpath, ", and so none will be used."
	if coloring == None:
		coloring = {}	

	if IsDir(path):
		L = [l for l in listdir(path) if l[-5:] == '.data' or l[-4:] == '.csv']
		coloring_names = path[len(rootpath):].split('.data/')[:-1]
	else:
		L = [path]
		coloring_names = []
	CSVList = []

	for l_idx, l in enumerate( L ):
		logging.info( 'Loading column %d of %d: %s' % ( l_idx, len( L ), l ) )
		parsed_filename = l.split('.')
		attribute_name = '.'.join(parsed_filename[:-2]).split('/')[-1]
		if parsed_filename[-1] == 'csv' and (ToLoad == None or attribute_name in ToLoad): 				
			CSVList += [attribute_name]
			if attribute_name not in attribute_names:
				attribute_string_pre = chomp(open_for_read(path + l if l != path else path)[0].read())
				attribute_string = attribute_string_pre.split('\n') if attribute_string_pre else []
				if Nrecords == None:
					Nrecords = len(attribute_string)
				if len(attribute_string) == Nrecords:
					try:
						type_data = eval(parsed_filename[-2])
						#for j, a in enumerate(attribute_string):
						#	try: float(a)
						#	except ValueError:
						#		print ' at line ', j, ' in ', l, ' bad value ', a, '.'

						if type_data == bool:
							attribute_typed = [attribute_string[j] == 'True' or attribute_string[j] != 'False' and bool(int(attribute_string[j])) for j in range(len(attribute_string))]
						else:
							attribute_typed = [type_data(attribute_string[j]) for j in range(len(attribute_string))]
						if len(rootheader) > 0 and attribute_name in rootheader:
							indvec = [attribute_names.index(j) for j in rootheader[:rootheader.index(attribute_name)] if j in attribute_names]
							insert_ind = max(indvec) + 1 if len(indvec) > 0 else 0
							data_list.insert(insert_ind, attribute_typed)
							attribute_names.insert(insert_ind, attribute_name)
						else:
							data_list += [attribute_typed]
							attribute_names += [attribute_name]		 		
					except:
						print "Warning: the data in the .csv file ", path + l if l != path else path, " does not match the given data type, ", parsed_filename[-2], ", and was not loaded"
						raise
				else:
					print "Warning: the column" + path + (l if l != path else path) + " has " + str(len(attribute_string)) + " records, which does not agree with the number of records in first column loaded, '" + attribute_names[0] + "', which has " + str(Nrecords) + " records -- only the first column loaded, as well as all the other columns which also have " + str(Nrecords) + " records, will be loaded."
					raise
							
		elif parsed_filename[-1] == 'data' and IsDir(path):
		#File is a .data subdirectory
			colorname = '.'.join(parsed_filename[:-1]).split('/')[-1]
			if ToLoad == None or not colorname in ToLoad:
				[data_list, attribute_names, coloring,rowdata] = DotDataListFromDirectory(path+l, data_list, attribute_names, rootpath, rootheader, coloring, ToLoad, Nrecords)
			elif colorname in ToLoad:
				[data_list, attribute_names, coloring,rowdata] = DotDataListFromDirectory(path+l, data_list, attribute_names, rootpath, rootheader, coloring, None, Nrecords)			
			
	if len(CSVList) > 0:
		for color in coloring_names:
			coloring[color] = list(set(coloring[color] + CSVList if color in coloring.keys() else CSVList))

	if path == rootpath and len(rootheader) > 1:
	#Use header in top directory to order attributes and colorings
		for color in coloring:
			coloring[color] = [name for name in attribute_names if name in coloring[color]]

	if len(data_list) > 0 and len(data_list[0]) > 0:
		logging.info( 'Loaded DotData from ' + path + ': ' + str(len( data_list[0] )) + ' rows, ' + str( len( data_list ) ) + ' columns.' )

	assert len( data_list ) == len( attribute_names )
	
	if  IsDir(path):
		if '__rowdata__.pickle' in listdir(path):
			try:
				rowdata = cPickle.load(open(Backslash(path) + '__rowdata__.pickle','r'))
				if len(rowdata) != len(data_list[0]):
					rowdata = None
			except:
				rowdata = None
		else:
			rowdata = None
	else:
		rowdata = None

	return [data_list, attribute_names, coloring,rowdata]
