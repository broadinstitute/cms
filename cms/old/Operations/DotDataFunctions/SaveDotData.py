from System.Utils import *
import logging, pickle

def SaveDotData(Data, TargetFolderName, HeaderOn = 1):

	logging.info( 'Saving DotData (%d rows, %d cols) to %s' % ( Data.numRows(), Data.numCols(), TargetFolderName ) )

	if TargetFolderName[-1] != '/':
		TargetFolderName = TargetFolderName + '/'
		
	MakeDir(TargetFolderName)

	KeyList = list(Data.coloring.keys())

	Nkeys = len(KeyList)		
	pairwise = [[set(Data.coloring[key1]) > set(Data.coloring[key2]) for key1 in KeyList] for key2 in KeyList]

	attribute_set = set(Data.dtype.names)
	
	for i in range(Nkeys):
		if sum(pairwise[i]) == 0:
			SaveDotData(Data[KeyList[i]], TargetFolderName + KeyList[i] + '.data/', HeaderOn)
			attribute_set = attribute_set - set(Data[KeyList[i]].dtype.names)
	attribute_list = list(attribute_set)		
	for attribute_num, attribute_name in enumerate(attribute_list):
		if Data.numRows() > 1000 and (attribute_num % 100) == 0:
			logging.info( 'Writing column %d of %d', ( attribute_num, len( attribute_list ) ) )
		typestr = Data.dtype[attribute_name].name.strip('0123456789').rstrip('ing')
		F = open_for_write(TargetFolderName + attribute_name + '.' + typestr + '.csv')[0]
		D = Data[attribute_name]
		if D.ndim > 1:
			D = D.flatten()	
		if typestr == 'str':
			F.write('\n'.join(D))
		else:
			F.write(str(D.tolist()).strip('[]').replace(', ','\n'))
			#F.write(str(D.tolist()).strip('[]'))
			F.close()
			
	if HeaderOn and Data.dtype.names > 1:
		G = open_for_write(TargetFolderName + TargetFolderName[:-6].split('/')[-1] + '.header.txt')[0]
		G.write('\n'.join(Data.dtype.names))
		G.close()
		
	if Data.rowdata != None:
		G = open(TargetFolderName + '__rowdata__.pickle','w')
		pickle.dump(Data.rowdata,G)
		G.close()

	logging.info( 'Saved DotData (%d rows, %d cols) to %s' % ( Data.numRows(), Data.numCols(), TargetFolderName ) )
