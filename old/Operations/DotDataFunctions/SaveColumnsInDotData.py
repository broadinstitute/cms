from System.Utils import *

def SaveColumnsInDotData(Data, TargetFolderName):
	if TargetFolderName[-1] != '/':
		TargetFolderName = TargetFolderName + '/'
		
	attribute_set = set(Data.dtype.names)

	for attribute_name in list(attribute_set):
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