import Classes.DotData 
import numpy

def datahstack(ListOfDatas):
	'''
	For "horizontal stacking" DotDatas, e.g. adding columns. 
	Requires all DotDatas to have same number of records. 
	If column appears to two DotDatas in ListOfDats, first is taken. 
	'''
#if ListAnd([isinstance(l,Classes.DotData.DotData) for l in ListOfDatas]):	
#else:
#	return numpy.hstack(ListOfDatas)
	dtype = []
	data = []
	coloring = {}
	for a in ListOfDatas:
		if 'coloring' in dir(a):
			coloring.update(a.coloring)
		dtype += a.dtype.descr
		for att in a.dtype.names:
			data += [a[att]]
	dtype = numpy.dtype(dtype)
	rowdata = SafeColumnStack([l.rowdata for l in ListOfDatas if l.rowdata != None])
	return Classes.DotData.DotData(Columns = data, dtype = dtype,coloring=coloring,rowdata=rowdata)



def SafeColumnStack(seq):
	'''
	'''
	seq = [l for l in seq if l != None]
	if len(seq) > 0:
		done = []
		Columns = []
		names = []
		for l in seq:
			Columns += [l[n] for n in l.dtype.names if n not in done]
			names += [n for n in l.dtype.names if n not in done]
			done += [n for n in l.dtype.names if n not in done]
		
		return numpy.rec.fromarrays(Columns,names=names)
	else:	
		return None