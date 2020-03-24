'''
"Vertical stacking" of DotDatas, e.g. adding rows. 
'''
import numpy
from System.Utils import uniqify, ListUnion, SimpleStack1, ListArrayTranspose
import Classes.DotData, pickle

def datavstack(ListOfDatas):
#if all([isinstance(l,DotData) for l in ListOfDatas]):	
#else:
#	return numpy.vstack(ListOfDatas)
	CommonAttributes = set(ListOfDatas[0].dtype.names)
	for l in ListOfDatas:
		CommonAttributes = CommonAttributes.intersection(set(l.dtype.names))
	CommonAttributes = list(CommonAttributes)
	CommonAttributes = [CommonAttributes[j] for j in [CommonAttributes.index(e) for e in ListOfDatas[0].dtype.names if e in CommonAttributes]]
	
	if len(CommonAttributes) == 0:
		try:
			return numpy.row_stack(ListOfDatas)
		except:
			print("The data arrays you tried to stack couldn't be stacked.")
	else:
		
		A = SimpleStack1([l[CommonAttributes] for l in ListOfDatas])
			
		if all(['coloring' in dir(l) for l in ListOfDatas]):
			restrictedcoloring = dict([(a,list(set(ListOfDatas[0].coloring[a]).intersection(set(CommonAttributes)))) for a in list(ListOfDatas[0].coloring.keys())])
			for l in ListOfDatas[1:]:
				restrictedcoloring.update(dict([(a,list(set(l.coloring[a]).intersection(set(CommonAttributes)))) for a in list(l.coloring.keys())]))
		else:
			restrictedcoloring = {}
		if not all ([l.rowdata == None for l in ListOfDatas]):
			rowdata = SafeSimpleStack([l.rowdata if l.rowdata != None else numpy.array(['']*len(l)) for l in ListOfDatas ])		
		else:
			rowdata = None

		return Classes.DotData.DotData(Array = A, dtype = A.dtype, coloring = restrictedcoloring,rowdata = rowdata)
		

def SafeSimpleStack(seq):
	'''
		Vertically stack sequences numpy record arrays. 
		Avoids some of the problems of numpy.v_stack
	'''

	names =  uniqify(ListUnion([list(s.dtype.names) for s in seq if s.dtype.names != None]))
	formats =  [max([s.dtype[att] for s in seq if s.dtype.names != None and att in s.dtype.names]).str for att in names]
	D = numpy.rec.fromarrays([ListUnion([s[att].tolist() if (s.dtype.names != None and att in s.dtype.names) else [nullvalue(format)] * len(s) for s in seq]) for (att,format) in zip(names,formats)], names = names)
	D = D.dumps()
	return numpy.loads(D)
	
def nullvalue(format):
	return 0 if format.startswith(('<i','|b')) else 0.0 if format.startswith('<f') else ''