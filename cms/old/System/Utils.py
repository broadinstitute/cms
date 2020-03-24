
"""Miscellaneous utilities.
"""

import sys, os, inspect, time, shutil, types, pickle, numpy, re
from string import Template
from traceback import print_exc
from numpy import nan, isnan
#import Tools.pp.pp as pp

def TypeInfer(column):
	"""Take a list of strings, and attempt to infer a datatype that fits them all.
	If the strings are all integers, returns a list of corresponding Python integers.
	If the strings are all floats, returns a list of corresponding Python floats.
	Otherwise, returns the original list of strings.

	Typically used to determine the datatype of a column read from a tab- 
	or comma-separated text file.
	"""
	try:
		return [int(x) if x != '' else numpy.nan for x in column]
	except:
		try:
			return [float(x) if x != '' else numpy.nan for x in column]
		except:
			return column


def PermInverse(s):
	'''
		Fast invert a numpy permutation.
	'''
	X = numpy.array(list(range(len(s))))
	X[s] = list(range(len(s)))
	return X


class BadCheckError: 
	'''
		Error class used for raising I/O errors (maybe should be moved to system_io_override?)
	'''
	def __init__(self,iofunc,readfiles,writefiles, Dependencies,Creates):
		print("\nCHECK_ERROR: An I/O exception occured in function", iofunc, ": either the files" , readfiles , "aren't in" , Dependencies, "or the files", writefiles, "aren't in ", Creates, '. \n')
		

def RecursiveFileList(ToList,Avoid=None):
	'''
	Given list of top-level directories, recursively gets a list of files in the directories. 
	
	ARGUMENTS:
	--ToList = list of top-level diretories as python list of path strings
	--Avoid = List of regexps of directory name patterns to NOT look in 
		e..g if a directory name matches at any level, the function will not 
		look further into that directory. 
	'''
	if isinstance(ToList,list):
		return ListUnion([RecursiveFileList(x,Avoid) for x in ToList])
	elif IsFile(ToList):
		return [ToList]
	elif IsDir(ToList):
		if ToList[-1] != '/':
			ToList += '/'
		L = []
		for l in listdir(ToList):
				if IsFile(ToList + l):
					L += [ToList + l]
				elif IsDir(ToList + l):
					if Avoid == None or not any([re.search(a,l) for a in Avoid]):
					#if Avoid == None or l not in Avoid :
						L += RecursiveFileList(ToList + l,Avoid)
		return L
	else:
		return []
	
def Max(x):
	if any(isnan(x)):
		return nan
	else:
		return max(x)

def enumeratefrom(i,A):
	assert i >= 0, 'index must be larger than 0 for enumeratefrom to work'
	return list(enumerate(['',]*i + list(A)))[i:]

def uniqify(seq, idfun=None): 
	'''
	Relatively fast pure python uniqification function that preservs ordering
	ARGUMENTS:
		seq = sequence object to uniqify
		idfun = optional collapse function to identify items as the same
	RETURNS:
		python list with first occurence of each item in seq, in order
	'''
	

	# order preserving
	if idfun is None:
		def idfun(x): return x
	seen = {}
	result = []
	for item in seq:
		marker = idfun(item)
		# in old Python versions:
		# if seen.has_key(marker)
		# but in new ones:
		if marker in seen: continue
		seen[marker] = 1
		result.append(item)
	return result
	

def FastArrayUniqify(X):
	'''
	Very fast uniqify routine for numpy array 
	ARGUMENT:
		X = a numpy array
	RETURNS:
		[D,s] where s is a permutation that will sort X, and D is the list of "first 
		differences" in the sorted verion of X
		This can be used to produce a uniqified version of X by simply taking:
			X[s][D]
		or
			X[s[D.nonzero()[0]]]
		But sometimes the information of D and s is useful. 
			
	'''
	s = X.argsort()
	X = X[s]
	return [numpy.append([True],X[1:] != X[:-1]),s]

	
def FastRecarrayUniqify(X):
	'''
	Record array version of FastArrayUniqify.
	ARGUMENT:
		X = numpy record array
	RETURNS:
		[D,s] where s is a permutation that will sort all the columsn of X in some order,
		and D is a list of "first differences" in the sorted version of X
		This can be used to produce a uniqified version of X by simply taking:
			X[s][D]
		or
			X[s[D.nonzero()[0]]]
		But sometimes the information of D and s is useful. 			
	'''
	N = X.dtype.names
	s = X.argsort(order=N)
	X = X[s]
	return [numpy.append([True],X[1:] != X[:-1]),s]


def DirName(path):
	'''
	utility that gets dir name;  sometimes this is the right thing to use intead of os.path.dirname itself
	'''
	if path[-1] == '/':
		path = path[:-1]
	return os.path.dirname(path)
				
				
def open_for_read(ToRead):
	return [open(ToRead,'r'),True]


def open_for_read_universal(ToRead):
	return [open(ToRead,'rU'),True]

def open_for_write(ToWrite):
	return [open(ToWrite,'w'),True]

def open_for_append(ToAddTo):
	return [open(ToAddTo,'a'),True]
		

def chkExists( path ):
    """If the given file or directory does not exist, raise an exception"""
    if not os.path.exists(path): raise IOError("Directory or file %s does not exist" % path)


def PathCompress(path):
	return path.replace('../','')


def redirect(x,To):
	'''
	utility that 'redirects' a path name from to a different directory. 
	ARGUMENTS:
		x = path to redirect
		To = location to redirect x to
	RETURNS:
		path created by taking file component of x and appending to To
	
	'''
	x = PathCompress(x)
	if '/' in x[:-1]:
		j = max([i for i in range(len(x)-1) if x[i] == '/'])
		if To[-1] != '/':
			To.append('/')
		return To + x[:j+1].replace('/','__') + x[j+1:]
	else:
		return To + x

def ListArrayTranspose(L):
	'''
	Tranposes the simple array presentation of a list of lists (of equal length). 
	Argument:
		L = [row1, row2, ...., rowN]
		where the rowi are python lists of equal length. 
	Returns:	
		LT, a list of python lists such that LT[j][i] = L[i][j]. 
	'''
	return [[row[i] for row in L] for i in range(len(L[0]))]


def GetFunctionsDefinedInModule(Module):
	'''
	Given live module, use inspect module to capture list of functions in
	module whose definition attribute is equal to the module name.
	'''
	Z = inspect.getmembers(Module)
	return dict([(a[0],a[1]) for a in Z if type(a[1]) == types.FunctionType and a[1].__module__ == Module.__name__ ])


def GetFunctionsMentionedInModule(Module):
	'''
	Given live module, use inspect module to capture list of functions that 
	appear somewhere live in the module (e.g. are loaded when the 
	module is imported)
	'''
	Z = inspect.getmembers(Module)
	return dict([(a[0],a[1]) for a in Z if type(a[1]) == types.FunctionType])


def RedirectList(ToRedirect,To):
	return tuple([redirect(x,To) for x in ToRedirect])
	

def FixedPath(path):
	'''
	Does some path compression, like os.path.normpath but proper for the 
	Data Environment (consider replacing with references to os.path.normpath)
	'''
	if path[:2] == './':
		path = path[2:]
	path = path.replace('/./','/')
	if path[:3] != '../':
		path = '../Temp/' + path
	return path


def GetDataEnvironmentDirectory():
	x = os.environ
	if 'DataEnvironmentDirectory' in list(x.keys()):
		return x['DataEnvironmentDirectory']
	else:
		print('DataEnvironmentDirectory not an environment variable, assuming it is ' , os.getcwd()[:os.getcwd().find('/')] + '/')
		return os.getcwd()[:os.getcwd().find('/')] + '/' 
	
def PathAlong(a,b):
	'''
	returns true when path a is inside the filetree under path b. 
	'''
	return PathStrictlyAlong(a,b) or (FixedPath(a) == FixedPath(b))
	
	
def PathStrictlyAlong(a,b):
	'''
	returns true when path a is strictly insider the filetree under path b. 
	'''
	a = FixedPath(a); b = FixedPath(b)
	return len(a) > len(b) and a[:len(b)] == b and (a[len(b)] == '/' or b[-1] == '/')


def funcname():
	'''
	returns name of function in call stack in which funcname() is being called 
	'''
	return sys._getframe(1).f_code.co_name 
	
def caller(level = 2):
	'''
	returns name of function in call stack which is calling the 
	function that  is calling caller()
	'''
	return sys._getframe(level).f_code.co_name 


def callermodule():
	'''
	returns name of module from which the function that  is calling caller() was imported
	'''
	return sys._getframe(2).f_code.co_filename
		
		
def TimeStamp(T = None):	
	'''
		deprecated TimeStamp function (should be replaced by time module formatters)
	'''
	return time.localtime(T).__str__().replace(',','_').replace(' ','').strip('()')
	
	
def Union(ListOfSets):
	'''
		takes python list of python sets [S1,S2, ..., SN] and returns their union
	'''
	u = set([])
	for s in ListOfSets:
		u.update(s)
	return u

def ListUnion(ListOfLists):
	'''
	takes python list of python lists
	
	[[l11,l12, ...], [l21,l22, ...], ... , [ln1, ln2, ...]] 
	
	and returns the aggregated list 
	
	[l11,l12, ..., l21, l22 , ...]
	'''
	u = []
	for s in ListOfLists:
		if s != None:
			u.extend(s)
	return u

def GetDefaultVal(func,varname,NoVal = None):
	'''
	given a live python function object "func", return the default value for 
	variable with name "varname" if it exists as a keyword variable, else 
	return NoVal
	'''
	V = inspect.getargspec(func)
	if varname in V[0]:
		varname_pos = min([i for i in range(len(V[0])) if V[0][i] == varname]) - len(V[0])
		return V[3][varname_pos]
	else:
		return NoVal


def MakeDir(DirName,creates = ()):
	'''
	is a "strong" directory maker -- if DirName already exists, this deletes it first 
	'''
	if os.path.exists(DirName):
			delete(DirName)
	os.mkdir(DirName)


def MakeDirs(DirName,creates = ()):
	'''
	"strong" version of os.makedirs
	'''
	if os.path.exists(DirName):
			delete(DirName)
	os.makedirs(DirName)

	
def strongcopy(tocopy,destination,use2 = False):
	'''
	"strong" version of copy -- if destination already exists, it removes 
	it first befire copying
	'''
	if os.path.isfile(tocopy):
		if use2:
			shutil.copy2(tocopy,destination)		
		else:
			shutil.copy(tocopy,destination)	
	elif os.path.isdir(tocopy):	
		if os.path.exists(destination):
			delete(destination)
		shutil.copytree(tocopy,destination)			
			
def delete(ToDelete):
	'''
	unified "strong" version of delete that uses os.remove for a file 
	and shutil.rmtree for a directory tree
	'''
	if os.path.isfile(ToDelete):
		os.remove(ToDelete)
	elif os.path.isdir(ToDelete):
		shutil.rmtree(ToDelete)		


def Log(s):
	"""Log a debug message"""
	# by default, do nothing.  


def TemplateInstance(templatepath,outpath, **kws):
	'''
	Apply python template at "templatepath" with substitutions from
	keyword arguments **kws passed in, and write out result to outpath
	Useful for html templating
	'''
	TemplateString = Template(open(templatepath,'r').read())
	OutFile = open(outpath,'w')
	NewString = TemplateString.substitute(kws)
	OutFile.write(NewString)
	OutFile.close()


def MakeT(r):
	'''
	If input 'r' is a comma-delimited string, return tuple split on 
	commas, else return tuple(r)
	'''
	return tuple(r.split(',')) if isinstance(r,str) else tuple(r)	


def getKalong(LL1,LL2,k):
	'''
	Fast version of "K-along" paths. 
	
	ARGUMENTS:
	--LL1 = numpy array of paths
	--LL2 = sorted numpy array of paths
	--k = nonnegative integer

	RETURNS:
	[A,B] where A and B are numpy arrays of indices in LL1 such that
	LL2[A[i]:B[i]] contains precisely those paths in B that are k 
	directory levels down from LL1[i] -- as path strings (no actual
	directory testing is done).  A[i] = B[i] = 0 if no paths in LL2 are k 
	directory levels down from LL1[i]
	
	E.g. if 
	
	LL1 = numpy.array(['../Data/Dan_Data/', '../Users/DanYamins/','../Users/SijiaWang/'])
	
	and
	
	LL2 = numpy.array(['../Data/Dan_Data/NPR_Puzzle_Solutions',
	'../Data/Dan_Data/RandomData','../Users/DanYamins/Finance/'])
	
	then
		
	getKalong(LL1,LL2,1) = [A,B] = [[0,2,0],[2,3,0]]
			
	'''

	SlashList1 = numpy.array([len(y.split('/')) - (1 if y[-1] == '/' else 0) for y in LL1])
	SlashList2 = numpy.array([len(z.split('/')) - (1 if z[-1] == '/' else 0) for z in LL2])
	
	M1 = numpy.rec.fromarrays([LL1,SlashList1],names=['Val','S'])
	M2 = numpy.rec.fromarrays([LL2,SlashList2],names=['Val','S'])
	s1 = M1.argsort(order=['S','Val']); M1 = M1[s1]; SlashList1 = SlashList1[s1]
	s2 = M2.argsort(order=['S','Val']); M2 = M2[s2]; SlashList2 = SlashList2[s2]
	
	Max = max(SlashList1)
	Min = min(SlashList1) 
	W = numpy.zeros(len(LL1),int)
	U = numpy.zeros(len(LL1),int)	
	
	for i in range(Min,Max+1):
		I1 = (SlashList1 == i); nz1 = I1.nonzero()[0]; 
		if len(nz1) > 0:
			st1 = min(nz1)
			I2 = (SlashList2 == i + k); nz2 = I2.nonzero()[0]; 
			if len(nz2) >0:
				st2 = min(nz2)
				[A,B] = getpathalong(M1['Val'][I1],M2['Val'][I2])
				W[st1:st1 + len(nz1)] = A + st2
				U[st1:st1 + len(nz1)] = B + st2

	return [s1,s2,W,U]
		

def maximalpathalong(YY,ZZ):
	'''
	Fast function for determining indices of elements of YY such that
	they are "path along" some element of ZZ, for numpy arrays YY 
	and ZZ.  When YY[i] is path along several element of ZZ, returns 
	index of the first occurence of the closest path.  If YY[i] is not 
	path-along any elements of Z, returns ''. 
		
	'''

	ZZ = ZZ.copy() ; ZZ.sort()
#	Y =  numpy.array([y + '/' if y[-1] != '/' else y for y in YY])
#	Z = numpy.array([y + '/' if y[-1] != '/' else y for y in ZZ])
	Y = YY
	Z = ZZ
	SlashList = numpy.array([len(y.split('/')) - (1 if y[-1] == '/' else 0) for y in Z])
	Max = max(SlashList) if len(SlashList) > 0 else 0
	Min = min(SlashList) if len(SlashList) > 0 else 0
	
	C = -1*numpy.ones((len(YY),),int)
	for i in range(Min,Max+1):
		T = numpy.array(['/'.join(z.split('/')[:i]) + ('/' if len(z.split('/')) > i else '')  for z in Y])
		[A,B] = fastequalspairs(T,Z)
		M = (B>A)
		C[M] = B[M] - 1
	
	z = numpy.append(ZZ,[''])
	return z[C]
	
	
def getpathalong(YY,ZZ):
	'''
	Fast version of path long for numpy arrays. 
	
	ARGUMENTS:
		LL1 = numpy array of paths
		LL2 = sorted numpy array of paths
	
	RETURNS:
	[A,B] where A and B are numpy arrays of indices in LL1 such 
	that LL2[A[i]:B[i]] contains precisely those paths in B that are
	in the directory tree of paths in LL1[i] ( as path strings -- no actual
	filesystem existence testing is done).  A[i] = B[i] = 0 if no paths in 
	LL2 are in the directory tree under LL1[i]
		
	E.g. if 
	
	LL1 = numpy.array(['../Data/Dan_Data/', '../Users/DanYamins/','../Users/SijiaWang/'])
	
	and
	
	LL2 = numpy.array(['../Data/Dan_Data/NPR_Puzzle_Solutions',	'../Data/Dan_Data/NPR_Puzzle_Solutions/AmericaPensacolaPuzzle/',
	'../Data/Dan_Data/RandomData','../Users/DanYamins/Finance/'])
	
	then
	
	getpathalong(LL1,LL2) = [A,B] = [[0,3,0],[3,4,0]]
		
		
	'''

#	Y = numpy.array([y[:-1] if y[-1] == '/' else y for y in YY])
#	Z = numpy.array([z[:-1] if z[-1] == '/' else z for z in ZZ])
	Y =  numpy.array([y + '/' if y[-1] != '/' else y for y in YY])
	Z = numpy.array([y + '/' if y[-1] != '/' else y for y in ZZ])
	SlashList = numpy.array([len(y.split('/')) - (1 if y[-1] == '/' else 0) for y in Y])
	Max = max(SlashList) if len(SlashList) > 0 else 0
	Min = min(SlashList) if len(SlashList) > 0 else 0
	W = numpy.zeros(len(Y),int)
	U = numpy.zeros(len(Y),int)

	for i in range(Min,Max+1):
		T = numpy.array(['/'.join(z.split('/')[:i]) + ('/' if len(z.split('/')) > i else '')  for z in Z ])    #get i-reduced slash list from Z, call it T
		R = (T[1:] != T[:-1]).nonzero()[0] 
		R = numpy.append(R,numpy.array([len(T)-1]))
		M = R[R.searchsorted(list(range(len(T))))]
		#get set of guys in Y with i slashes, call it L
		L = (SlashList == i)
		H = Y[L]
		D = T.searchsorted(H)
		T = numpy.append(T,numpy.array([0]))
		M = numpy.append(M,numpy.array([0]))
		W[L] = (T[D] == H) * D
		U[L] = (T[D] == H) * (M[D] + 1)
	
	return [W,U]
	
	
def getpathalongs(Y,Z):
	'''
		Returns numpy array of indices i in numpy array Y such that Y[i] is 
		path-along some path string in Z
	'''

	s = Z.argsort()
	Z = Z[s]
	[A,B] = getpathalong(Y,Z)
	L = ListUnion([list(range(A[i],B[i])) for i in range(len(A)) if A[i] < B[i]])
	return s[L]


def getpathstrictlyalong(YY,ZZ):
	'''
		Version of getpathalong that requires "strictly path along" 
	'''
	[A,B] = getpathalong(YY,ZZ)
	YY =  numpy.array([y + '/' if y[-1] != '/' else y for y in YY])
	ZZ = numpy.array([y + '/' if y[-1] != '/' else y for y in ZZ])	
	[C,D] = fastequalspairs(YY,ZZ)
	return [D,B]	
	

def fastequalspairs(Y,Z):
	'''
	
	ARGUMENTS:
		LL1 = numpy array of paths
		LL2 = sorted numpy array of paths
	
	RETURNS:
	
	[A,B] where A and B are numpy arrays of indices in LL1 such that:
		
		LL2[A[i]:B[i]] = LL1[i].   
		
	A[i] = B[i] = 0 if LL1[i] not in LL2
	'''	
	
	T = Z.copy()
	R = (T[1:] != T[:-1]).nonzero()[0] 
	R = numpy.append(R,numpy.array([len(T)-1]))
	M = R[R.searchsorted(list(range(len(T))))]
	D = T.searchsorted(Y)
	T = numpy.append(T,numpy.array([0]))
	M = numpy.append(M,numpy.array([0]))
	W = (T[D] == Y) * D
	U = (T[D] == Y) * (M[D] + 1)

	return [W,U]


def ModContents(obj,Cond = None):
	'''
	Modified version of Contents function, avoiding recursive inspection
	of objects that satisfy condition Cond 
	
	ARGUMENTS:
	--obj = BeautifulSoup object
	--Cond = two-place boolean function with arugments (o1,o2) where 
		o1,o2 are meant to be BeautifulSoup objects
	
	If Cond = None this is the same as Contents 
	'''
	if 'contents' not in dir(obj):
		return str(obj)
	else:
		if Cond == None:
			return ''.join([ModContents(o) for o in obj.contents])
		else:
			return ''.join([ModContents(o,Cond) for o in obj.contents if not Cond(o,obj)])

def Contents(obj):
	'''
	Convenience function for working with BeautifulSoup contents 
	objects (move this somewhere else?) 
	
	ARGUMENT:
	--obj = BeautifulSoup object. 
	
	Given BeautifulSoup object 'obj', extract the "string contents" recursively. 
	'''
	if 'contents' not in dir(obj):
		return str(obj)
	else:
		return ''.join([Contents(newobj) for newobj in obj.contents])


def fastisin(Y,Z):	
	'''
	fast routine for determining indices of elements in numpy array 
	Y that appear in numpy array Z
	
	returns boolean array of those indices
	'''
	if len(Z) > 0:
		T = Z.copy()
		T.sort()
		D = T.searchsorted(Y)
		T = numpy.append(T,numpy.array([0]))
		W = (T[D] == Y)
		if isinstance(W,bool):
			return numpy.zeros((len(Y),),bool)
		else:
			return (T[D] == Y) 
	else:
		return numpy.zeros((len(Y),),bool)


def FastRecarrayEquals(Y,Z):
	'''
	fast routine for determining whether numpy record array Y 
	equals record array Z
	'''

	if Y.dtype.names != Z.dtype.names or len(Y) != len(Z):
		return False
	else:
		NewY = numpy.array([str(l) for l in Y])
		NewZ = numpy.array([str(l) for l in Z])
		NewZ.sort(); NewY.sort()
		return all(NewY == NewZ)

def FastRecarrayEqualsPairs(Y,Z):

	NewY = numpy.array([str(l) for l in Y])
	NewZ = numpy.array([str(l) for l in Z])
	s = NewZ.argsort()  ; NewZ.sort()
	[A,B] = fastequalspairs(NewY,NewZ)
	return [A,B,s]
		
		
def IsDotPath(s,path=None):
	'''
	Determine whether s is possible valid dot path of a python module, and is 
	more accurante when the putative real (relative) file path is 
	given in path.   (If path argument is given this requires path to 
	be a DataEnvironment-relative path, starting with ../)
	'''
	step1 = set(s.replace('.','').lower()) <= set(['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','0','1','2','3','4','5','6','7','8','9','_'])
	if path == None:
		return step1
	else:
		return step1 and path.endswith('.py') and s == path[:-3].strip('./').replace('/','.')	


def FastRecarrayIsIn(Y,Z):
	'''
	Fast routine for determining which records in numpy record array
	Y appear in record array Z
	'''

	if Y.dtype.names != Z.dtype.names:
		return numpy.zeros((len(Y),),bool)
	else:
		NewY = numpy.array([str(l) for l in Y])
		NewZ = numpy.array([str(l) for l in Z])
		NewZ.sort()
		return fastisin(NewY,NewZ)
		
		
def FastRecarrayDifference(X,Y):
	'''
	fast routine for determining which records in numpy array X do
	not appear in numpy array Y
	'''
	
	if len(Y) > 0:
		Z = FastRecarrayIsIn(X,Y)
		return X[numpy.invert(Z)]
	else:
		return X

def fastarraymax(X,Y):
	'''
	fast way to achieve: 
	
	ARGUMENTS:
		X,Y numpy arrays of equal length
	RETURNS:
		Z where Z[i] = max(X[i],Y[i])
	'''
	Z = numpy.zeros((len(X),),int)
	A = X <= Y
	B = Y < X
	Z[A] = Y[A]
	Z[B] = X[B]
	return Z


def fastarraymin(X,Y):
	'''
	fast way to achieve: 
	
	ARGUMENTS:
		X,Y numpy arrays of equal length
	RETURNS:
		Z where Z[i] = min(X[i],Y[i])
	'''
	Z = numpy.zeros((len(X),),int)
	A = X <= Y
	B = Y < X
	Z[A] = X[A]
	Z[B] = Y[B]
	return Z
	

def SimpleStack(seq,UNIQIFY=False):
	'''
	Vertically stack sequences numpy record arrays. 
	Avoids some of the problems of numpy.v_stack
	'''
	newseq = [ss for ss in seq if len(ss) > 0]
	if len(newseq) > 1:
		seq = newseq
		names = seq[0].dtype.names
		formats = [max([a.dtype[att] for a in seq]).str for att in names]
		if UNIQIFY:
			X =  numpy.rec.fromarrays([ListUnion([a[att].tolist() for a in seq]) for att in names], names = names, formats = formats)
			[D,s] = FastRecarrayUniqify(X)
			return X[s][D]
		else:
			return numpy.rec.fromarrays([ListUnion([a[att].tolist() for a in seq]) for att in names], names = names, formats = formats)
	elif len(newseq) == 1:
		return newseq[0]
	else:
		return seq[0][0:0]
		
def SimpleStack1(seq,UNIQIFY=False):
	'''
	Vertically stack sequences numpy record arrays. 
	Avoids some of the problems of numpy.v_stack but is slower
	
	if UNIQIFY set to true, only retains unique records
	'''
	newseq = [ss for ss in seq if len(ss) > 0]
	if len(newseq) > 1:
		seq = newseq
		names = seq[0].dtype.names
		formats = [max([a.dtype[att] for a in seq]).str for att in names]
		if UNIQIFY:
			numpy.rec.fromrecords(uniqify(ListUnion([ar.tolist() for ar in newseq])), names = names, formats = formats)
		else:
			return numpy.rec.fromrecords(ListUnion([ar.tolist() for ar in newseq]), names = names, formats = formats)
	elif len(newseq) == 1:
		return newseq[0]
	else:
		return seq[0][0:0]		
	

		
def SimpleColumnStack(seq):
	'''
	Stack columns in sequences of numpy record arrays. 
	Avoids some of the problems of numpy.c_stack but is slower
	'''
	Columns = ListUnion([[a[l] for l in a.dtype.names] for a in seq])
	names = ListUnion([list(a.dtype.names) for a in seq])
	return numpy.rec.fromarrays(Columns,names=names)	
	
	
	
def RemoveColumns(recarray,ToRemove):
	'''
	Given numpy recarray and list of column names ToRemove,
	return recarray with columns whose names are not in ToRemove
	'''

	
	newdtype = numpy.dtype([x for x in recarray.dtype.descr if x[0] not in ToRemove])
	return numpy.rec.fromarrays([recarray[name] for name in recarray.dtype.names if name not in ToRemove],dtype = newdtype)
	
	
def MaximalCommonPath(PathList):
	'''
	Given list of paths, return common prefix.  Like 
	os.path.commonprefixbut proper for Data Environment purposes. 
	'''

	PathList = PathList[:]
	PathList = list(set(PathList))
	OriginalPathList = PathList
	if len(PathList) > 0:
		PathList = [p[:-1] if p[-1] ==  '/' else p for p in PathList]
		splitlist = [x.split('/') for x in PathList]
		minslash = min([len(x) for x in splitlist])
		done = False
		i = 0
		cpath = ''
		while not done and i <= minslash:
			L = ['/'.join(x[:i]) for x in splitlist]
			if len(set(L)) == 1:
				cpath = L[0]
				i = i+1
			else:
				done = True
		
		if (cpath in OriginalPathList):		# cpath is a file path
			if sum([p.startswith(cpath + '/') for p in OriginalPathList]) > 0:	 # cpath + '/' is also a directory
				cpath = '/'.join(cpath.split('/')[:-1])
		if cpath != '':
			if sum([p.startswith(cpath + '/') for p in OriginalPathList]) > 0:	 # cpath + '/' is a directory
				cpath += '/'
		return cpath
		
	else:
		return ''


def Backslash(Dir,Verbose=False):
	'''
	Adds '/' to end of a path (meant to make formatting of directory 
	Paths consistently have the slash)
	'''

	if Dir[-1] != '/':
		if Verbose:
			print("Warning: the directory name, ", Dir, ", was provided. The character '/' was appended to the end of the name.")
		return Dir + '/'
	else:
		return Dir


def MakeDirWithDummy(Dir):
	'''
		makes a directory with a empty file 'dummy' in it 
	'''
	Dir = Backslash(Dir)
	MakeDir(Dir)
	open_for_write(Dir + 'dummy')[0].write('')


def MakeDirWithInit(Dir):	
	'''
		makes a directory with an empty file '__init__.py' in it 
	'''

	Dir = Backslash(Dir)
	MakeDir(Dir)
	open_for_write(Dir + '__init__.py')[0].write('')
	

def GetTimeStampedArchiveName(toarchive):
	'''
	given path string, return corresponding name that it would have 
	in the Archive, with timestamp attached
	'''
	#this assumes path starts with '../'
	TS = TimeStamp()
	modifiedpath = toarchive[3:].replace('/','__')
	return 'Archive_' + TS + '_' + modifiedpath
	
	
def copy_to_archive(toarchive,depends_on=('../',),creates=('../Archive/',)):
	'''
		copy file or directory to archive with proper archive name
	'''
#archives the file or folder 'toarchive' to the archive with new name generated by timestamp

	ArchivedName = GetTimeStampedArchiveName(toarchive)

	if not PathExists('../Archive/'):
		print('Creating Archive ....')
		MakeDir('../Archive/')
		
	if PathExists(toarchive):
		strongcopy(toarchive,'../Archive/' + ArchivedName)
	else:
		print('ERROR: The path', toarchive, 'does not exist; nothing archived.')


def move_to_archive(toarchive,depends_on=('../',),creates=('../Archive/',)):
	'''
	move file or directory to archive with proper archive name
	'''
#archives the file or folder 'toarchive' to the archive with new name generated by timestamp
	
	ArchivedName = GetTimeStampedArchiveName(toarchive)
	
	if not PathExists('../Archive/'):
		print('Creating Archive ....')
		MakeDir('../Archive/')
		
	if PathExists(toarchive):
		Rename(toarchive,'../Archive/' + ArchivedName)
	else:
		print('ERROR: The path', toarchive, 'does not exist; nothing archived.')
		

def CompilerChecked(ToCheck):
	'''
	cleans a list of strings representing python regular expressions ToCheck, 
	returning only those that are not empty and properly compile
	'''
	X = []
	for L in ToCheck:
		LL = L[4:] if L.startswith('NOT ') else L
		try: 
			re.compile(LL)
		except:
			print("Error: the string, ", LL, " was found calling", funcname(), ". This string could not be compiled as a regular expression and will not be loaded.")
		else:
			X += [L if L != '' else '^$']
	return X		



def CheckInOutFormulae(ExpList,S):
	'''
	Given a list ExpList of Regular expression strings and "NOT '-prefixed
	regular expression strings, return list of all strings in list S that:
	-- match at least one of the expressions in ExpList that are _not_ prefixed by 'NOT '
	-- match none of the expressions in ExpList that _are_ prefixed by 'NOT '
	'''
	if isinstance(ExpList,str):
		ExpList = ExpList.split(',')
	InExpList = [exp for exp in ExpList if not exp.startswith('NOT ')]
	OutExpList = [exp[4:] for exp in ExpList if exp.startswith('NOT ')]
	F = lambda x,y,z : any([re.match(RegExp,z) != None for RegExp in x]) and not any([re.match(RegExp,z) != None for RegExp in y])
	return F(InExpList,OutExpList,S)


def AddInitsAbove(opfile):
	'''
	Given a python to a pytho module opfile, add an empty __init__.py file to the 
	directory containing opfile, if not such file exists.   
	Reset mod time of directory so it appears as if nothing as changed. 
	
	The intent of this is to allow python modules to be placed in directories
	in the Data Environment and then be accessed by package imports,
	without the user having to remember to put the '__init__.py' in the directory 
	that the module is in.   The timestamp of the containing directory is reset
	if the __init__.py is added to make sure that no stupid re-computations are
	done that make it appear as if things have changed when the havent. 
			
	'''
	DirList = opfile.split('/')[:-1]
	for ii in range(1,len(DirList)):
		DirName = '/'.join(DirList[:ii+1]) + '/'		
		oldatime = os.path.getatime(DirName)
		oldmtime = os.path.getmtime(DirName)
		if '__init__.py' not in listdir(DirName):
			F = open(DirName + '__init__.py','w')
			F.close()
		os.utime(DirName,(oldatime,oldmtime))
		
				
class multicaster():
	'''
	Class creating object that multicasts a string output stream to have both
	its original desired effect and also to print any output to a log file. 
	
	typical Usage:
		sys.stdout = multicaster(sys.__stdout__,'LogFile.txt')
	
	Then, whenever a 'print ' statement is made, output is directed both to
	original stdout as well as to the logfile "LogFile.txt"
	'''
	def __init__(self,filename,OldObject,New=False):
		'''
		ARGUMENTS:
			filename = name of file to write to 
			OldObject = original output stream to multicast
			NEW = boolean which overwrites log file if true; otherwise, 
			output of stream is _appended_ to 'filename'	
			
		'''
		self.file = filename
		self.old = OldObject
		
		if New:
			F = open(filename,'w')
			F.write('\n\n------------------------------------------------------------------------------------------------------------------------------------------------------\n')
			F.write('STARTING LOG: ' + time.strftime('%c %Z') + '\n')
			F.write('------------------------------------------------------------------------------------------------------------------------------------------------------\n\n')
			F.close()
			
	def __getattr__(self,name):
		'''
		This is intended to answer that whenever the stdout is asked to 
		do something other than write the  function is undisturbed.  If the
		stdout object were its own class (instead of it being merely a file 
		\object that is being used for the purpose of output), this would be 
		unneccesary, we'd simple subclass the stdout object. 
		'''
		if name != 'write':
			return self.old.__getattribute__(name)
					
	def write(self,s):
		F = open(self.file,'a')
		F.write(s)
		F.close()	
		return self.old.write(s)
		

def DictInvert(D):
	'''
		ARGUMENT:
			dictionary D
		OUTPUT:
		--dictionary whose keys are unique elements of values of D, and 
		whose values on key 'K' are lists of keys 'k' in D such that D[k] = K
	'''
	return dict([(v,set([j for j in list(D.keys()) if D[j] == v])) for v in set(D.values())])


def PathExists(ToCheck):
	'''
		convenient name for os function
		
		The reason it's done this way as opposed to merely setting 
			PathExists = os.path.exists
		in this module is that this will disturb the system i/o intercept because this module needs to be execfiled FIRST before system_io_override. 
		
	'''
	return os.path.exists(ToCheck)

def Rename(src,dest):
	'''
		convenient name for os function
		
		The reason it's done this way as opposed to merely setting 
			PathExists = os.path.exists
		in this module is that this will disturb the system i/o intercept because
		this module needs to be execfiled FIRST before system_io_override. 
		
	'''
	os.rename(src,dest)

def IsDir(ToCheck):
	'''
		convenient name for os function
		
		The reason it's done this way as opposed to merely setting 
			PathExists = os.path.exists
		in this module is that this will disturb the system i/o intercept because this
		module needs to be execfiled FIRST before system_io_override. 
		
	'''
	return os.path.isdir(ToCheck)
	
def IsFile(ToCheck):
	'''
		convenient name for os function
		
		The reason it's done this way as opposed to merely setting 
			PathExists = os.path.exists
		in this module is that this will disturb the system i/o intercept 
		because this module needs to be execfiled FIRST before system_io_override. 
		
	'''
	return os.path.isfile(ToCheck)
	
def FindAtime(ToAssay):
	'''
		convenient name for os function
		
		The reason it's done this way as opposed to merely setting 
			PathExists = os.path.exists
		in this module is that this will disturb the system i/o intercept
		because this module needs to be execfiled FIRST before system_io_override. 
		
	'''
	return os.path.getatime(ToAssay)
	
def listdir(ToList):
	'''
		convenient name for os function
		
		The reason it's done this way as opposed to merely setting 
			PathExists = os.path.exists
		in this module is that this will disturb the system i/o intercept
		because this module needs to be execfiled FIRST before system_io_override. 
		
	'''
	return os.listdir(ToList)
	
ListAnd = all
ListOr = any

