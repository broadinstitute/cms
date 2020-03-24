from numpy import *
from Classes.DotData import DotData

def xpopMerge(xpop_EastAsian_WestAfrican,xpop_European_EastAsian,xpop_European_WestAfrican):

	"""Reads Sweep xpop output for each pop pair and merges into a DotData with AllEHH logratio Deviation
	for both sides for each pop pair.
	"""

	dashFixer = lambda v: v if v != '-' else nan

	EA_WA = DotData(SVPath=xpop_EastAsian_WestAfrican,ToLoad=['Chrom','SNP pos (bases)','SNP pos (cM)','L AllEHH logratio Deviation','R AllEHH logratio Deviation'], SVValueFixer = dashFixer)
	EU_EA = DotData(SVPath=xpop_European_EastAsian,ToLoad=['Chrom','SNP pos (bases)','SNP pos (cM)','L AllEHH logratio Deviation','R AllEHH logratio Deviation'], SVValueFixer = dashFixer)
	EU_WA = DotData(SVPath=xpop_European_WestAfrican,ToLoad=['Chrom','SNP pos (bases)','SNP pos (cM)','L AllEHH logratio Deviation','R AllEHH logratio Deviation'], SVValueFixer = dashFixer)

	keyCol1 = array(EA_WA['SNP pos (bases)']) + array(EA_WA['Chrom'])*1000000
	keyCol2 = array(EU_EA['SNP pos (bases)']) + array(EU_EA['Chrom'])*1000000
	keyCol3 = array(EU_WA['SNP pos (bases)']) + array(EU_WA['Chrom'])*1000000

	keyData1 = DotData(Columns = [keyCol1,],names = ['key',])
	keyData2 = DotData(Columns = [keyCol2,],names = ['key',])
	keyData3 = DotData(Columns = [keyCol3,],names = ['key',])

	EA_WA = EA_WA.hstack(keyData1)
	EU_EA = EU_EA.hstack(keyData2)
	EU_WA = EU_WA.hstack(keyData3)

	blank = (nan,)*len(EA_WA.dtype.names)

	return DotData.mergeOnKeyCols((EA_WA,EU_EA,EU_WA),('key',)*3,(blank,)*3,
				      (' EastAsian_WestAfrican',' EastAsian_European',' European_WestAfrican'))
