from numpy import *
from Operations.Shari_Operations.localize.subs import normalize
import os
from Classes.DotData import DotData

def clean_hapmap2_region(Z, f = None, cleanIHS=True, cleanDer=True, cleanAnc=True, keepCausal = False,
                         returnInd = False):
    
    if f: f.write(str(len(Z)) + '\t')
    if cleanIHS:
        ind1 =  invert(isnan(Z.StdDiff))
        print 'SNPs with iHS scores: ',sum(ind1)
    else:
        ind1 = ones(len(Z))
    print 'ind1: ', sum(ind1)
    if f: f.write(str(sum(ind1)) + '\t')
    if cleanDer:
        ind2 = Z.derFreq > .2
    else:
        ind2 = ones(len(Z))
    print 'Derived > .2: ', sum(ind2)
    if f: f.write(str(sum(ind2)) + '\t')
    if cleanAnc:
        ind3 = Z.meanAnc > .4
    else:
        ind3 = ones(len(Z))
    print 'Ancestral > .4: ', sum(ind3)
    if f: f.write(str(sum(ind3)) + '\t')
    ind4 = zeros(len(Z))
    if keepCausal:
        ind4 = Z.Pos == 500000
    ind = all([ind1,ind2,ind3],axis=0)
    ind = any([ind,ind4],axis=0)
    print 'SNPs to keep: ', sum(ind)
    if f: f.write(str(sum(ind)) + '\t')

    if returnInd: return ind
    
    fZ = Z[ind]
    
    normlike = normalize(fZ.complike)[0]
    lik = exp(fZ.complike)
    Z = fZ.hstack(DotData(Columns = [normlike,lik],names = ['normLike','lik']))
    #Z=Z[['complike','derFreq','gdPos','iHS','max_xpop','meanAnc','meanFst','StdDiff','freqDiff','Chrom','Pos','normLike','lik']]

    return Z
    
