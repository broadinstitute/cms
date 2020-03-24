
from Classes.DotData import DotData
from Operations.MiscUtil import dbg, clamp
from Operations.Shari_Operations.localize.RunningStat import RunningStat 
from numpy import *
import numpy as np
import logging, itertools

## Normalizes an array by mean and std dev
def normalize(rawVals, ind = itertools.repeat( True ) ):
    notNan = []
    for val, useVal in zip( rawVals, ind ):
        if useVal and not isnan(val) and not isneginf(val):
            notNan.append(val)
#    print len(notNan)
    theMean = mean(notNan)
    theStd = std(notNan)
    
    normVals = []

    for val, useVal in zip( rawVals, ind ):
        normVals.append( (val - theMean)/theStd if useVal else nan )

    return normVals, theMean, theStd

## Bins values by frequencies and then normalizes within bins        
def normalizeByFreq(rawVals,ancfreq):

    Frequency = arange(0.05, 1.05, 0.05)
    #print Frequency

    ind = invert(isnan(rawVals))
    copyVals = rawVals

    rawVals = rawVals[ind]
    #print 'rawVals are', '\n'.join( map( str,  frozenset( map( type, rawVals ) ) ) )
    #print sum([isnan(val) for val in rawVals])
    der_freq = 1 - ancfreq[ind]
    temp = [ ]
    expectation = []
    normVal = []

    StdDev = std(rawVals)

    #print 'StdDev=', StdDev

    for i in range(len(Frequency)):
        for j in range(len(rawVals)):
            if ((Frequency[i] - der_freq[j]) < .05 and Frequency[i] - der_freq[j] > 0):
                temp.append(rawVals[j])
            
                
        theMean = mean(temp) if temp else 0
#        print ' theMean=', theMean
        expectation.append(theMean)
        temp = []

    # Bookkeeping
    found_freq = 0
    rawValCount = normValCount = 0
    normVal = zeros(len(copyVals))
    der_freq = 1 - ancfreq

    for j in range(len(copyVals)):
        found_freq = 0
        if ind[j]:
            for i in range(len(Frequency)):
                if ((Frequency[i] - der_freq[j]) < .05 and (Frequency[i] - der_freq[j]) >= 0): 
                    normVal[j] = (copyVals[j] - expectation[i])/StdDev
                    found_freq = 1
                    normValCount += 1
        else:
            normVal[j] = nan

        rawValCount +=1
#    print 'rawVal: ' + str(rawValCount) + ' normVal: ' + str(normValCount)
#    print 'zeros and ones: ' + str(sum(der_freq == 1) + sum(der_freq == 0))
    
    return normVal


## Calculates the "center of mass" of ihs scores and return its position
def CoM(iHS,gdPos):
    ihs = array(iHS)
    pos = array(gdPos)

    ind = ihs.argsort()
    sorted_pos = pos[ind]
    cm = mean(sorted_pos[0:19])
    return cm


## Calculates the "center of mass" of ihs scores and returns an array of the dist of each point from the center of mass
def CoM_dist(iHS,gdPos):
    ihs = array(iHS)
    pos = array(gdPos)

    ind = ihs.argsort()
    sorted_pos = pos[ind]
    cm = mean(sorted_pos[0:19])
    dist = abs(pos - cm)

    return dist

## Calculates the "center of mass" of ihs scores and returns an array of the dist of each point from the center of mass
def CoM_pos(iHS,gdPos,bpPos,selpos=500000):
    ihs = array(iHS)
    pos = array(gdPos)
    causal = where(bpPos == selpos)[0][0]
    causalPos = pos[causal]

    ind = ihs.argsort()
    sorted_pos = pos[ind]
    cm = mean(sorted_pos[0:19])

    return abs(cm - causalPos)

def scoreVal(stat,pos,nbin,start,end,distr_hit,distr_miss,selpos):
    """Bins vals by intervals in stat val, scores causal and non-causal.
    So, given a statistic for each SNP (e.g. a component of the normalized composite likelihood score),
    build a histogram giving for each bin the number of causal and non-causal
    SNPs with their statistic in that bin.

    Input params:

       stat - for each SNP, the statistic value (e.g. the normalized composite likelihood score)
       pos - position of this SNP on the chromosome, in bases (only so that we know
         whether it is the selected SNP, by comparing its position to selpos) 
       nbin, start, end - define the bins of this statistic
       selpos - position of the causal SNP on the chromosome, in bases

    Output params:

       distr_hit - for each bin, number of causal SNPs with the statistic in that bin
       distr_miss - for each bin, number of non-causal SNPs with the statistic in that bin
    
    """
    bin_size = float(end - start) / nbin

    for i in range(len(stat)):
        
        if not isfinite(stat[i]): continue

        bin = int((stat[i] - start)/bin_size)

        if bin >= nbin: bin = nbin - 1
        if bin < 0: bin = 0

        ( distr_hit if pos[i] == selpos else distr_miss )[bin] += 1

    logging.info( 'hits: %d miss: %d' % ( sum( distr_hit ), sum( distr_miss ) ) )

    return 

def likeVal(stat,hitsLikes,missLikes,stat_start,stat_end,nbin,lik,ind=None, defaultNumSnps = 10000,
            useDefaultNumSnps = True):
    """For eah SNP, calculate loglikelihood of that SNP being causal based on one statistic,
    given the SNP's bin for that statistic.

    Input params:

       stat - for each SNP, its statistic value (typically, normalized).
       hitsLikes/missLikes - for each bin, what fraction of causal/non-causal SNPs has
          the statistic in that bin?
       stat_start, stat_end, nbin - bin boundaries and number of bins for the statistic

       ind - optionally indicates which snps should be skipped.   if given,
          for snps for which ind[i] is false no likelihood value is computed
          (and lik[i] is set to nan).

    Output params:

       lik - for each SNP, log of the probability that the SNP is causal given what bin this statistic
          for this SNP falls into.  the lik array is allocated by the caller, and filled-in
          by this function.

    Use by: complike()      
       
    """

    # if the (optional) indicator of which SNPs to compute likelihood for is not given,
    # compute likelihood for all SNPs.
    if ind == None  or   len(ind) < 1:
        ind = ones(len(stat))

    bin_size = float(stat_end - stat_start) / nbin
    # var: p_causal - probability that a SNP chosen randomly from the replica is causal.
    #    within a replica (region), exactly one SNP is causal, so this is 1 / ( total number of SNPs ).
    #    For Bayes' rule, this is the prior probability that a SNP is causal, before we look at the value
    #    of the statistic.
    #p_causal = 1. / len(stat)

    numSnps = defaultNumSnps if useDefaultNumSnps else len( stat )
    dbg( 'numSnps' )
    p_causal = 1./numSnps if numSnps > 0 else np.nan
    p_noncausal = 1 - p_causal
    dbg( 'p_causal p_noncausal' )

    for i in range(len(stat)):

        if not ind[i]:
            p_causal_bin = 0

        elif not isfinite(stat[i]):
            p_causal_bin = 1. / nbin
        
        else:


            bin = int((stat[i] - stat_start)/bin_size)

            if bin >= nbin:
                bin = nbin - 1
            if bin < 0:
                bin = 0

                

            # P( causal | stat_in_this_bin ) = ( P( stat_in_this_bin | causal ) * P( causal ) ) / P( stat_in_this_bin )

            # here, hitLikes[bin] gives P( stat_in_this_bin )    
        
            num = hitsLikes[bin]*p_causal
            # var: denom - probability that a SNP, whether causal or not, has this statistic in this bin
            denom = hitsLikes[bin] * p_causal + missLikes[bin]*p_noncausal
            assert denom != 0

            #print num,denom

            p_causal_bin = num / denom

            if np.abs( stat[i] - 4.6387211 ) < 1e-5:
                dbg( '"here!" bin stat[i] hitsLikes[bin] p_causal num denom p_causal_bin' )

        # endif: whether we have the statistic value for this snp,
        #    and whether we are asked to compute the SNP's likelihood of being causal
        
        likscore = log(p_causal_bin) if p_causal_bin > 0.0 else nan
        lik[i] = likscore

    # endloop: for each SNP in the replica

# end: likeVal()
        
def calcRank(stat,nbin,distr_hit,distr_miss,pos,selpos):
    """Bins by rank instead of value: of the causal SNPs, how many
    are in the top 5%, 10%, 15% of the values of this statistic?  same for
    non-causal SNPs.

    Input params:

       stat - for each SNP, the value of a statistic
         (e.g. the composite likelihood score)
       nbin - number of bins

       pos - for each SNP, its position on the chromosome ( in bases );
          only needed so we can tell which SNP is the causal one.
       selpos - position of causal SNP on the chromosome, in bases.

    Output params:
       
      distr_hit/distr_miss - for each bin of ranks, how many causal/non-causal SNPs fall in that bin?
        NOTE: The counts are _added_ to the existing counts in these arrays.

     Returns:

       the bin number of the bin into which the causal SNP fell. 

    """
    bin_scale = 1. / nbin
    causalBin = -1

    sortInd = stat.argsort()

    for i in range(len(stat)):
        bin = int(float(i)/len(stat)/bin_scale)
        
        if bin == nbin:
            bin -= 1
        assert bin < nbin
        if pos[sortInd[i]] == selpos:
            distr_hit[bin] += 1
            causalBin = bin
        else:
            distr_miss[bin] += 1

    return causalBin

# end: calcRank()


# Returns top scoring SNPs for given stat
def highScorers(stat,nbin,ntopbins,pos):
    bin_scale = 1. / nbin
    
    highScorers = []
    CL = []
    
    sortInd = stat.argsort()

    for i in range(len(stat)):
        bin = int(float(i)/len(stat)/bin_scale)
        if bin > ntopbins:
            highScorers.append(pos[sortInd[i]])
            print(str(pos[sortInd[i]]) + '\t' + str( stat[sortInd[i]] ))
            CL.append(stat[i])
    
    output = array([highScorers,CL])
    return output

def calcCDF(ranksDotData, hitRanksName = 'hitsRanks'):
    hitRanks = array(ranksDotData[hitRanksName])
    totalHits = sum(hitRanks)
    print(totalHits)
    fractions = hitRanks / float(totalHits)
    print(fractions)

    return cumsum(fractions)

def calcFreqCDF(ranksDotData, hitRanksName = 'hitsRanks'):
    cdfs = {}
    for freq in ['high','low']:
        hitRanks = array(ranksDotData['hitsRanks_' + freq])
        totalHits = sum(hitRanks)
        print(totalHits)
        fractions = hitRanks / float(totalHits)
        print(fractions)

        cdfs[freq] =  cumsum(fractions)
    
    return cdfs

def nanMax(a, Axis = 1):
    ind = invert(isnan(a))
    nonNan = a[ind]
    if nonNan.any():
        return np.max(nonNan, axis = Axis)
    else:
        return nan



def meanStdDevOfTsvColumn( tsvFileName, columnName ):
    """Compute the mean and stddev of a column in a .tsv file, without loading
    the entire file into memory.  Useful when dealing with very large files.

    Returns: the tuple (mean, stddev).
    """

    stats = RunningStat()
    
    with open( tsvFileName ) as f:
        header = f.readline().strip().split()
        colNum = header.index( columnName )

        lineNum = 0
        for line in f:
            vals = line.strip().split()
            assert len( vals ) == len( header )
            v = float( vals[ colNum ] )
            if not isnan( v ) and not isneginf( v ):
                stats.add( v )
            lineNum += 1
            if ( lineNum % 1000000 ) == 0: logging.info( 'line %d' % lineNum )

    return stats.mean(), stats.std()


def meanStdDevOfTsvColumn2( tsvFileName, columnName ):
    """Compute the mean and stddev of a column in a .tsv file, without loading
    the entire file into memory.  Useful when dealing with very large files.

    Returns: the tuple (mean, stddev).
    """

    allVals = [ ]
    
    with open( tsvFileName ) as f:
        header = f.readline().strip().split()
        colNum = header.index( columnName )
        print('colnum=', colNum)

        lineNum = 0
        for line in f:
            vals = line.strip().split()
            assert len( vals ) == len( header )
            v = longdouble( float( vals[ colNum ] ) )
            if not isnan( v ) and not isneginf( v ):
                allVals.append( v )
            lineNum += 1
            #if lineNum > 100: break
            if ( lineNum % 100000 ) == 0: logging.info( 'line %d' % lineNum )

    #print 'allVals=', allVals           
            
    return mean( allVals ), std( allVals )

def likeRatio(stat,hitsLikes,missLikes,stat_start,stat_end,nbin,lik,ind=None):
    """For eah SNP, calculate loglikelihood of that SNP being causal based on one statistic,
    given the SNP's bin for that statistic.

    Input params:

       stat - for each SNP, its statistic value (typically, normalized).
       hitsLikes/missLikes - for each bin, what fraction of causal/non-causal SNPs has
          the statistic in that bin?
       stat_start, stat_end, nbin - bin boundaries and number of bins for the statistic

       ind - optionally indicates which snps should be skipped.   if given,
          for snps for which ind[i] is false no likelihood value is computed
          (and lik[i] is left unchanged).

    Output params:

       lik - for each SNP, log of the probability that the SNP is causal given what bin this statistic
          for this SNP falls into.  the lik array is allocated by the caller, and filled-in
          by this function.

    Use by: complike()      
       
    """

    # if the (optional) indicator of which SNPs to compute likelihood for is not given,
    # compute likelihood for all SNPs.
    if ind == None  or   len(ind) < 1:
        ind = ones(len(stat))

    bin_size = float(stat_end - stat_start) / nbin
    
    indNaN = hitsLikes != 1e-10
    missingVal = log(min(hitsLikes[indNaN]/missLikes[indNaN]))
    print(missingVal)

    for i in range(len(stat)):

        if not ind[i]:
            CLR = 0

        elif not isfinite(stat[i]):
            CLR = 0
        
        else:

            bin = int((stat[i] - stat_start)/bin_size)

            if bin >= nbin:
                bin = nbin - 1
            if bin < 0:
                bin = 0

            # P( causal | stat_in_this_bin ) = ( P( stat_in_this_bin | causal ) * P( causal ) ) / P( stat_in_this_bin )

            # here, hitLikes[bin] gives P( stat_in_this_bin )    
        
            num = hitsLikes[bin]
            # var: denom - probability that a SNP, whether causal or not, has this statistic in this bin
            denom = missLikes[bin]
            assert denom != 0

            #print num,denom

            CLR = num / denom

        # endif: whether we have the statistic value for this snp,
        #    and whether we are asked to compute the SNP's likelihood of being causal
        if CLR != 0.0:
            likratio = log( CLR ) if num != 1e-10 else missingVal
        else:
            likratio = nan
        lik[i] = likratio

    # endloop: for each SNP in the replica

# end: likeRatio()

            
def likeValBin(stat,hitsLikes,missLikes,stat_start,stat_end,nbin,lik,ind=None):
    """For eah SNP, calculate loglikelihood of that SNP being causal based on one statistic,
    given the SNP's bin for that statistic.

    Input params:

       stat - for each SNP, its statistic value (typically, normalized).
       hitsLikes/missLikes - for each bin, what fraction of causal/non-causal SNPs has
          the statistic in that bin?
       stat_start, stat_end, nbin - bin boundaries and number of bins for the statistic

       ind - optionally indicates which snps should be skipped.   if given,
          for snps for which ind[i] is false no likelihood value is computed
          (and lik[i] is left unchanged).

    Output params:

       lik - for each SNP, log of the probability that the SNP is causal given what bin this statistic
          for this SNP falls into.  the lik array is allocated by the caller, and filled-in
          by this function.

    Use by: complike()      
       
    """

    # if the (optional) indicator of which SNPs to compute likelihood for is not given,
    # compute likelihood for all SNPs.
    if ind == None  or   len(ind) < 1:
        ind = ones(len(stat))

    bin_size = float(stat_end - stat_start) / nbin
    
    indNaN = hitsLikes != 1e-10
    missingVal = log(min(hitsLikes[indNaN]/missLikes[indNaN]))
    print(missingVal)

    for i in range(len(stat)):

        if not ind[i]:
            bin = nan

        elif not isfinite(stat[i]):
            bin = nan
        
        else:

            bin = int((stat[i] - stat_start)/bin_size)

            if bin >= nbin:
                bin = nbin - 1
            if bin < 0:
                bin = 0

        lik[i] = bin


    # endloop: for each SNP in the replica

# end: likeVal()
            
        
        
    

            
        
        
    
