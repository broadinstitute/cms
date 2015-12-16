#!/usr/bin/python

import bisect
import util.file

class RecomMap(object):
    ''' Load in a recombination map which is a tab text file with four columns:
        chrom, pos (bp), recom rate (cM/Mb), and map pos (cM).  It must contain
        a header line but the contents of the header line are ignored.
        Provide methods to find the genetic map position and recombination rate
        of any physical coordinate in the genome, using linear interpolation to
        determine values for interior, unspecified positions, and truncating
        to the edge values for exterior, unspecified positions.
    '''
    def __init__(self, map_file, delim='\t'):
        self.header = None
        self.chr_to_gmap = {}
        self.chr_to_rate = {}
        with util.file.open_or_gzopen(map_file, 'rt') as inf:
            last_c = None
            last_p = 0
            last_mp = 0.0
            for line in inf:
                row = line.rstrip('\n').split(delim)
                assert len(row)==4, "there must be four columns, found %d.  line: %s" % (len(row), line.rstrip('\r\n'))
                if self.header is None:
                    self.header = row
                else:
                    c,p,rate,map_p = (row[0], int(row[1]), float(row[2]), float(row[3]))
                    assert 0.0 <= map_p #<= 100.01, "map position out of bounds, %s cM (%s:%d)" %(map_p, c, p)
                    # map_p = min(100.0, map_p) # commented out since the genetic distance in the map can be >100
                    assert c and p>0, "chr not specified or position non-positive %s:%d" % (c,p)
                    if last_c==c:
                        assert last_p<p, "non-monotonic increase in position at %s:%d" % (c,p)
                        assert last_mp<=map_p, "non-monotonic increase in map position at %s:%d" % (c,p)
                    else:
                        last_c = c
                        assert c not in self.chr_to_gmap and c not in self.chr_to_rate, "error: chromosome %s stopped, then started again later" % c
                        self.chr_to_gmap[c] = []
                        self.chr_to_rate[c] = []
                    last_p = p
                    last_mp = map_p
                    self.chr_to_gmap[c].append((p, map_p))
                    self.chr_to_rate[c].append((p, rate))
    def _interpolate(self, p, map_list, truncate=True):
        if p <= map_list[0][0]:
            # truncate to first seen value
            if truncate:
                val = map_list[0][1]
            else:
                val = None
        elif p >= map_list[-1][0]:
            # truncate to last seen value
            if truncate:
                val = map_list[-1][1]
            else:
                val = None
        else:
            # find i: the leftmost index where map_list[i][0]>=p  (using an efficient binary search)
            i = bisect.bisect_left(map_list, (p,0))
            if map_list[i][0]==p or map_list[i-1][1]==map_list[i][1]:
                # return exact value
                val = map_list[i][1]
            else:
                # interpolate
                l,r = map_list[i-1:i+1]
                val = l[1] + (r[1]-l[1]) * (p-l[0]) / (r[0]-l[0])
        return val
    def physToRate(self, c, pos, truncate=True):
        if c not in self.chr_to_rate:
            return None
        return self._interpolate(pos, self.chr_to_rate[c], truncate=truncate)
    def physToMap(self, c, pos, truncate=True):
        if c not in self.chr_to_gmap:
            return None
        return self._interpolate(pos, self.chr_to_gmap[c], truncate=truncate)
    def featuresBetween(self, chrom, mpos1, mpos2, features, fully=False):
        ''' Take an iterator of features and emit them as a filtered iterator
            where only the features that lie between two map positions are passed.
            fully=False allows for any partial overlap, fully=True requires that
            features are fully contained within these coordinates.
        '''
        assert mpos1 <= mpos2
        for c,start,stop in features:
            if c==chrom:
                start_cm = self.physToMap(c, start)
                stop_cm  = self.physToMap(c, stop)
                if fully:
                    if start_cm >= mpos1 and stop_cm <= mpos2:
                        yield (c,start,stop)
                else:
                    if stop_cm >= mpos1 and start_cm <= mpos2:
                        yield (c,start,stop)
