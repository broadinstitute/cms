#!/usr/bin/python

# built-ins
import os, re
import time

try:
    from itertools import izip as zip
except ImportError: # py3 zip is izip
    pass

try: 
    from itertools import count, zip_longest
except ImportError:
    from itertools import count, izip_longest as zip_longest
from collections import deque

# external
import pysam
from boltons.iterutils import chunked_iter

class VCFReader(object):

    ancestralAlleleRegex = re.compile(".*AA=(?P<allele>[^;,|\s]).*$")
    #ancestralAlleleRegex = re.compile(".*AA=(?P<allele>[^;,|\s])[,|\s].*")
    variantTypeRegex     = re.compile(".*VT=(?P<vt>[^;,|\s]+)[,|\s]*[;|]")
    contigLengthRegex    = re.compile('##contig=<ID=(?P<contig_id>[^,]*),.*length=(?P<length>\d*).*>')

    def __init__(self, file_path, parser=pysam.asVCF()):
        self.vcf_file_path = file_path
        self.tabix_file    = pysam.TabixFile(file_path, parser=parser)
        self.sample_names  = self.sample_names()
        self.clens         = self.contig_lengths()
        self.indexDelta    = -1 if tuple(map(int, pysam.__version__.split('.'))) > (0,5) else 0

    def contig_lengths(self):
        clens = []
        for line in self.tabix_file.header:
            if line.startswith('##contig=<ID=') and line.endswith('>'):
                matches = self.contigLengthRegex.match(line)
                c = matches.group("contig_id")
                clen = int(matches.group("length"))
                clens.append((c,clen))
        return dict(clens)

    @staticmethod
    def count_iter_items(iterable):
        """
        Consume an iterable not reading it into memory; return the number of items.
        """
        counter = count()
        deque(zip(iterable, counter), maxlen=0)  # (consume at C speed)
        return next(counter)

    @staticmethod
    def vcf_record_info(record):
        '''
            Given a tabix record iterator, this returns a dict() of the values 
            described in the INFO column, where "key=val" becomes key:val, and 'val' becomes val:val
        '''
        recordVariantDetails = record[7].split(";")
        recordVariantDetailsDict = {}

        for item in recordVariantDetails:
            try:
                key,val = item.split("=")
                recordVariantDetailsDict[key] = val
            except ValueError:
                recordVariantDetailsDict[item] = item

        return recordVariantDetailsDict

    def sample_header_line(self):
        '''
            Returns a string of the VCF file line containing the CHROM line.
        '''

        elem = None
        headerIter = self.tabix_file.header
        while True:
            try:
                elem = headerIter.next()
                if str(elem).find("CHROM", 1) > 0:
                    break
            except StopIteration:
                raise ValueError("Header not found in VCF file: {}".format(vcfFilePath))
        return str(elem)

    def sample_names(self):
        '''
            Returns a list of the sample names described in the VCF file
        '''
        samples = None
        headerRowString = self.sample_header_line()
        headerRowList = headerRowString.rstrip('\r\n').split('\t')
        samples = headerRowList[9:]

        return samples

    def records(self, chromosome_num, start_pos_bp, end_pos_bp, parser=pysam.asVCF()):
        '''
            Returns an iterator for the file records.
        '''

        start_pos_bp = 1 if start_pos_bp == None else start_pos_bp
        end_pos_bp = self.clens[str(chromosome_num)] if end_pos_bp == None else end_pos_bp

        # subtract one since pysam uses incorrect 0-indexed positions
        records = self.tabix_file.fetch( str(chromosome_num), start_pos_bp+self.indexDelta, end_pos_bp+self.indexDelta, parser)

        return records

    @classmethod
    def parse_ancestral(cls, infoString):
        match = cls.ancestralAlleleRegex.match(infoString)
        if match:
            ancestralAllele = match.group("allele")
            return ancestralAllele
        return ""

    @classmethod
    def variant_is_type(cls, infoString, variantType="SNP"):
        match = cls.variantTypeRegex.match(infoString)
        if match:
            return match.group("vt") == variantType
            #return match.group("vt").lower() == variantType.lower()
        return False

    @classmethod
    def allele_is_snp(cls, record):
        VALID_BASES = ["A","C","G","T","N","a","c","g","t","n"]
        if (len(record.ref) == 1 and len(record.alt) == 1) or ( all(variant in VALID_BASES for variant in record.ref.split(",")) and 
             all(variant in VALID_BASES for variant in record.alt.split(",")) ):
            return True
        return False
        

# def test():
#     rootDir = os.path.dirname(os.path.realpath(__file__))
#     # actual number of records in chromosome 1: 6468094
#     fileName = "/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

#     filePath = rootDir + fileName
    
#     tb = pysam.TabixFile(filePath)
#     print "size of header: {}".format(count_iter_items(tb.header))
    
#     print getVCFSampleNames(filePath)

#     colHeaderSeen = False
#     #for item in tb.header:
#     #    print item
#     #records = tb.fetch("1", 0, 250240543, parser=pysam.asTuple())
#     records = tb.fetch("1", 0, 1066000, parser=pysam.asTuple())

#     #print sum(1 for _ in records)

#     chunkIter = chunked_iter(records, 1000)

#     #print count_iter_items(chunkIter)
#     #print chunkIter.next()

    

#     numRecords = 0
#     currentIter = None
#     try:
#         while chunkIter.next():
#             currentIter = chunkIter
#             numRecords += 1
#     except:
#         print count_iter_items(currentIter)
#         print "end reached, iterator stopped, numRecords: {}".format(numRecords)
    


    # print "end reached, numRecords: {}".format(numRecords)
    # totalSnps = 0
    # for record in records:
    #     totalSnps += 1
    #     if totalSnps > 1:
    #         break

    #     print record

    #     #print parseStats(record)

    # print totalSnps

# class TabixReader(pysam.Tabixfile):
#     ''' A wrapper around pysam.Tabixfile that provides a context and
#         allows us to query using 1-based coordinates.
        
#         We should request upstream to the Pysam people to implement
#         methods __reversed__() and __len__() in pysam.TabixIterator.
#         __getitem__ could be a bonus, but prob unnecessary.
#     '''
#     def __init__(self, inFile, parser=pysam.asTuple()):
#         # because of odd Cython weirdness, we don't actually want to call super.__init__ here..
#         #super(TabixReader, self).__init__(inFile=inFile, parser=parser, mode=mode)
#         self.parser = parser
#     def __enter__(self):
#         return self
#     def __exit__(self, exc_type, exc_val, exc_tb):
#         self.close()
#         return 0
#     def close(self):
#         if hasattr(super(TabixReader, self), 'close'):
#             super(TabixReader, self).close()
#         else:
#             log.warn("warning: pysam-0.5 lacks a pysam.Tabixfile.close() method.  The input file may not be closed.")
#             # NOTE: this does not exist in pysam-0.5
#             # I'm not sure if that means it actually doesn't close the file ever.
#     def chroms(self):
#         return self.contigs
#     def get(self, chrom=None, start=None, stop=None, region=None):
#         if start!=None:
#             start -= 1
#         return self.fetch(reference=chrom, start=start, end=stop,
#             region=region, parser=self.parser)


# if tuple(map(int, pysam.__version__.split('.'))) > (0,5):
#     # new versions of pysam return a zero-based position here
#     def get_pos_from_vcf_record(vcfrec):
#         return vcfrec.pos + 1
# else:
#     # old versions of pysam return the correct one-based position (same as the VCF file) here
#     def get_pos_from_vcf_record(vcfrec):
#         return vcfrec.pos


# def bytes_to_string(o):
#     if type(o) == bytes:
#         o = o.decode('utf-8')
#     return o

# class VcfReader(TabixReader):
#     ''' Same as TabixReader with a few more perks for VCF files:
#         - emit results parsed as pysam VCF rows
#         - provide self.chrlens(), a dict mapping chrom name to chrom length
#         - provide self.samples(), a list of sample names in order of appearance
#         - provide get_range(c,start,stop) and get_snp_genos(c,pos)
#     '''

#     ancestralAlleleRegex = re.compile(".*AA=(?P<allele>[^;,|\s])[,|\s].*")
#     variantTypeRegex     = re.compile(".*VT=(?P<vt>[^;,|\s]+)[,|\s]*")

#     def __init__(self, inFile, ploidy=1, parser=pysam.asVCF()):
#         super(VcfReader, self).__init__(inFile=inFile, parser=parser)
#         assert ploidy in (1,2)
#         self.ploidy = ploidy
#         self.clens = []
#         self.sample_names = None
#         for line in self.header:
#             line = bytes_to_string(line)
#             if line.startswith('##contig=<ID=') and line.endswith('>'):
#                 exp = re.compile('##contig=<ID=(?P<contig_id>[^,]*),.*length=(?P<length>\d*).*>')
#                 matches = exp.match(line)
#                 c = matches.group("contig_id")
#                 clen = int(matches.group("length"))
#                 self.clens.append((c,clen))
#             elif line.startswith('#CHROM'):
#                 row = line.split('\t')
#                 assert len(row) > 8, "No samples are present in the VCF file: {}".format(inFile)
#                 self.sample_names = row[9:]
#         self.clens = dict(self.clens)
#         assert self.sample_names
#     def samples(self):
#         return self.sample_names
#     def chrlens(self):
#         return self.clens
#     def get_positions(self, c=None, start=None, stop=None, region=None):
#         for snp in self.get(c,start,stop,region):
#             yield (bytes_to_string(snp.contig), get_pos_from_vcf_record(snp), get_pos_from_vcf_record(snp)+len(snp.ref)-1)
    
#     @classmethod
#     def parse_ancestral(cls, infoString):
#         match = cls.ancestralAlleleRegex.match(infoString)
#         if match:
#             ancestralAllele = match.group("allele")
#             return ancestralAllele
#         return ""

#     @classmethod
#     def variant_is_type(cls, infoString, variantType="SNP"):
#         match = cls.variantTypeRegex.match(infoString)
#         if match:
#             return match.group("vt") == variantType
#             #return match.group("vt").lower() == variantType.lower()
#         return False

#     def get_range(self, c=None, start=None, stop=None, region=None, as_strings=True, more=False, omit_genotypes=False):
#         ''' Read a VCF file (optionally just a piece of it) and return contents
#             as an iterator with parsed contents.  Each row is returned as
#             a 4-tuple:
#                 1: chr
#                 2: pos (int)
#                 3: list of allele strings (in order)
#                 4: list of genotypes as 2-tuples:
#                     haploid: (sample, allele)
#                     diploid: (sample, [allele, allele])
#             If as_strings, the alleles will be actual alleles.  Otherwise,
#             alleles will be integers.
#             If more is true, a fifth column will be emitted with the pysam VCF object.
#         '''
#         for snp in self.get(c,start,stop,region):
#             alleles = [bytes_to_string(snp.ref)] + bytes_to_string(snp.alt).split(',')
#             alleles = [a for a in alleles if a != '.']
#             if not omit_genotypes:
#                 if self.ploidy==1:
#                     genos = [(self.sample_names[i], int(bytes_to_string(snp[i])[0]))
#                         for i in range(len(self.sample_names))
#                         if bytes_to_string(snp[i])[0] != '.']
#                     if as_strings:
#                         genos = [(s,alleles[a]) for s,a in genos]
#                 else:
#                     #genos = [i for i in range(len(self.sample_names)) if bytes_to_string(snp[i])[0] != '.']
#                     genos = [(self.sample_names[i], [int(bytes_to_string(snp[i])[j*2]) for j in range(self.ploidy)])
#                         for i in range(len(self.sample_names))
#                         if bytes_to_string(snp[i])[0] != '.']
#                     #print genos
#                     if as_strings:
#                         genos = [(s,[alleles[a] for a in g]) for s,g in genos]
#             else:
#                 genos = []
#             if more:
#                 yield (bytes_to_string(snp.contig), get_pos_from_vcf_record(snp), alleles, genos, snp)
#             else:
#                 yield (bytes_to_string(snp.contig), get_pos_from_vcf_record(snp), alleles, genos)
#     def get_snp_genos(self, c, p, as_strings=True):
#         ''' Read a single position from a VCF file and return the genotypes
#             as a sample -> allele map.  If there is not exactly one matching
#             row in the VCF file at this position (if there are none or multiple)
#             then we return an empty map: {}.
#         '''
#         snps = [x for x in self.get_range(c,p,p,as_strings=as_strings)]
#         return len(snps)==1 and dict(snps[0][3]) or {}
#     def getFullSequences(self, c, start, stop, samples,
#         na='-', refSeq=None, refInvariantsOnly=False, ignoreIndels=False):
#         ''' chr - chromosome name
#             start - start position
#             stop - default = start
#             samples is a list of samples.  None can be used for the ref sample.
#             if refSeq is a string with length = stop-start+1, then use this as the
#                 base sequence for missing data, otherwise, use the na string.
#                 if refInvariantsOnly is True, refSeq is only used for invariant
#                 sites or sites with no entries for any samples.  if False, refSeq
#                 is used for all missing data.
#         '''
#         assert 1<=start<=stop
#         assert len(na)==1
        
#         # get all the VCF records
#         vcf_records = [(p-start,alleles,dict(genos)) for chrom,p,alleles,genos
#             in self.get_range(c, start, stop, as_strings=True)
#             if not ignoreIndels or set(map(len, alleles)) == set([1])]
            
#         # Construct a list called "seq" into which we will replace alleles as
#         # we discover them.  This is a list, not a string, because each position
#         # might be replaced with an arbitrary length string (deletions will be
#         # empty strings, insertions will be strings > length 1).  This all gets
#         # converted to a string at the end.
#         if refSeq:
#             # assume refSeq alleles as a baseline for missing data
#             assert len(refSeq)==(stop-start+1)
#             seq = list(refSeq)
#             if refInvariantsOnly:
#                 # but don't assume refSeq alleles for any known variant sites
#                 for i,alleles,genos in vcf_records:
#                     if len(set(genos[s] for s in samples if s in genos))>1:
#                         for j in range(len(alleles[0])):
#                             if i+j < len(seq):
#                                 seq[i+j] = na
#         else:
#             # assume nothing when data is missing
#             seq = list(na * (stop-start+1))
        
#         for sample in samples:
#             assert sample==None or sample in self.samples()
            
#             # Make copy of reference sequence
#             newseq = [s for s in seq]
    
#             # Replace alleles, one site at a time.
#             replaceAlleles(sample, newseq, vcf_records)
            
#             # Convert list back to string.  This string may or may not be the same
#             # length as the original (if ignoreIndels is True, it must be the same
#             # length as the original).
#             newseq = ''.join(newseq)
#             assert len(newseq)==(stop-start+1) or not ignoreIndels
#             yield (sample, newseq)

# if __name__ == "__main__":
#     startTime = time.time()
#     print __file__
#     rootDir = os.path.dirname(os.path.realpath(__file__))
#     # actual number of records in chromosome 1: 6468094
#     #for i in range(1,22):

#     chromosomeNum = 1
#     filePath = rootDir + "/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz".format(chromosomeNum)
#     #print filePath

#     #vcf = VCFProcessor(filePath)
#     #records = vcf.records(chromosomeNum, 0, 250240543)
#     #print vcf.sample_names
#     #records = vcf.records(chromosomeNum, 0, 240543)
#     #for record in records:
#     #    print vcf.vcf_record_info(record[-1])
#     #    break

#     reader = VcfReader(filePath)
#     reader.ploidy = 2
#     #records = reader.get_range(chromosomeNum, 0, 5000000, more=True)
#     records = reader.get_range(chromosomeNum, 0, 250240543, more=True)
#     #print list(records)

#     snpCount = 0

#     highQualityAACalls = 0
#     lowQualityAACalls  = 0
#     unknownAACalls     = 0
#     missingAACalls     = 0

#     for record in records:
#         #print record[-2]
        
#         infoString = record[-1].info

#         if reader.variant_is_type(infoString, "SNP"):
#             ancestralAllele = reader.parse_ancestral(infoString)

#             if len(ancestralAllele):
#                 # if the AA is high-quality call
#                 if ord(ancestralAllele) in range(65,91):
#                     highQualityAACalls += 1
#                 if ord(ancestralAllele) in range(97,123):
#                     lowQualityAACalls += 1
#                 if ord(ancestralAllele) == 46:
#                     unknownAACalls += 1
#             else:
#                 missingAACalls += 1

#             #print ""
#             #print "reference allele:  {}".format(record[-1].ref)
#             #print "variant allele:    {}".format(record[-1].alt)
#             #print "ancestral allele:  {}".format(ancestralAllele)
            
#             if snpCount % 500 == 0:
#                 print "  "
#                 print "SNP count: {}".format(snpCount)
#                 print 'Elapsed runtime: {:.2} minutes.'.format((float(time.time()-startTime)/60.0))

#             snpCount += 1
#         if snpCount == 10:
#             break

#     totalAACalls = highQualityAACalls + lowQualityAACalls
#     highQualityCallPercentage = (float(highQualityAACalls) / float(totalAACalls))
#     missingAACallPercentage = (float(missingAACalls) / float(totalAACalls))

#     print "========"
#     print "highQualityAACalls: {}".format(highQualityAACalls)
#     print "lowQualityAACalls:  {}".format(lowQualityAACalls)
#     print "unknownAACalls:     {}".format(unknownAACalls)
#     print "missingAACalls:     {}".format(missingAACalls)
#     print "---"
#     print "highQualityCallPercentage: {:.2%}".format(highQualityCallPercentage)
#     print "missingAACallPercentage:   {:.2%}".format(missingAACallPercentage)
    #print vcf.sample_names

    

    #print vcf.count_iter_items(records)

    #test()
