#!/usr/bin/python
''' This script contains command-line utilities for calculating EHH-based scans
    for positive selection in genomes, including EHH, iHS, and XP-EHH.
    '''

__author__   = "tomkinsc@broadinstitute.org"
__commands__ = []

# built-ins
import array, argparse, logging, os, sys, array, bisect

# from external modules
from Bio import SeqIO

# from this package
import tools.selscan
import util.cmd, util.file, util.vcf

log = logging.getLogger(__name__)


class Selection(object):
    def __init__(self):
        pass

    def main(self):
        return 0

# === Selscan common parser ===

def parser_selscan_common(parser=argparse.ArgumentParser()):
    parser.add_argument("inVcf", help="Input VCF file")
    parser.add_argument("outPrefix", help="Base name for all output files")
    parser.add_argument("--ancestralVcf", type=str,
        help="""A one-sample VCF file describing the ancestral allele for all
        SNPs.  File must be on the same reference and coordinate space as inVcf.  All
        het genotypes will be ignored.  This will allow us to code ancestral alleles as
        1 and derived alleles as 0.  If this file is not specified, we will use major
        alleles as 1 (ancestral) and minor alleles as 0 (derived).""",
        default=None)
    parser.add_argument('--mapFile', type=str,
        help="""Recombination map file; tab delimited with four columns:
        chrom, pos (bp), recom rate (cM/Mb), and map pos (cM)""",
        default=None)
    parser.add_argument('--gap-scale', default=20000, type=int,
        help="""Gap scale parameter in bp. If a gap is encountered between
        two snps > GAP_SCALE and < MAX_GAP, then the genetic distance is
        scaled by GAP_SCALE/GA (default: %(default)s).""")
    parser.add_argument('--chromosome', type=int,
        help="""Chromosome number.""")
    parser.add_argument('--start-bp', default=0, type=int,
        help="""Coordinate in bp of start position. (default: %(default)s).""")
    parser.add_argument('--end-bp', type=int,
        help="""Coordinate in bp of end position.""")
    parser.add_argument('--maf', default=0.05, type=float,
        help="""Minor allele frequency. If a site has a MAF below this value, the program will not use
        it as a core snp. (default: %(default)s).""")
    parser.add_argument('--threads', default=1, type=int,
        help="""The number of threads to spawn during the calculation.
        Partitions loci across threads. (default: %(default)s).""")

    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser("filter")
    subparser.add_argument("inCallSampleFile", help="The file containing four columns: sample, pop, super_pop, gender")
    subparser.add_argument('pops', nargs='+',
        help='Populations to include in the calculation (ex. "FIN")')

    return parser

# === Selscan EHH ===

def parser_selscan_ehh(parser=argparse.ArgumentParser()):
    parser = parser_selscan_common(parser)

    parser.add_argument('--window', default=100000, type=int,
        help="""When calculating EHH, this is the length of the window in bp in 
        each direction from the query locus (default: %(default)s).""")
    parser.add_argument('--cutoff', default=0.05, type=float,
        help="""The EHH decay cutoff (default: %(default)s).""")
    parser.add_argument('--max-extend', default=1000000, type=int,
        help="""The maximum distance an EHH decay curve is allowed to extend from the core.
        Set <= 0 for no restriction. (default: %(default)s).""")    

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_ehh)
    return parser

def main_selscan_ehh(args):
    # selscan expects the map file to be in the format: 
    # <chr#> <locusID> <genetic pos> <physical pos>.
    #print "ARGS: ", args
    pass
__commands__.append(('selscan_ehh', parser_selscan_ehh))

# === Selscan IHS ===

def parser_selscan_ihs(parser=argparse.ArgumentParser()):
    parser = parser_selscan_common(parser)

    parser.add_argument('--skip-low-freq', default=False, action='store_true',
        help=' Do not include low frequency variants in the construction of haplotypes (default: %(default)s).')
    parser.add_argument('--trunc-ok', default=False, action='store_true',
        help="""If an EHH decay reaches the end of a sequence before reaching the cutoff,
        integrate the curve anyway.
        Normal function is to disregard the score for that core. (default: %(default)s).""")

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_ihs)
    return parser

def main_selscan_ihs(args):
    pass
__commands__.append(('selscan_ihs', parser_selscan_ihs))

# === Selscan XPEHH ===

def parser_selscan_xpehh(parser=argparse.ArgumentParser()):
    parser = parser_selscan_common(parser)

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_xpehh)
    return parser

def main_selscan_xpehh(args):
    parser.add_argument('--trunc-ok', default=False, action='store_true',
        help="""If an EHH decay reaches the end of a sequence before reaching the cutoff,
        integrate the curve anyway.
        Normal function is to disregard the score for that core. (default: %(default)s).""")
    
    pass
__commands__.append(('selscan_xpehh', parser_selscan_xpehh))

# === Parser setup ===

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)