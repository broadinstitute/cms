#!/usr/bin/env python
''' This script contains command-line utilities for calculating EHH-based scans
    for positive selection in genomes, including EHH, iHS, and XP-EHH.
    '''

__author__   = "tomkinsc@broadinstitute.org"
__commands__ = []

# built-ins
import argparse, logging, os
import imp # for programmatic/dynamic import

# from external modules
# from Bio import SeqIO

# from this package
import tools.selscan
import util.cmd, util.file, util.vcf_reader, util.call_sample_reader, util.json_helpers, util.db_helper

log = logging.getLogger(__name__)


class Selection(object):
    def __init__(self):
        pass

    def main(self):
        return 0

# File conversion parser to create selscan-compatible input, including population filtering

def parser_selscan_file_conversion(parser=argparse.ArgumentParser()):
    """
            Process a bgzipped-VCF (such as those included in the Phase 3 1000 Genomes release) into a gzip-compressed
            tped file of the sort expected by selscan. 
    """
    parser.help = """This function converts an input VCF into a TPED file usable by selscan.
                     It can optionally filter samples by population given a callsample file, and by ancestral call quality."""

    parser.add_argument("inputVCF", help="Input VCF file")
    parser.add_argument("genMap", help="Genetic recombination map tsv file with four columns: (Chromosome, Position(bp), Rate(cM/Mb), Map(cM))")
    parser.add_argument("outPrefix", help="Output file prefix")
    parser.add_argument("outLocation", help="Output location")        

    parser.add_argument('chromosomeNum', type=str, help="""Chromosome number.""")
    parser.add_argument('--startBp', default=0, type=int,
        help="""Coordinate in bp of start position. (default: %(default)s).""")
    parser.add_argument('--endBp', type=int,
        help="""Coordinate in bp of end position.""")
    parser.add_argument('--ploidy', default=2, type=int,
        help="""Number of chromosomes expected for each genotype. (default: %(default)s).""")
    #parser.add_argument("--ancestralVcf", type=str,
    #    help="""A one-sample VCF file describing the ancestral allele for all
    #    SNPs.  File must be built on the same reference and coordinate space as inVcf.  All
    #    het genotypes will be ignored.  This will allow us to code ancestral alleles as
    #    1 and derived alleles as 0.  If this file is not specified, we will use major
    #    alleles as 1 (ancestral) and minor alleles as 0 (derived).""",
    #    default=None)

    parser.add_argument('--considerMultiAllelic', default=False, action='store_true',
        help='Include multi-allelic variants in the output as separate records')
    parser.add_argument('--rescaleGeneticDistance', default=False, action='store_true',
        help='Genetic distance is rescaled to be out of 100.0 cM')
    parser.add_argument('--includeLowQualAncestral', default=False, action='store_true',
        help='Include variants where the ancestral information is low-quality ' + 
        '(as indicated by lower-case x for AA=x in the VCF info column) (default: %(default)s).')
    


    #subparsers = parser.add_subparsers()
    #subparser = subparsers.add_parser("filter")
    parser.add_argument("--codingFunctionClassFile", help="A python class file containing a function used to code each genotype as '1' and '0'. coding_function(current_value, reference_allele, alternate_allele, ancestral_allele)")
    parser.add_argument("--multiAllelicMergeFunction", help=argparse.SUPPRESS, default="OR", choices=["OR","AND","XOR"])
    parser.add_argument("--sampleMembershipFile", help="The call sample file containing four columns: sample, pop, super_pop, gender")
    parser.add_argument('--filterPops', nargs='+',
        help='Populations to include in the calculation (ex. "FIN")')
    parser.add_argument('--filterSuperPops', nargs='+',
        help='Super populations to include in the calculation (ex. "EUR")')

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_file_conversion)
    return parser

def main_selscan_file_conversion(args):
    
    # define coding functions here, in accordance with spec arg preferences

    if args.ploidy < 1:
        raise argparse.ArgumentError('Argument "--ploidy" must be one or greater.')

    if (args.filterPops or args.filterSuperPops) and not args.sampleMembershipFile:
        raise argparse.ArgumentTypeError('Argument "--sampleMembershipFile" must be specifed if --filterPops or --filterSuperPops are used.')
    
    samples_to_include = None
    if (args.filterPops or args.filterSuperPops) and args.sampleMembershipFile:
        csr = util.call_sample_reader.CallSampleReader(args.sampleMembershipFile)
        samples_to_include = list(csr.filterSamples(pop=args.filterPops, super_pop=args.filterSuperPops))

    # Write JSON file containing metadata describing this TPED file

    if args.codingFunctionClassFile is not None:
        genoCoder = importPythonClassFromFile(args.codingFunctionClassFile)
        coding_function = genoCoder.coding_function
    else:
        #def coding_function(current_value, val_representing_alt_allele, reference_allele, ancestral_allele, alt_allele):
        def coding_function(current_value, reference_allele, alternate_allele, ancestral_allele, value_of_current_allele):
            """for a given genotype call, returns the proper tped/selscan encoding, such that the ancestral allele is 
            always 1 and any derived allele is coded as 0. Collapses multiallelic sites to biallelic, whether on one
            or two lines (??) together with bitwise row merge.
            current_value: char, 0-3; the genotype call as printed in the VCF
            reference_allele: char, ACGT; the allele that is coded as "0"
            alternate_allele: char, ACGT; the allele 
            """

            if ancestral_allele == reference_allele: 
                if current_value == "0":
                    return "1"

            if ancestral_allele == alternate_allele:
                if current_value == value_of_current_allele:
                    return "1"
            
            return "0"

    tools.selscan.SelscanFormatter().process_vcf_into_selscan_tped(  vcf_file = args.inputVCF, 
                                    gen_map_file                              = args.genMap, 
                                    samples_to_include                        = samples_to_include, 
                                    outfile_location                          = args.outLocation, 
                                    outfile_prefix                            = args.outPrefix, 
                                    chromosome_num                            = args.chromosomeNum, 
                                    start_pos_bp                              = args.startBp, 
                                    end_pos_bp                                = args.endBp, 
                                    ploidy                                    = args.ploidy, 
                                    consider_multi_allelic                    = args.considerMultiAllelic, 
                                    rescale_genetic_distance                  = args.rescaleGeneticDistance, 
                                    include_variants_with_low_qual_ancestral  = args.includeLowQualAncestral,
                                    coding_function                           = coding_function,
                                    multi_alleli_merge_function               = args.multiAllelicMergeFunction)

    outTpedFile = args.outLocation  + "/" + args.outPrefix + ".tped.gz"
    outMetadataFile = args.outLocation  + "/" + args.outPrefix + ".metadata.json"
    outTpedAlleleMetaFile = args.outLocation  + "/" + args.outPrefix + ".tped.allele_metadata.gz"

    metaDataDict = {}
    metaDataDict["vcf_file"]                  = args.inputVCF
    metaDataDict["chromosome_num"]            = args.chromosomeNum,
    metaDataDict["call_membership_file"]      = args.sampleMembershipFile
    metaDataDict["genetic_map"]               = args.genMap
    metaDataDict["tped_file"]                 = outTpedFile
    metaDataDict["tped_allele_metadata_file"] = outTpedAlleleMetaFile
    metaDataDict["samples_included"]          = samples_to_include
    metaDataDict["populations"]               = args.filterPops if args.filterPops is not None else []
    metaDataDict["super_populations"]         = args.filterSuperPops if args.filterSuperPops is not None else []
    util.json_helpers.JSONHelper.annotate_json_file(outMetadataFile, metaDataDict)

    return 0
__commands__.append(('selscan_file_conversion', parser_selscan_file_conversion))

# === Selscan common parser ===

def parser_selscan_common(parser=argparse.ArgumentParser()):
    """
        Build a parser to which arguments are added which are common to 
        several selscan functions.
    """
    parser.add_argument("inputTped", help="Input tped file")
    parser.add_argument("outFile", help="Output filepath")
    
    #parser.add_argument('--mapFile', type=str,
    #    help="""Recombination map file; tab delimited with four columns:
    #    chrom, pos (bp), recom rate (cM/Mb), and map pos (cM)""",
    #    default=None)
    parser.add_argument('--gapScale', default=20000, type=int,
        help="""Gap scale parameter in bp. If a gap is encountered between
        two snps > GAP_SCALE and < MAX_GAP, then the genetic distance is
        scaled by GAP_SCALE/GA (default: %(default)s).""")
    
    parser.add_argument('--maf', default=0.05, type=float,
        help="""Minor allele frequency. If a site has a MAF below this value, the program will not use
        it as a core snp. (default: %(default)s).""")
    parser.add_argument('--threads', default=1, type=int,
        help="""The number of threads to spawn during the calculation.
        Partitions loci across threads. (default: %(default)s).""")

    return parser

# === Selscan EHH ===

def parser_selscan_ehh(parser=argparse.ArgumentParser()):
    """
        Perform selscan's calculation of EHH.
    """
    parser = parser_selscan_common(parser)

    parser.add_argument("locusID", help="The locus ID")

    parser.add_argument('--window', default=100000, type=int,
        help="""When calculating EHH, this is the length of the window in bp in 
        each direction from the query locus (default: %(default)s).""")
    parser.add_argument('--cutoff', default=0.05, type=float,
        help="""The EHH decay cutoff (default: %(default)s).""")
    parser.add_argument('--maxExtend', default=1000000, type=int,
        help="""The maximum distance an EHH decay curve is allowed to extend from the core.
        Set <= 0 for no restriction. (default: %(default)s).""")

    parser.epilog = """Output format: <physicalPosDelta_bp> <geneticPosDelta_cm> <'1' EHH> <'0' EHH><EHH>"""

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_ehh)
    return parser

def main_selscan_ehh(args):
    if args.threads < 1:
        raise argparse.ArgumentTypeError("You must specify more than 1 thread. %s threads given." % args.threads)

    # Check to see if json file exists for given tped, and if so annotate it with the output filenames

    tools.selscan.SelscanTool().execute_ehh(
        locus_id        = args.locusID,
        tped_file       = args.inputTped,
        out_file        = args.outFile,
        window          = args.window,
        cutoff          = args.cutoff,
        max_extend      = args.maxExtend,
        threads         = args.threads,
        maf             = args.maf,
        gap_scale       = args.gapScale
    )

    locusMetadata = {}
    locusMetadata["locus"]   = args.locusID
    locusMetadata["ehh"]     = args.outFile + ".ehh." + args.locusID + ".out"
    locusMetadata["ehh_anc"] = args.outFile + ".ehh." + args.locusID + ".out" + ".anc.colormap"
    locusMetadata["ehh_der"] = args.outFile + ".ehh." + args.locusID + ".out" + ".der.colormap"
    locusMetadata["ehh_log"] = args.outFile + ".ehh." + args.locusID + ".log"

    jsonFilePath = os.path.abspath(args.inputTped).replace(".tped.gz",".metadata.json")
    util.json_helpers.JSONHelper.annotate_json_file(jsonFilePath, locusMetadata, key_to_act_on="ehh", append=True)

    return 0
__commands__.append(('selscan_ehh', parser_selscan_ehh))

# === Selscan IHS ===

def parser_selscan_ihs(parser=argparse.ArgumentParser()):
    """
        Perform selscan's calculation of iHS.
    """
    parser = parser_selscan_common(parser)

    parser.add_argument('--skipLowFreq', default=False, action='store_true',
        help=' Do not include low frequency variants in the construction of haplotypes (default: %(default)s).')
    parser.add_argument('--dontWriteLeftRightiHH', default=False, action='store_true',
        help=' When writing out iHS, do not write out the constituent left and right ancestral and derived iHH scores for each locus.(default: %(default)s).')
    parser.add_argument('--truncOk', default=False, action='store_true',
        help="""If an EHH decay reaches the end of a sequence before reaching the cutoff,
        integrate the curve anyway.
        Normal function is to disregard the score for that core. (default: %(default)s).""")

    parser.epilog = """Output format: <locusID> <physicalPos_bp> <'1' freq> <ihh1> <ihh0> <unstandardized iHS> <derived_ihh_left> <derived_ihh_right> <ancestral_ihh_left> <ancestral_ihh_right>"""

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_ihs)
    return parser

def main_selscan_ihs(args):
    if args.threads < 1:
        raise argparse.ArgumentTypeError("You must specify more than 1 thread. %s threads given." % args.threads)

    tools.selscan.SelscanTool().execute_ihs(
        tped_file          = args.inputTped,
        out_file           = args.outFile,
        write_detailed_ihh = True if not args.dontWriteLeftRightiHH else False,
        skip_low_freq      = args.skipLowFreq,
        trunc_ok           = args.truncOk,
        threads            = args.threads,
        maf                = args.maf,
        gap_scale          = args.gapScale   
    )

    metaData = {}
    metaData["ihs"] = args.outFile+".ihs.out"
    metaData["ihs_log"] = args.outFile+".ihs.log"
    jsonFilePath = os.path.abspath(args.inputTped).replace(".tped.gz",".metadata.json")
    util.json_helpers.JSONHelper.annotate_json_file(jsonFilePath, metaData)

    return 0
__commands__.append(('selscan_ihs', parser_selscan_ihs))

# === Selscan nSL

def parser_selscan_nsl(parser=argparse.ArgumentParser()):
    """
        Perform selscan's calculation of nSL.
    """
    parser = parser_selscan_common(parser)

    parser.add_argument('--truncOk', default=False, action='store_true',
        help="""If an EHH decay reaches the end of a sequence before reaching the cutoff,
        integrate the curve anyway.
        Normal function is to disregard the score for that core. (default: %(default)s).""")
    parser.add_argument('--maxExtendNsl', default=100, type=int,
        help="""The maximum distance an nSL haplotype is allowed to extend from the core.
        Set <= 0 for no restriction. (default: %(default)s).""")

    #parser.epilog = """Output format: <locusID> <physicalPos_bp> <'1' freq> <ihh1> <ihh0> <unstandardized iHS> <derived_ihh_left> <derived_ihh_right> <ancestral_ihh_left> <ancestral_ihh_right>"""

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_nsl)
    return parser

def main_selscan_nsl(args):
    if args.threads < 1:
        raise argparse.ArgumentTypeError("You must specify more than 1 thread. %s threads given." % args.threads)

    tools.selscan.SelscanTool().execute_nsl(
        tped_file          = args.inputTped,
        out_file           = args.outFile,
        trunc_ok           = args.truncOk,
        threads            = args.threads,
        maf                = args.maf,
        max_extend_nsl     = args.maxExtendNsl,
        gap_scale          = args.gapScale   
    )

    metaData = {}
    metaData["nsl"] = args.outFile+".nsl.out"
    metaData["nsl_log"] = args.outFile+".nsl.log"
    jsonFilePath = os.path.abspath(args.inputTped).replace(".tped.gz",".metadata.json")
    util.json_helpers.JSONHelper.annotate_json_file(jsonFilePath, metaData)

    return 0
__commands__.append(('selscan_nsl', parser_selscan_nsl))

# === Selscan XPEHH ===

def parser_selscan_xpehh(parser=argparse.ArgumentParser()):
    """
        Perform selscan's calculation of XPEHH.
    """
    parser = parser_selscan_common(parser)

    parser.add_argument("inputRefTped", help="Input tped for the reference population to which the first is compared")

    parser.add_argument('--truncOk', default=False, action='store_true',
        help="""If an EHH decay reaches the end of a sequence before reaching the cutoff,
        integrate the curve anyway.
        Normal function is to disregard the score for that core. (default: %(default)s).""")

    parser.epilog = """Output format: <locusID> <physicalPos_bp> <geneticPos_cm> <popA '1' freq> <ihhA> <popB '1' freq> <ihhB> <unstandardized XPEHH>"""

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_xpehh)
    return parser

def main_selscan_xpehh(args):
    if args.threads < 1:
        raise argparse.ArgumentTypeError("You must specify more than 1 thread. %s threads given." % args.threads)
         
    tools.selscan.SelscanTool().execute_xpehh(
        tped_file       = args.inputTped,
        tped_ref_file   = args.inputRefTped,
        out_file        = args.outFile,
        trunc_ok        = args.truncOk,
        threads         = args.threads,
        maf             = args.maf,
        gap_scale       = args.gapScale
    )

    metaData = {}
    metaData["xpehh"] = args.outFile+".xpehh.out"
    metaData["xpehh_log"] = args.outFile+".xpehh.log"
    metaData["xpehh_pop_A_tped"] = args.inputTped
    metaData["xpehh_pop_B_tped"] = args.inputRefTped
    jsonFilePath = os.path.abspath(args.inputTped).replace(".tped.gz",".metadata.json")

    util.json_helpers.JSONHelper.annotate_json_file(jsonFilePath, metaData, key_to_act_on="xpehh", append=True)

    return 0
__commands__.append(('selscan_xpehh', parser_selscan_xpehh))

# === Selscan normalization ===

def parser_selscan_norm_common(parser=argparse.ArgumentParser()):
    """
        Build a parser to which arguments are added which are common to 
        several selscan functions.
    """

    input_file_help_string = """A list of files delimited by whitespace for joint normalization.
    Expected format for iHS/nSL files (no header):\n
    <locus name> <physical pos> <freq> <ihh1/sL1> <ihh2/sL0> <ihs/nsl>
    Expected format for XP-EHH files (one line header):
    <locus name> <physical pos> <genetic pos> <freq1> <ihh1> <freq2> <ihh2> <xpehh>"""

    parser.add_argument("inputFiles", help=input_file_help_string, nargs='+')

    parser.add_argument('--bins', default=100, type=int,
        help="""The number of frequency bins in [0,1] for score normalization (default: %(default)s)""")

    parser.add_argument('--critPercent', default=-1.00, type=float,
        help="""Set the critical value such that a SNP with iHS in the most extreme CRIT_PERCENT tails 
        (two-tailed) is marked as an extreme SNP. Not used by default (default: %(default)s)""")

    parser.add_argument('--critValue', default=2.00, type=float,
        help="""Set the critical value such that a SNP with |iHS| > CRIT_VAL is marked as an extreme SNP.  
        Default as in Voight et al. (default: %(default)s)""")

    parser.add_argument('--minSNPs', default=10, type=int,
        help="""Only consider a bp window if it has at least this many SNPs (default: %(default)s)""")

    parser.add_argument('--qbins', default=20, type=int,
        help="""Outlying windows are binned by number of sites within each
        window.  This is the number of quantile bins to use. (default: %(default)s)""")

    parser.add_argument('--winSize', default=100000, type=int,
        help="""GThe non-overlapping window size for calculating the percentage
        of extreme SNPs (default: %(default)s)""")

    parser.add_argument('--bpWin', default=False, action='store_true',
        help='If set, will use windows of a constant bp size with varying number of SNPs')

    return parser

def parser_selscan_norm_nsl(parser=argparse.ArgumentParser()):
    parser = parser_selscan_norm_common(parser)
    
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_norm_nsl)

    return parser

def main_selscan_norm_nsl(args):
    """
        Normalize Selscan's nSL output
    """
    # impose sanity checks on values here, perhaps

    tools.selscan.SelscanNormTool().execute_nsl_norm(
        input_file_list = args.inputFiles,
        bins            = args.bins,
        crit_percent    = args.critPercent,
        crit_val        = args.critValue,
        min_snps        = args.minSNPs,
        qbins           = args.qbins,
        winsize         = args.winSize,
        bp_win          = args.bpWin
    )
    return 0
__commands__.append(('selscan_norm_nsl', parser_selscan_norm_nsl))

def parser_selscan_norm_ihs(parser=argparse.ArgumentParser()):
    parser = parser_selscan_norm_common(parser)
    
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_norm_ihs)

    return parser

def main_selscan_norm_ihs(args):
    """
        Normalize Selscan's iHS output
    """
    # impose sanity checks on values here, perhaps

    tools.selscan.SelscanNormTool().execute_ihs_norm(
        input_file_list = args.inputFiles,
        bins            = args.bins,
        crit_percent    = args.critPercent,
        crit_val        = args.critValue,
        min_snps        = args.minSNPs,
        qbins           = args.qbins,
        winsize         = args.winSize,
        bp_win          = args.bpWin
    )
    return 0
__commands__.append(('selscan_norm_ihs', parser_selscan_norm_ihs))


def parser_selscan_norm_xpehh(parser=argparse.ArgumentParser()):
    parser = parser_selscan_norm_common(parser)
    
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_norm_xpehh)

    return parser

def main_selscan_norm_xpehh(args):
    """
        Normalize Selscan's XPEHH output
    """
    # impose sanity checks on values here, perhaps

    tools.selscan.SelscanNormTool().execute_xpehh_norm(
        input_file_list = args.inputFiles,
        bins            = args.bins,
        crit_percent    = args.critPercent,
        crit_val        = args.critValue,
        min_snps        = args.minSNPs,
        qbins           = args.qbins,
        winsize         = args.winSize,
        bp_win          = args.bpWin
    )
    return 0
__commands__.append(('selscan_norm_xpehh', parser_selscan_norm_xpehh))

# === For storage and retrieval of data from SQLite ===

def parser_selscan_store_results_in_db(parser=argparse.ArgumentParser()):
    """
        Aggregate results from selscan in to a SQLite database via helper JSON metadata file.
    """
    parser.add_argument("inputFile", help="Input *.metadata.json file")
    parser.add_argument("outFile", help="Output SQLite filepath")

    util.cmd.common_args(parser, (('loglevel', 'INFO'), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_selscan_store_results_in_db)
    return parser

def main_selscan_store_results_in_db(args):
    
    inFilePath = os.path.abspath(args.inputFile)
    if os.path.isfile(inFilePath):
        metaData = util.json_helpers.JSONHelper.read_data(inFilePath)

        cr = util.db_helper.ScanStatStorer(inFilePath, args.outFile)

        if "locus_info_stored_in_db" not in metaData:
            log.info("Storing locus info in DB...")
            cr.store_locus_info()
            util.json_helpers.JSONHelper.annotate_json_file(inFilePath, {"locus_info_stored_in_db":str(args.outFile)})

        if "ihs_info_stored_in_db" not in metaData:
            log.info("Storing IHS info in DB...")
            cr.store_ihs()
            util.json_helpers.JSONHelper.annotate_json_file(inFilePath, {"ihs_info_stored_in_db":str(args.outFile)})

        if "ehh_info_stored_in_db" not in metaData:
            log.info("Storing EHH info in DB...")
            cr.store_ehh()
            util.json_helpers.JSONHelper.annotate_json_file(inFilePath, {"ehh_info_stored_in_db":str(args.outFile)})

        if "xpehh_info_stored_in_db" not in metaData:
            log.info("Storing XPEHH info in DB...")
            cr.store_xpehh()
            util.json_helpers.JSONHelper.annotate_json_file(inFilePath, {"xpehh_info_stored_in_db":str(args.outFile)})

        log.info("DB storage complete: %s", args.outFile)
    else:
        raise IOError('Input metadata file specified does not exist: %s', inFilePath )
    return 0

__commands__.append(('store_selscan_results_in_db', parser_selscan_store_results_in_db))

# ====== util ========

def importPythonClassFromFile(uri, absl=False):
    if not absl:
        uri = os.path.normpath(os.path.join(os.path.dirname(__file__), uri))
    filePath, fileName = os.path.split(uri)
    fileNameSansExtension, ext = os.path.splitext(fileName)
        
    filePathSansExtension = os.path.join(filePath, fileNameSansExtension)
        
    if os.path.exists(filePathSansExtension + '.pyc'):
        try:
            return imp.load_compiled(fileNameSansExtension, filePathSansExtension + '.pyc')
        except:
            pass
    if os.path.exists(filePathSansExtension + '.py'):
        try:
            return imp.load_source(fileNameSansExtension, filePathSansExtension + '.py')
        except:
            pass

# === Parser setup ===

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

