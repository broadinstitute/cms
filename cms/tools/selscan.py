'''
    SelScan: A program to calculate EHH-based scans for positive selection in genomes
    https://github.com/szpiech/selscan
'''

__author__ = "tomkinsc@broadinstitute.org"

# built-ins
import operator
import os, os.path, subprocess, re
import logging
from datetime import datetime, timedelta
import argparse
import operator

try:
    from itertools import izip as zip # pylint:disable=redefined-builtin
except ImportError: # py3 zip is izip
    pass

# intra-module dependencies
import tools, util.file
from util.vcf_reader import VCFReader
from util.recom_map import RecomMap

#third-party dependencies
import pysam
from boltons.timeutils import relative_time
import numpy as np

tool_version = '1.1.0'
url = 'https://github.com/szpiech/selscan/archive/{ver}.zip'

log = logging.getLogger(__name__)

class SelscanFormatter(object):
    genoRegex = re.compile("^(\S+\s+){9}(?P<genos>.*)$")

    @staticmethod
    def _moving_avg(x, prevAvg, N):
        return float( sum([int(prevAvg)] * (N-1)) + x) / float(N)

    @staticmethod
    def remove_last_line_of_file(fileObj, close=False):
        # Move the file pointer to the end of the file. 
        fileObj.seek(0, os.SEEK_END)

        # offset to skip the null character at the end of the file
        pos = fileObj.tell() - 1

        # Read each character in the file one at a time from the last 
        # character going backwards, searching for a newline character
        # If we find a new line, exit the search
        while pos > 0 and fileObj.read(1) != "\n":
            pos -= 1
            fileObj.seek(pos, os.SEEK_SET)

        # As we're not at the start of the file, 
        # delete all the characters ahead of this position
        if pos > 0:
            #fileObj.seek(pos, os.SEEK_SET)
            fileObj.truncate()

        if close:
            fileObj.close()

    @staticmethod
    def _build_variant_output_strings(chrm, idx, pos_bp, map_pos_cm, genos, 
                                                                ref_allele, alt_allele, ancestral_call, allele_freq):
        outputStringDict = dict()
        # outputStringDict["tpedString"]     = "{chr} {pos_bp}-{idx} {map_pos_cm} {pos_bp} {genos}\n".format(chr=chrm, 
        #     idx=idx, pos_bp=pos_bp, map_pos_cm=map_pos_cm, genos=" ".join(genos))
        outputStringDict["tpedString"]     = "{chr} {pos_bp} {map_pos_cm} {pos_bp} {genos}\n".format(chr=chrm,
                                                        pos_bp=pos_bp, map_pos_cm=map_pos_cm, genos=" ".join(genos))
        outputStringDict["metadataString"] = "{chr} {pos_bp}-{idx} {pos_bp} {map_pos_cm} {ref_allele}".format( chr=chrm, 
            idx        = idx, 
            pos_bp     = pos_bp, 
            map_pos_cm = map_pos_cm, 
            ref_allele = ref_allele) + "{alt_allele} {ancestral_call} {allele_freq}\n".format( alt_allele = alt_allele, 
                                                                                        ancestral_call = ancestral_call, 
                                                                                        allele_freq    = allele_freq)

        return outputStringDict

    @staticmethod
    def _build_map_output_string(chrm, pos_bp, map_pos_cm):
        return "{chr} {pos_bp} {map_pos_cm} {pos_bp}\n".format(chr=chrm, pos_bp=pos_bp, map_pos_cm=map_pos_cm)

    @classmethod
    def process_vcf_into_selscan_tped(cls, vcf_file, gen_map_file, outfile_location,
        outfile_prefix, chromosome_num, samples_to_include=None, start_pos_bp=None, end_pos_bp=None, ploidy=2, 
        consider_multi_allelic=True, rescale_genetic_distance=False, include_variants_with_low_qual_ancestral=False, coding_function=None, 
        multi_alleli_merge_function="AND"):
        """
            Process a bgzipped-VCF (such as those included in the Phase 3 1000 Genomes release) into a gzip-compressed
            tped file of the sort expected by selscan. 
        """
        assert ploidy > 0

        processor = VCFReader(vcf_file)

        end_pos = processor.clens[str(chromosome_num)] if end_pos_bp is None else end_pos_bp

        records = processor.records( str(chromosome_num), start_pos_bp, end_pos, pysam.asVCF())

        outTpedFile = outfile_location + "/" + outfile_prefix + ".tped.gz"
        outTpedMetaFile = outfile_location + "/" + outfile_prefix + ".tped.allele_metadata.gz"

        if samples_to_include is not None and len(samples_to_include) > 0:
            indices_of_matching_samples = sorted([processor.sample_names.index(x) for x in samples_to_include])
        else:
            indices_of_matching_samples = range(0,len(processor.sample_names))

        indices_of_matching_genotypes = [(x*2, (x*2)+1) for x in indices_of_matching_samples]
        indices_of_matching_genotypes = list(np.ravel(np.array(indices_of_matching_genotypes)))

        rm = RecomMap(gen_map_file)

        for filePath in [outTpedFile, outTpedMetaFile]:
            assert not os.path.exists(filePath), "File {} already exists. Consider removing this file or specifying a different output prefix. Processing aborted.".format(filePath)

        mergeOperatorString = ""
        if multi_alleli_merge_function == "OR":
            mergeOperatorString = "|"
        if multi_alleli_merge_function == "AND":
            mergeOperatorString = "&"
        if multi_alleli_merge_function == "XOR":
            mergeOperatorString = "^"

        startTime = datetime.now()
        sec_remaining_avg = 0
        current_pos_bp = 1

        with util.file.open_or_gzopen(outTpedFile, 'w') as of1, util.file.open_or_gzopen(outTpedMetaFile, 'w') as of2:
            # WRITE header for metadata file here with selected subset of sample_names
            headerString = "CHROM VARIANT_ID POS_BP MAP_POS_CM REF_ALLELE ALT_ALLELE ANCESTRAL_CALL ALLELE_FREQ_IN_POP\n".replace(" ","\t")
            of2.write(headerString)

            of1linesToWrite = []
            of2linesToWrite = []

            recordLengths = set() #

            recordCount = 0
            mostRecentRecordPosSeen = -1
            positionHasBeenSeenBefore = False
            previouslyCodedGenotypes = GenoRecord([])
            lineToWrite1 = None
            lineToWrite2 = None
            previousAncestral = None
            ancestralDiffersFromPrevious = True
            for record in records:
                # in some cases, there may be records with duplicate positions in the VCF file
                # to account for that we collapse rows that pass our filter and then write out the rows when we 
                # encounter a record with a new position 
                if record.pos != mostRecentRecordPosSeen:

                    if positionHasBeenSeenBefore and not consider_multi_allelic:
                        lineToWrite1 = None
                        lineToWrite2 = None

                    if lineToWrite1 is not None and lineToWrite2 is not None:
                        if len(previouslyCodedGenotypes) == ploidy*len(indices_of_matching_samples):
                            # write the output line
                            of1.write(lineToWrite1)
                            of2.write(lineToWrite2)
                               
                        else:
                            genotypesCount            = len(previouslyCodedGenotypes)
                            countSpecificTpedName     = outTpedFile.replace(outfile_prefix, outfile_prefix + "_" + str(genotypesCount))
                            countSpecificMetafileName = outTpedMetaFile.replace(outfile_prefix, outfile_prefix + "_" + str(genotypesCount))
                            with util.file.open_or_gzopen(countSpecificTpedName, 'a') as of1l, util.file.open_or_gzopen(countSpecificMetafileName, 'a') as of2l:
                                of1l.write(lineToWrite1)
                                of2l.write(lineToWrite2)

                        lineToWrite1 = None
                        lineToWrite2 = None
                    mostRecentRecordPosSeen = record.pos
                    positionHasBeenSeenBefore = False
                else:
                    positionHasBeenSeenBefore = True

                # if the variant is a SNP
                # OLD style looking at INFO VT value: 
                # processor.variant_is_type(record.info, "SNP"):
                VALID_BASES = ["A","C","G","T","N","a","c","g","t","n"]
                if (len(record.ref) == 1 and len(record.alt) == 1) or ( all(variant in VALID_BASES for variant in record.ref.split(",")) and 
                     all(variant in VALID_BASES for variant in record.alt.split(",")) ):

                    alternateAlleles = [record.alt]
                    if record.alt not in ['A','T','C','G']:
                        #print record.alt
                        if consider_multi_allelic:
                            pass
                            alternateAlleles = record.alt.split(",")
                        else:
                            # continue on to next variant record
                            continue

                    ancestral_allele = processor.parse_ancestral(record.info)
                    chromStr = "chr{}".format(record.contig)

                    # if the AA is populated, and the call meets the specified criteria
                    if (ancestral_allele in ['A','T','C','G']) or (include_variants_with_low_qual_ancestral and ancestral_allele in ['a','t','c','g']):
                        
                        if previousAncestral != ancestral_allele:
                            previousAncestral = ancestral_allele
                            ancestralDiffersFromPrevious = True
                        else:
                            ancestralDiffersFromPrevious = False

                        recordString = record.__str__()

                        match = cls.genoRegex.match(recordString)
                        if match:
                            rawGenos = match.group("genos")
                            genos = rawGenos[::2]
                            recordPosition = record.pos+1
                            if chromosome_num.upper() != "X":
                                try:
                                    genotypes_for_selected_samples = operator.itemgetter(*indices_of_matching_genotypes)(genos)
                                except Exception: # if this is a record of mixed ploidy, that is to say the X chromosome
                                    raise
                            else:
                                matching_genotypes = np.array(record[:len(record)])[indices_of_matching_samples]
                                genotypes_for_selected_samples_split = [x.split("|") for x in matching_genotypes]
                                genotypes_for_selected_samples = [y for x in genotypes_for_selected_samples_split for y in x]

                            recordLengths.add(len(genotypes_for_selected_samples))

                            map_pos_cm = rm.physToMap(chromStr, record.pos, rescale=rescale_genetic_distance)

                            numberOfHaplotypes = float(len(genotypes_for_selected_samples))
                            
                            codingFunc = np.vectorize(coding_function)

                            coded_genotypes_for_selected_samples = GenoRecord(["0"] * len(genotypes_for_selected_samples))
                            if consider_multi_allelic:
                                #coded_genotypes_for_selected_samples = GenoRecord(["1"] * len(genotypes_for_selected_samples))
                                for idx, altAllele in enumerate(alternateAlleles):
                                    value_of_current_allele = str(idx+1)
                                    coded_genotypes_for_selected_samples_for_allele = GenoRecord(codingFunc(genotypes_for_selected_samples, record.ref, altAllele, ancestral_allele, value_of_current_allele))
                                    #coded_genotypes_for_selected_samples |= coded_genotypes_for_selected_samples_for_allele
                                    if idx==0:
                                        coded_genotypes_for_selected_samples = coded_genotypes_for_selected_samples_for_allele
                                    else:
                                        coded_genotypes_for_selected_samples = coded_genotypes_for_selected_samples.f[mergeOperatorString](coded_genotypes_for_selected_samples_for_allele)
                                    #coded_genotypes_for_selected_samples = np.array(list(str(bin(int("".join(coded_genotypes_for_selected_samples),2) | int("".join(coded_genotypes_for_selected_samples_for_allele),2)))[2:].zfill(len(coded_genotypes_for_selected_samples))))
                            else:
                                coded_genotypes_for_selected_samples = GenoRecord(codingFunc(genotypes_for_selected_samples, record.ref, record.alt, ancestral_allele))

                            # if this is the first record in the file, create an array filled with zeros for the previously coded alleles
                            if recordCount == 0:
                                previouslyCodedGenotypes = GenoRecord(["1"] * len(genotypes_for_selected_samples))
                                
                            # bitwise OR coded genotypes for duplicate records, merge variants
                            # ...except selscan logic is inverted, so bitwise AND
                            # TODO: invert?
                            # record @ pos1 = 001001
                            # record @ pos1 = 100001
                            #                -------
                            # coded result  = 101001
                            #log.debug(genotypes_for_selected_samples)
                            #log.debug(coded_genotypes_for_selected_samples)

                            if positionHasBeenSeenBefore:
                                #coded_genotypes_for_selected_samples |= previouslyCodedGenotypes
                                coded_genotypes_for_selected_samples = coded_genotypes_for_selected_samples.f[mergeOperatorString](previouslyCodedGenotypes)
                                #coded_genotypes_for_selected_samples = np.array(list(str(bin(int("".join(coded_genotypes_for_selected_samples),2) & int("".join(previouslyCodedGenotypes),2)))[2:].zfill(len(coded_genotypes_for_selected_samples))))

                            previouslyCodedGenotypes = coded_genotypes_for_selected_samples

                            allele_freq_for_pop = float(list(coded_genotypes_for_selected_samples).count("1")) / numberOfHaplotypes

                            outStrDict = cls._build_variant_output_strings(record.contig, str(1), recordPosition, 
                                map_pos_cm, coded_genotypes_for_selected_samples, record.ref, record.alt, 
                                ancestral_allele, allele_freq_for_pop)

                            lineToWrite1 = outStrDict["tpedString"]
                            lineToWrite2 = outStrDict["metadataString"].replace(" ","\t")

                            recordCount += 1
                            current_pos_bp = int(recordPosition)

                            if recordCount % 1000 == 0:
                                number_of_seconds_elapsed = (datetime.now() - startTime).total_seconds()
                                bp_per_sec = float(current_pos_bp) / float(number_of_seconds_elapsed)
                                bp_remaining = end_pos - current_pos_bp
                                sec_remaining = bp_remaining / bp_per_sec
                                sec_remaining_avg = cls._moving_avg(sec_remaining, sec_remaining_avg, 10)
                                time_left = timedelta(seconds=sec_remaining_avg)
                            
                                if sec_remaining > 10:
                                    human_time_remaining = relative_time(datetime.utcnow()+time_left)
                                    print("")
                                    print("Completed: {:.2%}".format(float(current_pos_bp)/float(end_pos)))
                                    print("Estimated time of completion: {}".format(human_time_remaining))
                                    #log.info("Genotype counts found: %s", str(list(recordLengths)))

            if positionHasBeenSeenBefore and not consider_multi_allelic:
                lineToWrite1 = None
                lineToWrite2 = None

            if lineToWrite1 is not None and lineToWrite2 is not None:
                if len(previouslyCodedGenotypes) == ploidy*len(indices_of_matching_samples):
                    # write the output lines
                    of1.write(lineToWrite1)
                    of2.write(lineToWrite2)
                       
                else:
                    genotypesCount            = len(previouslyCodedGenotypes)
                    countSpecificTpedName     = outTpedFile.replace(outfile_prefix, outfile_prefix + "_" + str(genotypesCount))
                    countSpecificMetafileName = outTpedMetaFile.replace(outfile_prefix, outfile_prefix + "_" + str(genotypesCount))
                    with util.file.open_or_gzopen(countSpecificTpedName, 'a') as of1l, util.file.open_or_gzopen(countSpecificMetafileName, 'a') as of2l:
                        of1l.write(lineToWrite1)
                        of2l.write(lineToWrite2)

            log.info("Genotype counts found: %s", str(list(recordLengths)))

class SelscanBaseTool(tools.Tool):
    def __init__(self, install_methods = None):
        if install_methods is None:
            install_methods = []
            os_type                 = self.get_os_type()
            binaryPath              = self.get_selscan_binary_path( os_type )
            
            target_rel_path = 'selscan-{ver}/{binPath}'.format(ver=tool_version, binPath=binaryPath)
            verify_command  = '{dir}/selscan-{ver}/{binPath} --help > /dev/null 2>&1'.format(dir=util.file.get_build_path(), 
                                                                                ver=tool_version, binPath=binaryPath) 

            install_methods.append(
                    DownloadAndBuildSelscan(  url.format( ver=tool_version ),
                                            target_rel_path = target_rel_path,
                                            verifycmd       = verify_command,
                                            verifycode      = 1 # selscan returns exit code of 1 for the help text...
                    )
            )

        tools.Tool.__init__(self, install_methods = install_methods)

    def version(self):
        return tool_version

    @classmethod
    def uncomment_line(cls, line):
        line = line.strip()
        while line.startswith( ("#", " ") ):
            line = line[1:]
        return line.strip()

    @classmethod
    def comment_line(cls, line):
        return "#" + cls.uncomment_line(line)

    @classmethod
    def change_line_comment(cls, line, comment_action):
        if comment_action == 1:
            return cls.comment_line(line)
        if comment_action == -1:
            return cls.uncomment_line(line)
        return line

    @classmethod
    def reform_makefile(cls, inputFileName, outputFileName):
        OSX_DELIMITER = "#for osx systems\n"
        LINUX_DELIMITER = "#for linux systems\n"

        inputFilePath = os.path.abspath(inputFileName)

        outputFilePath = os.path.dirname(os.path.realpath(inputFilePath))
        outputFilePath = os.path.join("/{}".format(outputFileName))

        current_os = cls.get_os_type()
        comment_action = 0 #-1 = remove comments until blank line, 0 = do not change line, 1 = comment line

        with open(inputFilePath, "r") as inFile:
            with open(outputFilePath, "w") as outFile:
                for line in inFile:
                    if line == OSX_DELIMITER:
                        outFile.write(line)
                        if current_os == "osx":
                            comment_action = -1
                        else: 
                            comment_action = 1
                        continue
                    if line == LINUX_DELIMITER:
                        outFile.write(line)
                        if current_os == "linux":
                            comment_action = -1
                        else: 
                            comment_action = 1
                        continue
                    if line.replace(" ","") == '\n' and comment_action != 0:
                        comment_action = 0

                    line = cls.change_line_comment(line, comment_action)
                    outFile.write(line.replace("\n","")+"\n")

    @staticmethod
    def get_os_type():
        ''' inspects the system uname and returns a string representing the OS '''

        uname = os.uname()
        if uname[0] == "Darwin":
            return "osx"
        if uname[0] == "Linux":
            return "linux"

    @staticmethod
    def get_selscan_binary_path(os_type, binaryName="selscan", full=True):
        ''' returns the location of the binary relative to the extracted archive, for the given os '''

        selscanPath = "bin/"

        if os_type == "osx":
            selscanPath += "osx/"
        elif os_type == "linux":
            selscanPath += "linux/"
        elif os_type == "win":
            selscanPath += "win/"

        if full:
            selscanPath += binaryName

            if os_type == "win":
                selscanPath += ".exe"

        return selscanPath

class SelscanTool(SelscanBaseTool):
    def __init__(self, install_methods = None):
        super(SelscanTool, self).__init__()

    def execute_ehh(self, locus_id, tped_file, out_file, window, cutoff, max_extend, threads, maf, gap_scale):
        out_file = os.path.abspath(out_file)

        toolCmd = [self.install_and_get_path()]
        toolCmd.append("--ehh")
        toolCmd.append(locus_id)
        toolCmd.append("--tped")
        toolCmd.append(tped_file)
        toolCmd.append("--out")
        toolCmd.append(out_file)
        if window:
            toolCmd.append("--ehh-win")
            toolCmd.append(window)
        if cutoff:
            toolCmd.append("--cutoff")
            toolCmd.append("{:.6}".format(cutoff))
        if max_extend:
            toolCmd.append("--max-extend")
            toolCmd.append((max_extend))
        if threads > 0:
            toolCmd.append("--threads")
            toolCmd.append((threads))
        else:
            raise argparse.ArgumentTypeError("You must specify more than 1 thread. %s threads given." % threads)
        if maf:
            toolCmd.append("--maf")
            toolCmd.append("{:.6}".format(maf))
        if gap_scale:
            toolCmd.append("--gap-scale")
            toolCmd.append((gap_scale))

        toolCmd = [str(x) for x in toolCmd]
        log.debug(' '.join(toolCmd))
        subprocess.check_call( toolCmd )

    def execute_ihs(self, tped_file, out_file, threads, maf, gap_scale, write_detailed_ihh=True, skip_low_freq=True, trunc_ok=False):
        toolCmd = [self.install_and_get_path()]
        toolCmd.append("--ihs")
        toolCmd.append("--tped")
        toolCmd.append(tped_file)
        toolCmd.append("--out")
        toolCmd.append(out_file)
        if skip_low_freq:
            toolCmd.append("--skip-low-freq")
        if write_detailed_ihh:
            toolCmd.append("--ihs-detail")
        if trunc_ok:
            toolCmd.append("--trunc-ok")
        if threads > 0:
            toolCmd.append("--threads")
            toolCmd.append((threads))
        else:
            raise argparse.ArgumentTypeError("You must specify more than 1 thread. %s threads given." % threads)
        if maf:
            toolCmd.append("--maf")
            toolCmd.append("{:.6}".format(maf))
        if gap_scale:
            toolCmd.append("--gap-scale")
            toolCmd.append((gap_scale))

        toolCmd = [str(x) for x in toolCmd]        
        log.debug(' '.join(toolCmd))
        subprocess.check_call( toolCmd )

    def execute_xpehh(self, tped_file, tped_ref_file, out_file, threads, maf, gap_scale, trunc_ok=False):
        toolCmd = [self.install_and_get_path()]
        toolCmd.append("--xpehh")
        toolCmd.append("--tped")
        toolCmd.append(tped_file)
        toolCmd.append("--tped-ref")
        toolCmd.append(tped_ref_file)
        toolCmd.append("--out")
        toolCmd.append(out_file)
        if trunc_ok:
            toolCmd.append("--trunc-ok")
        if threads > 0:
            toolCmd.append("--threads")
            toolCmd.append((threads))
        else:
            raise argparse.ArgumentTypeError("You must specify more than 1 thread. %s threads given." % threads)
        if maf:
            toolCmd.append("--maf")
            toolCmd.append("{:.6}".format(maf))
        if gap_scale:
            toolCmd.append("--gap-scale")
            toolCmd.append((gap_scale))

        toolCmd = [str(x) for x in toolCmd]
        log.debug(' '.join(toolCmd))
        subprocess.check_call( toolCmd )

    

class SelscanNormTool(SelscanBaseTool):
    def __init__(self, install_methods = None):
        super(SelscanNormTool, self).__init__()

        if install_methods is None:
            install_methods = []
            os_type                 = self.get_os_type()
            normBinaryPath          = self.get_selscan_binary_path( os_type, binaryName="norm" )
            
            target_norm_rel_path = 'selscan-{ver}/{binPath}'.format(ver=tool_version, binPath=normBinaryPath)
            verify_norm_command  = '{dir}/selscan-{ver}/{binPath} --help > /dev/null 2>&1'.format(dir=util.file.get_build_path(), 
                                                                                ver=tool_version, binPath=normBinaryPath) 

            install_methods.extend([
                    tools.DownloadPackage( url.format( ver=tool_version ),
                                            target_rel_path = target_norm_rel_path,
                                            verifycmd       = verify_norm_command,
                                            verifycode      = 1 ) # norm --help returns exit code of 1
            ])

        tools.Tool.__init__(self, install_methods = install_methods)

    def execute_ihs_norm(self, input_file_list, bins, crit_percent, crit_val, min_snps, qbins, winsize, bp_win=False):
        toolCmd = [self.install_and_get_path()]
        toolCmd.append("--ihs")

        if len(input_file_list):
            toolCmd.append("--files")
            toolCmd.append(" ".join(input_file_list))
        else:
            raise argparse.ArgumentTypeError("You must specify at least one file to normalize.")

        if bins:
            toolCmd.append("--bins")
            toolCmd.append(bins)
        if crit_percent:
            toolCmd.append("--crit-percent")
            toolCmd.append("{:.6}".format(crit_percent))
        if crit_val:
            toolCmd.append("--crit-val")
            toolCmd.append("{:.6}".format(crit_val))
        if min_snps:
            toolCmd.append("--min-snps")
            toolCmd.append(min_snps)
        if qbins:
            toolCmd.append("--qbins")
            toolCmd.append(qbins)
        if winsize:
            toolCmd.append("--winsize")
            toolCmd.append(winsize)
        if bp_win:
            toolCmd.append("--bp-win")

        toolCmd = [str(x) for x in toolCmd]        
        log.debug(' '.join(toolCmd))
        subprocess.check_call( toolCmd )

    def execute_xpehh_norm(self, input_file_list, bins, crit_percent, crit_val, min_snps, qbins, winsize, bp_win=False):
        toolCmd = [self.install_and_get_path()]
        toolCmd.append("--xpehh")

        
        if len(input_file_list):
            toolCmd.append("--files")
            toolCmd.append(" ".join(input_file_list))
        else:
            raise argparse.ArgumentTypeError("You must specify at least one file to normalize.")

        if bins:
            toolCmd.append("--bins")
            toolCmd.append(bins)
        if crit_percent:
            toolCmd.append("--crit-percent")
            toolCmd.append("{:.6}".format(crit_percent))
        if crit_val:
            toolCmd.append("--crit-val")
            toolCmd.append("{:.6}".format(crit_val))
        if min_snps:
            toolCmd.append("--min-snps")
            toolCmd.append(min_snps)
        if qbins:
            toolCmd.append("--qbins")
            toolCmd.append(qbins)
        if winsize:
            toolCmd.append("--winsize")
            toolCmd.append(winsize)
        if bp_win:
            toolCmd.append("--bp-win")

        toolCmd = [str(x) for x in toolCmd]        
        log.debug(' '.join(toolCmd))
        subprocess.check_call( toolCmd )

class GenoRecord(list):
    def __init__(self, genotypes):
        """
            the input parameter for instantiation, "genotypes", is expected
            to be a list of genotype values encoded as "0" or "1" (strings)
        """
        #self.genotypes = genotypes
        self.extend(genotypes)

        self.f =        { 
                            '~': self.__invert__,
                            '^': self.__xor__,
                            '|': self.__or__,
                            '&': self.__and__,
                            '' : self.__self__
                        }

    def __self__(self):
        return self

    #def __str__(self):
    #    return "".join(self)

    def __invert__(self):
        return GenoRecord(["0" if x=="1" else "1" for x in self])

    def __and__(self, other):
        return self._perform_bitwise(other, operator.and_)

    def __or__(self, other):
        return self._perform_bitwise(other, operator.or_)

    def __xor__(self, other):
        return self._perform_bitwise(other, operator.xor)

    def _perform_bitwise(self, other, operation):
        return GenoRecord(list( str( bin( operation(int("".join(self),2), int("".join(other),2)) ) )[2:].zfill(len(other)) ))

class DownloadAndBuildSelscan(tools.DownloadPackage) :
    def post_download(self) :
        selscanDir = os.path.join(self.destination_dir, 'selscan-{ver}'.format(ver=tool_version))
        selscanSrcDir = os.path.join(selscanDir,'src')
        # In this version, comment/uncomment the appropriate build flags
        makeFilePath = os.path.join(selscanSrcDir, 'Makefile')

        os.rename(makeFilePath, makeFilePath + '.orig')

        SelscanTool.reform_makefile(makeFilePath + '.orig', makeFilePath)

        # Now we can make:
        os.system('cd "{}" && make -s && mv ./selscan ./norm {}'.format(selscanSrcDir, 
            os.path.join(selscanDir, SelscanTool.get_selscan_binary_path(SelscanTool.get_os_type(), full=False))))






