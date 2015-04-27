'''
    SelScan: A program to calculate EHH-based scans for positive selection in genomes
    https://github.com/szpiech/selscan
'''

__author__ = "tomkinsc@broadinstitute.org"


from Bio import SeqIO
import logging, tools, util.file
import os, os.path, subprocess

tool_version = '1.0.4'
url = 'https://github.com/szpiech/selscan/archive/{ver}.zip'

log = logging.getLogger(__name__)

class SelScanTool(tools.Tool):
    def __init__(self, install_methods = None):
        if install_methods == None:
            install_methods = []
            os_type                 = get_os_type()
            binaryPath              = get_selscan_binary_path( os_type    )
            binaryDir               = get_selscan_binary_path( os_type, full=False )
            
            target_rel_path = 'selscan-{ver}/{binPath}'.format(ver=tool_version, binPath=binaryPath)
            verify_command  = '{dir}/selscan-{ver}/{binPath} --help > /dev/null 2>&1'.format(dir=util.file.get_build_path(), ver=tool_version, binPath=binaryPath) 

            install_methods.append(
                    tools.DownloadPackage(  url.format( ver=tool_version ),
                                            target_rel_path = target_rel_path,
                                            verifycmd       = verify_command,
                                            verifycode      = 1 # selscan returns exit code of 1 for the help text...
                    )
            )

        tools.Tool.__init__(self, install_methods = install_methods)

    def version(self):
        return tool_version


    def execute_ehh():
        toolCmd = [self.install_and_get_path()]
        toolCmd.append("--ehh")

    def execute_ihs():
        toolCmd = [self.install_and_get_path()]
        toolCmd.append("--ihs")

    def execute_xpehh():
        toolCmd = [self.install_and_get_path()]
        toolCmd.append("--xpehh")


    # OLD: to use for fleshing out other execution functions
    def execute(self, inFastas, outFile, localpair, globalpair, preservecase, reorder, 
                outputAsClustal, maxiters, gapOpeningPenalty=None, offset=None, threads=-1, verbose=True):

        inputFileName         = ""
        tempCombinedInputFile = ""

        # get the full paths of input and output files in case the user has specified relative paths
        inputFiles = []
        for f in inFastas:
            inputFiles.append(os.path.abspath(f))
        outFile = os.path.abspath(outFile)

        # ensure that all sequence IDs in each input file are unique 
        # (otherwise the alignment result makes little sense)
        # we can check before combining to localize duplications to a specific file
        for filePath in inputFiles:
            self.__seqIdsAreAllUnique(filePath)

        # if multiple fasta files are specified for input
        if len(inputFiles)>1:
            # combined specified input files into a single temp FASTA file so MAFFT can read them
            tempFileSuffix = ""
            for filePath in inputFiles:
                tempFileSuffix += "__" + os.path.basename(filePath)
            tempCombinedInputFile = util.file.mkstempfname('__combined.{}'.format(tempFileSuffix))
            with open(tempCombinedInputFile, "w") as outfile:
                for f in inputFiles:
                    with open(f, "r") as infile:
                        outfile.write(infile.read())
                #outFile.close()
            inputFileName = tempCombinedInputFile
        # if there is only once file specified, just use it
        else:
            inputFileName = inputFiles[0]

        # check that all sequence IDs in a file are unique
        self.__seqIdsAreAllUnique(inputFileName)

        # change the pwd, since the shell script that comes with mafft depends on the pwd
        # being correct
        pwdBeforeMafft = os.getcwd()
        os.chdir(os.path.dirname(self.install_and_get_path()))

        # build the MAFFT command
        toolCmd = [self.install_and_get_path()]
        toolCmd.append("--auto")
        toolCmd.append("--thread {}".format( max( int(threads), 1 )) )

        if localpair and globalpair:
            raise Exception("Alignment type must be either local or global, not both.")

        if localpair:
            toolCmd.append("--localpair")
            if not maxiters:
                maxiters = 1000
        if globalpair:
            toolCmd.append("--globalpair")
            if not maxiters:
                maxiters = 1000
        if preservecase:
            toolCmd.append("--preservecase")
        if reorder:
            toolCmd.append("--reorder")
        if gapOpeningPenalty:
            toolCmd.append("--op {penalty}".format(penalty=gapOpeningPenalty))
        if offset:
            toolCmd.append("--ep {num}".format(num=offset))
        if not verbose:
            toolCmd.append("--quiet")
        if outputAsClustal:
            toolCmd.append("--clustalout")
        if maxiters:
            toolCmd.append("--maxiterate {iters}".format(iters=maxiters))
        
        toolCmd.append(inputFileName)

        log.debug(' '.join(toolCmd))

        # run the MAFFT alignment
        with open(outFile, 'w') as outf:
            subprocess.check_call(toolCmd, stdout=outf)

        if len(tempCombinedInputFile):
            # remove temp FASTA file
            os.unlink(tempCombinedInputFile)

        # restore pwd
        os.chdir(pwdBeforeMafft)

def get_os_type():
    ''' inspects the system uname and returns a string representing the OS '''

    uname = os.uname()
    if uname[0] == "Darwin":
        return "osx"
    if uname[0] == "Linux":
        return "linux"

def get_selscan_binary_path(os_type, full=True):
    ''' returns the location of the binary relative to the extracted archive, for the given os '''

    selscanPath = "bin/"

    if os_type == "osx":
        selscanPath += "osx/"
    elif os_type == "linux":
        selscanPath += "linux/"
    elif os_type == "win":
        selscanPath += "win/"

    if full:
        selscanPath += "selscan"

        if os_type == "win":
            selscanPath += ".exe"

    return selscanPath




