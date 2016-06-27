"""This is a test; foobarxyz"""

## experimental. //cd model_target_stats.c? last updated: 06.27.16 vitti@broadinstitute.org
## usage: python cms_modeller.py {bootstrap/grid} -tped {infile} ....

import sys
import argparse
import util.cmd

__author__   = "vitti@broadinstitute.org"
__commands__ = []



def parser_bootstrap(parser):
	if parser == None:
		parser=argparse.ArgumentParser(description="exploratory fitting of demographic models to population genetic data as part of CMS 2.0 package.\nSource code and binaries can be found at https://github.com/broadinstitute/cms>")
	#parser.help = '''test test'''
	#parser.add_argument('subcommand', type=str, default="bootstrap", help="also a test")
	parser.add_argument('-tped',type=str,default=None)
	parser.add_argument('-recom',type=str,default=None)
	parser.add_argument('-regions',type=str,default=None)
	parser.add_argument('-out',type=str,default=None)
	#subparsers = parser.add_subparsers(title='subcommands', dest='command')
	return parser
__commands__.append(('bootstrap', parser_bootstrap))

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)

def main():

	util.cmd.main_argparse(__commands__, __doc__)

	#parser = cms_modeller_parser()
	#parsed=parser.parse_args()

	#if parsed.subcommand == "bootstrap":
	#	print "bootstrap!"

	#	if not parsed.tped:
	#	    print "Need input tped filename"
	#	    sys.exit(1)
	#	if not parsed.out:
	#	    print "Need output filename"
	#	    sys.exit(1)

	#elif parsed.grid:
	#	print "grid!"
	#else:
	#	print "foobar!"

main()