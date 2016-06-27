## experimental. //cd model_target_stats.c? last updated: 06.27.16 vitti@broadinstitute.org
## usage: python cms_modeller.py {bootstrap/grid} -tped {infile} ....

import sys
import argparse
import util.cmd

__author__   = "vitti@broadinstitute.org"
__commands__ = []



def parser_bootstrap():
	parser=argparse.ArgumentParser(description="exploratory fitting of demographic models to population genetic data")
	#parser.help = '''test test'''
	parser.add_argument('subcommand', type=str, default="bootstrap", help="this is required")
	parser.add_argument('-tped',type=str,default=None)

	subparsers = parser.add_subparsers()
	subparser.add_argument('ref', type=str, help='foo1 help')

	parser.add_argument('-recom',type=str,default=None)
	parser.add_argument('-regions',type=str,default=None)
	parser.add_argument('-out',type=str,default=None)
	subparsers = parser.add_subparsers(title='subcommands', dest='command')
	return parser