## experimental; one of CMS 2.0 suite top-level scripts. last updated: 06.27.16 vitti@broadinstitute.org
## usage: python likes_from_model.py {}

## experimental. //cd model_target_stats.c? last updated: 06.27.16 vitti@broadinstitute.org
## usage: python cms_modeller.py {bootstrap/grid} -tped {infile} ....

import sys
import argparse
import util.cmd

__author__   = "vitti@broadinstitute.org"
__commands__ = []



def parser_likes():
	parser=argparse.ArgumentParser(description="generation of likelihood tables")
	parser.help = '''test test wahoo'''
	parser.add_argument('subcommand', type=str, default="bootstrap", help="this is required")
	#parser.add_argument('-tped',type=str,default=None)
	#parser.add_argument('-recom',type=str,default=None)
	#parser.add_argument('-regions',type=str,default=None)
	#parser.add_argument('-out',type=str,default=None)
	subparsers = parser.add_subparsers(title='subcommands', dest='command')
	return parser