## experimental. last updated: 06.27.16 vitti@broadinstitute.org
## usage: python cms_modeller.py {subcommand: bootstrap / grid / }

import sys
import argparse
import util.cmd

def full_parser_cms_modeller():
	parser=argparse.ArgumentParser(description="exploratory fitting of demographic models to population genetic data")
	subparsers = parser.add_subparsers(help="sub-commands")

	bootstrap_parser = subparsers.add_parser('bootstrap', help='perform bootstrap estimates of population summary statistics for model target values')
	bootstrap_parser.add_argument('tped', action='store', help='input tped file')

	grid_parser = subparsers.add_parser('grid', help='run grid search')
	grid_parser.add_argument('inputparamfile', action='store', help='file with specifications of grid search')

	return parser

def parser_bootstrap():
	
	parser.add_argument('subcommand', type=str, default="bootstrap", help="subcommand: bootstrap / grid / gradient ...")
	parser.add_argument('-tped',type=str,default=None)

	parser.add_argument("--verbose", help="increase output verbosity",
                    action="store_true")

	return parser


#args = parser.parse_args()
#if args.verbose:
#    print("verbosity turned on")