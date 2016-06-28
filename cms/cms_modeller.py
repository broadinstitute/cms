## top-level script for demographic modeling as part of CMS 2.0.
## last updated: 06.28.16 vitti@broadinstitute.org

import argparse
import sys

def full_parser_cms_modeller():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for exploratory fitting of demographic models to population genetic data.")
	subparsers = parser.add_subparsers(help="sub-commands")

	######################
	## CALCULATE TARGET ##
	#######################

	target_stats_parser = subparsers.add_parser('target_stats', help='perform per-site(/per-site-pair) calculations of population summary statistics for model target values')
	target_stats_parser.add_argument('tped', action='store', help='input tped file')

	bootstrap_parser = subparsers.add_parser('bootstrap', help='perform bootstrap estimates of population summary statistics in order to finalize model target values')

	######################
	## VISUALIZE MODEL  ##
	######################

	point_parser = subparsers.add_parser('point', help='run simulates of a point in parameter-space')

	#########################
	## FIT MODEL TO TARGET ##
	#########################

	grid_parser = subparsers.add_parser('grid', help='run grid search')
	grid_parser.add_argument('inputparamfile', action='store', help='file with specifications of grid search')

	optimize_parser = subparsers.add_parser('optimize', help='run optimization algorithm to fit model parameters')
	optimize_parser.add_argument('inputparamfile', action='store', help='file with specifications of optimization search to run')

	return parser

#parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
# ...
#args = parser.parse_args()
#if args.verbose:
#    print("verbosity turned on")