## top-level script for combining scores into composite statistics as part of CMS 2.0.
## last updated: 06.28.16 vitti@broadinstitute.org

import argparse
import sys

def full_parser_composite():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for combining component statistics -- i.e., the final step of the CMS 2.0 pipeline.")
	subparsers = parser.add_subparsers(help="sub-commands")

	#################
	## GENOME-WIDE ##
	#################

	bayesian_gw_parser = subparsers.add_parser('bayesian_gw', help='default algorithm and weighting, genome-wide')
	bayesian_gw_parser.add_argument('inputparamfile', action='store', help='file with specifications for input')

	###################
	## WITHIN REGION ##
	###################

	bayesian_region_parser = subparsers.add_parser('bayesian_region', help='default algorithm and weighting, within-region')

	ml_region_parser = subparsers.add_parser('ml_region', help='machine learning algorithm (within-region)')

	return parser