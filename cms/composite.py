## top-level script for combining scores into composite statistics as part of CMS 2.0.
## last updated: 06.28.16 vitti@broadinstitute.org

import argparse
import sys

def full_parser_composite():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for combining component statistics -- i.e., the final step of the CMS 2.0 pipeline.")
	subparsers = parser.add_subparsers(help="sub-commands")

	bayesian_parser = subparsers.add_parser('bayesian', help='default algorithm and weighting')
	bayesian_parser.add_argument('inputparamfile', action='store', help='file with specifications for input')

	return parser