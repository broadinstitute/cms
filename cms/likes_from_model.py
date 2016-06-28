## top-level script for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 06.28.16 vitti@broadinstitute.org

import argparse
import sys

def full_parser_likes_from_model():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for generating probability distributions for component scores from pre-specified demographic model(s).")
	subparsers = parser.add_subparsers(help="sub-commands")

	run_sims_parser = subparsers.add_parser('run_sims', help='run simulations')
	run_sims_parser.add_argument('inputparamfile', action='store', help='file with specifications for input')

	return parser