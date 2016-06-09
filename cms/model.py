## top-level script for CMS 2.0 demographic modeling (sub-)pipeline
## last updated 06.07.16	vitti@broadinstitute.org

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("subcommand", help="")
args = parser.parse_args()
subcommand = args.subcommand
print subcommand + "?"