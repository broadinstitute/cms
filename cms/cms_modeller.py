## experimental. //cd model_target_stats.c? last updated: 06.21.16 vitti@broadinstitute.org
## usage: python cms_modeller.py {bootstrap/grid} -tped {infile} ....

import sys
import argparse

def cms_modeller_parser():
	parser=argparse.ArgumentParser(description="exploratory fitting of demographic models to population genetic data as part of CMS 2.0 package.\nSource code and binaries can be found at https://github.com/broadinstitute/cms>")
	parser.add_argument('subcommand', type=str, default="bootstrap")
	parser.add_argument('-tped',type=str,default=None)
	parser.add_argument('-recom',type=str,default=None)
	parser.add_argument('-regions',type=str,default=None)
	parser.add_argument('-out',type=str,default=None)
	return parser


def main():
	parser = cms_modeller_parser()
	parsed=parser.parse_args()

	if parsed.subcommand == "bootstrap":
		print "bootstrap!"

		if not parsed.tped:
		    print "Need input tped filename"
		    sys.exit(1)
		if not parsed.out:
		    print "Need output filename"
		    sys.exit(1)

	elif parsed.grid:
		print "grid!"
	else:
		print "foobar!"

main()