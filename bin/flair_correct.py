#! /usr/bin/env python3

import sys
import argparse
import subprocess
import os
import tempfile
import glob
from ssCorrect import *
from bed_to_psl import bed_to_psl
from utils import handle_prog_errors

def parseCommandLine():
	parser = argparse.ArgumentParser(description='flair-correct parse options',
		usage='python flair.py correct -q query.bed12 [-f annotation.gtf]v[-j introns.tab] -g genome.fa [options]')
	parser.add_argument('correct')
	required = parser.add_argument_group('required named arguments')
	atleastone = parser.add_argument_group('at least one of the following arguments is required')
	required.add_argument('-q', '--query', type=str, default='', required=True,
			action='store', dest='q', help='uncorrected bed12 file')
	required.add_argument('-g', '--genome', action='store', dest='g',
		type=str, required=True, help='FastA of reference genome')
	atleastone.add_argument('-j', '--shortread', action='store', dest='j', type=str, default=None,
		help='bed format splice junctions from short-read sequencing')
	atleastone.add_argument('-f', '--gtf', default=None,
		action='store', dest='f', help='GTF annotation file')
	parser.add_argument('-c', '--chromsizes', type=str, action='store',
		dest='c', default='', help='chromosome sizes tab-separated file, specify to get .psl output (in addition to .bed)')
	parser.add_argument('--nvrna', action='store_true', dest='n', default=False,
		help='''specify this flag to keep the strand of a read consistent after correction''')
	parser.add_argument('-t', '--threads', type=int, action='store', dest='t', default=4,
		help='splice site correction script number of threads (4)')
	parser.add_argument('-w', '--ss_window', type=int, action='store', dest='w', default=10,
		help='window size for correcting splice sites (W=10)')
	parser.add_argument('-o', '--output',
		action='store', dest='o', default='flair', help='output name base (default: flair)')
	parser.add_argument('--print_check',
		action='store_true', dest='p', default=False, help='Print err.txt with step checking.')
	args, unknown = parser.parse_known_args()
	if unknown:
		sys.stderr.write('Correct unrecognized arguments: {}\n'.format(' '.join(unknown)))
		if not aligned_reads:
			return 1

	if not args.j and not args.f:
		sys.stderr.write('Please specify at least one of the -f or -j arguments for correction\n')
		return 1
	return args

def correct(aligned_reads=False):
	args = parseCommandLine()
	if aligned_reads:
		args.q = aligned.reads
	elif not args.q:
		sys.stderr.write('Please specify --query\n')
		return 1
	resolveStrand = False if args.n else True
	try:
		ssCorrect(bed=args.q, gtf=args.f, otherJuncs=args.j, wiggle=args.w, threads=args.t, outFile=args.o, 
			keepTemp=False, resolveStrand=resolveStrand, tempDirName=None, genomeFasta=args.g, verbose=True, printErr=args.p)
	except:
		sys.stderr.write('Correction command did not exit with success status\n')
		return 1

	if args.c:
		try:
			bed_to_psl(args.c, args.o+'_all_corrected.bed', args.o+'_all_corrected.psl')
		except Exception as ex:
			sys.stderr.write('bed_to_psl command did not exit with success status\n')
			handle_prog_errors(ex, True)

#		except:
#			return 1
#	if args.c and subprocess.call([sys.executable, path+'bin/bed_to_psl.py', args.c, args.o+'_all_corrected.bed',
#		args.o+'_all_corrected.psl']):
#		return 1

	return args.o+'_all_corrected.bed'

if __name__ == "__main__":
	correct()


