#! /usr/bin/env python3

import sys
import argparse
import subprocess
import os
from sam_to_psl import sam_to_psl
from bam2Bed12 import bam2Bed12


def parseCommandLine():
	parser = argparse.ArgumentParser(description='flair-align parse options',
		usage='python flair.py -g genome.fa -r <reads.fq>|<reads.fa> [options]')
	parser.add_argument('align')
	required = parser.add_argument_group('required named arguments')
	atleastone = parser.add_argument_group('Either one of the following arguments is required')
	required.add_argument('-r', '--reads', action='store', dest='r',
		nargs='+', type=str, required=True, help='FastA/FastQ files of raw reads')
	atleastone.add_argument('-g', '--genome', action='store', dest='g',
		type=str, help='FastA of reference genome, can be minimap2 indexed')
	atleastone.add_argument('--mm_index', action='store', dest='mm_index', type=str, default='',
		help='minimap2 index .mmi file')
	parser.add_argument('-o', '--output', action='store', dest='o', default='flair.aligned',
		help='output file name base (default: flair.aligned)')
	parser.add_argument('-t', '--threads', type=str,
		action='store', dest='t', default='4', help='minimap2 number of threads (4)')
	parser.add_argument('--junction_bed', action='store', dest='junction_bed', default='',
		help='annotated isoforms/junctions bed file for splice site-guided minimap2 genomic alignment')
	parser.add_argument('--pychopper', type=str, default='', action='store', dest='pychopper',
		help='specify cdna_classifier.py here to trim reads prior to aligning')
	parser.add_argument('-m', '--minimap2', type=str, default='minimap2',
		action='store', dest='m', help='path to minimap2 if not in $PATH')
	parser.add_argument('--nvrna', action='store_true', dest='n', default=False,
		help='specify this flag to use native-RNA specific alignment parameters for minimap2')
	parser.add_argument('-sam', '--samtools', action='store', dest='sam', default='samtools',
		help='samtools executable path if not in $PATH')
	parser.add_argument('-c', '--chromsizes', type=str, action='store', dest='c', default='',
		help='''chromosome sizes tab-separated file, used for converting sam to genome-browser
		compatible psl file''')
	parser.add_argument('--psl', action='store_true', dest='p',
		help='also output sam-converted psl')
	parser.add_argument('-v1.3', '--version1.3', action='store_true', dest='v',
		help='specify if samtools version 1.3+')
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1,
		help='minimum MAPQ of read alignment to the genome (1)')
	parser.add_argument('-N', type=int, action='store', dest='N', default=0,
		help='retain at most INT secondary alignments from minimap2 alignment (0)')
	parser.add_argument('--quiet', default=False, action='store_true', dest='quiet',
			help='''Suppress progress statements from being printed''')

	args, unknown = parser.parse_known_args()
	if unknown and not args.quiet:
		sys.stderr.write('Align unrecognized arguments: {}\n'.format(' '.join(unknown)))

	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'
	for i in range(len(args.r)):
		if not os.path.exists(args.r[i]):
			sys.stderr.write('Check that read file {} exists\n'.format(args.r[i]))
			return 1
		if args.pychopper:
			subprocess.check_call([args.pychopper, '-t', args.t, '-r', args.o+'.'+args.r[i]+'.pychopper_report.pdf',
				'-u', args.o+'.'+args.r[i]+'.unclassified.fastq', '-w', args.o+'.'+args.r[i]+'.rescued.fastq',
				args.r, args.o+'.'+args.r[i]+'.trimmed.fastq'])
			args.r[i] = args.o+'.'+args.r[i]+'.trimmed.fastq'
	return args


def align():
	args = parseCommandLine()
	mm2_command = [args.m, '-ax', 'splice', '-t', args.t, args.g]+args.r
	if args.mm_index:
		mm2_command[5] = args.mm_index
	if args.n:
		mm2_command[3:3] = ['-uf', '-k14']
	if args.junction_bed:
		mm2_command[3:3] = ['--junc-bed', args.junction_bed]
	if str(args.N) != '0':
		mm2_command[3:3] = ['-N', str(args.N)]
	else:
		mm2_command[3:3] = ['--secondary=no']

	try:
		if args.quiet:
			if subprocess.call(mm2_command, stdout=open(args.o+'.sam', 'w'),
			stderr=open(args.o+'.mm2_stderr', 'w')):
				return 1
		elif subprocess.call(mm2_command, stdout=open(args.o+'.sam', 'w')):
			return 1
	except:
		sys.stderr.write('Possible minimap2 error, specify executable path with -m?\n')
		if args.quiet:
			sys.stderr.write('Check {}\n'.format(args.o+'.mm2_stderr'))
		return 1

	if args.quality != 0:
		if subprocess.call([args.sam, 'view', '-q', str(args.quality), '-h', '-S', args.o+'.sam'],
		stdout=open(args.o+'.q.sam', 'w'), stderr=open(args.o+'.samtools_stderr', 'w')):
			sys.stderr.write('Possible issue with samtools, see {}\n'.format(args.o+'.samtools_stderr'))
			return 1

		subprocess.check_call(['mv', args.o+'.q.sam', args.o+'.sam'])
		subprocess.check_call(['rm', args.o+'.samtools_stderr'])

	if args.p:
		try:
			sam_to_psl(args.o+'.sam', args.o+'.psl', args.c)
		except:
			return 1

	if not args.v:  # samtools version is < 1.3 or unspecified --> detect version; some versions don't have --version
		ver = subprocess.Popen([args.sam], stderr=subprocess.PIPE, universal_newlines=True)
		for line in ver.stderr:
			if 'Version:' in line:
				v = line.rstrip()[line.find('Version:')+9:line.find('Version:')+13]
				v = v[:-1] if v[-1] == '.' else v
				try:
					if float(v) >= 1.3:
						if not args.quiet:
							sys.stderr.write('Samtools version >= 1.3 detected\n')
						args.v = True
						break
				except:
					if not args.quiet:
						sys.stderr.write('Could not detect samtools version, assuming < 1.3\n')

	if args.v:  # samtools verison 1.3+
		if subprocess.call([args.sam, 'sort', '-@', args.t, args.o+'.sam', '-o', args.o+'.bam'],
			stderr=open(args.o+'.bam.stderr', 'w')):
			sys.stderr.write('Samtools issue with sorting minimap2 sam\n')
			return 1
		subprocess.check_call(['rm', args.o+'.bam.stderr'])
	else:
		if subprocess.call([args.sam, 'view', '-h', '-Sb', '-@', args.t, args.o+'.sam'],
				stdout=open(args.o+'.unsorted.bam', 'w')):
			sys.stderr.write('Possible issue with samtools executable\n')
			return 1
		if subprocess.call([args.sam, 'sort', '-@', args.t, args.o+'.unsorted.bam', args.o],
				stderr=open(args.o+'.unsorted.bam.stderr', 'w')):
			sys.stderr.write('If using samtools v1.3+, please specify -v1.3 argument\n')
			return 1
		subprocess.check_call(['rm', args.o+'.unsorted.bam', args.o+'.unsorted.bam.stderr'])

	subprocess.check_call([args.sam, 'index', args.o+'.bam'])

	keep_supplementary = True if args.N != 0 else False
	bam2Bed12(args.o+'.bam', args.o+'.bed', keep_supplementary)
	return args.o+'.bed'


if __name__ == "__main__":
	align()
