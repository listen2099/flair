#! /usr/bin/env python3

import sys
import argparse
import subprocess
import os
import tempfile
import glob
import re
import shutil
import uuid
import pybedtools
from multiprocessing import Pool
from tqdm import *
from bed_to_psl import bed_to_psl
from utils import handle_prog_errors
from ssPrep import ssPrep

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

def addOtherJuncs(juncs, bedJuncs, chromosomes, fa, known, printErr, printErrFname, verbose):

    lineNum = 0
    if verbose: sys.stderr.write("Step 2/5: Processing additional junction file  %s ..." % (bedJuncs))
    cols = None

    with open(bedJuncs) as l:
        for num,ll in enumerate(l,0):
            cols = ll.rstrip().split()
            if num >10:
                break


    # guess what kind of bedFile
    if cols is None:
        raise Exception("Empty junctions BED file, not supported")

    if cols[-1] == "+" or cols[-1] == "-":
        # normal bed
        reverseSS = "-"
        strandCol = -1
        starOffset = 0

    elif len(cols) == 12:
        # bed12
        bedType   = "bed12"
        raise Exception("Bed12 not currently supported for other_juncs.bed. Please convert to bed6.")

    elif cols[3] == "0" or cols[3] == "1" or cols[3] == "2":
        # star junc.tab
        reverseSS = "2"
        strandCol = 3
        starOffset = 1

    else:
        raise Exception("Cannot find strand info for %s. Is this bed6 or STAR_juncs.tab file?" % bedJuncs)

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Adding other juncs, assuming file is %s" % "bed6" if strandCol == -1 else "STAR", file=fo)

    tempJuncs = list()
    with open(bedJuncs,'r') as bedLines:
        for line in bedLines:
            cols = line.rstrip().split()
            chrom, c1, c2, strand = cols[0], int(cols[1])-starOffset, int(cols[2]), cols[strandCol]

            if chrom not in juncs:
                juncs[chrom] = dict()

            if c2-c1 < 5:
                continue

            if starOffset:
                if strand == "1": strand = "+"
                elif strand == "2": strand = "-"
                else: continue

            chromosomes.add(chrom)
            key = (c1, c2, strand)
            if key in juncs[chrom]:
                juncs[chrom][key] = "both"
                continue
            tempJuncs.append((chrom,c1,c2,"%s,%s,%s,%s" % (chrom,c1,c2,strand),0,strand))

    try:
        btJuncs = pybedtools.BedTool(tempJuncs)
        dinucSeq = btJuncs.sequence(fi=fa, s=True, tab=True, name=True)
        with open(dinucSeq.seqfn) as fileObj:
            for i in fileObj:
                header,seq = i.rstrip().split()
                chrom,c1,c2,strand = header.split(",")
                c1,c2 = int(c1),int(c2)
                if "+" in strand:
                    strand = strand[strand.rfind('+')]
                elif "-" in strand:
                    strand = strand[strand.rfind('-')]
                key = (c1,c2, strand)
                known1,known2 = known.get((chrom,c1),None),known.get((chrom,c2),None)

                if known1 != None:
                    if known1 != strand:
                        continue
                    else:
                        pass
                else:
                    pass

                if known2 != None:
                    if known2 != strand:
                        continue
                    else:
                        pass
                else:
                    pass

                fivePrime = seq[:2]
                if key not in juncs[chrom]:
                    juncs[chrom][key] = "sr"
                elif fivePrime == "GT":
                    juncs[chrom][key] = "sr"

    except Exception as e:
        print(e,"Splice site motif filtering failed. Check pybedtools and bedtools is properly install and in $PATH",file=sys.stderr)
        sys.exit(1)

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** GTF Juncs + other juncs now total %s juncs from %s chromosomes." % (sum([len(x)for x in juncs.values()]), len(list(juncs.keys()))), file=fo)

    return juncs, chromosomes

def gtfToSSBed(infile, knownSS, printErr, printErrFname, verbose):
    ''' Convenience function, reformats GTF to bed'''

    # First: get all exons per transcript.
    exons = dict()
    chromosomes = set()
    with open(infile,'r') as lines:
        for l in lines:
            if l[0] == "#": # skip header lines
                continue

            cols = l.split("\t")

            if "exon" == cols[2]:

                # -1 for 1 to 0 based conversion
                chrom, c1, c2, strand =  cols[0], int(cols[3])-1, int(cols[4]), cols[6]
                chromosomes.add(chrom)
                #txn info is in the SECOND position of the shoutout column
                txn = re.search('transcript_id "([^\"]+)"', l).group(1)#
                #cols[-1].split(";")[1].split()[-1].replace('"','')

                key = (chrom, txn, strand)

                if key not in exons:
                    exons[key] = list()
                exons[key].append(c1)
                exons[key].append(c2)
    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Read GTF. Got %s transcripts" % len(list(exons.keys())), file=fo)
            print("** Getting introns...Read GTF", file=fo)

    # Second: get junction and splice sites from transcript exons.
    txnList = list(exons.keys())
    juncs = dict()

    for exonInfo in tqdm(txnList, total=len(txnList), desc="Step 1/5: Splitting junctions from GTF by chromosome", dynamic_ncols=True, position=1) if verbose else txnList:
        chrom, txn, strand = exonInfo

        if chrom not in juncs:
            juncs[chrom] = dict()

        coords = list(exons[exonInfo])

        # assume lowest and highest as TSS and TES, and remove them
        coords.sort()
        coords = coords[1:-1]

        # Coords is list of exons, so a list less than 2 is a single exon gene.
        if len(coords)<2: continue

        for pos in range(0,len(coords)-1,2):
            c1 = coords[pos]
            c2 = coords[pos+1]

            if abs(c2 - c1) <= 5:
                continue

            juncs[chrom][(c1,c2,strand)] = "gtf"
            knownSS[(chrom, c1)] = strand
            knownSS[(chrom, c2)] = strand

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Created %s juncs from %s chromosomes." % (sum([len(x)for x in juncs.values()]), len(list(juncs.keys()))), file=fo)
    return juncs, chromosomes, knownSS


def runCMD(x):

    tDir, prefix,juncs,reads, rs, f, err, errFname = x
    resolveStrand = True if rs else False
    checkFname = errFname if err else False
    try:
        ssPrep(bed=reads, knownJuncs=juncs, fa=f, wiggle=15, out=prefix, resolveStrand=resolveStrand, 
		workingDir=tDir, checkFname=checkFname)
    except Exception as ex:
        handle_prog_errors(ex, True)
        sys.stderr.write('ssPrep command did not exit with success status\n')
#    except:
        return 1
    return 0


def ssCorrect(bed, gtf, otherJuncs, wiggle, threads, outFile, keepTemp, resolveStrand, tempDirName, genomeFasta, verbose, printErr):
    globals()['printErr'] = printErr
    if os.path.isfile(genomeFasta+".fai"):
        pass
    else:
        testString =  """
            chrX 1    100   feature1  0 +
        """
        test = pybedtools.BedTool(testString, from_string=True)
        a = test.sequence(fi=genomeFasta)

    # make temp dir for dumping
    if tempDirName == None:
        tempDirName = "tmp_%s" % str(uuid.uuid4())
    try:
        current_directory = os.getcwd()
        tempDir = os.path.join(current_directory, tempDirName)
        os.mkdir(tempDir)
    except OSError:
        print ("Creation of the directory %s failed" % tempDirName)
        sys.exit(1)

    printErrFname = "err_%s.txt" % tempDirName

    # Convert gtf to bed and split by cromosome.
    juncs, chromosomes, knownSS  = dict(), set(), dict() # initialize juncs for adding to db
    if gtf != None: juncs, chromosomes, knownSS = gtfToSSBed(gtf, knownSS, printErr, printErrFname, verbose)

    # Do the same for the other juncs file.
    if otherJuncs != None: juncs, chromosomes = addOtherJuncs(juncs, otherJuncs, chromosomes, genomeFasta, knownSS, printErr, printErrFname, verbose)
    knownSS = dict()

    # added to allow annotations not to be used.
    if len(list(juncs.keys()))<1:
        print("No junctions from GTF or junctionsBed to correct with. Exiting...", file=sys.stderr)
        sys.exit(1)

    annotations = dict()
    for chrom, data in tqdm(juncs.items(), desc="Step 3/5: Preparing annotated junctions to use for correction", total=len(list(juncs.keys())), dynamic_ncols=True, position=1) if verbose else juncs.items():
        annotations[chrom] = os.path.join(tempDir,"%s_known_juncs.bed" % chrom)
        with open(os.path.join(tempDir,"%s_known_juncs.bed" % chrom),"w") as out:
            sortedData = sorted(list(data.keys()), key=lambda item: item[0])
            for k in sortedData:
                annotation = data[k]
                c1, c2, strand = k
                print(chrom,c1,c2,annotation,".",strand, sep="\t", file=out)

    sortedData = None
    skippedChroms = set()
    readDict = dict()
    with open(bed) as lines:
        outDict = dict()
        for line in tqdm(lines, desc="Step 4/5: Preparing reads for correction", dynamic_ncols=True, position=1) if verbose else lines:
            cols  = line.rstrip().split()
            chrom = cols[0]
            if chrom not in chromosomes:
                if chrom not in skippedChroms:
                    skippedChroms.add(chrom)
                    #if verbose: tqdm.write("Reference sequence not found in annotations, skipping: %s" % (chrom), file=sys.stderr)
                    continue
            else:
                if chrom not in outDict:
                    readDict[chrom] = os.path.join(tempDir,"%s_temp_reads.bed" % chrom)
                    outDict[chrom] = open(os.path.join(tempDir,"%s_temp_reads.bed" % chrom),'w')
                print(line.rstrip(),file=outDict[chrom])

    cmds = list()
    for chrom in readDict:
        juncs = annotations[chrom]
        reads = readDict[chrom]

        outDict[chrom].close()

        cmds.append((tempDir, chrom, juncs,reads, resolveStrand, genomeFasta, printErr, printErrFname))

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Prepared correct commands for %s read files" % len(cmds), file=fo)


    juncs = None
    annotations = None
    p = Pool(threads)
    childErrs = set()
    for i in tqdm(p.imap(runCMD, cmds), total=len(cmds), desc="Step 5/5: Correcting Splice Sites", dynamic_ncols=True,position=1) if verbose else p.imap(runCMD,cmds):
        childErrs.add(i)
    if len(childErrs)>1:
        print(childErrs,file=sys.stderr)
        sys.exit(1)


    with open("%s_all_inconsistent.bed" % outFile,'wb') as inconsistent:
        for chrom in readDict:
            with open(os.path.join(tempDir, "%s_inconsistent.bed" % chrom),'rb') as fd:
                shutil.copyfileobj(fd, inconsistent, 1024*1024*10)


    with open("%s_all_corrected.bed" % outFile,'wb') as corrected:
        for chrom in readDict:
            with open(os.path.join(tempDir, "%s_corrected.bed" % chrom),'rb') as fd:
                shutil.copyfileobj(fd, corrected, 1024*1024*10)
    if keepTemp:
        pass
    else:
        try:
            shutil.rmtree(tempDir)
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror), file=sys.stderr)

def correct(aligned_reads=False):
	args = parseCommandLine()
	if aligned_reads:
		args.q = aligned.reads
	elif not args.q:
		sys.stderr.write('Please specify --query\n')
		return 1
	resolveStrand = False if args.n else True
	verbose = False if args.q else True
	try:
		ssCorrect(bed=args.q, gtf=args.f, otherJuncs=args.j, wiggle=args.w, threads=args.t, outFile=args.o, 
			keepTemp=False, resolveStrand=resolveStrand, tempDirName=None, genomeFasta=args.g, verbose=verbose, printErr=args.p)
	except Exception as ex:
		handle_prog_errors(ex, True)
		sys.stderr.write('Correction command did not exit with success status\n')
	#except:
	#	return 1

	if args.c:
		try:
			bed_to_psl(args.c, args.o+'_all_corrected.bed', args.o+'_all_corrected.psl')
		except Exception as ex:
			sys.stderr.write('bed_to_psl command did not exit with success status\n')
			handle_prog_errors(ex, True)

#		except:
#			return 1

	return args.o+'_all_corrected.bed'

if __name__ == "__main__":
	correct()

