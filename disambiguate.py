#!/usr/bin/env python
"""
This is the main function to call for disambiguating between a human and 
mouse BAM files that have alignments from the same source of fastq files.
It is part of the explant RNA/DNA-Seq workflow where an informatics 
approach is used to distinguish between human and mouse RNA/DNA reads.
	
For reads that have aligned to both organisms, the functionality is based on 
comparing quality scores from either Tophat of BWA (under development). Read
name is used to collect all alignments for both mates (_1 and _2) and 
compared between human and mouse alignments.

For tophat (default, can be changed using option -a), the sum of the flags XO,
NM and NH is evaluated and the lowest sum wins the paired end reads. For equal
scores, the reads are assigned as ambiguous.

The alternative algorithm (bwa) disambiguates (for aligned reads) by tags AS 
(alignment score, higher better), NM (edit distance, lower better) and XS 
(suboptimal alignment score, higher better), by first looking at AS, then 
NM and finally XS.


Usage: 
disambiguate.py [-h <human.bam>] [-m <mouse.bam>] [-o <outputdir>] [-d] [-i <intermediatedir>] [-s <samplenameprefix>] [-a ALGORITHM]

If <samplenameprefix> is provided it will be used in the output filenames, 
otherwise the output filenames will be named after <mouse.bam> and <human.bam>

N.B. that Tophat2 always renames its output as accepted_hits.bam 
regardless of the fastq filenames and therefore it is a very good 
idea to provide a <samplenameprefix>

Options:
  -h <human.bam>           Input human BAM file location (REQUIRED)
  -m <mouse.bam>           Input mouse BAM file location (REQUIRED)
  -o <outputdir>           Output directory. If the directory does not exist it
                           will be generated. If not provided, 'disambres/' 
                           will be used.
  -d                       Disable BAM file sorting by name. Using this option 
                           it is possible to serve BAM files that have already
                           been name sorted using samtools
  -i <intermediatedir>     Location where to store intermediate name sorted BAM
                           files (if name sorting was not disabled). If not 
                           provided 'intermfiles/' will be used
  -s <samplenameprefix>    A prefix (e.g. sample name) to use for the output 
                           BAM files. If not provided, the human BAM file 
                           prefix will be used.
  -a ALGORITHM             tophat (default) or bwa. Note that for bwa 
                           alignments the tag provided in the BAM files are 
                           different to those from tophat.
		
Code by Miika Ahdesmaki July-August 2013.
"""


from __future__ import print_function
import sys, getopt, re, time, pysam
from array import array
from os import path, makedirs
from docopt import docopt

# "natural comparison" for strings
def nat_cmp(a, b):
	convert = lambda text: int(text) if text.isdigit() else text # lambda function to convert text to int if number present
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] # split string to piecewise strings and string numbers
	#return cmp(alphanum_key(a), alphanum_key(b)) # use internal cmp to compare piecewise strings and numbers
	return (alphanum_key(a) > alphanum_key(b))-(alphanum_key(a) < alphanum_key(b))

# read reads into a list object for as long as the read qname is constant (sorted file). Return the first read with new qname or None
def read_next_reads(fileobject, listobject):
	qnamediff = False
	while not qnamediff:
		try:
			myRead=fileobject.next()
		except StopIteration:
			#print("5")
			return None # return None as the name of the new reads (i.e. no more new reads)
		if nat_cmp(myRead.qname, listobject[0].qname)==0:
			listobject.append(myRead)
		else:
			qnamediff = True
	return myRead # this is the first read with a new qname

# disambiguate between two lists of reads
def disambiguate(humanlist, mouselist, disambalgo):
	if disambalgo == 'tophat':
		dv = 2**13 # a high quality score to replace missing quality scores (no real quality score should be this high)
		sa = array('i',(dv for i in range(0,4))) # score array, with [human_1_QS, human_2_QS, mouse_1_QS, mouse_2_QS]
		for read in humanlist:
			if 0x4&read.flag: # flag 0x4 means unaligned
				continue
			QScore = read.opt('XO') + read.opt('NM') + read.opt('NH')
		   # directionality (_1 or _2)
			d12 = 0 if 0x40&read.flag else 1
			if sa[d12]>QScore:
				sa[d12]=QScore # update to lowest (i.e. 'best') quality score
		for read in mouselist:
			if 0x4&read.flag: # flag 0x4 means unaligned
				continue
			QScore = read.opt('XO') + read.opt('NM') + read.opt('NH')
		   # directionality (_1 or _2)
			d12 = 2 if 0x40&read.flag else 3
			if sa[d12]>QScore:
				sa[d12]=QScore # update to lowest (i.e. 'best') quality score
		if min(sa[0:2])==min(sa[2:4]) and max(sa[0:2])==max(sa[2:4]): # ambiguous
			return 0
		elif min(sa[0:2]) < min(sa[2:4]) or min(sa[0:2]) == min(sa[2:4]) and max(sa[0:2]) < max(sa[2:4]):
			# assign to human
			return 1
		else:
			# assign to mouse
			return -1
	elif disambalgo == 'bwa':
		dv = -2^13 # default value, low
		bwatags = ['AS','NM','XS'] # in order of importance (compared sequentially, not as a sum as for tophat)
		bwatagsigns = [1,-1,1] # for AS and XS higher is better. for NM lower is better, thus multiply by -1
		AS = list()
		for x in range(0, len(bwatagsigns)):
			AS.append(array('i',(dv for i in range(0,4)))) # alignment score array, with [human_1_Score, human_2_Score, mouse_1_Score, mouse_2_Score]
		#		
		for read in humanlist:
			if 0x4&read.flag: # flag 0x4 means unaligned
				continue
			# directionality (_1 or _2)
			d12 = 0 if 0x40&read.flag else 1
			for x in range(0, len(bwatagsigns)):
				QScore = bwatagsigns[x]*read.opt(bwatags[x])
				if AS[x][d12]<QScore:
					AS[x][d12]=QScore # update to highest (i.e. 'best') quality score
		#
		for read in mouselist:
			if 0x4&read.flag: # flag 0x4 means unaligned
				continue
		   # directionality (_1 or _2)
			d12 = 2 if 0x40&read.flag else 3
			for x in range(0, len(bwatagsigns)):
				QScore = bwatagsigns[x]*read.opt(bwatags[x])
				if AS[x][d12]<QScore:
					AS[x][d12]=QScore # update to highest (i.e. 'best') quality score
		#
		for x in range(0, len(bwatagsigns)):
			if max(AS[x][0:2]) > max(AS[x][2:4]) or max(AS[x][0:2]) == max(AS[x][2:4]) and min(AS[x][0:2]) > min(AS[x][2:4]):
				# assign to human
				return 1
			elif max(AS[x][0:2]) < max(AS[x][2:4]) or max(AS[x][0:2]) == max(AS[x][2:4]) and min(AS[x][0:2]) < min(AS[x][2:4]):
				# assign to mouse
				return -1
		return 0 # ambiguous
	else:
		print("Not implemented yet")
		sys.exit(2)		


#code
def main(argv):
	numhum = nummou = numamb = 0	
	starttime = time.clock()
	# parse inputs
	humanfile = argv['-h'] if argv['-h'] is not None else ''
	mousefile = argv['-m'] if argv['-h'] is not None else ''
	samplenameprefix = argv['-s'] if argv['-s'] is not None else ''
	outputdir = argv['-o'] if argv['-o'] is not None else 'disambres/'
	intermdir = argv['-i'] if argv['-i'] is not None else 'intermfiles/'
	disablesort = True if argv['-d'] else False
	disambalgo = argv['-a'] if argv['-a'] is not None else 'tophat'
	supportedalgorithms = set()
	supportedalgorithms.add('tophat')
	supportedalgorithms.add('bwa')
	
	# check existence of input BAM files
	if len(humanfile) < 1 or len(mousefile) < 1 or not path.isfile(humanfile) or not not path.isfile(mousefile):
		print(__doc__)
		sys.stderr.write("\nERROR in disambiguate.py: Two existing input BAM files must be specified using options -h and -m\n")
		sys.exit(2)
	if len(samplenameprefix) < 1:
		humanprefix = path.basename(humanfile.replace(".bam",""))
		mouseprefix = path.basename(mousefile.replace(".bam",""))
	else:
		humanprefix = samplenameprefix
		mouseprefix = samplenameprefix
	
	samplenameprefix = None # clear variable
	if disambalgo not in supportedalgorithms:
		print(disambalgo+" is not a supported disambiguation scheme at the moment.")
		sys.exit(2)
	
	if disablesort:
		humanfilesorted = humanfile # assumed to be sorted externally...
		mousefilesorted = mousefile # assumed to be sorted externally...
	else:
		if not path.isdir(intermdir):
			makedirs(intermdir)
		humanfilesorted = path.join(intermdir,humanprefix+".human.namesorted.bam")
		mousefilesorted = path.join(intermdir,mouseprefix+".mouse.namesorted.bam")
		#print("Name sorting human and mouse BAM files using samtools")
		if not path.isfile(humanfilesorted):
			pysam.sort("-n","-m","2000000000",humanfile,humanfilesorted.replace(".bam",""))
		if not path.isfile(mousefilesorted):
			pysam.sort("-n","-m","2000000000",mousefile,mousefilesorted.replace(".bam",""))
		#print("Intermediate name sorted BAM files stored under " + intermdir(
	
	#print("Processing human and mouse files for ambiguous reads")
   # read in human reads and form a dictionary
	myHumanFile = pysam.Samfile(humanfilesorted, "rb" )
	myMouseFile = pysam.Samfile(mousefilesorted, "rb" )
	if not path.isdir(outputdir):
		makedirs(outputdir)
	myHumanUniqueFile = pysam.Samfile(path.join(outputdir, humanprefix+".human.bam"), "wb", template=myHumanFile) 
	myHumanAmbiguousFile = pysam.Samfile(path.join(outputdir, humanprefix+".ambiguousHuman.bam"), "wb", template=myHumanFile)
	myMouseUniqueFile = pysam.Samfile(path.join(outputdir, mouseprefix+".mouse.bam"), "wb", template=myMouseFile)
	myMouseAmbiguousFile = pysam.Samfile(path.join(outputdir, mouseprefix+".ambiguousMouse.bam"), "wb", template=myMouseFile)
	summaryFile = open(path.join(outputdir,humanprefix+'_summary.txt'),'w')
	
	#initialise
	try: 
		nexthumread=myHumanFile.next()
		nextmouread=myMouseFile.next()
	except StopIteration:
		print("No reads in one or either of the input files")
		sys.exit(2)
	
	EOFmouse = EOFhuman = False
	prevHumID = '-+=RANDOMSTRING=+-'
	prevMouID = '-+=RANDOMSTRING=+-'
	while not EOFmouse&EOFhuman:
		while not (nat_cmp(nexthumread.qname,nextmouread.qname) == 0):
			# check order between current human and mouse qname (find a point where they're identical, i.e. in sync)		
			while nat_cmp(nexthumread.qname,nextmouread.qname) > 0 and not EOFmouse: # mouse is "behind" human, output to mouse disambiguous
				myMouseUniqueFile.write(nextmouread)
				if not nextmouread.qname == prevMouID:
					nummou+=1 # increment mouse counter for unique only
				prevMouID = nextmouread.qname
				try:
					nextmouread=myMouseFile.next()
				except StopIteration:
					#print("1")
					EOFmouse=True
			while nat_cmp(nexthumread.qname,nextmouread.qname) < 0 and not EOFhuman: # human is "behind" mouse, output to human disambiguous
				myHumanUniqueFile.write(nexthumread)
				#print(nexthumread.qname + " " + prevHumID +" " +str(not nexthumread.qname == prevHumID))
				#temp = raw_input("press enter: ")
				if not nexthumread.qname == prevHumID:
					numhum+=1 # increment human counter for unique only
				prevHumID = nexthumread.qname
				try:
					nexthumread=myHumanFile.next()
				except StopIteration:
					#print("2")
					EOFhuman=True
			if EOFhuman or EOFmouse:	
				break
		# at this point the read qnames are identical and/or we've reached EOF
		humlist = list()
		moulist = list()
		#print(nexthumread.qname + " " + nextmouread.qname + " " + str(nat_cmp(nexthumread.qname,nextmouread.qname)))
		if nat_cmp(nexthumread.qname,nextmouread.qname) == 0:
			humlist.append(nexthumread)
			nexthumread = read_next_reads(myHumanFile, humlist) # read more reads with same qname (the function modifies humlist directly)
			if nexthumread == None:
				EOFhuman = True
			moulist.append(nextmouread)
			nextmouread = read_next_reads(myMouseFile, moulist) # read more reads with same qname (the function modifies moulist directly)
			if nextmouread == None:
				EOFmouse = True
		
		# perform comparison to check mouse, human or ambiguous
		if len(moulist) > 0 and len(humlist) > 0:
			myAmbiguousness = disambiguate(humlist, moulist, disambalgo)
			#print(myAmbiguousness)
			if myAmbiguousness < 0: # mouse
				nummou+=1 # increment mouse counter
				for myRead in moulist:
					myMouseUniqueFile.write(myRead)
			elif myAmbiguousness > 0: # human
				numhum+=1 # increment human counter
				for myRead in humlist:
					myHumanUniqueFile.write(myRead)
			else: # ambiguous
				numamb+=1 # increment ambiguous counter
				for myRead in moulist:
					myMouseAmbiguousFile.write(myRead)
				for myRead in humlist:
					myHumanAmbiguousFile.write(myRead)
		if EOFhuman:
			#flush the rest of the mouse reads
			#prevMouID = None
			while not EOFmouse:
				myMouseUniqueFile.write(nextmouread)
				if not nextmouread.qname == prevMouID:
					nummou+=1 # increment mouse counter for unique only
				prevMouID = nextmouread.qname
				try:
					nextmouread=myMouseFile.next()
				except StopIteration:
					#print("3")
					EOFmouse=True
		if EOFmouse:
			#flush the rest of the human reads
			#prevHumID = None
			while not EOFhuman:
				myHumanUniqueFile.write(nexthumread)
				if not nexthumread.qname == prevHumID:
					numhum+=1 # increment human counter for unique only
				prevHumID = nexthumread.qname
				try:
					nexthumread=myHumanFile.next()
				except StopIteration:
					#print("4")
					EOFhuman=True

		#end while not

	summaryFile.write("sample\tunique human pairs\tunique mouse pairs\tambiguous pairs\n")
	summaryFile.write(humanprefix+"\t"+str(numhum)+"\t"+str(nummou)+"\t"+str(numamb)+"\n")
	summaryFile.close()
	myHumanFile.close()
	myMouseFile.close()
	myHumanUniqueFile.close()
	myHumanAmbiguousFile.close()
	myMouseUniqueFile.close()
	myMouseAmbiguousFile.close()
   
	#print("Time taken in minutes " + str((time.clock() - starttime)/60))

if __name__ == "__main__":
	opts = docopt(__doc__,  help=False)  
	main(opts)
	#main(sys.argv[1:])
