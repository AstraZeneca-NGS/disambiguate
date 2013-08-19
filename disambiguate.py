# usage: disambiguate.py -h <human.sam> -m <mouse.sam>

#imports
import sys, getopt, re, time
from array import array
from os import path, makedirs

#code
def main(argv):
	humanfile = ''
	mousefile = ''
	outputdir = 'disambres/'
	starttime = time.clock()
	# parse input arguments
	try:
		opts, args = getopt.getopt(argv,"h:m:o:",["help"])
		if len(opts) < 1:
			print 'Usage: disambiguate.py -h <human.sam> -m <mouse.sam> -o <outputdir>'
			sys.exit()
	except getopt.GetoptError:
		print 'Usage: disambiguate.py -h <human.sam> -m <mouse.sam> -o <outputdir>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '--help':
			print 'disambiguate.py -h <human.sam> -m <mouse.sam> -o <outputdir>'
			sys.exit(2)
		elif opt in ("-h"):
			humanfile = arg
		elif opt in ("-m"):
			mousefile = arg
		elif opt in ("-o"):
			outputdir = arg

   #set some defaults
	dv = 2**13 # a high quality score to replace missing quality scores (no real quality score should be this high)
   # form 'defaultArray' to fill new entries in the dictionry.
   # the fields are: human_1, human_2, mouse_1, mouse_2, human (bool), mouse (bool), ambiguous (bool)
	combineDict = dict()

   # read in human reads and form a dictionary
	myHumanFile = open(humanfile,'r')
	myMouseFile = open(mousefile,'r')
	if not path.isdir(outputdir):
		makedirs(outputdir)
	myHumanUniqueFile = open(path.join(outputdir,path.basename(humanfile.replace(".sam",".disambiguousHuman.sam"))),'w')
	myHumanAmbiguousFile = open(path.join(outputdir,path.basename(humanfile.replace(".sam",".ambiguousHuman.sam"))),'w')
	myMouseUniqueFile = open(path.join(outputdir,path.basename(mousefile.replace(".sam",".disambiguousMouse.sam"))),'w')
	myMouseAmbiguousFile = open(path.join(outputdir,path.basename(mousefile.replace(".sam",".ambiguousMouse.sam"))),'w')
	for line in myHumanFile.readlines():
		if re.compile('@').match(line): # check if line begins with @. If so, ignore
			continue
		splitLine=line.replace('\n','').split('\t')
      # generate quality score as sum of XO, NM and NH (integers)
		temp = re.compile('XO:i:(\d+).*NM:i:(\d+).*NH:i:(\d+)').search(line)
		if not temp:
			print "ERROR: Human score can not be determined for read "+splitLine[0]
			sys.exit(2)
		QScore = int(temp.group(1)) + int(temp.group(2)) + int(temp.group(3))
      # directionality (_1 or _2)
		d12 = 0 if 0x40&int(splitLine[1]) else 1

		if splitLine[0] in combineDict: # check if dictionary entry exists
      	# check if current entry is higher..
			if combineDict[splitLine[0]][d12]>QScore:
				combineDict[splitLine[0]][d12]=QScore # update to lowest (i.e. 'best') quality score
		else:
      	# generate new dictionary entry and update field for quality score
			#combineDict[splitLine[0]] = np.array([dv, dv, dv, dv], dtype=np.uint16)
			combineDict[splitLine[0]] = array('i',(dv for i in range(0,4)))
			combineDict[splitLine[0]][d12] = QScore
		QScore = temp = d12 = None
		#continue	# redundant continue
	
	myHumanFile.close()
	line = None
	# read in mouse reads and append dictionary
	for line in myMouseFile.readlines():
		if re.compile('@').match(line): # check if line begins with @. If so, ignore
			continue
		splitLine=line.replace('\n','').split('\t')
      # generate quality score as sum of XO, NM and NH (integers)
		temp = re.compile('XO:i:(\d+).*NM:i:(\d+).*NH:i:(\d+)').search(line)
		if not temp:
			print "ERROR: Mouse score can not be determined for read "+splitLine[0]
			sys.exit(2)
		QScore = int(temp.group(1)) + int(temp.group(2)) + int(temp.group(3))
      # directionality (_1 or _2)
		d12 = 2 if 0x40&int(splitLine[1]) else 3
		if splitLine[0] in combineDict: # check if dictionary entry exists
      	# check if current entry is higher..
			if combineDict[splitLine[0]][d12]>QScore:
				combineDict[splitLine[0]][d12]=QScore # update to lowest (i.e. 'best') quality score
		else:
      	# generate new dictionary entry and update field for quality score
			#combineDict[splitLine[0]] = np.array([dv, dv, dv, dv], dtype=np.uint16)
			combineDict[splitLine[0]] = array('i',(dv for i in range(0,4)))
			combineDict[splitLine[0]][d12] = QScore
		QScore = temp = d12 = None
		#continue # redundant continue
	   
	myMouseFile.close()

	# read in human and mouse reads one more time and divide into unique and ambiguous
	myHumanFile = open(humanfile,'r')
	myMouseFile = open(mousefile,'r')
	for line in myHumanFile.readlines():
		if re.compile('@').match(line):
			myHumanUniqueFile.write(line)
			myHumanAmbiguousFile.write(line)
			continue
		splitLine=line.replace('\n','').split('\t')
		h1,h2,m1,m2=combineDict[splitLine[0]][0:4]
		myUnique = min(h1,h2) < min(m1,m2) or min(h1,h2) == min(m1,m2) and max(h1,h2) < max(m1,m2)
		myAmbiguous = min(h1,h2) == min(m1,m2) and max(h1,h2) == max(m1,m2)
		if myUnique:
			myHumanUniqueFile.write(line)
		elif myAmbiguous:
			myHumanAmbiguousFile.write(line)
		
	for line in myMouseFile.readlines():
		if re.compile('@').match(line):
			myMouseUniqueFile.write(line)
			myMouseAmbiguousFile.write(line)
			continue
		splitLine=line.replace('\n','').split('\t')
		h1,h2,m1,m2=combineDict[splitLine[0]][0:4]
		myUnique = min(h1,h2) > min(m1,m2) or min(h1,h2) == min(m1,m2) and max(h1,h2) > max(m1,m2)
		myAmbiguous = min(h1,h2) == min(m1,m2) and max(h1,h2) == max(m1,m2)
		if myUnique:
			myMouseUniqueFile.write(line)
		elif myAmbiguous:
			myMouseAmbiguousFile.write(line)

   # close the rest of the files
   
	myHumanUniqueFile.close()
	myHumanAmbiguousFile.close()
	myMouseUniqueFile.close()
	myMouseAmbiguousFile.close()
   
	print "Time taken in minutes " + str((time.clock() - starttime)/60)
	print "Memory taken by dictionary in Mbytes " + str(sys.getsizeof(combineDict)/(2**30))

if __name__ == "__main__":
	main(sys.argv[1:])


