/*
 * dismain.cpp
 *
 *  Created on: Sep 19, 2013
 *  Author: kknj005
 *
*       The MIT License (MIT)

		Copyright (c) <2013> <AstraZeneca PLC 'and Simon Gray, Miika Ahdesmaki, Zhongwu Lai>

		Permission is hereby granted, free of charge, to any person obtaining a copy
		of this software and associated documentation files (the "Software"), to deal
		in the Software without restriction, including without limitation the rights
		to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
		copies of the Software, and to permit persons to whom the Software is
		furnished to do so, subject to the following conditions:

		The above copyright notice and this permission notice shall be included in
		all copies or substantial portions of the Software.

		THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
		IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
		FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
		AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
		LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
		OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
		THE SOFTWARE.
 */
//you don't really need to pay any attention to this part, sweep it under the rug.
//they just include various libraries which are needed for the code to function.
//leave them alone.
#include <iostream>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <list>
#include <vector>
#include <algorithm>
//this thing let's us find the size of an array as c++ isn't smart enough for that
#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))
//namespaces allow the use of commands that come with them
using namespace std;
using namespace BamTools;

//much to my dismay these had to be global due to complications with method calls later on
//they represent the ends of the input files
bool EOFmouse,EOFhuman;

//this method is used to compare strings and returns a number
//it is taken from the samtools c api and performs a natural comparison
//it is very complex.
static int strnum_cmp(const char *_a, const char *_b)
{
	const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
	const unsigned char *pa = a, *pb = b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			while (*pa == '0') ++pa;
			while (*pb == '0') ++pb;
			while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
			if (isdigit(*pa) && isdigit(*pb)) {
				int i = 0;
				while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
				return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
			} else if (isdigit(*pa)) return 1;
			else if (isdigit(*pb)) return -1;
			else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
		} else {
			if (*pa != *pb) return (int)*pa - (int)*pb;
			++pa; ++pb;
		}
	}
	return *pa? 1 : *pb? -1 : 0;
}

//this function takes a Bamreader object, a list of BamAlignments and a boolean to determine if its human or mouse
//you don't even know what some of those words mean do you?
//A Bamreader is the .bam iterator for the input files and an Alignment is the read from the bam file outputted by the reader
BamAlignment read_next_reads(BamReader& temp,list<BamAlignment>& alist, bool isHuman )
{
	bool qnamediff; //this boolean represents whether the read name is the same or not
	qnamediff = false;
	BamAlignment myread; //just a temporary BamAlignment object used to store the output from the bam reader
	while(!qnamediff)//while qnamediff is false
	{

		if((temp.GetNextAlignment(myread))== true)//as long as there is a next alignment in the file this will be true
		{
		if(strnum_cmp((myread.Name).c_str(),((alist.front()).Name).c_str())==0){ //comparing the names 0 means the same name
			alist.push_back(myread); //inserting the read into the list
		}
		else{
			qnamediff = true;//setting qname
		}
		}
		else{
			qnamediff = true;// i think you know by now
			if(isHuman)//test if it's a human
			{
				EOFhuman = true;//setting the global EOF (end of file) variable
			}
			else
			{
				EOFmouse = true;//this is does the same thing for mice
			}
		}
	}
	return myread;
}

//are you sitting down?
//good.
//you're in for the long haul with this function, go grab a drink and a snack first.
//right now you are finally ready let's get started
//this is the disambiguate function it is probably the most important function
//this is used to differentiate between mouse aligned reads, human aligned reads
//and ones that just can't be clearly distinguished.
//it takes two lists one of human reads and one of mouse reads and a string identifying the algorithm
int disambiguate(list<BamAlignment>& hlist,list<BamAlignment>& mlist, string disambalgo)
{
	int result; //this is returned at the end and is used to tell where the read should be sent
	//i.e it's the result.....
	if (disambalgo == "tophat")//this checks for the algorithm
		//for tophat the lower the score the better
	{
		int dv; //this is used to fill the array with a impossibly high number
		dv = 8192;//it is then checked against with the scores from the reads
		int sa [4];//creating the array
		for(int i = 0; i < 4; i++) //filling the array
		{
			sa[i] = dv;
		}
		for(std::list<BamAlignment>::iterator it = hlist.begin(); it!= hlist.end(); it++)//for starting at the beginning
			//of the human list until the end
		{
			BamAlignment temp = *it;
			if(!temp.IsMapped())//this checks if the read has aligned properly
			{
				//cout << "this is not mapped" << endl;
				continue;//will skip over the next part if it hasn't aligned properly
			}
			int QScore;

			if(temp.IsMapped())//double checking whether it has aligned properly
			{
				uint32_t NM; //declaring ints to assign the value returned below
				uint32_t XO;//i literally just explained this
				uint32_t NH;//ask google.
				temp.GetTag<uint32_t>("NM", NM);//calling the gettag method on a bamalignment object
				temp.GetTag<uint32_t>("XO", XO);//checking for the " " tag using a string and using
				temp.GetTag<uint32_t>("NH", NH);//the previously created ints to give the return value to
				QScore = NM + XO + NH; //this is simple maths.

			}

			int d12;//this is used to reference positions in the array
			if(temp.IsFirstMate())//this returns true if it is _1

			{
				d12 = 0;//if it's _1 it is in slot 0
			}
			else {d12 = 1;}//else it is _2

			if(sa[d12] > QScore){ //next it checks whether sa[d12] is greater than the QScore of the read
				sa[d12] = QScore;//if so it is replaced by QScore
			}

		}//human for loop
		//the for loop below follows the same structure as the human for loop
		//if you replace the hs with ms it'll be fine i promise.
		//i.e mlist -> hlist
		for(std::list<BamAlignment>::iterator it = mlist.begin(); it!= mlist.end(); it++)
				{
					BamAlignment temp = *it;
					if(!temp.IsMapped())
					{
						//cout << "this is not mapped" << endl;
						continue;
					}
					int QScore;

					if(temp.IsMapped())
					{
						uint32_t NM;
						uint32_t XO;
						uint32_t NH;
						temp.GetTag<uint32_t>("NM", NM);
						temp.GetTag<uint32_t>("XO", XO);
						temp.GetTag<uint32_t>("NH", NH);
						QScore = NM + XO + NH;

					}

					int d12;
					if(temp.IsFirstMate())
					{
						d12 = 2;
					}
					else {d12 = 3;}

					if(sa[d12] > QScore){
						sa[d12] = QScore;
					}

				}//mouse for loop
		//To explain the array sa has four slots 0,1 are human _1 and _2 respectively
		//slots 2,3 are mouse _1 and _2 respectively
		//this is the complex mathy bit
		//the first if here checks first, whether the min of the human reads is equal to the min of the mouse reads
		// AND then it checks if this max of the human reads is equal to the max of the mouse reads
		//if it is then the read is ambiguous
		if(((min(sa[0],sa[1]))==(min(sa[2],sa[3])))&((max(sa[0],sa[1]))==(max(sa[2],sa[3]))))
		{
			result = 0;

		}
		//this compares the same as before apart it checks if the human minimum is less than the mouse minimum
		//OR if the minimums are equal
		//AND if the human max is less than the mouse max
		else if((min(sa[0],sa[1]))<(min(sa[2],sa[3]))|| (((min(sa[0],sa[1]))==(min(sa[2],sa[3])))&((max(sa[0],sa[1]))<(max(sa[2],sa[3])))))
		{
			result = 1;

		}
		else //otherwise mouse i.e -1
		{
			result = -1;

		}
	}//if tophat
	//the next algorithm is BWA, told you it was a long one didn't i
	//have a drink, rest your eyes
	//watch some cat videos on youtube
	//right you must be ready by now!
	//let's go
	else if (disambalgo == ("BWA")||disambalgo == ("bwa"))
		//with bwa the higher the score the better
	{

		int dv;//you'll remember this from above right?
		dv = -8192;//you should re read it then...
		//this is where bwa gets a little different
		string bwatagname[3] = {"AS","NM","XS"};//so we have an array of names to be used later
		int bwatagsigns[3] = {1,-1,1};//because of the scale between the different numbers we have to multiply
		//the tag NM by -1 so the scale fits
		int size;

		size = ARRAY_SIZE(bwatagsigns);//this is used to make the size of certain variabls change based
		//on how many tags you want


		int sa [size][4];//this array is 2 dimensional the columns represent the _1 and _2 of human and mouse reads
		//the size is used to make the row count dynamic as rows represent the different tags

		for (int j = 0; j < size; j++ )//as the same as in the tophat this loops through all the slots and replaces
			//with dv
		{

		for(int i = 0; i < 4; i++)
			{

			sa[j][i] = dv;
			}
		}

		for(std::list<BamAlignment>::iterator it = hlist.begin(); it!= hlist.end(); it++)
			{
				BamAlignment temp = *it;

				if(!temp.IsMapped())//checking whether it aligned properly
				{

					continue;
				}
				int QScore;
				int d12;
				if(temp.IsFirstMate())//checking if it's _1
				{
					d12 = 0;
				}
				else {d12 = 1;}//setting d12
				for (int x = 0; x < size; x++)//this loops through all the sizes of x
								{
									//creating a temp value to store the gettag results
									uint32_t *bwatagvaluetemp = new uint32_t(1);
									//we use x to get the right tagname
									temp.GetTag<uint32_t>(bwatagname[x],*bwatagvaluetemp);
									//and then get the corresponding tag sign
									//bwa is done sequentially over grouped so the Qscores are calculated per Tag
									//maths tagsign * tagvalue from gettag
									QScore = ((bwatagsigns[x]) * (int)(*bwatagvaluetemp));

								if(sa[x][d12] < QScore){
									sa[x][d12] = QScore;
								}

								//empty the value otherwise gettag doesn't write!
								delete bwatagvaluetemp;
								}

			}//human for loop

		//again replace the ms with hs!
		for(std::list<BamAlignment>::iterator it = mlist.begin(); it!= mlist.end(); it++)
			{
				BamAlignment temp = *it;
				if(!temp.IsMapped())
				{
					//cout << "this is not mapped" << endl;
					continue;
				}
				int QScore;
				int d12;
				if(temp.IsFirstMate())
				{
					d12 = 2;
				}
				else {d12 = 3;}
				for (int x = 0; x < size; x++)
				{
					uint32_t *bwatagvaluetemp = new uint32_t(1);

					temp.GetTag<uint32_t>(bwatagname[x],*bwatagvaluetemp);

					QScore = ((bwatagsigns[x]) * (int)(*bwatagvaluetemp));

				if(sa[x][d12] < QScore){
					sa[x][d12] = QScore;
				}

				delete bwatagvaluetemp;
				}
			}



				//the for and ifs are slightly different but you should be able to follow based on earlier comments
				for (int x = 0; x < size; x++)
				{
						//the main difference being returning over changing result and returning after
						//this is due to the sequential methodology of BWA
						if((max((sa[x][0]),(sa[x][1])))>(max((sa[x][2]),(sa[x][3])))|| (((max((sa[x][0]),(sa[x][1])))==(max((sa[x][2]),(sa[x][3]))))&((min((sa[x][0]),(sa[x][1])))>(min((sa[x][2]),(sa[x][3]))))))
						{
							return 1;

						}
						else if((max(sa[x][0],sa[x][1]))<(max(sa[x][2],sa[x][3]))|| (((max(sa[x][0],sa[x][1]))==(max(sa[x][2],sa[x][3])))&((min(sa[x][0],sa[x][1]))<(min(sa[x][2],sa[x][3])))))
						{
							return -1;

						}
						else {
							result = 0;
						}

						}

	}
	else {
		//just in case someone puts in something stupid for the algorithm value
		//users huh?
		cerr << "Selected alignment algorithm option not yet implemented" << endl;
		exit(1);
		}
	return result;

}

int main(int argc, char **argv) {

	string outputPath; //created now set by user later using ovalue
	string suffix; //created now set by user later using svalue
	string algorithm;//are you detecting a pattern yet? avalue...
	//input files are set here
	string human; //this is set by hvalue
	string mouse;//this is set by mvalue



			char *hvalue = NULL; //human file
			char *mvalue = NULL; //mouse file
			char *ovalue = NULL; //output path
			char *svalue = NULL; //suffix value actually a prefix but shhh
			char *avalue = NULL; //algorithm name

			//the following lines of code may or may not set up the user input
			//i am not at liberty to discuss this
			//just don't tell anyone and you won't die
			//it just sets up the options for user input via command line

	       int index;//this is used later for error checking
	       int c; //this is used to accept arguments
	       opterr = 0;
	       while ((c = getopt (argc, argv, "h:m:o:s:a:di")) != -1)

	         switch (c)
	           {
	           case 'h':
	             hvalue = optarg; //optarg means that -h will accept an argument i.e -h human.bam
	             break;
	           case 'm':
	             mvalue = optarg;
	             break;
	           case 'o':
	             ovalue = optarg;
	             break;
	           case 's':
	        	 svalue = optarg;
	        	 break;
	           case 'a':
	        	   avalue = optarg;
	        	   break;
	           case 'd':
				  cerr << "Disablesort option is not needed." << endl; //not yet implemented
				   break;
	           case 'i':
				  cerr << "intermediate directory not needed." << endl; //not yet implemented
				   break;
	           case '?':
	        	   if (optopt == 'c')//HELP COMMENTING MIIKA
	        	                   fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	        	               else if (isprint (optopt))
	        	                   fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	        	               else
	        	                   fprintf (stderr,
	        	                            "Unknown option character `\\x%x'.\n",
	        	                            optopt);
	             return 1;
	           default: abort ();
	           }

	       //printing values this probably won't be here
	       printf ("hvalue = %s\n mvalue = %s\n ovalue = %s\n svalue = %s\n avalue = %s\n",
	               hvalue, mvalue, ovalue, svalue, avalue);

	       //checking for too many arguments?
	       for (index = optind; index < argc; index++)
	       printf ("Non-option argument %s\n", argv[index]);

	       //this is just a help file basically
	       //argc is equal to 1 if the only thing you type is the ./programname etc
	       //argc is argument count
	       if(argc<2)
	       	       {
	    	   	   	   cout <<"   "<< endl;
	       	    	   cout <<"Usage:"<< endl;
	       	    	   cout <<"./disambiguate [-h <human.bam>] [-m <mouse.bam>] [-o <outputdir>] [-s <samplenameprefix>] [-a ALGORITHM] "<< endl;
	       	    	   cout <<"   "<< endl;


	       	    	   cout <<"   "<< endl;
	       	    	   cout <<"   "<< endl;
	       	    	   cout <<"   "<< endl;

	       	    	   cout <<"Options:   "<< endl;

	       	    	   cout <<"-h <human.bam>            Input human BAM file location (REQUIRED)  "<< endl;

	       	    	   cout <<"-m <mouse.bam>            Input mouse BAM file location (REQUIRED)   "<< endl;

	       	    	   cout <<"-o <outputdir>            Output directory. If the directory does not exist it   "<< endl;
	       	    	   cout <<"                          will be generated. If not provided, 'disambres/'    "<< endl;
	       	    	   cout <<"                          will be used.   "<< endl;

	       	    	   cout <<"-s <samplenameprefix>     A prefix (e.g. sample name) to use for the output   "<< endl;
	       	    	   cout <<"                          BAM files. If not provided, the human BAM file   "<< endl;
	       	    	   cout <<"                          prefix will be used.  "<< endl;

	       	    	   cout <<"-a ALGORITHM              tophat (default) or bwa. Note that for bwa    "<< endl;
	       	    	   cout <<"                          alignments the tags provided in the BAM files are    "<< endl;
	       	    	   cout <<"                          different to those from tophat.   "<< endl;
	       	    	   cout <<"   "<< endl;
	       	    	   exit(1);
	       	       }
	       //these next set of ifs give miniature help based on certain nulls
	       if(hvalue == NULL)
	       {
	    	   cerr << "***********************************************" << endl;
	    	   cerr << "You must set all arguments properly." << endl;
	    	   cerr << "-h takes a file path + filename with extension." << endl;
	    	   cerr << "run the program without any commands for help. " << endl;
	    	   cerr << "***********************************************" << endl;
	    	   exit(1);
	       }
	       if(mvalue == NULL)
	       {
			   cerr << "***********************************************" << endl;
			   cerr << "You must set all arguments properly." << endl;
			   cerr << "-m takes a file path + filename with extension." << endl;
			   cerr << "run the program without any commands for help. " << endl;
			   cerr << "***********************************************" << endl;
			   exit(1);
		   }
	       if(ovalue == NULL)
		   {
			   cerr << "**********************************************" << endl;
			   cerr << "You must set all arguments properly." << endl;
			   cerr << "-o takes a directory path for output files." << endl;
			   cerr << "run the program without any commands for help." << endl;
			   cerr << "**********************************************" << endl;
			   exit(1);
		   }
	       if(svalue == NULL)
		   {
			   cerr << "**********************************************" << endl;
			   cerr << "You must set all arguments properly." << endl;
			   cerr << "-s takes a string that will be the prefix." << endl;
			   cerr << "run the program without any commands for help." << endl;
			   cerr << "**********************************************" << endl;
			   exit(1);
		   }
	       if(avalue == NULL)
		   {
			   cerr << "**********************************************" << endl;
			   cerr << "You must set all arguments properly." << endl;
			   cerr << "-a takes a string that is the algorithm name.[BWA,tophat]" << endl;
			   cerr << "run the program without any commands for help. " << endl;
			   cerr << "**********************************************" << endl;
			   exit(1);
		   }
	       //this assigns the user values to the program variables
	       suffix = svalue;
	       outputPath = ovalue;
	       algorithm = avalue;
	       human = hvalue;
	       mouse = mvalue;
   //output files are set here
	//out file names are locked however can be customised with a prefix
	       //prefix is called suffix for some unknown reason i may change this
	       //output path is set by user and then the file location is built using insertions
	string humanDisambiguous = "human.bam"; //this stays the same/becomes new default
	humanDisambiguous.insert(0,"."); //same
	humanDisambiguous.insert(0,suffix); //set by svalue
	humanDisambiguous.insert(0,"/"); //same
	humanDisambiguous.insert(0,outputPath); //set by ovalue
	string humanAmbiguous = "humanAmbiguous.bam";
	humanAmbiguous.insert(0,"."); //same
	humanAmbiguous.insert(0,suffix); //set by svalue
	humanAmbiguous.insert(0,"/"); // same
	humanAmbiguous.insert(0,outputPath); //set by ovalue
	string mouseDisambiguous = "mouse.bam";
	mouseDisambiguous.insert(0,"."); //same
	mouseDisambiguous.insert(0,suffix); //set by svalue
	mouseDisambiguous.insert(0,"/");//same
	mouseDisambiguous.insert(0,outputPath); //set by ovalue
	string mouseAmbiguous = "mouseAmbiguous.bam";
	mouseAmbiguous.insert(0,"."); //same
	mouseAmbiguous.insert(0,suffix); //set by svalue
	mouseAmbiguous.insert(0,"/"); // same
	mouseAmbiguous.insert(0,outputPath); //set by ovalue



	//File input Readers
	BamReader humanReader; //makes new read object
	if( !humanReader.Open(human) ){ //checks whether it can open the file i.e if it is a bam file
		cerr << "Could not open input BAM file." << endl;
		}
	BamReader mouseReader;
	if( !mouseReader.Open(mouse) ){
		cerr << "Could not open input BAM file." << endl;
		}

	//Get headers and reference data from input files
	//these are needed for the output files!
	const SamHeader humanHeader = humanReader.GetHeader();
	const RefVector humanReferences = humanReader.GetReferenceData();

	const SamHeader mouseHeader = mouseReader.GetHeader();
	const RefVector mouseReferences = mouseReader.GetReferenceData();

	//File output writers
	BamWriter HDWriter; //creating writer object
	if ( !HDWriter.Open(humanDisambiguous, humanHeader, humanReferences) ) { //checking whether the file can be
		//and written to.
	cerr << "Could not open output BAM file" << endl;
	}
	BamWriter HAWriter;
	if ( !HAWriter.Open(humanAmbiguous, humanHeader, humanReferences) ) {
	cerr << "Could not open output BAM file" << endl;
	}
	BamWriter MDWriter;
	if ( !MDWriter.Open(mouseDisambiguous, mouseHeader, mouseReferences) ) {
	cerr << "Could not open output BAM file" << endl;
	}
	BamWriter MAWriter;
	if ( !MAWriter.Open(mouseAmbiguous, mouseHeader, mouseReferences) ) {
	cerr << "Could not open output BAM file" << endl;
	}


	BamAlignment ah; //creating alignment objects
	BamAlignment am;

	EOFhuman = false; //setting end of file global variables
	EOFmouse = false;
	string prevHumID = "RANDOMSTRING"; //setting to random string for comparison reasons
	string prevMouID = "RANDOMSTRING";
	int nummou = 0; //setting up unique read counters
	int numhum = 0;
	int numamb = 0;
	//if there are no next alignments the program will exit
	if((humanReader.GetNextAlignment(ah)==false)||(mouseReader.GetNextAlignment(am)==false))
	{
		exit(1);
	}

	while (!EOFmouse & !EOFhuman){
		//while both EOF globals are false
		while(!strnum_cmp((ah.Name).c_str(),(am.Name).c_str())==0){
			//while the names of a human and mouse reads are NOT the same
			while((strnum_cmp((ah.Name).c_str(),(am.Name).c_str())>0) && !EOFmouse){
				//while the names of the human and mouse read comparison is greater than 0 and not EOFmouse
				MDWriter.SaveAlignment(am);
				//Write alignment to file

					if (am.Name != prevMouID){ //if mouse name is unique
						nummou+=1;	}
						prevMouID == am.Name;

					if((mouseReader.GetNextAlignment(am)==false)) //getting next alignment and testing for EOF
						{
							EOFmouse = true;
						}

			}//strnum_cmp>0 & !EOFmouse loop

			while((strnum_cmp((ah.Name).c_str(),(am.Name).c_str())<0) && !EOFhuman){
				//while the names of the human and mouse read comparison is less than 0 and not EOFhuman
				HDWriter.SaveAlignment(ah);
				//i think you can get the jist from before right?
						if (ah.Name != prevHumID){
									numhum+=1;	}
									prevHumID == ah.Name;
					if((humanReader.GetNextAlignment(ah)==false))
						{
							EOFhuman = true;
						}

			}//strnum_cmp<0 & !EOFhuman loop
			if (EOFhuman || EOFmouse){ //testing for EOFs being true
					break;}

		}//strnum_cmp == 0 loop
	//Now the reads should be identical or reached EOF
	//creating the lists
	list<BamAlignment> humlist;
	list<BamAlignment> moulist;
	//if the names are the same
	if (strnum_cmp((ah.Name).c_str(),(am.Name).c_str())==0){

		bool isHuman;
		isHuman = true; //telling between human and mouse
		humlist.push_back(ah); //putting alignment into list
		ah = read_next_reads(humanReader,humlist,isHuman); //using read_next_reads function
		moulist.push_back(am);
		isHuman = false;
		am = read_next_reads(mouseReader,moulist,isHuman);

	}//if strnum == 0

	if (((moulist.size()) > 0) & ((humlist.size()) > 0)) //if the lists both have reads in them
	{

		int myAmbiguousness; //temp int to store the result of the disambiguate function

		myAmbiguousness = disambiguate(humlist,moulist,algorithm);//calling disambiguate function
		if(myAmbiguousness < 0)//testing for results
			//then processing the results
		{
			nummou++; //increment unique read count
			for(std::list<BamAlignment>::iterator it = moulist.begin(); it!= moulist.end(); it++)
					{
						BamAlignment temp = *it;
						MDWriter.SaveAlignment(temp); //interating through and writing the reads to file
					}
		}
		else if(myAmbiguousness > 0)
		{
			numhum++;
			for(std::list<BamAlignment>::iterator it = humlist.begin(); it!= humlist.end(); it++)
					{
						BamAlignment temp = *it;
						HDWriter.SaveAlignment(temp);
					}
		}
		else
		{
			numamb++; //incrementing the ambiguous read count and then writing to files
			for(std::list<BamAlignment>::iterator it = moulist.begin(); it!= moulist.end(); it++)
					{
									BamAlignment temp = *it;
									MAWriter.SaveAlignment(temp);
					}
			for(std::list<BamAlignment>::iterator it = humlist.begin(); it!= humlist.end(); it++)
					{
									BamAlignment temp = *it;
									HAWriter.SaveAlignment(temp);
					}
		}
	}//if mou/humlist > 0

		if(EOFhuman) //testing for EOF human
		{

			while(!EOFmouse) //this empties the mouse file so all reads are processed
			{

				MDWriter.SaveAlignment(am);
				if(am.Name != prevMouID)
					{

					nummou++;
					}
				prevMouID = am.Name;
				if((mouseReader.GetNextAlignment(am)==false))
					{

						EOFmouse = true;
					}
			}
		}
		if(EOFmouse) //save as above but flushing human files
				{
					while(!EOFhuman)
					{
						HDWriter.SaveAlignment(ah);
						if(ah.Name != prevHumID)
						{
							numhum++;
						}
						prevHumID = ah.Name;
						if((humanReader.GetNextAlignment(ah)==false))
						{
							EOFhuman = true;
						}

					}
				}


	}//EOFmouse & EOFhuman loop

	//closing the humanReaders and writers
	humanReader.Close();
	mouseReader.Close();
	HDWriter.Close();
	HAWriter.Close();
	MDWriter.Close();
	MAWriter.Close();



	return 0;

}


