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
#include <fstream>
#include <stdio.h>
//#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <list>
#include <vector>
#include <algorithm>
#include <tclap/CmdLine.h>

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
  string prefix; //created now set by user later using svalue
  string algorithm;//are you detecting a pattern yet? avalue...
  //input files are set here
  string human; //this is set by hvalue
  string mouse;//this is set by mvalue

  // use TCLAP to parse input command line arguments:
  // Wrap everything in a try block.  Do this every time,
  // because exceptions will be thrown for problems.
  try {
    // Define the command line object, and insert a message
    // that describes the program. The "Command description message"
    // is printed last in the help text. The second argument is the
    // delimiter (usually space) and the last one is the version number.
    // The CmdLine object parses the argv array based on the Arg objects
    // that it contains.
    TCLAP::CmdLine cmd("Command description message", ' ', "1.0");
    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-n Bishop".
    TCLAP::ValueArg<string> alignerArg("a","aligner","Aligner option {tophat(default),bwa}",false,"tophat","string"\
					    );
    TCLAP::ValueArg<string> outpudirArg("o","output-dir","Output directory",true,"disambres/","string");
    TCLAP::ValueArg<string> prefixArg("s","prefix","Sample ID or name used as prefix. Do not include .bam",true,"my\
Sample","string");
    // Add the argument nameArg to the CmdLine object. The CmdLine object
    // uses this Arg to parse the command line.
    cmd.add( alignerArg );
    cmd.add( outpudirArg );
    cmd.add( prefixArg );
    // Define a switch and add it to the command line.
    // A switch arg is a boolean argument and only defines a flag that
    // indicates true or false.  In this example the SwitchArg adds itself
    // to the CmdLine object as part of the constructor.  This eliminates
    // the need to call the cmd.add() method.  All args have support in
    // their constructors to add themselves directly to the CmdLine object.
    // It doesn't matter which idiom you choose, they accomplish the same thing.
    TCLAP::SwitchArg reverseSwitch("d","no-sort","(Deprecated option for backwards compatibility)", cmd, false);
    // Positional arguments A and B
    TCLAP::UnlabeledValueArg<string>  speciesA( "A", "Species A BAM file", true,
						     "A.bam", "A"  );
    TCLAP::UnlabeledValueArg<string>  speciesB( "B", "Species B BAM file", true,
						     "B.bam", "B"  );
    cmd.add( speciesA );
    cmd.add( speciesB );
    // Parse the argv array.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    outputPath = outpudirArg.getValue();
    prefix = prefixArg.getValue(); 
    algorithm = alignerArg.getValue();
    // input files are set here
    human = speciesA.getValue(); 
    mouse = speciesB.getValue();
  } catch (TCLAP::ArgException &e)  // catch any exceptions
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  
  //output files are set here
  //out file names are locked however can be customised with a prefix
  //output path is set by user and then the file location is built using insertions
  string humanDisambiguous = "disambiguatedSpeciesA.bam"; //this stays the same/becomes new default
  humanDisambiguous.insert(0,"."); //same
  humanDisambiguous.insert(0,prefix); //set by svalue
  humanDisambiguous.insert(0,"/"); //same
  humanDisambiguous.insert(0,outputPath); //set by ovalue
  string humanAmbiguous = "ambiguousSpeciesA.bam";
  humanAmbiguous.insert(0,"."); //same
  humanAmbiguous.insert(0,prefix); //set by svalue
  humanAmbiguous.insert(0,"/"); // same
  humanAmbiguous.insert(0,outputPath); //set by ovalue
  string mouseDisambiguous = "disambiguatedSpeciesB.bam";
  mouseDisambiguous.insert(0,"."); //same
  mouseDisambiguous.insert(0,prefix); //set by svalue
  mouseDisambiguous.insert(0,"/");//same
  mouseDisambiguous.insert(0,outputPath); //set by ovalue
  string mouseAmbiguous = "ambiguousSpeciesB.bam";
  mouseAmbiguous.insert(0,"."); //same
  mouseAmbiguous.insert(0,prefix); //set by svalue
  mouseAmbiguous.insert(0,"/"); // same
  mouseAmbiguous.insert(0,outputPath); //set by ovalue
  string summaryFile = "_summary.txt";
  summaryFile.insert(0,prefix);
  summaryFile.insert(0,"/");
  summaryFile.insert(0,outputPath);

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
	prevMouID = am.Name;

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
	prevHumID = ah.Name;
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
  // write summary stats to file
  ofstream myfile(summaryFile.c_str());
  myfile << "sample\tunique species A pairs\tunique species B pairs\tambiguous pairs" <<endl;
  myfile << prefix <<"\t" << numhum<<"\t" << nummou<<"\t" << numamb << endl;
  myfile.close();
  return 0;

}


