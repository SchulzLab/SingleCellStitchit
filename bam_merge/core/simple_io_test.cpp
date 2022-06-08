#include "../bam/SamFile.h"
#include "../bam/SamTags.h"
#include "../bam/SamRecord.h"
#include "../bam/SamValidation.h"
#include "../bam/BaseUtilities.h"
#include <stdio.h>
#include <stdlib.h>
#include <queue>

int ReadIndexedBam(const char* inputFilename,
                   const char* outputFilename,
                   const char* indexFilename)
{

   SamFileHeader samHeader;
   // Open the input file for reading.
   SamFile samIn(inputFilename, SamFile::READ, &samHeader);
   std::cout << "opened input file" << std::endl;

   // Open the bam index file for reading.
   samIn.ReadBamIndex(indexFilename);
   std::cout << "opened index file" << std::endl;

   // Open the output file for writing.
   SamFile samOut(outputFilename, SamFile::WRITE, &samHeader);
   std::cout << "opened output file" << std::endl;
   

   SamRecord samRecord;

   // Loop through each Reference.
   for(int i = -1; i < 23; i++)
   {
      int numSectionRecords = 0;
      samIn.SetReadSection(i);
      std::cout << "set read section" << i << std::endl;
      // Keep reading records until they aren't more.
      while(samIn.ReadRecord(samHeader, samRecord))
      {
         numSectionRecords++;
         // Successfully read a record from the file, so write it.
         samOut.WriteRecord(samHeader, samRecord);
      }

      std::cout << "Reference ID " << i << " has " << numSectionRecords 
                << " records" << std::endl;
   }
   
   std::cout << "Number of records = " << samIn.GetCurrentRecordCount() << std::endl;
   
   return(0);
}

int main(int argc, char ** argv)
{   /*
    printf("num of args: %d\n", argc);
    if (argc != 4) {
	printf("exit 1\n");
	return 1;
    }
    printf("%s, %s, %s, %s\n", argv[0], argv[1], argv[2], argv[3]);
    */

    ReadIndexedBam("../NA12878.bam", "../output.bam", "../NA12878.bai");
    printf("finishing compilation\n");
    return 0;
}

/*
int main(int argc, char ** argv)
{

   // Open the input file for reading.
   // SamFile samIn;
   SamFileHeader samHeader;
   SamFile samIn("../NA12878.bam", SamFile::READ, &samHeader);

   // Open the output file for writing.
   // SamFile samOut;
   SamFile samOut("../output.bam", SamFile::WRITE, &samHeader);
   // samOut.OpenForWrite("../output.bam");

   // Read the sam header.

   // Write the sam header.
   // samOut.WriteHeader(samHeader);

   SamRecord samRecord;

    // Set returnStatus to success.  It will be changed
    // to the failure reason if any of the writes fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;

    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Successfully read a record from the file, so write it.
        samOut.WriteRecord(samHeader, samRecord);
    }

    std::cout << std::endl << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl;
    std::cout << "Number of records written = " << 
        samOut.GetCurrentRecordCount() << std::endl;

    // Return success since a failure would have thrown
    // an exception.
    return(returnStatus);
 }
*/
