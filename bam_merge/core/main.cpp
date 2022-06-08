/*
 * author: Laura Rumpf
 * distribute bam file entries to metacell bam files*
 *
 */

#include "../bam/SamFile.h"
#include "../bam/SamTags.h"
#include "../bam/SamRecord.h"
#include "../bam/SamValidation.h" 
#include "../bam/BaseUtilities.h"
#include <stdio.h>
#include <stdlib.h>
#include <functional>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <queue>
#include <string.h>
#include <filesystem>
#include <array>
#include <unordered_map>

#define MAXLINE 8000

/** @brief returns true if the TAG of record specified in @param tag is 
 * identical to @param value
 */
bool getBarcode(SamRecord& record, const char *tag, std::string &value) {
  char vtype = 'Z';
  char *tmp_tag = new char();
  String *tmp_val = new String;

 
  value.erase(remove(value.begin(), value.end(), '"'), value.end());
  

  while (record.getNextSamTag(tmp_tag, vtype, (void**)(&tmp_val))) {


    if (*tmp_tag == tag[0] && *(tmp_tag + 1) == tag[1]) {
   
        //std::cout << "CB tag found!" << std::endl;
      if ((*(tmp_val)).Compare(value.c_str()) == 0) {
     	          
        return true;
      }
    }
  }
  return false;
}

/** @brief returns true if the cellID of the record is identical to the value
 * specified in @param value.
 */
std::string getCellID(SamRecord& record, int c_len, int pos){
  const char *tmp_name = record.getReadName();
  std::string name(tmp_name);
  
  name = name.substr(pos, c_len);
    
  return(name);
  
}

bool process_10xdata(const char* input_filename, const char *output_dir, const char *bam_filename) {
  bool successRead = false;
  std::ifstream input; 
  input.open(input_filename); 
  if (!input.is_open()) {
    fprintf(stderr, "input file can't be opened or invalid\n");
    std::cout << "format: ./core/bam_merge <flag> <input bam file> <input txt file>" << std::endl;
    std::cout << "flag -X for 10x genomic sequences, for non-10x sequences: -D DNA sequence, -R for RNA sequence" << std::endl;
    exit(1);
  }
  std::string line;
  if (!std::getline(input, line)) {
    std::cout << "invalid input" << std::endl;
    exit(1);
  }
   
  // initialize read
  // write output file name
  // there could be multiple so need to specify which one
  
  SamRecord read;

  char outputFileName[128];
  
  std::vector<std::string> metacells;
  std::vector<std::string> outfilenames;
  std::unordered_map<std::string, ptrdiff_t> mc_hash;
  // int csvCount = 0;
  // get all possible output sam file names from the metacell ids in the csv file
   
   while (std::getline(input, line)) {
	    
	      std::string delim = ",";
	      int d1 = line.find(delim);                         
	      std::string metacell = line.substr(d1 + 1);
	      std::string cell_id = line.substr(0, d1);


	       metacell.erase(remove(metacell.begin(), metacell.end(), '"'), metacell.end());	
	       cell_id.erase(remove(cell_id.begin(), cell_id.end(), '"'), cell_id.end());    

	      if (std::find(metacells.begin(), metacells.end(), metacell.c_str()) == metacells.end()){
		 metacells.push_back(metacell.c_str());
	         std::cout << metacells.back() << std::endl;	 
                 snprintf (outputFileName, sizeof(outputFileName), "%s/output_%s.bam", output_dir, metacell.c_str());
                 puts (outputFileName);
 		 outfilenames.push_back(outputFileName);		 
	      }

	       ptrdiff_t pos = distance(metacells.begin(), find(metacells.begin(), metacells.end(), metacell));

	       std::string hashkey = cell_id;
	       mc_hash[hashkey] = pos;
		      
   }
  //input.clear();
  //input.seekg(0, input.beg);
   input.close();
  // check hash entries
  for( const auto& [key, value] : mc_hash ) {
	          std::cout << "Key:[" << key << "] Value:[" << value << "]\n";
  }

 
  //int s = size(mc_hash); 
  //std::cout << "#elements in hashmap: " << s << std::endl;
  //std::cout << "#entries csv: " << count << std::endl;
  std::cout << "Hashmap is created!" << std::endl;
  std::cout << "Output filenames are stored!" << std::endl;

  std::vector<std::shared_ptr<SamFile>> openedFiles;
  
  
   // SamFile f;
   SamFileHeader fileHeader;
   SamFile fileIn(bam_filename, SamFile::READ, &fileHeader);

   // keep all sam files open for writing
   // store pointer to output bam files in vector 
   for(std::size_t i = 0; i < outfilenames.size(); ++i) {
     	   
     openedFiles.reserve(outfilenames.size());
   
      auto f_ptr = std::shared_ptr<SamFile>(new SamFile(outfilenames[i].c_str(), SamFile::WRITE, &fileHeader));
      openedFiles.push_back(f_ptr);
     
   }
   //  fileIn.Close(); 
  
   std::cout << "Output files are opened for writing!" << std::endl; 
   
    int write_cnt = 0;
    
    // loop through all the records to look for ones that matches cell ID
    do {
      successRead = fileIn.ReadRecord(fileHeader, read);
      if (!read.checkTag("CB",'Z')) { 
       		 const char* val = read.getReadName();
		 std::cout << "Entry for Read Name: " << val << "has no CB tag!" << std::endl;
		 continue;
      }
     
       String cell_barcode = read.getString("CB");
       std::string cb(cell_barcode);
       	
       // check if cell barcode is in hashmap
       if (mc_hash.find(cb) == mc_hash.end()) {	
	       //std::cout << cb << " is not in the hashmap!" << std::endl;
	       continue;
       }
	
       //std::string mc = mc_hash[cb];
       
       std::ptrdiff_t mc = mc_hash[cb];
       std::shared_ptr<SamFile>* fileOutPtr(&openedFiles[mc]);       

   
      (*fileOutPtr)->WriteRecord(fileHeader, read);
       write_cnt++;	
        
    } while (successRead);
    
    fileIn.Close();
  
  std::cout << "input bam is closed" << std::endl;

  for(std::size_t i = 0; i < openedFiles.size(); ++i) {
	openedFiles[i] -> Close();
  }
  std::cout << "all output bam files are closed" << std::endl;

  return true;	
  
}

bool process_non10xdata(const char* input_filename, const char *output_dir, const char *bam_filename, int pos_cell_id) {
  	
  bool successRead = false;

  std::ifstream input;
  input.open(input_filename); 

/*  if(input.good()){ 
	  std::cout << "No Reading Errors" << std::endl;
  }*/

  if (!input.is_open()) {
    fprintf(stderr, "input file can't be opened or invalid\n");
    std::cout << "format: ./core/bam_merge <flag> <input bam file> <input txt file>" << std::endl;
    std::cout << "flag -X for 10x genomic sequences, for non-10x sequences: -D DNA sequence, -R for RNA sequence" << std::endl;
    exit(1);
  }
  std::string line;
  if (!std::getline(input, line)) {
    std::cout << "invalid input" << std::endl;
    exit(1);
  }
  SamRecord read; // initialize read

  // write output file name
  // there could be multiple so need to specify which one
  char outputFileName[128];
   
  std::vector<std::string> metacells;
  std::vector<std::string> outfilenames;
  std::unordered_map<std::string, ptrdiff_t> mc_hash;
  int cell_id_len;

                                                                                                                                                                                                                     int csv_count = 0;
																										    // get all possible output sam file names from the metacell ids in the csv file
   while (std::getline(input, line)) {
	    
	      std::string delim = ",";
	      int d1 = line.find(delim);                         
	      std::string metacell = line.substr(d1 + 1);
	      std::string cell_id = line.substr(0, d1);

	      metacell.erase(remove(metacell.begin(), metacell.end(), '"'), metacell.end());
	      cell_id.erase(remove(cell_id.begin(), cell_id.end(), '"'), cell_id.end());

	      if (csv_count < 1){
		      cell_id_len = cell_id.length();
	      }

	      if (std::find(metacells.begin(), metacells.end(), metacell.c_str()) == metacells.end()){
		 metacells.push_back(metacell.c_str());
	         std::cout <<  metacells.back() << std::endl;	
		// metacell.erase(remove(metacell.begin(), metacell.end(), '"'), metacell.end());
                 snprintf (outputFileName, sizeof(outputFileName), "%s/output_%s.bam", output_dir, metacell.c_str());
                 puts (outputFileName);
 		 outfilenames.push_back(outputFileName);		 
	      }
	
	      //value: position of outputfile 
      	      ptrdiff_t pos = distance(metacells.begin(), find(metacells.begin(), metacells.end(), metacell));
	      
	      //std::cout << "cell_id: " << cell_id << std::endl;
	      //std::cout << "metacell: " << metacell << std::endl;
	      //std::cout << "mc_index: " << pos << std::endl;
	       
	      //key: cell id
	      std::string hashkey = cell_id;

	      mc_hash[hashkey] = pos;

	      csv_count++;
   }
 // input.clear();
 // input.seekg(0, input.beg);
    input.close();
 // check hash entries
    for( const auto& [key, value] : mc_hash ) {
    	std::cout << "Key:[" << key << "] Value:[" << value << "]\n";
    }

  std::cout << "Output filenames are stored!" << std::endl;
  std::cout << "Hashmap is created!" << std::endl;
    
  std::vector<std::shared_ptr<SamFile>> openedFiles;
   // SamFile f;
  SamFileHeader fileHeader;
  SamFile fileIn(bam_filename, SamFile::READ, &fileHeader);

   // keep all sam files open for writing
   // store pointer to output bam files in vector 
   for(std::size_t i = 0; i < outfilenames.size(); ++i) {
     	   
      openedFiles.reserve(outfilenames.size());
   
      auto f_ptr = std::shared_ptr<SamFile>(new SamFile(outfilenames[i].c_str(), SamFile::WRITE, &fileHeader));
      openedFiles.push_back(f_ptr);
     
   }
 //  fileIn.Close(); 
  
   std::cout << "Output files are opened for writing!" << std::endl; 
   
    int write_cnt = 0;

    // loop through all the records to look for ones that matches cell ID
    do {
      successRead = fileIn.ReadRecord(fileHeader, read);
    
      const char *tmp_name = read.getReadName();
      std::string cb(tmp_name);

      if (cb.compare("UNKNOWN") == 0) {
      		continue;
      }	 
	    
      cb = cb.substr(pos_cell_id, cell_id_len);
	
      //check if cell_id is in hashmap
      if (mc_hash.find(cb) == mc_hash.end()) {
      		std::cout << cb << " is not in the hashmap!" << std::endl;
		continue;
      }

      std::ptrdiff_t mc = mc_hash[cb];
     
      std::shared_ptr<SamFile>* fileOutPtr(&openedFiles[mc]);

      (*fileOutPtr)->WriteRecord(fileHeader, read);
      write_cnt++;
    } while (successRead);
  
    fileIn.Close();
  
  std::cout << "input bam is closed" << std::endl;

  for(std::size_t i = 0; i < openedFiles.size(); ++i) {
	openedFiles[i] -> Close();
  }
  std::cout << "all output bam files are closed" << std::endl; 
 
  
  return true;
}



int main(int argc, char **argv) {

  if (argc < 3) {
    std::cout << "format: ./core/bam_merge <flag> <input bam file> <input txt file> <position cell id in read name for non-10x input> <optional: output directory>" << std::endl;
    std::cout << "flag: -X for 10x genomic sequences, for non-10x sequences: -D DNA sequence" << std::endl;
    exit(1);
  }

  // checking flags
  bool f_none = false;
  bool f_10x = false;
  std::string flag = std::string(argv[1]);
 
  int pos_int;
  const char *output_dir;
  if (flag.compare("-X") == 0) {
	  f_10x = true;
          //optional 4th parameter output directory
	  if (argc < 4 || argv[4] == NULL){
	 	 output_dir = "outputs";
	  }
	  else {
	  	std::string out = std::string(argv[4]);
	 	output_dir = out.c_str();
	  }
  }
  else if(flag.compare("-D") == 0) {
	  f_none = true;
	 
	  if(argc < 4 || (argv[4] == NULL)){
	  	std::cout << "for non-10x sequences the start position of the cell id substring in the read name is required" << std::endl;
		exit(1);
	  }


	  std::string pos = std::string(argv[4]);

	  if (std::all_of(pos.begin(), pos.end(), isdigit)){
		std::cout << pos << std::endl;
	        pos_int = std::stoi(pos);
		if (pos_int < 0){
     			std::cout << "for non-10x sequences the start position of the cell id substring in readname must be a positve integer" << std::endl;	
			exit(1);
	    	}
		 //optional fifth parameter output directory
		if (argc < 5 || argv[5] == NULL){
			output_dir = "outputs";
		}
		else {
			std::string out = std::string(argv[5]);
			output_dir = out.c_str();
		}
	   }
           else {
	   	std::cout << "for non-10x sequences the start position of the cell id substring in readname must be a positve integer" << std::endl;
                exit(1);
	   }
	 }
  else {
	std::cout << "format: ./core/bam_merge <flag> <input bam file> <input txt file> <position cell id in read name for non-10x input> <optional: output directory>" << std::endl;
	std::cout << "flag: -X for 10x genomic sequences, for non-10x sequences: -D DNA sequence" << std::endl;
        std::cout << "an invalid flag was entered" << std::endl;
	exit(1);
  }
  // input text file w metacell number and cell barcode
  // taken in as the second parameter
  // path& p = std::filesystem::current_path();

  
  if (std::filesystem::create_directory(output_dir)) {
  	std::cerr<< "output directory created" << std::endl;
  }
  else  {
	std::cerr << strerror(errno) << std::endl;
  }




 const char* bam_filename = argv[2];
 const char* input_filename = argv[3];
 //const int* pos_cell_id = argv[4];
 
   std::cout << "input csv: " << input_filename << std::endl;
   std::cout << "input bam: " << bam_filename << std::endl;
 
 
  //###########################################################################
  // Set returnStatus to success.  It will be changed
  // to the failure reason if any of the writes fail.
  SamStatus::Status returnStatus = SamStatus::SUCCESS;
 
  //###########################################################################
  // std::cout << f_10x << std::endl;
  // std::cout << f_none << std::endl;
 
  if (! f_10x) {
 	 std::cout << "non-10xATAC" << std::endl;
 	 process_non10xdata(input_filename, output_dir, bam_filename, pos_int);
  }
  else {
  	std::cout << "10x" << std::endl;
 	process_10xdata(input_filename, output_dir, bam_filename);
 }
 
 
 for (const auto & entry : std::filesystem::directory_iterator(output_dir)){
 	std::filesystem::path p = entry.path();
 	std::string path_string_old = p.string();
 	std::string path_string_new = p.string();
 
 // path_string_new.erase(remove_if(path_string_new.begin(), path_string_new.end(), isspace), path_string_new.end());
 
	for(std::string::iterator it = path_string_new.begin(); it != path_string_new.end(); ++it) {
 		if(*it == ' ' || *it == '-') {
 			*it = '_';
		}
		 }

	 int result;
	 result= rename(path_string_old.c_str() ,  path_string_new.c_str());
	 if ( result != 0 )
		 perror( "Error renaming file" );
 }

    return(returnStatus);


}
