#include "ExpressionReader.h"
#include <iostream>
#include <math.h> 

ExpressionReader::ExpressionReader(const std::string& expFileName)
	:expFileName_(expFileName)
{
	expressionMap_.clear();
}
          
/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @param
 * @return
 * @throw
 */
void ExpressionReader::loadExpressionData(const std::string& targetGeneID,bool log2Transform=false){
          std::ifstream expressionFile;
          expressionFile.open(expFileName_);
          if (!expressionFile) throw std::invalid_argument("Expression file "+expFileName_+" could not be opened");

	expressionMap_.clear();
          //Reading sample Names
          std::string temp, buf;
          std::getline(expressionFile,temp);
          std::vector<std::string> sampleNames;
          std::stringstream sN(temp);
          while(sN >> buf){
                    sampleNames.push_back(buf);
          }         
	unsigned int numberSamples = sampleNames.size();
          //Generating expression map
          double value;
          unsigned int counter=0;
	bool flag=false;
          while (!expressionFile.eof()){
                    std::getline(expressionFile,temp);
                    std::stringstream sE(temp);
                    if (temp != ""){
                              if (sE >> buf){				      
			          if (buf == targetGeneID){
					  std::cout<< "target found!" <<std::endl;
                                                  while (sE >> value){
						// std::cout << value << std::endl;	  
						if (log2Transform)
						 	expressionMap_[sampleNames[counter]]=log2(value+1.0);
						else	    
                                                        expressionMap_[sampleNames[counter]]=value;
							//std::cout<< sampleNames[counter] << std::endl;
							//std::cout<< expressionMap_[sampleNames[counter]] << std::endl;
                                                        counter+=1;
                                                  }
 						//std::cout << "Expression Map: " << std::endl;
	         				/*for(std::map<std::string, double>::const_iterator it = expressionMap_.begin(); it != expressionMap_.end(); ++it)
	       			 		{
		                     			std::cout << it->first << " " << it->second << "\n";
        					}*/
						/*std::cout << expressionMap_.size() << std::endl;
						std::map<std::string,double>::iterator it = expressionMap_.begin();
						for (it=expressionMap_.begin(); it!=expressionMap_.end(); ++it)
							std::cout << it->first << " => " << it->second << '\n';*/		
						flag=true;
                                                expressionFile.close();
					if (counter != numberSamples){
						throw std::invalid_argument("The number of gene expression values for gene "+targetGeneID+" does not match the number of expected entries.");
                                        	}
                              	}
			}
                              else{
				expressionFile.close();
                                        throw std::invalid_argument("Expression file "+expFileName_+" is not properly formatted");
                              	}
                    	}
          	}
	if (!flag){
	          expressionFile.close();
	          throw std::invalid_argument("Expression data for the gene "+targetGeneID+" could not be found");
	}
}

void ExpressionReader::checkDiversity(){
	std::map<int,unsigned int> occurences;
	for (const auto& element : expressionMap_){
		if (occurences.find(element.second)==occurences.end())
			occurences[element.second]=1;
		else
			occurences[element.second]+=1;
	}
	unsigned int mapSize=occurences.size();
	if (mapSize==1)
		throw std::runtime_error("Expression file "+expFileName_+" contains only one value for the specified gene");
	else{
		double minSize=expressionMap_.size()*0.03;
		for (const auto& element : occurences){
			if (element.second < minSize)
				throw std::runtime_error("At least three percent of the data should have the discrete expression value of "+std::to_string(element.first));
		}
	}
}

std::map<std::string, double>& ExpressionReader::getExpressionMap(){
	return expressionMap_;
}

const std::string& ExpressionReader::getFilename(){
	return expFileName_;

}

std::ostream& operator<<(std::ostream& os, const ExpressionReader& r){
          for (const auto& element : r.expressionMap_){
                    os << element.first << " " << element.second << std::endl;
          }
	return os;
}   
