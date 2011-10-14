//BioFilterReader.h

#ifndef __BIOFILTERREADER_H__
#define __BIOFILTERREADER_H__

#include <vector>
#include <string>
#include <fstream>
#include "BioReader.h"


/// Reads models from file 
class BioFilterReader: public BioReader{

  public:
  
    BioFilterReader();
    
    /// Fills vector with models from file
    int GetModels(std::vector<BioModel>& models, std::string filename, unsigned int max_read);
  
  private:
  
    void initialize();
  
    std::string filename;
    unsigned int maximum_reads;
    
    std::ifstream reader;
    
};

#endif


