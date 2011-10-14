//BioReader.h

#ifndef __BIOREADER_H__
#define __BIOREADER_H__

#include <vector>
#include <string>
#include <fstream>

/// Contains information related to Biofilter models provided
struct BioModel{
  std::vector<std::string> idstring;
  float implication_index;
};


///
/// Abstract base class for bio filter readers
///
class BioReader{

  public:
    virtual ~BioReader(){}

    /// Fills vector with models from file
    virtual int GetModels(std::vector<BioModel>& models, std::string filename, unsigned int max_read)=0;

};

#endif
