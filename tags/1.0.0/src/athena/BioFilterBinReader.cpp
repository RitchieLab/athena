//BioFilterBinReader.cpp

#include "BioFilterBinReader.h"
#include "AthenaExcept.h"
#include <sstream>
#include <iostream>

using namespace std;

///
/// Constructor
///
BioFilterBinReader::BioFilterBinReader(){
  initialize();
}

///
/// Initializes object
///
void BioFilterBinReader::initialize(){
  maximum_reads = 500;
}


///
/// Fills vector with models from file
/// @param models Vector will contain models in order read from the file
/// @param filename File to read from
/// @param max_read Maximum number of models to read
/// @return Number of models actually read
///
int BioFilterBinReader::GetModels(std::vector<BioModel>& models, string filename, unsigned int max_read){
  
  
  maximum_reads = max_read;
  if(!reader.is_open()){
    reader.open(filename.c_str(), ios::binary);
    if(!reader.is_open()){
      throw AthenaExcept("Unable to open bio filter file " + filename);
    }
  }
  
  models.clear();
  
  unsigned int modelSize, groupCount, ddCount, ddRedundant, geneCount;
  unsigned int locus, group, ddGroup, gene;
  unsigned int impScore;
  
  while(!reader.eof() && (models.size() < maximum_reads) && reader.good()){

    BioModel mod;
  
    // read in all data -- but only using ID of snps and implication index score
    reader.read((char*)&modelSize, sizeof(unsigned int));
    reader.read((char*)&groupCount, sizeof(unsigned int));
    reader.read((char*)&ddCount, sizeof(unsigned int));
    reader.read((char*)&ddRedundant, sizeof(unsigned int));
    reader.read((char*)&geneCount, sizeof(unsigned int));
    reader.read((char*)&impScore, sizeof(float));
    mod.implication_index = impScore;
    
    // get model id numbers ( rs numbers with 'rs' removed)
    // put 'rs' back on to ID before setting ID string in model
    for(size_t i=0; i<modelSize; i++){
      reader.read((char*)&locus, sizeof(unsigned int));
      stringstream ss;
      ss << locus;
      mod.idstring.push_back("rs" + ss.str());
    }

    // get group IDs -- not currently used in ATHENA
    vector<unsigned int> groups(groupCount,0);
    for(size_t i=0; i<groupCount; i++){
      reader.read((char*)&group, sizeof(unsigned int));
      groups[i] = group;
    }
    
    vector<unsigned int> ddGroups(ddCount, 0);
    for(size_t i=0; i<ddCount; i++){
      reader.read((char*)&ddGroup, sizeof(unsigned int));
      ddGroups[i] = group;
    }
    
    vector<unsigned int> genes(geneCount, 0);
    for(size_t i=0; i<geneCount; i++){
      reader.read((char*)&gene, sizeof(unsigned int));
      genes[i] = gene;
    }
    models.push_back(mod);
  }
  if(models.size() < maximum_reads || reader.eof() || !reader.good())
    reader.close();  
  return models.size();
}

