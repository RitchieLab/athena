//BioFilterReader.cpp

#include "BioFilterReader.h"
#include "HemannExcept.h"
#include <sstream>
#include <iostream>

using namespace std;

///
/// Constructor
///
BioFilterReader::BioFilterReader(){
  initialize();
}

///
/// Initializes object
///
void BioFilterReader::initialize(){
  maximum_reads = 500;
}


///
/// Fills vector with models from file
/// @param models Vector will contain models in order read from the file
/// @param filename File to read from
/// @param max_read Maximum number of models to read
/// @return Number of models actually read
///
int BioFilterReader::GetModels(std::vector<BioModel>& models, string filename, unsigned int max_read){
  
  maximum_reads = max_read;

  if(!reader.is_open()){
    reader.open(filename.c_str(), ios::in);
    if(!reader.is_open())
      throw HemannExcept("Unable to open bio filter file " + filename);
  }
  
  models.clear();
  
  // assume 2 locus models for now
  string id, line;
  while(!reader.eof() && (models.size() < maximum_reads)){
    getline(reader, line);
    if(line.find_first_of("0123456789") == string::npos)
      continue;
    stringstream ss(line);
    BioModel mod;
    ss >> id;
    mod.idstring.push_back(id);
    ss >> id;
    mod.idstring.push_back(id);
    ss >> mod.implication_index;
    models.push_back(mod);
  }

  // when all models done reading close the stream
  if(models.size() < maximum_reads)
    reader.close();

  return models.size();

}
