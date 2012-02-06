//BioFilterModelCollection.cpp

#include "BioFilterModelCollection.h"
#include "BioFilterReader.h"
#include "BioFilterBinReader.h"
#include <vector>
#include "genegenemodelreader.h"

using namespace std;

/// 
/// Default constructor
///
BioFilterModelCollection::BioFilterModelCollection(){

}


/// 
/// Alternative constructor
/// @param filename Filename containing biofilter models
/// @param nModels Number of models to retrieve
/// @param bioFileType
///
BioFilterModelCollection::BioFilterModelCollection(std::string filename, int nModels, string bioFileType){
  initialize(nModels);
  fill_tree(filename, nModels, bioFileType);
}


///
/// Alternative constructor that uses the gene-gene files and archive file
/// type to extract snp-snp models.
/// @param genegeneFile Gene-gene file
/// @param archiveFile 
/// @param Number of models to retrieve
///
BioFilterModelCollection::BioFilterModelCollection(std::string genegeneFile, std::string archiveFile,
  int nModels){
  initialize(nModels);
  fill_tree_archive(genegeneFile, archiveFile, nModels);
  
}


///
/// Initializes biofile type map
/// @param nModels Number of models that will be placed in the tree
///
void BioFilterModelCollection::initialize(int nModels){
  biomodels.SetMaxSize(nModels);
  randmodels.SetMaxSize(nModels);

  bioFileMap["TEXT"] = textFile;
  bioFileMap["BINARY"] = binaryFile;  
}

///
/// Returns a BioReader object based on string passed
///
BioReader* BioFilterModelCollection::getBioReader(string bioFileType){
  BioReader* reader;
  
  switch(bioFileMap[bioFileType]){
    case textFile:
      reader = new BioFilterReader();
      break;
    case binaryFile:
      reader = new BioFilterBinReader();
      break;
    case NoMatch:
    default:
      throw AthenaExcept(bioFileType + " is not a valid parameter for specifying biofilter file type");
  }
  return reader;
  
}



///
/// Fills collection tree by getting all models from reader.  The tree will only hold up to
/// the maximum set in the tree.
///
void BioFilterModelCollection::fill_tree(std::string filename, int nModels, string bioFileType){

//   BioFilterReader reader;
  BioReader* reader = getBioReader(bioFileType);
  
  vector<BioModel> models;
  vector<BioModel>::iterator modIter;
  
  
  do{
    reader->GetModels(models, filename, nModels);
    for(modIter = models.begin(); modIter != models.end(); modIter++){
      biomodels.Insert(modIter->implication_index, *modIter);
    }
  }while(int(models.size()) == nModels);

  // tree will now hold as many models as requested or as many as were in the file
  // determine the percentage increments for each 
  float impTotal = getImplicationTotal();

  // need to fill the random tree now with percentages
  float currPercent = 0.0;

 // creating tree where the keys are the fractional minimum for selecting the model
 // in the node.  A randomly generated number from 0 to 1 will be able to select
 // the appropriate model by using the find nearest min method of the tree and
 // returning those loci.
 // -- changed so that probability is simply implication / total implication -- smd 11/2/09
 for(BioModelTreeNode* currNode = biomodels.GetFirst();currNode != NULL; currNode = currNode->GetNext()){
    currPercent += currNode->GetKey() / impTotal;
    randmodels.Insert(currPercent, currNode->GetData());
  }

  delete reader;
}


///
/// Returns total of all implication scores in the tree
/// @returnn implication total
///
float BioFilterModelCollection::getImplicationTotal(){
   float impTotal=0.0;
   for(BioModelTreeNode* currNode = biomodels.GetFirst();currNode != NULL; currNode = currNode->GetNext()){ 
    impTotal += currNode->GetKey();
  }
  return impTotal;
}



///
/// Fills tree with samples from genegene file from biofilter
///
void BioFilterModelCollection::fill_tree_archive(std::string genegeneFile, std::string archiveFile, 
  int nModels){
  bool isbinary=false;

  // if it is binary the first unsigned int will equal 0
  // otherwise it will be some other number that forms part of
  // the count in the text file
  ifstream reader;
  unsigned int binaryCheck;
  reader.open(archiveFile.c_str(), ios::binary);
  if(!reader.is_open()){
    throw AthenaExcept("Unable to open bio filter file " + archiveFile);
  }
  reader.read((char*)&binaryCheck, sizeof(unsigned int));
  if(binaryCheck==0){
    isbinary = true;
  }
  reader.close(); 
  
  Biofilter::GeneGeneModelReader modelArchive(genegeneFile.c_str(), archiveFile.c_str(), isbinary);

  Biofilter::GeneGeneModelReader::iterator iter = modelArchive.begin();
  
  Biofilter::GeneGeneModelReader::iterator endIter = modelArchive.end();
  
  int currModels = 0;
  
  Biofilter::SnpSnpModel* snp_model;
  set<uint>::iterator snpiter;
  float imp_totals = 0.0;

  while(currModels < nModels){
    // when end is reached of gene gene models
    // loop through again 
    if(iter == endIter){
      iter = modelArchive.begin();
    }

    snp_model = iter->RandomModel();
    
    BioModel biomod;
    for(snpiter = snp_model->snps.begin(); snpiter != snp_model->snps.end(); snpiter++){
      stringstream ss;
      ss << *snpiter;
      biomod.idstring.push_back("rs" + ss.str());
    }
    biomod.implication_index = iter->ImplicationIndex();
    imp_totals += biomod.implication_index;
    // place in tree
    biomodels.Insert(biomod.implication_index, biomod);
    currModels++;
    iter++;
  }
  
  float currPercent = 0.0;
  for(BioModelTreeNode* currNode = biomodels.GetFirst();currNode != NULL; currNode = currNode->GetNext()){
    currPercent += currNode->GetKey()/imp_totals;
    randmodels.Insert(currPercent, currNode->GetData());
  }
 
}


///
/// Calculates and returns a vector containing the fractional
/// percentage increases for the models of each index
/// @param totals Map with key as index and total as the number of models with that implication score
/// @return vector of fractional 
///
map<int, float> BioFilterModelCollection::GetPercentages(map<int, int>& totals){

  float total_size = 0;
  map<int, int>::iterator iter;
  
  for(iter = totals.begin(); iter != totals.end(); iter++){
    total_size += iter->first; 
  }
  
  map<int, float> each_div;
  
  for(iter = totals.begin(); iter != totals.end(); iter++){
    each_div[iter->first] = iter->first / total_size;
  }
  
  map<int, float> model_fract;
  
  for(iter = totals.begin(); iter != totals.end(); iter++){
    model_fract[iter->first] = each_div[iter->first] / iter->second;
  }

  return model_fract;
}



///
/// need to determine how many implication indexes are present 
// in the tree
/// @return map with unique indexes in it from tree and totals for each
///
map<float, int> BioFilterModelCollection::GetImplicationTotals(){
  BioModelTreeNode* currNode;

  map<float, int> implication_totals;

  for(currNode = biomodels.GetFirst(); currNode != NULL;
    currNode = currNode->GetNext()){   
    if(implication_totals.find(currNode->GetKey()) == implication_totals.end()){
      implication_totals[currNode->GetKey()] = 1;
    }
    else{
      implication_totals[currNode->GetKey()]++;
    }
  }
  
  return implication_totals;
}


///
/// Returns random model from tree.
///
BioModel BioFilterModelCollection::GetRandomModel(){
  int randomint = rand();
  float randomNum = float(randomint) / (RAND_MAX);

  RandomTreeNode* randNode = randmodels.FindNearestMax(randomNum);
  return randNode->GetData();
}

