//BioFilterModelCollection.h

#ifndef __BIOFILTERMODELCOLLECTION_H__
#define __BIOFILTERMODELCOLLECTION_H__

#include "rbtree.h"
#include "BioReader.h"
#include "HemannExcept.h"
#include <map>

///
/// Holds biomodels and provides models either by position or
/// by random.  Each model is grouped with like models and then
/// every model within a grouping has an equal chance to be 
/// returned on every call for a random model.
///

// has static variable that allows for sorting either in descending (0) or ascending order (2)
struct intLT{
   int operator()(const int l, const int r) const {
      if (l>r) return 1;
      if (l<r) return -1;
      return 0;
   }
};

struct impLT{
   int operator()(const float l, const float r) const {
      if (l>r) return 1;
      if (l<r) return -1;
      return 0;
   }  
};


struct randFloatLT{
   int operator()(const float l, const float r) const {
      if (l<r) return 1;
      if (l>r) return -1;
      return 0;
   }
};

typedef Utility::RBTree<float, BioModel, impLT> BioModelTree;
typedef Utility::RBTreeNode<float, BioModel, impLT> BioModelTreeNode;
typedef Utility::RBTree<float, BioModel, randFloatLT> RandomTree;
typedef Utility::RBTreeNode<float, BioModel, randFloatLT> RandomTreeNode;

#include<iostream>
using namespace std;

class BioFilterModelCollection{

  public:

    BioFilterModelCollection();

    BioFilterModelCollection(std::string filename, int nModels, std::string bioFileType);
    
    BioFilterModelCollection(std::string genegeneFile, std::string archiveFile,
      int nModels);
    
    BioModel GetRandomModel();
    
    std::map<float, int> GetImplicationTotals();
    
    std::map<int,float> GetPercentages(std::map<int, int>& totals);
    
    BioModel GetNextModel(){       
      BioModelTreeNode* returnNode = currBioModelNode;
      currBioModelNode = currBioModelNode->GetPrev();
      if(currBioModelNode == NULL)
        currBioModelNode = biomodels.GetLast();
      return returnNode->GetData();
    }
    
    void SetStartModel(int startIndex){
      int count=0;
      startBioModelNode = biomodels.GetLast();
      while(count++ < startIndex){
        startBioModelNode = startBioModelNode->GetPrev();
        if(startBioModelNode == NULL) // wrap around
          startBioModelNode = biomodels.GetLast();
      }
      currBioModelNode = startBioModelNode;
    }
    
  private:
    
    BioReader* getBioReader(string bioFileType);
  
    enum BioModelFileType{
      NoMatch,
      textFile,
      binaryFile
    };
  
    map<std::string, BioModelFileType> bioFileMap;
  
    void initialize(int nModels);
    
    void fill_tree(std::string filename, int nModels, std::string bioFileType);
    
    void fill_tree_archive(std::string genegeneFile, std::string archiveFile, 
      int nModels);
    
    float getImplicationTotal();
    
    BioModelTree biomodels;
    RandomTree randmodels;
    
    int start_index;
    BioModelTreeNode* currBioModelNode, *startBioModelNode;

};

#endif
