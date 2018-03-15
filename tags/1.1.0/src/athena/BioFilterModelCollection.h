/*
Copyright Marylyn Ritchie 2011

This file is part of ATHENA.

ATHENA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ATHENA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ATHENA.  If not, see <http://www.gnu.org/licenses/>.
*/
//BioFilterModelCollection.h

#ifndef __BIOFILTERMODELCOLLECTION_H__
#define __BIOFILTERMODELCOLLECTION_H__

#include "rbtree.h"
#include "BioReader.h"
#include "AthenaExcept.h"
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

class BioFilterModelCollection{

	public:

		BioFilterModelCollection();

		BioFilterModelCollection(std::string filename, int nModels, std::string bioFileType);
		
		BioFilterModelCollection(std::string genegeneFile, std::string archiveFile,
			int nModels);
		
		BioModel getRandomModel();
		
		std::map<float, int> getImplicationTotals();
		
		std::map<int,float> getPercentages(std::map<int, int>& totals);
		
		BioModel getNextModel(){       
			BioModelTreeNode* returnNode = currBioModelNode;
			currBioModelNode = currBioModelNode->GetPrev();
			if(currBioModelNode == NULL)
				currBioModelNode = bioModels.GetLast();
			return returnNode->GetData();
		}
		
		void setStartModel(int startIndex){
			int count=0;
			startBioModelNode = bioModels.GetLast();
			while(count++ < startIndex){
				startBioModelNode = startBioModelNode->GetPrev();
				if(startBioModelNode == NULL) // wrap around
					startBioModelNode = bioModels.GetLast();
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
		
		void fillTree(std::string filename, int nModels, std::string bioFileType);
		
		void fillTreeArchive(std::string genegeneFile, std::string archiveFile, 
			int nModels);
		
		float getImplicationTotal();
		
		BioModelTree bioModels;
		RandomTree randModels;
		
		int startIndex;
		BioModelTreeNode* currBioModelNode, *startBioModelNode;

};

#endif
