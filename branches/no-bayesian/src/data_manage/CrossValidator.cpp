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
#include "CrossValidator.h"

#include <algorithm>
#include <stdlib.h>
#include <fstream>

namespace data_manage
{

CrossValidator::CrossValidator()
{
}

CrossValidator::~CrossValidator()
{
}


///
///  Splits data into testing and training sets for the number of
///  CV intervals selected
///  @param numCrossVal Number of crossvalidation intervals
///  @param holder Dataholder containing data for this analysis
///	 @return CVSet
///
CVSet CrossValidator::splitData(unsigned int numCrossVal, Dataholder* holder){

	if(holder->getTestSplit() > 0)
		return splitByNum(holder);

	vector<Individual*> shuffled, affected, unaffected;
	statusBin(holder, affected, unaffected);

	if(numCrossVal > 1){
		shuffleInds(affected);
		shuffleInds(unaffected);
	}
	else{
		statusBinAffOnly(holder, affected, unaffected);
	}

	vector<Individual*> temp;
	
	splits.assign(numCrossVal, temp);
	distributeInds(numCrossVal, affected, splits);
	distributeInds(numCrossVal, unaffected, splits);

	return createSet(numCrossVal, holder);
}



/// 
/// Creates CVSet based on splits and number of cross-validations
///
CVSet CrossValidator::createSet(unsigned int numCrossVal, Dataholder* holder){

	CVSet set;
	if(numCrossVal > 1){
		unsigned int group;
		// using the splits stored in vector construct the CV Intervals and fill the set
		for(unsigned int currCV=0; currCV < numCrossVal; currCV++){
			Dataset training(holder->getMissingCoValue(), holder->getMissingGenotype(), holder->isCaseControl()), 
				testing(holder->getMissingCoValue(), holder->getMissingGenotype(), holder->isCaseControl());
			for(group=0; group < numCrossVal; group++){
				if(group != currCV)
					training.addInds(splits[group]);
				else
					testing.addInds(splits[group]);
			}
			CVInterval interval;
			training.calcSSTotal();
			testing.calcSSTotal();
			training.setHolder(holder);
			testing.setHolder(holder);
// 			training.setAllLevels(holder->getAllLevels());
// 			testing.setAllLevels(holder->getAllLevels());
			interval.addSet(training);
			interval.addSet(testing);
			set.addInterval(interval);
		}
	}
	else{ // only one interval so don't split data
		CVInterval interval;
		Dataset training(holder->getMissingCoValue(), holder->getMissingGenotype(), 
			holder->isCaseControl());
		training.addInds(splits[0]);
		training.calcSSTotal();
// 		training.setAllLevels(holder->getAllLevels());
		training.setHolder(holder);
		interval.addSet(training);
		set.addInterval(interval);
	}	
	return set;
}



///
/// Splits data according to index in the holder class.
/// For use with user-specified splits.
/// @param holder Dataholder
///
CVSet CrossValidator::splitByNum(Dataholder* holder){
	int splitIndex = holder->getTestSplit();
	CVSet set;

	Dataset training(holder->getMissingCoValue(), holder->getMissingGenotype(), holder->isCaseControl()), 
		testing(holder->getMissingCoValue(), holder->getMissingGenotype(), holder->isCaseControl());
	
	vector<Individual*> testSet, trainSet;
	
	int i=0;
	for(i=0; i<splitIndex; i++){
		trainSet.push_back(holder->getInd(i));
	}  
	
	int nInds = int(holder->numInds());
	for(; i<nInds; i++){
		testSet.push_back(holder->getInd(i));
	}

	training.addInds(trainSet);
	testing.addInds(testSet);
	
	CVInterval interval;
	training.calcSSTotal();
	testing.calcSSTotal();
	interval.addSet(training);
	interval.addSet(testing);
	set.addInterval(interval);
	
	return set;
}

///
/// Splits vector into specified set
/// @param numSplits
/// @param inds vector of Individual pointers
/// @param splits two-dimensional vector of Individual pointers that will contain
/// individual splits
///
void CrossValidator::distributeInds(unsigned int numSplits, vector<Individual*>& inds,
		vector<vector<Individual*> >& splits){

	for(unsigned int i=0; i<inds.size(); i++){
		splits[i%numSplits].push_back(inds[i]);
	}

}


///
/// Distributes individuals between affected and unaffected arrays
/// @param holder Dataholder
/// @param affected Vector will contain pointers to affected inds
/// @param unaffected Vector will contain pointers to unaffected inds
///
void CrossValidator::statusBin(Dataholder* holder, vector<Individual*>& affected,
		vector<Individual*>& unaffected){

	affected.clear();
	unaffected.clear();

	unsigned int n_inds = holder->numInds();
	Individual* ind;
	for(unsigned int i=0; i < n_inds; i++){
		ind = holder->getInd(i);
		if(ind->getStatus())
			affected.push_back(ind);
		else
			unaffected.push_back(ind);
	}
}



///
/// Places all individuals in affected 
/// @param holder Dataholder
/// @param affected Vector will contain pointers to affected inds
/// @param unaffected Vector will contain pointers to unaffected inds
///
void CrossValidator::statusBinAffOnly(Dataholder* holder, vector<Individual*>& affected,
		vector<Individual*>& unaffected){

	affected.clear();
	unaffected.clear();

	unsigned int n_inds = holder->numInds();
	Individual* ind;
	for(unsigned int i=0; i < n_inds; i++){
		ind = holder->getInd(i);
		affected.push_back(ind);
	}
}


///
/// Shuffles individuals in the dataset
/// so that the splitting of the data will be independent of the
/// original order in the data file
/// @param indexes vector will contain pointers to now shuffled individuals
///
void CrossValidator::shuffleInds(vector<Individual*> & inds){
	TRandom myrand;
	random_shuffle(inds.begin(), inds.end(), myrand);
}


///
/// Loads the splits from a file.
/// @param filename
/// @param holder Dataholder
///	 @return CVSet
///
CVSet CrossValidator::loadSplits(std::string filename, Dataholder* holder){
	ifstream is(filename.c_str());
		if(!is.is_open()){
			throw DataExcept("ERROR: Unable to open CV split file " + filename + "\n");
		}
	string id;
	splits.clear();
	vector<Individual*> inds;
	// remove first 'split' line
	is >> id;
	while(is >> id){
		if(id.compare("split")==0){
			splits.push_back(inds);
			inds.clear();
			continue;
		}
		inds.push_back(holder->getIndByID(id));
	}
	splits.push_back(inds);
	is.close();
	return createSet(splits.size(), holder);
}

///
/// Saves splits to file as ID from the individuals
/// @param filename
///
void CrossValidator::saveSplits(std::string filename){
	ofstream os;
	os.open(filename.c_str(), ios::out);
	for(vector<vector<Individual*> >::iterator splitIter=splits.begin(); splitIter != splits.end();
		splitIter++){
		os << "split" << endl;
		for(vector<Individual*>::iterator indIter=splitIter->begin(); indIter != splitIter->end();
			indIter++){
			os << (*indIter)->getID() << endl;
		}
	}
	os.close();
}


}
