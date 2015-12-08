/*
Copyright Marylyn Ritchie 2015

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

#include "GAFunct.h"

float	GAFunct::connProb =0.1;
data_manage::Dataset* GAFunct::caseDataset = NULL;
data_manage::Dataset* GAFunct::controlDataset = NULL;
vector<Variable*> GAFunct::varList;
GABayesSolutionCreator GAFunct::caseBayesCreator;
GABayesSolutionCreator GAFunct::controlBayesCreator;
// double GAFunct::fitnessTime = 0.0;
// double GAFunct::loopTime = 0.0;
// double GAFunct::maxCheckTime = 0.0;

///
/// Used to assign fitness values in GA
/// @param g Genome to evaluate
///
float GAFunct::GACaseObjective(GAGenome& g){

	time_t startTime, endTime;

  GA2DArrayGenome<int> & genome = (GA2DArrayGenome<int> &)g;
  removeSelfAndDup(genome);
// time (&startTime);
//   caseBayesCreator.checkNodeLimits(genome);
// time (&endTime);
// double dif = difftime (endTime,startTime);
// maxCheckTime += dif;

// time (&startTime);
//   caseBayesCreator.fixLoops(genome);
	caseBayesCreator.limitChildren(genome);
	caseBayesCreator.breakLoops(genome);

// time(&endTime);
// dif = difftime (endTime,startTime);
// loopTime += dif;

// cout << "calculating case fitness" << endl;
// time(&startTime);
	float score=caseBayesCreator.calcScore(genome, varList, caseDataset);
// time(&endTime);
// dif = difftime (endTime,startTime);
// fitnessTime += dif;
// cout << "final score=" << score << endl;


// cout << "output test case set" << endl;
// for(unsigned int i=0; i<caseDataset->numInds(); i++){
// 	for(unsigned int j=0; j<caseDataset->numGenos(); j++){
// 	  cout << caseDataset->getInd(i)->getGenotype(j) << " ";
// 	}
// 	cout << endl;
// }
// cout << "==================" << endl;
// exit(1);

	return score;
// 	return caseBayesCreator.calcScore(genome, varList, caseDataset);
}

///
/// Used to assign fitness values in GA
/// @param g Genome to evaluate
///
float GAFunct::GAControlObjective(GAGenome& g){
  GA2DArrayGenome<int> & genome = (GA2DArrayGenome<int> &)g;
  removeSelfAndDup(genome);
  controlBayesCreator.checkNodeLimits(genome);
//   controlBayesCreator.fixLoops(genome);
  controlBayesCreator.breakLoops(genome);
// cout << "calculating control fitness" << endl;
	return controlBayesCreator.calcScore(genome, varList, controlDataset);
}


vector<std::string> GAFunct::getAdditionalFinalOutput(float score){
	vector<std::string> addOutputValues;
	if(score > caseDataset->getConstant() ){
		addOutputValues.push_back("+");
	}
	else{
		addOutputValues.push_back("-");
	}
	return addOutputValues;
}

///
/// Return additional column names for output
///
vector<std::string>  GAFunct::getAdditionalOutputNames(){
	vector<std::string> names(1, "not-improved");
	return names;
}

/// conducts random initialization of genomes
void GAFunct::initCase(GAGenome &g){
	init(g, caseBayesCreator);
}

/// conducts random initialization of genomes
void GAFunct::initControl(GAGenome &g){
	init(g, controlBayesCreator);
}


///
/// conducts random initialization of genomes
/// @param g GAGenome to initialize
/// @param gaBayesCreator GABayesSolutionCreator to use in initialization
///
void GAFunct::init(GAGenome &g, GABayesSolutionCreator& gaBayesCreator){
	GA2DArrayGenome<int> & genome = (GA2DArrayGenome<int> &)g;
	for(int y=0; y<genome.height(); y++){
		vector<int> conns;
		for(int x=0; x<genome.height(); x++){
			if(GARandomFloat() <= connProb)
				conns.push_back(x);
		}
		if(conns.size() > genome.width()){
			set<int> keep=gaBayesCreator.limitConnections(y,conns,genome.width());
			int x=0;
			for(set<int>::iterator iter=keep.begin(); iter != keep.end(); ++iter){
				genome.gene(x,y,*iter);
				x++;
			}
		}
		else{
			size_t x=0;
			for(; x<conns.size(); x++){
				genome.gene(x,y,conns[x]);
			}
			for(; x<genome.width(); x++){
				genome.gene(x,y,-1);
			}
		}
	}
}

///
/// unset any connections to self in matrix
/// @param g GAGenome to initialize
///
void GAFunct::removeSelfAndDup(GA2DArrayGenome<int>& genome){
	set<int> dups;
	int parentIdx;

	for(int y=0;y<genome.height();y++){
		dups.clear();
		dups.insert(y);
		for(int x=0; x<genome.width(); x++){
			parentIdx = genome.gene(x,y);
			if(parentIdx == -1)
				continue;
			else{
				if(dups.find(parentIdx) != dups.end()){
					genome.gene(x,y,-1);
				}
				else{
					dups.insert(parentIdx);
				}
			}
		}
	}
}


/// conducts random initialization of genomes
int GAFunct::mutateCase(GAGenome & c, float pmut){
	return customMutator(c, pmut, caseBayesCreator);
}

/// conducts random initialization of genomes
int GAFunct::mutateControl(GAGenome & c, float pmut){
	return customMutator(c, pmut, controlBayesCreator);
}

///
/// Custom mutator has same effect as flipmutator for binary string genome
///
int GAFunct::customMutator(GAGenome & c, float pmut, GABayesSolutionCreator& gaBayesCreator){
  GA2DArrayGenome<int> &genome=DYN_CAST(GA2DArrayGenome<int> &, c);
  register int n, m, i, j, k;
  if(pmut <= 0.0) return(0);

  float nMut = pmut * genome.size();
  set<int> muts;
  if(nMut < 1.0){		// we have to do a flip test on each bit
    nMut = 0;
    for(i=0; i<genome.height(); i++){
			muts.clear();
    	for(j=0; j<genome.height(); j++){
    		if(GAFlipCoin(pmut)){
    			muts.insert(j);
    		}
    	}
    	// determine which will be parents for this variable
			vector<int> parents;
			for(k=0; k<genome.width(); k++){
				int par = genome.gene(i,k);
				if(par == -1){
					continue;
				}
				if(muts.find(par) != muts.end()){
					muts.erase(par);
				}
				else{
					parents.push_back(par);
				}
			}
			for(set<int>::iterator iter=muts.begin(); iter != muts.end(); ++iter){
				parents.push_back(*iter);
			}

			if(parents.size() > genome.width()){
				set<int> keep=gaBayesCreator.limitConnections(i,parents,genome.width());
				int j=0;
				for(set<int>::iterator iter=keep.begin(); iter != keep.end(); ++iter){
					genome.gene(i,j,*iter);
					j++;
				}
			}
			else{
				size_t j=0;
				for(; j<parents.size(); j++){
					genome.gene(i,j,parents[j]);
				}
				for(; j<genome.width(); j++){
					genome.gene(i,j,-1);
				}
			}
    }
  }
  else{				// only flip the number of bits we need to flip
  	map<int, set<int> > muts;
    for(n=0; n<nMut; n++){
      m = GARandomInt(0, genome.size()-1);
      i = m % genome.height();
      j = m / genome.height();
			muts[j].insert(i);
    }
// 		for(i=0; i<genome.height(); i++){
		for(map<int, set<int> >::iterator mutIter = muts.begin(); mutIter != muts.end(); ++mutIter){
			vector<int> parents;
			i = mutIter->first;
			for(k=0; k<genome.width(); k++){
				int par = genome.gene(i,k);
				if(par == -1){
					continue;
				}
				if(muts[i].find(par) != muts[i].end()){
					muts[i].erase(par);
				}
				else{
					parents.push_back(par);
				}
			}
			for(set<int>::iterator iter=muts[i].begin(); iter != muts[i].end(); ++iter){
				parents.push_back(*iter);
			}

			if(parents.size() > genome.width()){
				set<int> keep=gaBayesCreator.limitConnections(i,parents,genome.width());
				int j=0;
				for(set<int>::iterator iter=keep.begin(); iter != keep.end(); ++iter){
					genome.gene(i,j,*iter);
					j++;
				}
			}
			else{
				size_t j=0;
				for(; j<parents.size(); j++){
					genome.gene(i,j,parents[j]);
				}
				for(; j<genome.width(); j++){
					genome.gene(i,j,-1);
				}
			}
		}

  }
  return(STA_CAST(int,nMut));
}

///
/// sets the Dataset for objective function to work with
///
void GAFunct::setDatasets(data_manage::Dataset* caseDS, data_manage::Dataset* controlDS,
	std::vector<Variable*> vList){
	caseDataset = caseDS;
	controlDataset = controlDS;
	varList=vList;

	caseBayesCreator.setMIScores(caseDataset,vList);
	caseBayesCreator.setNoParentScores(caseDataset,vList);

	controlBayesCreator.setMIScores(controlDataset,vList);
	controlBayesCreator.setNoParentScores(controlDataset, vList);

// cout << "output test case set" << endl;
// for(unsigned int i=0; i<caseDS->numInds(); i++){
// 	for(unsigned int j=0; j<caseDS->numGenos(); j++){
// 	  cout << caseDS->getInd(i)->getGenotype(j) << " ";
// 	}
// 	for(unsigned int j=0; j<caseDS->numCovariates(); j++){
// 	  cout << caseDS->getInd(i)->getCovariate(j) << " ";
// 	}
// 	cout << endl;
// }
// cout << "==================" << endl;

// 	if(!set->isCaseControl() && solCreator->getCalculator()->requiresCaseControl()){
// 		throw AthenaExcept(solCreator->getCalculator()->getName() + " requires a case-control dataset");
// 	}
// 	solCreator->setCalculatorConstant(ds);
}

///
/// sets the Dataset for objective function to work with
///
void GAFunct::setDataset(data_manage::Dataset* caseDS, std::vector<Variable*> vList){
	caseDataset = caseDS;
	varList=vList;

	caseBayesCreator.setMIScores(caseDataset,vList);
	caseBayesCreator.setNoParentScores(caseDataset,vList);
}





