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
// data_manage::Dataset* GAFunct::controlDataset = NULL;
vector<data_manage::Dataset*> GAFunct::datasets;
vector<Variable*> GAFunct::varList;
GABayesSolutionCreator GAFunct::caseBayesCreator;
// GABayesSolutionCreator GAFunct::controlBayesCreator;
vector<GABayesSolutionCreator> GAFunct::bayesCreators;
//double GAFunct::fitnessTime = 0.0;
//double GAFunct::loopTime = 0.0;
//double GAFunct::maxCheckTime = 0.0;
//double GAFunct::calcK2Time = 0.0;
// int GAFunct::initConnections=0;
// int GAFunct::dupConnections=0;
// int GAFunct::limitChildConnections=0;
// int GAFunct::brokenLoopConnections=0;

///
/// Used to assign fitness values in GA
/// @param g Genome to evaluate
///
float GAFunct::GACaseObjective(GAGenome& g){

// 	time_t startTime, endTime;

  Athena2DArrayGenome<int> & genome = (Athena2DArrayGenome<int> &)g;

// initConnections+=countConnections(genome);

  removeSelfAndDup(genome);
//time (&startTime);
//   caseBayesCreator.checkNodeLimits(genome);
  caseBayesCreator.limitChildren(genome);
//time (&endTime);
// double dif = difftime (endTime,startTime);
//maxCheckTime += dif;
// dupConnections+=countConnections(genome);
///time (&startTime);
//   caseBayesCreator.fixLoops(genome);
// limitChildConnections+=countConnections(genome);
	caseBayesCreator.breakLoops(genome);
// brokenLoopConnections+=countConnections(genome);
//time(&endTime);
//dif = difftime (endTime,startTime);
//loopTime += dif;

//time(&endTime);
//dif = difftime (endTime,startTime);
//loopTime += dif;

// cout << "calculating case fitness" << endl;
//time(&startTime);
	float score=caseBayesCreator.calcScore(genome, varList, caseDataset);
//time(&endTime);
//dif = difftime (endTime,startTime);
//fitnessTime += dif;
//calcK2Time += GABayesSolutionCreator::calcK2Time;
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
/// counts total connections
///
int GAFunct::countConnections(Athena2DArrayGenome<int>& genome){
	int conns=0;
	for(int y=0; y<genome.height(); y++){
		for(int x=0; x<genome.width(); x++){
			if(genome.gene(x,y)!=-1){
				conns++;
			}
		}
	}
	return conns;
}


///
/// Used to assign fitness values in GA
/// @param g Genome to evaluate
///
// float GAFunct::GAControlObjective(GAGenome& g){
//   Athena2DArrayGenome<int> & genome = (Athena2DArrayGenome<int> &)g;
//   removeSelfAndDup(genome);
//   controlBayesCreator.checkNodeLimits(genome);
// //   controlBayesCreator.fixLoops(genome);
//   controlBayesCreator.breakLoops(genome);
// // cout << "calculating control fitness" << endl;
// 	return controlBayesCreator.calcScore(genome, varList, controlDataset);
// }


///
/// Used to assign fitness values in GA
/// @param g Genome to evaluate
///
float GAFunct::GAObjective(GAGenome &g){
	Athena2DArrayGenome<int> & genome = (Athena2DArrayGenome<int> &)g;
  removeSelfAndDup(genome);
	int category = genome.category();
	bayesCreators[category].checkNodeLimits(genome);
	bayesCreators[category].breakLoops(genome);
	return bayesCreators[category].calcScore(genome, varList, datasets[category]);
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
// void GAFunct::initControl(GAGenome &g){
// 	init(g, controlBayesCreator);
// }


/// conducts initialization using appropriate bayes creator
void GAFunct::init(GAGenome &g){
	Athena2DArrayGenome<int> & genome = (Athena2DArrayGenome<int> &)g;
	init(g, bayesCreators[genome.category()]);
}


///
/// conducts random initialization of genomes
/// @param g GAGenome to initialize
/// @param gaBayesCreator GABayesSolutionCreator to use in initialization
///
void GAFunct::init(GAGenome &g, GABayesSolutionCreator& gaBayesCreator){
	Athena2DArrayGenome<int> & genome = (Athena2DArrayGenome<int> &)g;
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
void GAFunct::removeSelfAndDup(Athena2DArrayGenome<int>& genome){
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


/// conducts mutation of genomes
int GAFunct::mutateCase(GAGenome & c, float pmut){
	return customMutator(c, pmut, caseBayesCreator);
}

/// conducts mutation of genomes
// int GAFunct::mutateControl(GAGenome & c, float pmut){
// 	return customMutator(c, pmut, controlBayesCreator);
// }

/// conducts mutation of genomes
int GAFunct::mutate(GAGenome& c, float pmut){
	Athena2DArrayGenome<int> & genome = (Athena2DArrayGenome<int> &)c;
	return customMutator(c, pmut, bayesCreators[genome.category()]);
}

///
/// Custom mutator has same effect as flipmutator for binary string genome
///
int GAFunct::customMutator(GAGenome & c, float pmut, GABayesSolutionCreator& gaBayesCreator){
  Athena2DArrayGenome<int> &genome=DYN_CAST(Athena2DArrayGenome<int> &, c);
  register int n, m, i, j, k;
  if(pmut <= 0.0) return(0);

	int matrixSize=genome.height()*genome.height();
  float nMut = pmut * matrixSize;


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
				int par = genome.gene(k,i);
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
					genome.gene(j,i,*iter);
					j++;
				}
			}
			else{
				size_t j=0;
				for(; j<parents.size(); j++){
					genome.gene(j,i,parents[j]);
				}
				for(; j<genome.width(); j++){
					genome.gene(j,i,-1);
				}
			}
    }
  }
  else{				// only flip the number of bits we need to flip
  	map<int, set<int> > muts;
  	int matrixSize=genome.height()*genome.height();
    for(n=0; n<nMut; n++){
      m = GARandomInt(0, matrixSize-1);
      i = m % genome.height();
      j = m / genome.height();
			muts[i].insert(j);
    }

		for(map<int, set<int> >::iterator mutIter = muts.begin(); mutIter != muts.end(); ++mutIter){
			vector<int> parents;
			i = mutIter->first;
			for(k=0; k<genome.width(); k++){
				int par = genome.gene(k,i);
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
					genome.gene(j,i,*iter);
					j++;
				}
			}
			else{
				size_t j=0;
				for(; j<parents.size(); j++){
					genome.gene(j,i,parents[j]);
				}
				for(; j<genome.width(); j++){
					genome.gene(j,i,-1);
				}
			}
		}

  }
  return(STA_CAST(int,nMut));
}

///
/// Prunes network connections for the appropriate dataset
/// @returns updated score for network
///
float GAFunct::prune(vector<vector<int> >& conns, int category){
	return bayesCreators[category].pruneNetwork(conns, varList, datasets[category]);
}

///
/// Prunes network connections for a case dataset network
/// @returns updated score for network
///
float GAFunct::pruneCase(vector<vector<int> >& conns){
	return caseBayesCreator.pruneNetwork(conns, varList, caseDataset);
}

///
/// Prunes network connections for a control dataset network
/// @returns updated score for network
///
// float GAFunct::pruneControls(vector<vector<int> >& conns){
// 	return controlBayesCreator.pruneNetwork(conns, varList, controlDataset);
// }


///
/// sets the Dataset for objective function to work with
///
// void GAFunct::setDatasets(data_manage::Dataset* caseDS, data_manage::Dataset* controlDS,
// 	std::vector<Variable*> vList, bool needMI){
// 	caseDataset = caseDS;
// 	controlDataset = controlDS;
// 	varList=vList;
//
// 	if(needMI)
// 		caseBayesCreator.setMIScores(caseDataset,vList);
// 	caseBayesCreator.setNoParentScores(caseDataset,vList);
//
// 	if(needMI)
// 		controlBayesCreator.setMIScores(controlDataset,vList);
// 	controlBayesCreator.setNoParentScores(controlDataset, vList);

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
// }


///
/// sets the Dataset for objective function to work with
///
void GAFunct::setDatasets(vector<data_manage::Dataset*> categorySets,
	std::vector<Variable*> vList, bool needMI){

// 	GABayesSolutionCreator newCreator;
	datasets = categorySets;
	varList=vList;

	bayesCreators.assign(datasets.size(),caseBayesCreator);

	for(size_t i=0; i<datasets.size(); i++){
		if(needMI)
			bayesCreators[i].setMIScores(datasets[i],vList);
		bayesCreators[i].setNoParentScores(datasets[i], vList);
	}
// exit(1);
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





