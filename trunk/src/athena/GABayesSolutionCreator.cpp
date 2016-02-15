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

#include "GABayesSolutionCreator.h"
#include "Tarjans.h"
#include <cmath>
#include <algorithm>
#include "CycleBreaker.h"
//double GABayesSolutionCreator::calcK2Time = 0.0;

GABayesSolutionCreator::GABayesSolutionCreator(){
	calculator=NULL;
	maxParents=maxChildren=0;
	limitPtr=&GABayesSolutionCreator::removeLowMI;
	childLimitPtr=&GABayesSolutionCreator::removeLowChildMI;
}

GABayesSolutionCreator::~GABayesSolutionCreator(){
	if(calculator != NULL){
		delete calculator;
	}
}


void GABayesSolutionCreator::setNodeLimitMethod(string method){
	if(method.compare("BESTMI")==0){
		limitPtr=&GABayesSolutionCreator::removeLowMI;
		childLimitPtr=&GABayesSolutionCreator::removeLowChildMI;
	}
	else if(method.compare("RANDOM")==0){
		limitPtr=&GABayesSolutionCreator::removeRandom;
		childLimitPtr=&GABayesSolutionCreator::removeRandom;
	}
	else{
		throw AthenaExcept("No match for node limit method " + method);
	}
}


void GABayesSolutionCreator::breakLoops(Athena2DArrayGenome<int>& genome){
	CycleBreaker breaker;

// for(int i=0; i<genome.height(); i++){
// 	for(int j=0; j<genome.width(); j++){
// 		cout << " " << genome.gene(i,j);
// 	}
// 	cout << endl;
// }

// 	int cycBroken=breaker.breakCycles(genome);
	breaker.breakCycles(genome);
// 	cout << "broke " << cycBroken << endl;
// if(cycBroken > 0){
// for(int i=0; i<genome.height(); i++){
// 	for(int j=0; j<genome.width(); j++){
// 		cout << " " << genome.gene(i,j);
// 	}
// 	cout << endl;
// }
// 	cout << "============================" << endl;
// }

}

///
/// Fix loops in bayesian network by breaking the lowest information ones
/// @param g GA2DVinStrGenome g
///
void GABayesSolutionCreator::fixLoops(Athena2DArrayGenome<int>& genome){

	// genome is a matrix and width=height
	Tarjans tarj(genome.width());

// cout << endl;
// for(int i=0; i<genome.height(); i++){
// // 	cout << "i=" << i << " ";
// 	for(int j=0; j<genome.width(); j++){
// // 	cout <<  genome.gene(i,j) << " ";
// 		if(genome.gene(i,j))
// 		cout << "adjMatrix[" << i << "][" << j << "] = true;\n";
// 	}
// // 	cout << "\n";
// }
// cout << "\n";
// exit(1);
	// add edges
	for(int i=0; i<genome.height(); i++){
		for(int j=0; j<genome.width(); j++){
			// add edge for each connection
			if(genome.gene(i,j)){
				tarj.addEdge(i,j);
// cout << "add edge " << i << " to " << j << endl;
			}
		}
	}
// 	int fromVar, toVar;
// cout << "start SCC" << endl;
	while(tarj.SCC()){
		 std::vector<std::vector<int> > loops =tarj.getLoops();

// cout << "loops size=" << loops.size() << endl;
// for(size_t i=0; i<loops.size(); i++){
// cout << "loop " << i+1 << " => ";
// for(size_t j=0; j<loops[i].size(); j++){
//  cout << loops[i][j] << " ";
// }
// cout << endl;
// }

			int parIdx, childIdx;

			// find lowest MI for any parent->child in loop
		 for(std::vector<std::vector<int> >::iterator loopIter=loops.begin();
			loopIter != loops.end(); ++loopIter){
				double worstMI = 1e10;
				for(std::vector<int>::iterator varIter=(*loopIter).begin();
					varIter != (*loopIter).end(); ++varIter){
					for(std::vector<int>::iterator childIter=(*loopIter).begin(); childIter != (*loopIter).end();
						++childIter){
// cout << "look for " << *varIter << " connect to " << *childIter << endl;
						if(!genome.gene(*varIter, *childIter)){
							continue;
// cout << "no connection" << endl;
						}
// cout << "MI=" << miScores[*varIter][*childIter] << endl;
						if(miScores[*varIter][*childIter] < worstMI){
							worstMI=miScores[*varIter][*childIter];
							parIdx=*varIter;
							childIdx=*childIter;
						}
					}
				}
			}
// cout << "break " << parIdx << "->" << childIdx << endl;
		// break lowest MI in genome
		genome.gene(parIdx, childIdx, 0);
		tarj.removeEdge(parIdx, childIdx);
	}
// cout << "done fixing loops" << endl;

// for(int i=0; i<genome.height(); i++){
// 	cout << "i=" << i << " ";
// 	for(int j=0; j<genome.width(); j++){
// 	cout <<  genome.gene(i,j) << " ";
// 	}
// 	cout << "\n";
// }
}


void GABayesSolutionCreator::limitChildren(Athena2DArrayGenome<int>& genome){
	if(maxChildren <= 0)
		return;

	int genomeWidth=genome.width();
	int genomeHeight=genome.height();

// genome.write(cout);
	vector<int> empty;
	vector<vector<int> > children(genomeHeight, empty);
	for(int y=0; y<genomeHeight; y++){
		for(int x=0; x<genomeWidth; x++){
			if(genome.gene(x,y) != -1){
				children[genome.gene(x,y)].push_back(y);
			}
		}
	}

	for(size_t parIdx=0; parIdx < children.size(); parIdx++){
		if(children[parIdx].size() <= maxChildren)
			continue;
// cout << "parIdx=" << parIdx << "\t";
		set<int> keep=(this->*childLimitPtr)(parIdx,children[parIdx],maxChildren);
// cout << "keep=" << *(keep.begin()) << endl;
		for(vector<int>::iterator iter=children[parIdx].begin(); iter!=children[parIdx].end();
			++iter){
				int childIdx=*iter;
			if(keep.find(*iter)==keep.end()){
				for(int x=0; x<genomeWidth; x++){
					if(genome.gene(x,childIdx) == parIdx){
						genome.gene(x,childIdx,-1);
					}
				}
			}
		}
	}

}

void GABayesSolutionCreator::checkNodeLimits(Athena2DArrayGenome<int>& genome){
	vector<int> parents;
	int genomeWidth=genome.width();
	int genomeHeight=genome.height();


// for(int i=0; i<genomeHeight; i++){
// 	cout << i+1;
// 	for(int j=0; j<1; j++){
// 		cout << " " << genome.gene(i,j);
// 	}
// 	cout << endl;
// }

	// check for limits on number of parents
	if(maxParents > 0){
		for(int j=0; j<genomeWidth;j++){
			parents.clear();
			for(int i=0; i<genomeHeight; i++){
				if(genome.gene(i,j)){
					parents.push_back(i);
				}
			}
// cout << parents.size() << endl;
			if(parents.size() > maxParents){
// 				vector<int> keep=removeLowMI(j,parents,maxParents);
				set<int> keep = (this->*limitPtr)(j,parents,maxParents);
				for(vector<int>::iterator iter=parents.begin(); iter!=parents.end();
					++iter){
					if(keep.find(*iter)==keep.end()){
						genome.gene(*iter,j,0);
// cout << "removing " << *iter << " now=" << genome.gene(*iter,j) << endl;
					}
				}
// exit(1);
			}
		}
	}


// for(int i=0; i<1; i++){
// 	cout << i+1;
// 	for(int j=0; j<genomeWidth; j++){
// 		cout << " " << genome.gene(i,j);
// 	}
// 	cout << endl;
// }
	vector<int> children;
	if(maxChildren > 0){
		for(int i=0; i<genomeHeight; i++){
// cout << "parent=" << i << endl;
			children.clear();
			for(int j=0; j<genomeWidth; j++){
				if(genome.gene(i,j)){
					children.push_back(j);
				}
			}
// cout << "children=" << children.size() << endl;
			if(children.size() > maxChildren){
				set<int> keep=(this->*childLimitPtr)(i,children,maxChildren);
// cout << "keep size=" << keep.size() << endl;
				for(vector<int>::iterator iter=children.begin(); iter!=children.end();
					++iter){
					if(keep.find(*iter)==keep.end()){
						genome.gene(i,*iter,0);
// cout << "removing " << *iter << " now=" << genome.gene(i,*iter) << endl;
					}
				}
// exit(1);
			}
		}
	}
// exit(1);

}

///
/// Returns set with indexes of connections to keep
///
set<int> GABayesSolutionCreator::limitConnections(int childIndex, vector<int>& parents,
	int maxConn){
	return (this->*limitPtr)(childIndex,parents,maxConn);
}


set<int> GABayesSolutionCreator::removeRandom(int childIndex,vector<int>& parents,int maxConn){
	std::random_shuffle(parents.begin(), parents.end());
	std::set<int> keep;
	for(size_t i=0; i<maxConn; i++){
		keep.insert(parents[i]);
	}
	return keep;
}

set<int> GABayesSolutionCreator::removeLowMI(int childIndex,vector<int>& parents,int maxConn){
	// need to map the parents by MI score
	set<int> keep;
// 		int worstScore;

	map<float,int> parentsToKeep;
// cout << "removeLowMI" << endl;
	for(vector<int>::iterator iter=parents.begin(); iter != parents.end(); ++iter){
// cout << "parent=" << *iter << endl;
		if(parentsToKeep.size() < maxConn){
// cout << "added " << *iter << endl;
			parentsToKeep[miScores[*iter][childIndex]]=*iter;
		}
		// need to track scores with elements (need a map)
		// have to redo this part
		else{
		map<float,int>::iterator worstIter= parentsToKeep.begin();
// cout << "worst score=" << worstIter->first << " best score=" << parentsToKeep.rbegin()->first << endl;
// cout << "this score=" << miScores[*iter][childIndex] << endl;
			if(miScores[*iter][childIndex]>=worstIter->first){
// cout << "keep" << endl;
				parentsToKeep.erase(worstIter);
				parentsToKeep[miScores[*iter][childIndex]]=*iter;
			}
		}
	}
	for(map<float,int>::iterator mapIter=parentsToKeep.begin(); mapIter != parentsToKeep.end();
		++mapIter){
// cout << "final keep " << mapIter->second << endl;
		keep.insert(mapIter->second);
	}

	return keep;
}


set<int> GABayesSolutionCreator::removeLowChildMI(int parentIndex,vector<int>& children,int maxConn){
	set<int> keep;
// 	int worstScore;

	map<float,int> childrenToKeep;

	for(vector<int>::iterator iter=children.begin(); iter != children.end(); ++iter){
		if(childrenToKeep.size() < maxConn){
			childrenToKeep[miScores[parentIndex][*iter]]=*iter;
		}
		// need to track scores with elements (need a map)
		// have to redo this part
		else{
			map<float,int>::iterator worstIter = childrenToKeep.begin();
			if(miScores[parentIndex][*iter] >= worstIter->first){
				childrenToKeep.erase(worstIter);
				childrenToKeep[miScores[parentIndex][*iter]]=*iter;
			}
		}
	}
	for(map<float,int>::iterator mapIter=childrenToKeep.begin(); mapIter != childrenToKeep.end();
		++mapIter){
		keep.insert(mapIter->second);
	}
	return keep;
}


/// write out equation
void GABayesSolutionCreator::writeGenoNet(vector<vector<int> >& eq){
	for(size_t i=0; i<eq.size(); i++){
		cout << "[G" << i+1;
		if(!eq[i].empty()){
			cout << "|G" << (eq[i][0]+1);
		}
		for(size_t j=1; j<eq[i].size(); j++){
			cout << ":G" << (eq[i][j]+1);
		}
		cout << "]";
	}
	cout << "\n";
}

///
/// Construct equation string from genome
/// @genome Athena2DArrayGenome<int>
///
vector<vector<int> > GABayesSolutionCreator::constructEquation(Athena2DArrayGenome<int>& genome,
	 std::vector<Variable*> varList){
	vector<int> empty;
	vector<vector<int> > conns(varList.size(), empty);
	for(int i=0; i<genome.height(); i++){
		for(int j=0; j<genome.width(); j++){
			if(genome.gene(i,j)==1){
				conns[j].push_back(i);
			}
		}
	}
	return conns;
}

///
/// prunes network by removing connections that do not improve the overall network score
/// @param conns contains connections for network -- changes to network are done in place
/// @param varList
/// @param dSet
/// @returns score of pruned network
///
float GABayesSolutionCreator::pruneNetwork(vector<vector<int> >& conns, vector<Variable*>& varList,
		data_manage::Dataset* dSet){


// map<int, string> convert;
// convert[0]="a";
// convert[1]="b";
// convert[2]="c";
// data_manage::Dataholder* holder = dSet->getHolder();
// cout << "G" + varList[0]->getName(holder);
// for(unsigned int j=1; j<varList.size(); j++){
// 	cout << "\tG" + varList[j]->getName(holder);
// }
// cout << endl;
//
// for(unsigned int i=0; i<dSet->numInds(); i++){
// 	data_manage::Individual * ind=dSet->getInd(i);
// 	cout << convert[ind->getGenotype(varList[0]->getIndex())];
// 	for(unsigned int j=1; j<varList.size(); j++){
// 		cout << "\t" << convert[ind->getGenotype(varList[j]->getIndex())];
// 	}
// 	cout << endl;
// }
// exit(1);

// cout << "Connections:" << endl;
// for(size_t i=0; i<conns.size(); i++){
// 	if(conns[i].empty())
// 		continue;
// 	cout << i << ": ";
// 	for(size_t j=0; j<conns[i].size(); j++){
// 		cout << conns[i][j] << " ";
// 	}
// 	cout << endl;
// }

// 	float score=0.0;
	float nodeScore, connScore;
	calculator->reset();
	int nParams,thisParams;
	vector<int> parents;

	for(size_t childIdx=0; childIdx < conns.size();childIdx++){
		// when no connections use no parent score and go to next variable
		if(conns[childIdx].empty()){
// cout << "childIdx=" << childIdx << " empty" << endl;
			calculator->addIndScore(noParentScores[childIdx],0);
			continue;
		}
		parents.clear();
		set<int> nodeParents;
		for(size_t parentIdx=0; parentIdx < conns[childIdx].size(); parentIdx++){
// 			parents.push_back(conns[childIdx][parentIdx]);
			nodeParents.insert(conns[childIdx][parentIdx]);
		}

		for(size_t parentIdx=0; parentIdx<conns[childIdx].size(); parentIdx++){
// 			parents.clear();
			// score with all variables
			parents.assign(nodeParents.begin(), nodeParents.end());
			nodeScore = k2Calc(childIdx,parents,varList,dSet,nParams);
// cout << childIdx << " nodeScore=" << nodeScore;
			nodeParents.erase(conns[childIdx][parentIdx]);
			// score without this connection
			if(nodeParents.empty()){
				connScore = noParentScores[childIdx];
				thisParams = 0;
			}
			else{
				parents.assign(nodeParents.begin(), nodeParents.end());
				connScore = k2Calc(childIdx,parents,varList,dSet,thisParams);
			}
// cout << " connScore=" << connScore;
			if(connScore > nodeScore){ // if the score is better without the conection
// cout << " replace ";
				nodeScore = connScore;
				nParams = thisParams;
			}
			else{
// cout << " keep ";
				// put the connection back into set
				nodeParents.insert(conns[childIdx][parentIdx]);
			}
		}
// cout << endl;
		// set connections in original structure to final set
		conns[childIdx].assign(nodeParents.begin(), nodeParents.end());
		calculator->addIndScore(nodeScore,nParams);
	}
// cout << " final score=" << -calculator->getScore() << endl;

// cout << "Connections:" << endl;
// for(size_t i=0; i<conns.size(); i++){
// 	if(conns[i].empty())
// 		continue;
// 	cout << i << ": ";
// 	for(size_t j=0; j<conns[i].size(); j++){
// 		cout << conns[i][j] << " ";
// 	}
// 	cout << endl;
// }
// cout << "=====================================" << endl;
// 	return score;
	return -calculator->getScore();
}


///
/// Calculate and return network score
/// @param genome
/// @param varList
/// @param dSet
/// @return network score
///
double GABayesSolutionCreator::calcScore(Athena2DArrayGenome<int>& genome, vector<Variable*> varList,
	data_manage::Dataset* dSet){

// 	double totalScore=0.0;
	double score;
	int nParams;
	calculator->reset();


//time_t startTime, endTime;
//calcK2Time=0.0;
// cout << "\n";
// cout << "   ";
// for(int i=0; i<genome.height(); i++){
// 	cout << " " << i;
// }
// cout << "\n-----------------------\n";
// for(int i=0; i<genome.height(); i++){
// 	cout << "i=" << i << " ";
// 	for(int j=0; j<genome.width(); j++){
// 	cout <<  genome.gene(i,j) << " ";
// 	}
// 	cout << "\n";
// }
//
// vector<vector<int> > eq = constructEquation(genome,varList);
// writeGenoNet(eq);


// int j=1;
// vector<int>pars;
// pars.push_back(2);
// cout << k2Calc(j,pars,varList,dSet,nParams) << endl;
// exit(1);
	for(int y=0; y<genome.height(); y++){
		vector<int> parents;
// cout << "j=" << j << endl;
		for(int x=0; x<genome.width(); x++){
			// connection
			if(genome.gene(x,y) != -1){
// cout << "parent  " << i << " for " << j << endl;
				parents.push_back(genome.gene(x,y));
			}
		}

// cout << "ch=G" << j+1 << " parents: ";
// for(size_t p=0; p<parents.size(); p++){
// cout << "G" << parents[p] << " ";
// }

		if(parents.empty()){
			// the additional parameter term (if any) is already included in noParentScores
			calculator->addIndScore(noParentScores[y],0);
// cout << "score=" << noParentScores[j] << "\n";
		}
		else{
			score = k2Calc(y,parents,varList,dSet,nParams);
//time(&endTime);
//calcK2Time += difftime (endTime,startTime);
			calculator->addIndScore(score, nParams);
		}
	}
//exit(1);
// cout << "final score=" << -calculator->getScore() << "\n";
	return -calculator->getScore();
}


///
/// Calculate k2 score when parents connected to node
/// @param childIdx
/// @param parIndexes
/// @param varList
/// @param dSet
///
double GABayesSolutionCreator::k2Calc(int childIdx, vector<int>& parIndexes,
	vector<Variable*> varList, data_manage::Dataset* dSet, int& nP){

	nP=0;
	Variable * node = varList[childIdx];
	vector<Variable*> parents;
	for(size_t i=0; i<parIndexes.size(); i++){
		parents.push_back(varList[parIndexes[i]]);
	}

	// construct table with node values and the values for the parent combinations
	vector<int> parentValues;
// cout << " start configPar "; cout.flush();
	int parentLevels = configParentData(parentValues, parents, dSet);
// cout << " end configPar "; cout.flush();
	int nodeLevels = node->getNumLevels();
	int imaginary = 0;
	float alpha = 1.0;
	unsigned int numInds = dSet->numInds();
	vector<int> inner(parentLevels,0);
	vector<vector<int> > totals(nodeLevels, inner);
	vector<int> parentTotals(parentLevels, 0);
	vector<int> nodeTotals(nodeLevels, 0);
	data_manage::Individual* ind;
	for(unsigned int i=0; i < numInds; i++){
		ind = (*dSet)[i];
		totals[node->getValue(ind)][parentValues[i]]++;
		parentTotals[parentValues[i]]++;
		nodeTotals[node->getValue(ind)]++;
	}
// cout << " start condit "; cout.flush();
	// calculate conditional posterior probability
	double score = 0.0;
	for(int i=0; i<nodeLevels; i++){
		for(int j=0; j<parentLevels; j++){
			score += lgamma(totals[i][j] + alpha); // - lgamma(alpha) /* always zero for k2 */
		}
	}
	int finalParentLevels = parentLevels;
	for(int j=0; j < parentLevels; j++){
		if(parentTotals[j] == 0){
			finalParentLevels--;
		}
	}
	int finalNodeLevels=0;
	for(int j=0; j < nodeLevels; j++){
		if(nodeTotals[j] > 0)
			finalNodeLevels++;
	}
	imaginary = finalNodeLevels * finalParentLevels;
	nP = (finalNodeLevels-1) * finalParentLevels;
	for(int j=0; j < parentLevels; j++){
		if(parentTotals[j] > 0){
			score += lgamma(double(imaginary)/finalParentLevels) -
				lgamma(parentTotals[j] + double(imaginary)/finalParentLevels);
		}
	}
	return score;
}

///
/// Create parent data combination
/// @param parentValues [out]
/// @param parents
/// @returns number of different levels(factors) in the parent combined values
///
int GABayesSolutionCreator::configParentData(vector<int>& parentValues,
	vector<Variable*> &parents, data_manage::Dataset* dSet){

	// set number of levels for each parent
	vector<int> nLevels(parents.size(), 0);
	int nl = 1;
	for(size_t i=0; i<parents.size(); i++){
		nLevels[i] = parents[i]->getNumLevels();
		nl *= nLevels[i];
	}

	vector<int> cumulativeLevels(parents.size(), 1);
	for(unsigned int i=1; i<cumulativeLevels.size(); i++){
// 		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels;
		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels[i-1];
	}
// 	int nl = cumulativeLevels.back() * nLevels;
	deque<float> args;
	unsigned int nParents = parents.size();
	data_manage::Individual* ind;

	for(unsigned int i=0; i < dSet->numInds(); i++){
		ind = (*dSet)[i];
		int value = 0;
		for(unsigned int j=0; j < nParents; j++){
			value += parents[j]->getValue(ind) * cumulativeLevels[j];
		}
		parentValues.push_back(value);
	}
	return nl;
}

///
/// Calculates and stores Mutual Information scores for use in breaking loops.
/// @param ds
/// @param varList
///
void GABayesSolutionCreator::setMIScores(data_manage::Dataset* ds, std::vector<Variable*>& varList){
	size_t gridSize=varList.size();
	vector<double> row(gridSize, 0.0);
	miScores.assign(gridSize, row);
	for(size_t i=0; i<gridSize; i++){
		for(size_t j=0; j<gridSize; j++){
			if(i==j){
				continue;
			}
			miScores[i][j]=calcMI(varList[i], varList[j], ds);
// cout << "parent=" << i << " child=" << j << " MI=" << miScores[i][j] << endl;
		}
	}
}

///
/// Calculate mutual information for the two variables passed on the dataset
///
double	GABayesSolutionCreator::calcMI(Variable* parentVar, Variable* childVar, data_manage::Dataset* ds){
		double mutualInfo=0.0;

// 	vector<IndividualTerm*> terms(2,NULL);
// 	terms[0] = static_cast<IndividualTerm*>(v1);
// 	terms[1] = static_cast<IndividualTerm*>(v2);
// 	int parentLevels = configParentData(combinedValues, terms);

	// need total of each individual
	// need total of combination
	// loop through each individuals and add the individual totals
	// add the combination total
	// assume three levels (to hold SNP data)
// 	int nLevels = 3; // this has to be changed also

	vector<Variable*> terms(2,NULL);
	terms[0] = parentVar;
	terms[1] = childVar;

	// set number of levels for each parent
	vector<int> nLevels(terms.size(), 0);
	int nl = 1;
	for(size_t i=0; i<terms.size(); i++){
		nLevels[i] = terms[i]->getNumLevels();
		nl *= nLevels[i];
	}

	vector<int> combinedTotals(nl,0);
	size_t v1Levels = terms[0]->getNumLevels();
	size_t v2Levels = terms[1]->getNumLevels();
	vector<int> v1Totals(v1Levels, 0);
	vector<int> v2Totals(v2Levels, 0);


	vector<int> cumulativeLevels(terms.size(), 1);
	for(unsigned int i=1; i<cumulativeLevels.size(); i++){
		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels[i-1];
	}

// 	unsigned int nTerms = terms.size();
	data_manage::Individual* ind;
	unsigned int N = ds->numInds();

	for(unsigned int i=0; i < ds->numInds(); i++){
		ind = (*ds)[i];

		int value = 0;
		value += terms[0]->getValue(ind) * cumulativeLevels[0];
		v1Totals[terms[0]->getValue(ind)]++;

		value += terms[1]->getValue(ind) * cumulativeLevels[1];
		v2Totals[terms[1]->getValue(ind)]++;

		combinedTotals[value]++;
	}

	// calculate mutual information
	// for each combination
	int index;
	for(unsigned int i=0; i<v1Levels; i++){
		for(unsigned int j=0; j<v2Levels; j++){
			index = i * cumulativeLevels[0] + j * cumulativeLevels[1];
// cout << "v1 " << i << "=" << v1Totals[i] << endl;
// cout << "v2 " << j << "=" << v1Totals[j] << endl;
			if(combinedTotals[index] > 0){
				mutualInfo += combinedTotals[index] / double(N) * log(N*combinedTotals[index]/double(v1Totals[i]*v2Totals[j]));
			}
// cout << "i=" << i << " j=" << j << " tot=" << combinedTotals[index] << " mi=" << mutualInfo << endl;
// cout << "first part=" << combinedTotals[index] / double(N) << endl;
// cout <<  "second part without log=" << N*combinedTotals[index]/double(v1Totals[i]*v2Totals[j]) << endl;
// cout << "second part=" << log(N*combinedTotals[index]/double(v1Totals[i]*v2Totals[j])) << endl;
		}
	}
	return mutualInfo;
}

///
/// store scores for each variable in case where it has no parents
///
void GABayesSolutionCreator::setNoParentScores(data_manage::Dataset* ds, vector<Variable*>&  vList){

	noParentScores.clear();
	varParams.clear();
	double total=0.0;

	int nP;
	for(vector<Variable*>::iterator varIter=vList.begin(); varIter != vList.end();
		++varIter){
		calculator->reset();
		noParentScores.push_back(k2CalcNoParent(*varIter,ds,nP));
		varParams.push_back(nP);
		calculator->addIndScore(noParentScores.back(), nP);
		noParentScores.back() = -calculator->getScore();
		total += noParentScores.back();
// cout  << "var " << (*varIter)->getIndex() << " K2 score=" << noParentScores.back()<< endl;
	}
	ds->setConstant(total);
}


///
/// Calculate k2 score
/// @param var Variable to analyze
/// @ds Dataset to use
/// @param nP (out) number of parameters in calculation
/// @returns k2 score for node
///
double GABayesSolutionCreator::k2CalcNoParent(Variable* var, data_manage::Dataset* ds,
	int& nP){
		// create empty vector to hold totals
	vector<unsigned int> totals(var->getNumLevels(),0);
	unsigned int nInds = ds->numInds();
	double result = 0.0;

	for(unsigned int i=0; i < nInds; i++){
// 		totals[int(ds->getInd(i)->getCovariate(cIndex))]++;
		totals[int(var->getValue(ds->getInd(i)))]++;
	}

	nP =-1;
	// alpha = 1
	// imaginary is simply number of levels of totals for K2
	double imaginary = double(totals.size());
	double alpha = 1.0;
	for(unsigned int i=0; i<totals.size(); i++){
		result += lgamma(totals[i]+alpha); // no need to get gamma of alpha (it is 0)
		// reduce imaginary when a cell is empty
		if(totals[i]==0)
			imaginary--;
		else
			nP++;
	}
	result += lgamma(imaginary) - lgamma(imaginary + nInds);
	return result;
}



