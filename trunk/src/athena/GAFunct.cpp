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

float	GAFunct::connProb =0.5;
data_manage::Dataset* GAFunct::caseDataset = NULL;
data_manage::Dataset* GAFunct::controlDataset = NULL;
vector<Variable*> GAFunct::varList;
GABayesSolutionCreator GAFunct::caseBayesCreator;
GABayesSolutionCreator GAFunct::controlBayesCreator;

///
/// Used to assign fitness values in GA
/// @param g Genome to evaluate
///
float GAFunct::GACaseObjective(GAGenome& g){
  GA2DBinaryStringGenome & genome = (GA2DBinaryStringGenome &)g;
  removeSelfConns(genome);
  caseBayesCreator.fixLoops(genome);
// cout << "calculating case fitness" << endl;
// float score=caseBayesCreator.calcScore(genome, varList, caseDataset);
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

// return score;
	return caseBayesCreator.calcScore(genome, varList, caseDataset);
}

///
/// Used to assign fitness values in GA
/// @param g Genome to evaluate
///
float GAFunct::GAControlObjective(GAGenome& g){
  GA2DBinaryStringGenome & genome = (GA2DBinaryStringGenome &)g;
  removeSelfConns(genome);
  controlBayesCreator.fixLoops(genome);
// cout << "calculating control fitness" << endl;
	return controlBayesCreator.calcScore(genome, varList, controlDataset);
}


vector<std::string> GAFunct::getAdditionalFinalOutput(GAGenome& g){
	vector<std::string> addOutputValues;
	if(g.score() > caseDataset->getConstant() ){
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
	GA2DBinaryStringGenome & genome = (GA2DBinaryStringGenome &)g;

// cout << "genome height=" << genome.height() << endl;
// cout << "genome width=" << genome.width() << endl;
// exit(1);

	for(int i=0; i<genome.height(); i++){
		for(int j=0;j< genome.width(); j++){
float r=GARandomFloat();
// cout << r << endl;
// 			if(GARandomFloat() > connProb)
			if(r > connProb)
				genome.gene(i,j,0);
			else{
				genome.gene(i,j,1);
// 				cout << "CONNECTED" << endl;
			}
		}
	}

// 	removeSelfConns(genome);

// cout << "\n";
// for(int i=0; i<genome.height(); i++){
// 	cout << "i=" << i << " ";
// 	for(int j=0; j<genome.width(); j++){
// 	cout <<  genome.gene(i,j) << " ";
// 	}
// 	cout << "\n";
// }

// 	gaBayesCreator.fixLoops(genome);
}

///
/// unset any connections to self in matrix
/// @param g GAGenome to initialize
///
void GAFunct::removeSelfConns(GA2DBinaryStringGenome& genome){
	for(int i=0; i<genome.width(); i++){
		genome.gene(i,i,0);
	}
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





