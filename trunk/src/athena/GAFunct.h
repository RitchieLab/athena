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

#ifndef _GAFUNCT_H
#define	_GAFUNCT_H

#include <ga/GA2DBinStrGenome.h> // and the 2D binary string genome
#include "GABayesSolutionCreator.h"
#include <Dataset.h>

///
/// Provides objective function for use with GALIb and GE library
///
class GAFunct{

public:
	/// Used to assign fitness values in GA
	static float GACaseObjective(GAGenome& g);

	/// Used to assign fitness values in GA
	static float GAControlObjective(GAGenome& g);

	/// conducts random initialization of genomes
	static void initCase(GAGenome &g);

	/// conducts random initialization of genomes
	static void initControl(GAGenome &g);

	/// sets rank
	static void setInitConP(float p){connProb=p;}

	/// sets the Dataset for objective function to work with
	static void setDatasets(data_manage::Dataset* caseDS, data_manage::Dataset* controlDS,
		std::vector<Variable*> vList);

	/// sets Solution type for objective function
	static void setSolutionType(std::string calculatorName){
// 				solCreator = SolutionFactory::createSolution(solutionName);
		caseBayesCreator.setCalculator(calculatorName);
		controlBayesCreator.setCalculator(calculatorName);
	}

	static double getWorstScore(){return caseBayesCreator.getWorst();}

	static Solution* getBlankSolution(){return caseBayesCreator.createNewSolution();}

	static vector<std::string> getAdditionalFinalOutput(float score);

	static vector<std::string>  getAdditionalOutputNames();

	static void setDataset(data_manage::Dataset* caseDS, std::vector<Variable*> vList);

	static void setNodeMaximums(int maxP, int maxC){
		caseBayesCreator.setNodeMax(maxP,maxC);
		controlBayesCreator.setNodeMax(maxP,maxC);
	}

	static void setNodeLimitMethod(std::string method){
		caseBayesCreator.setNodeLimitMethod(method);
		controlBayesCreator.setNodeLimitMethod(method);
	}
static double fitnessTime, loopTime, maxCheckTime;

private:
	static void init(GAGenome &g, GABayesSolutionCreator& gaBayesCreator);
	static void removeSelfConns(GA2DBinaryStringGenome& genome);

	static data_manage::Dataset* caseDataset, *controlDataset;
	static float connProb;
	static vector<Variable*> varList;
	static GABayesSolutionCreator caseBayesCreator, controlBayesCreator;


};
#endif
