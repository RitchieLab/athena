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

#include "Athena2DArrayGenome.h" // and the 2D genom
#include "GABayesSolutionCreator.h"
#include <Dataset.h>

///
/// Provides objective function for use with GALIb and GE library
///
class GAFunct{

public:
	/// Used to assign fitness values in GA
	static float GACaseObjective(GAGenome& g);

	static float GAObjective(GAGenome &g);

	/// Used to assign fitness values in GA
// 	static float GAControlObjective(GAGenome& g);

	/// conducts random init
	static void init(GAGenome& g);

	/// conducts random initialization of genomes
	static void initCase(GAGenome &g);

	/// conducts random initialization of genomes
// 	static void initControl(GAGenome &g);

	/// sets rank
	static void setInitConP(float p){connProb=p;}

	/// sets the Dataset for objective function to work with
// 	static void setDatasets(data_manage::Dataset* caseDS, data_manage::Dataset* controlDS,
// 		std::vector<Variable*> vList, bool needMI=true);

	/// sets Solution type for objective function
	static void setSolutionType(std::string calculatorName){
// 				solCreator = SolutionFactory::createSolution(solutionName);
		caseBayesCreator.setCalculator(calculatorName);
		for(size_t i=0; i<bayesCreators.size(); i++)
			bayesCreators[i].setCalculator(calculatorName);
// 		controlBayesCreator.setCalculator(calculatorName);
	}

	static double getWorstScore(){return caseBayesCreator.getWorst();}

	static Solution* getBlankSolution(){return caseBayesCreator.createNewSolution();}

	static vector<std::string> getAdditionalFinalOutput(float score);

	static vector<std::string>  getAdditionalOutputNames();

	static void setDataset(data_manage::Dataset* caseDS, std::vector<Variable*> vList);

	static void setDatasets(vector<data_manage::Dataset*> categorySets,
		std::vector<Variable*> vList, bool needMI);

	static void setNodeMaximums(int maxP, int maxC){
		caseBayesCreator.setNodeMax(maxP,maxC);
// 		controlBayesCreator.setNodeMax(maxP,maxC);
	}

	static void setNodeLimitMethod(std::string method){
		caseBayesCreator.setNodeLimitMethod(method);
// 		controlBayesCreator.setNodeLimitMethod(method);
	}

	static int mutate(GAGenome & c, float pmut);
	static int mutateCase(GAGenome & c, float pmut);
// 	static int mutateControl(GAGenome & c, float pmut);

	static float prune(vector<vector<int> >& conns, int category);
	static float pruneCase(vector<vector<int> >& conns);
// 	static float pruneControls(vector<vector<int> >& conns);

// static double fitnessTime, loopTime, maxCheckTime;
// 	static int initConnections, dupConnections, limitChildConnections, brokenLoopConnections;

private:
	static void init(GAGenome &g, GABayesSolutionCreator& gaBayesCreator);
	static void removeSelfAndDup(Athena2DArrayGenome<int>& genome);
	static int customMutator(GAGenome & c, float pmut, GABayesSolutionCreator& gaBayesCreator);
	static int countConnections(Athena2DArrayGenome<int>& genome);

	static data_manage::Dataset* caseDataset;//, *controlDataset;
	static vector<data_manage::Dataset*> datasets;
	static float connProb;
	static vector<Variable*> varList;
	static GABayesSolutionCreator caseBayesCreator; //controlBayesCreator;
	static vector<GABayesSolutionCreator> bayesCreators;


};
#endif
