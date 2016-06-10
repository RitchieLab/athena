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

#ifndef _GABAYESSOLUTIONCREATOR_H
#define	_GABAYESSOLUTIONCREATOR_H

#include "Athena2DArrayGenome.h"
#include <Dataset.h>
#include "Variable.h"
#include "SolutionCalculator.h"
#include "GABayesSolution.h"
#include <set>

///
/// Provides objective function for use with GALIb and GE library
///
class GABayesSolutionCreator{

public:

	GABayesSolutionCreator();

	~GABayesSolutionCreator();

	GABayesSolutionCreator(const GABayesSolutionCreator& orig);

	GABayesSolutionCreator& operator=(const GABayesSolutionCreator& orig);

	void fixLoops(Athena2DArrayGenome<int>& g);

	void checkNodeLimits(Athena2DArrayGenome<int>& g);

	void setMIScores(data_manage::Dataset* ds, std::vector<Variable*>&  vList);

	/// store scores
	void setNoParentScores(data_manage::Dataset* ds, std::vector<Variable*>&  vList);

	/// adds a solution calculator for determining fitness
	virtual void setCalculator(std::string calc_type){
		if(calculator != NULL)
			delete calculator;
		calculator = CalculatorFactory::getFactory().create(calc_type);
	}

	/// Calculate and return network score
	double calcScore(Athena2DArrayGenome<int>& genome, std::vector<Variable*> varList,
		data_manage::Dataset* dSet);

		/// Returns worst score
		inline float getWorst(){return calculator->getWorst();}

		virtual Solution* createNewSolution(){
				GABayesSolution* sol = new GABayesSolution;
				return sol;}
	void setNodeMax(int maxP, int maxC){
		maxParents=maxP;
		maxChildren=maxC;
	}

	void setNodeLimitMethod(std::string method);

	set<int> limitConnections(int childIndex, vector<int>& parents,
		int maxConn);
//static double calcK2Time;

	void breakLoops(Athena2DArrayGenome<int>& genome);

	void limitChildren(Athena2DArrayGenome<int>& genome);

	/// prunes network by removing connections that do not improve the overall network score
	float pruneNetwork(vector<vector<int> >& conns, vector<Variable*>& varList,
		data_manage::Dataset* dSet);

private:
	double calcMI(Variable* parentVar, Variable* childVar, data_manage::Dataset* ds);

	int configParentData(std::vector<int>& parentValues, std::vector<Variable*> &parents,
		data_manage::Dataset* dSet,std::set<int>& missingVals);

	double k2Calc(int childIdx, vector<int>& parIndexes,
		vector<Variable*> varList, data_manage::Dataset* dSet, int& nP);

	double k2CalcNoParent(Variable* var, data_manage::Dataset* ds, int& nP);

	vector<vector<int> > constructEquation(Athena2DArrayGenome<int>& genome, std::vector<Variable*> varList);
	void writeGenoNet(vector<vector<int> >& eq);

	std::set<int> removeLowMI(int childIndex,vector<int>& parents,int maxConn);
	std::set<int> removeLowChildMI(int childIndex,vector<int>& parents,int maxConn);
	std::set<int> removeRandom(int childIndex,vector<int>& parents,int maxConn);

	void	setMissingIndexes(std::vector<Variable*> &parents, std::set<int>& missingVals,
		std::vector<int>& nLevels, vector<int>& cumulativeLevels);

	std::vector<std::vector<double> > miScores;
	std::vector<double> noParentScores;
	std::vector<int> varParams;
	SolutionCalculator * calculator;
	int maxParents, maxChildren;
	std::set<int> (GABayesSolutionCreator::*limitPtr)(int,vector<int>&,int);
	std::set<int> (GABayesSolutionCreator::*childLimitPtr)(int,vector<int>&,int);


};

#endif
