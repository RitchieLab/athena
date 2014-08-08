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
/* 
 * File:   BayesSolutionCreator.h
 * Author: dudeksm
 *
 * Created on March 4, 2014, 5:18 PM
 */

#ifndef _BAYESSOLUTIONCREATOR_H
#define	_BAYESSOLUTIONCREATOR_H

#include "SolutionCreator.h"
#include "SolutionCalculator.h"
#include "BayesSolution.h"
#include <set>
#include "DAGraph.h"

class BayesSolutionCreator: public SolutionCreator{
		
public:
	 
		/// Constructor
		BayesSolutionCreator();
		
		/// Alternative constructor
		BayesSolutionCreator(vector<string>& symbols);
		
		void initialize();

		/// creates solution from vector of strings
		virtual void establishSolution(vector<string>& symbols, Dataset* set);
 
		 /// creates solution from vector of strings
		virtual void establishSolution(vector<string>& symbols);
		
		virtual void establishSolutionEquation(std::vector<std::string>& symbols);
 
		void establishSolutionOrig(vector<string>& symbols);
 
		/// returns fitness score through evaluation of solution
		virtual float evaluate(Dataset* set);
		
		virtual vector<std::string> getAdditionalOutputNames(){
			vector<std::string> calcNames = calculator->getAdditionalOutputNames();
			return calcNames;
		}
		
		virtual vector<std::string> getAdditionalFinalOutput(){
			vector<std::string> returnValues;
			if(notImproved){
				returnValues.push_back("*");
			}
			else{
				returnValues.push_back(" ");
			}
			return returnValues;
		}
		
		/// optimize solution by running back propagation
		int optimizeSolution(std::vector<std::string>& symbols, Dataset* set);
		
		/// returns optimized score
		float getOptimizedScore(){return optimizedScore;}
			 
		 virtual Solution* createNewSolution(){
				BayesSolution* sol = new BayesSolution;
				return sol;}
		
		/// frees memory associated with constants
		void freeSolution(){
// 			network.clearNodes();
		}
		
		virtual void restrict(vector<string>& vars){}
		
		/// outputs each individual with score of the network 
		virtual float evaluateWithOutput(Dataset* set, ostream& os);
		
		/// evaluate for further output information (such as AUC)
		virtual void evaluateForOutput(Dataset* set);
		
		/// writes a dot compatible text file representing the network
		virtual void graphicalOutput(ostream& os, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed);
 
 		/// writes solution out as an equation
 		virtual void equationOutput(ostream& os, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed);
 
		vector<int> getGeneIndexes();
		vector<int> getCovarIndexes();
		
		/// Returns worst score
		inline float getWorst(){return calculator->getWorst();}
 
		inline unsigned int getNumGenes(){return getGeneIndexes().size();}
		
		inline unsigned int getNumCovars(){return getCovarIndexes().size();}
		
		inline unsigned int getNumNodes(){return network.numNodes();}
		
		inline unsigned int getComplexity(){return getNumNodes();};

		inline int getNumIndsEvaluated(){return nIndsEvaluated;}

		/// Returns symbol that corresponds to start of optimization symbols  
		string getStartOptSymbol(){return startOpt;}
		
		std::set<string> getOptIncluded(){return optSymbols;}

		char getLeftOptBound(){return leftOptBound;}
		char getRightOptBound(){return rightOptBound;}
		std::set<string> getOptArgSymbols(){return optArgSymbols;}

		void detailedLogging();
		unsigned int getDetailedLog();
		
		virtual void setMapper(AthenaGrammarSI* m);

protected:
		
		/// evaluates single individual and returns value for that individual
		virtual float evaluateInd(Individual* ind);
		
		/// resets number of genotypes and covariates when necessary
		void setVariables(Dataset* set);
		
		/// returns true when ind has complete data for solution
		bool useInd(Individual* ind, Dataset* set);
		
		void setParentScores(Dataset* set);
		
		int configParentData(std::vector<int>& parentValues, 
			std::vector<IndividualTerm*> &parents);
		
		double k2calcNoParent(unsigned int gIndex, double& nP);
		double k2calcPhenoNoParent(double& nP);
		double k2calcWithParent(IndividualTerm* node, std::vector<IndividualTerm*> &parents,
			double& nP);
		
		std::string getLabel(GraphNodeIter& node, Dataholder* holder);
		
		TerminalSymbCreator termHolder;
		
		DAGraph network;
		
		float optimizedScore;
		bool terminalsSet, parentScoresSet, notImproved;
		unsigned int nnTerminalSize, nnDepth, numNodes;
		std::map<TerminalSymbol*, double> parentScore, parentParams;
		Dataset* currentSet;
// 		vector<TerminalSymbol *> postFixStack;
		int nIndsEvaluated;
		string startOpt;
		std::set<string> optSymbols, optArgSymbols;
		
		char leftOptBound, rightOptBound;
		
		// vectors for holding list of covariates and genotypes in current solution
		// these can be checked against an ind and the individual can be skipped 
		// when needed
		map<TerminalSymbol*, int> covars, genos;
		
// vector holds pointers to constants that can be destroyed after all evaluations
// 		vector<TerminalSymbol* > constants;
};


#endif	/* _BAYESSOLUTIONCREATOR_H */