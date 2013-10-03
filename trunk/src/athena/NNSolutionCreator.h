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
 * File:   NNSolutionCreator.h
 * Author: dudeksm
 *
 * Created on December 1, 2008, 4:18 PM
 */

#ifndef _NNSOLUTIONCREATOR_H
#define	_NNSOLUTIONCREATOR_H

#include "SolutionCreator.h"
#include "TerminalSymbCreator.h"
#include "SolutionCalculator.h"
#include "NNSolution.h"
#include "NNLog.h"
#include <set>

class NNSolutionCreator: public SolutionCreator{
		
public:
	 
		/// Constructor
		NNSolutionCreator();
		
		/// Alternative constructor
		NNSolutionCreator(vector<string>& symbols);
		
		void initialize();
		
		/// creates solution from vector of strings
		virtual void establishSolution(vector<string>& symbols, Dataset* set);
 
		 /// creates solution from vector of strings
		virtual void establishSolution(vector<string>& symbols);
 
		/// returns fitness score through evaluation of solution
		virtual float evaluate(Dataset* set);
		
		/// optimize solution by running back propagation
		int optimizeSolution(std::vector<std::string>& symbols, Dataset* set);
		
		/// returns optimized score
		float getOptimizedScore(){return optimizedScore;}
			 
		 virtual Solution* createNewSolution(){
				 NNSolution* sol = new NNSolution;
				return sol;}
		
		/// frees memory associated with constants
		void freeSolution(){
			for(unsigned int i=0; i<constants.size(); i++){
				delete constants[i];
			}
			constants.clear();
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

		inline int getNumIndsEvaluated(){return nIndsEvaluated;}

		/// Returns symbol that corresponds to start of optimization symbols  
		string getStartOptSymbol(){return startOpt;}
		
		std::set<string> getOptIncluded(){return optSymbols;}

		char getLeftOptBound(){return leftOptBound;}
		char getRightOptBound(){return rightOptBound;}
		std::set<string> getOptArgSymbols(){return optArgSymbols;}

		void detailedLogging();
		unsigned int getDetailedLog();

protected:
		
		void compressOperator(vector<TerminalSymbol*> & postFixStack,
			vector<TerminalSymbol*>& newStack);
		
		/// evaluates single individual and returns value for that individual
		virtual float evaluateInd(Individual* ind);
		
		/// resets number of genotypes and covariates when necessary
		void setVariables(Dataset* set);
		
		/// returns true when ind has complete data for solution
		bool useInd(Individual* ind, Dataset* set);
		
		TerminalSymbCreator termHolder;
		
		float optimizedScore;
		bool terminalsSet;
		unsigned int nnTerminalSize, nnDepth;
		vector<TerminalSymbol *> postFixStack;
		int nIndsEvaluated;
		string startOpt;
		std::set<string> optSymbols, optArgSymbols;
		
		char leftOptBound, rightOptBound;
		
		// vectors for holding list of covariates and genotypes in current solution
		// these can be checked against an ind and the individual can be skipped 
		// when needed
		map<TerminalSymbol*, int> covars, genos;
		
		// vector holds pointers to constants that can be destroyed after all evaluations
		vector<TerminalSymbol* > constants;
};


#endif	/* _NNSOLUTIONCREATOR_H */
