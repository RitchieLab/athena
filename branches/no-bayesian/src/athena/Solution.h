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
 * File:   Solution.h
 * Author: dudeksm
 *
 * Created on November 13, 2008, 3:35 PM
 */

#ifndef _SOLUTION_H
#define	_SOLUTION_H

#include <Dataset.h>
#include <string>
#include <vector>
#include <iostream>
#include <Dataholder.h>

using namespace data_manage;
///
/// Solution is the abstract base class for solutions in the HEMANN
/// system.  Solutions are the structure created by the algorithms.
/// An example is NeuralNetwork.  Future possible solutions would be
/// support vector machines, etc.
///
class Solution{

public:

		/// default constructor
		Solution(){solutionName = "";}

		/// Named constructor
		Solution(std::string name){setName(name);}

		/// Destructor
		virtual ~Solution(){}

		/// Clones the current solution
		virtual Solution* clone();

		/// return number of symbols
		unsigned int getNumSymbols(){return symbols.size();}

		/// operator [] for accessing symbols
		inline std::string& operator[](unsigned int index){return symbols[index];}

		/// sets solution namem
		void setName(std::string name){solutionName = name;}

		/// returns solution name (type)
		std::string getName(){return solutionName;}

		/// set symbols
		void setSymbols(std::vector<std::string>& sym){ symbols=sym;}

		/// returns symbols
		std::vector<std::string>& getSymbols(){return symbols;}

		/// returns fitness
		float fitness(){return solFitness;}

		/// set fitness
		void fitness(float fit){solFitness = fit;}

		/// returns test score
		float testVal(){return testScore;}

		/// sets test score
		void testVal(float val){testScore = val;}

		/// adjusts output so that
		virtual void adjustDummyEncoding(){};

		/// outputs as a text file
		virtual void outputSolution(std::ostream& os);

		/// outputs as a graphviz compatible file
		virtual void outputGraph(std::ostream& os);

		/// outputs a more human-readable version of the network
		virtual void outputClean(std::ostream& os, data_manage::Dataholder& data,
			bool mapUsed, bool ottDummy, bool continMapUsed);

		void copy(Solution* other);

		/// returns the genotypes present in the solution
		virtual std::vector<int> getGenotypes(bool dummyEncoded = true){
				std::vector<int> blank;
				return blank;}

		/// returns covariates present in the solution
		virtual std::vector<int> getCovariates(){
				std::vector<int> blank;
				return blank;}

		virtual void setCovariates(std::vector<int> c){}
		virtual void setGenotypes(std::vector<int> g){}

		/// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
		virtual void adjustScoreOut(Dataset* trainSet, Dataset* testSet, std::string calcName){}

		/// Adjusts output of the scores when needed
		virtual void adjustScoreOut(Dataset* trainSet, std::string calcName){}

		/// Adjusts score passed and returns value
		virtual float adjustScoreOut(float score, int nIndsTested, float constant,
			std::string calcName){return score;}

		/// returns complexity
		int getComplexity(){return complexity;}

		/// set complexity
		void setComplexity(int c){complexity=c;}

		/// sets additional output
		inline void setAdditionalOutput(std::vector<std::string> output){
			addOut.insert(addOut.end(), output.begin(), output.end());
		}

		/// returns additional output
		inline std::vector<std::string>& getAdditionalOutput(){return addOut;}


protected:

		std::vector<std::string> symbols, addOut;
		float solFitness, testScore;
		int complexity;

private:

		std::string solutionName;
};


#endif	/* _SOLUTION_H */

