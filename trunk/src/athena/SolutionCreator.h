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
 * File:   SolutionCreator.h
 * Author: dudeksm
 *
 * Created on December 1, 2008, 4:10 PM
 */

#ifndef _SOLUTIONCREATOR_H
#define	_SOLUTIONCREATOR_H

#include "Solution.h"
#include "CalculatorFactory.h"
#include "AlgorithmLog.h"
#include "Structs.h"
#include <Dataset.h>
#include <set>
#include <AUCCalc.h>

#include <iostream>

using namespace data_manage;
///
/// Solution is the abstract base class for solutions in the HEMANN
/// system.  Solutions are the structure created by the algorithms.
/// An example is NeuralNetwork.  Future possible solutions would be
/// support vector machines, etc.
///
class SolutionCreator{
		
public:

		virtual ~SolutionCreator(){}
		
		/// creates solution from a vector of strings
		virtual void establishSolution(std::vector<std::string>& symbols, Dataset* set)=0;

		/// creates solution from a vector of strings
		virtual void establishSolution(std::vector<std::string>& symbols)=0;
		
		/// returns fitness score through evaluation of solution
		virtual float evaluate(Dataset* set)=0;
		
		/// adds a solution calculator for determining fitness
		virtual void setCalculator(std::string calc_type){
				if(calculator != NULL)
					delete calculator;
				calculator = CalculatorFactory::createCalculator(calc_type);
		}
		
		virtual void freeSolution()=0;
		
		/// returns solution
		Solution* getSolution(){return sol;}
		
		/// returns fitness
		float fitness(){return solFitness;}
		
		/// set fitness
		void fitness(float fit){solFitness = fit;}
		
		/// optimize solution
		virtual int optimizeSolution(std::vector<std::string>& symbols, Dataset* set)=0;
		
		/// returns vector of optimized values for use in adapting the original
		vector<float> getOptimizedValues(){return optValues;}
		
		/// returns optimized score
		virtual float getOptimizedScore()=0;
		
		/// returns vector of structs containing optimized values stored as strings
		vector<symbVector> getOptimizedSymbols(){return optValSymbols;}
		
		/// returns a blank solution of appropriate type
		virtual Solution* createNewSolution()=0;
		
		virtual void setCalculatorConstant(float constant){calculator->setConstant(constant);}
		virtual float getCalculatorConstant(){return calculator->getConstant();}
		
		virtual bool maxBest(){return calculator->maxBest();}
		
		virtual float getWorst(){return calculator->getWorst();}
		
		virtual void restrict(vector<string>& restrictions)=0;
		
		virtual float evaluateWithOutput(Dataset* set, ostream& os)=0;
		
		virtual void captureEvaluation(Dataset* set, vector<stat::TestResult>& results)=0;
		
		/// writes a graphical or file that can be converted to a graphic representation of the solution
		virtual void graphicalOutput(ostream& os, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed)=0;
		
		virtual std::string graphicExt(){return graphicExtension;}
		
		virtual unsigned int getNumGenes()=0;
		
		virtual unsigned int getNumCovars()=0;
		
		virtual vector<int> getGeneIndexes()=0;
		virtual vector<int> getCovarIndexes()=0;
		
		virtual string getStartOptSymbol()=0;
		virtual std::set<string> getOptIncluded()=0;

		virtual char getLeftOptBound()=0;
		virtual char getRightOptBound()=0;
		virtual std::set<string> getOptArgSymbols()=0;
		
		virtual int getNumIndsEvaluated()=0;
		
		virtual void detailedLogging()=0;
		virtual unsigned int getDetailedLog()=0;
		
		std::string calculatorName(){return calculator->getName();}
		
protected:
		Solution* sol;
		SolutionCalculator * calculator;
		std::string graphicExtension;
		float solFitness;
		vector<float> optValues;
		vector<symbVector> optValSymbols;
};



#endif	/* _SOLUTIONCREATOR_H */

