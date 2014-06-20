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
 * File:   NNSolution.h
 * Author: dudeksm
 *
 * Created on November 13, 2008, 4:38 PM
 */

#ifndef _NNSOLUTION_H
#define	_NNSOLUTION_H

#include "Solution.h"
#include "TerminalSymbCreator.h"
#include "SolutionCalculator.h"

class NNSolution: public Solution{
		
public:
	 
		/// Constructor
		NNSolution(){setName("Neural Network");}
 
		virtual Solution* clone();
		
		/// returns the genotypes present in the solution
		virtual vector<int> getGenotypes(bool dummyEncoded=true);
		
		/// returns covariates present in the solution
		virtual vector<int> getCovariates();
		
		/// outputs a more human-readable version of the network
		virtual void outputClean(std::ostream& os, data_manage::Dataholder& data,
			bool mapUsed,  bool ottDummy, bool continMapUsed);
		
		/// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
		virtual void adjustScoreOut(Dataset* trainSet, Dataset* testSet, std::string calcName);
		
		/// Adjusts output of the scores when needed for training set only
		virtual void adjustScoreOut(Dataset* trainSet, std::string calcName);
		
		/// Adjusts score passed and returns value
		virtual float adjustScoreOut(float score, int nIndsTested, float ssTotal, std::string calcName);
		
		void setGramDepth(int dep){gramDepth = dep;}
		int getGramDepth(){return gramDepth;}
		
		void setNNDepth(int dep){nnDepth = dep;}
		int getNNDepth(){return nnDepth;}
		
private:

		/// adjusts indexes back to original values if dummy-encoded
		int adjustDummyEncoding(int genotype);
		
		float alterScore(float mse, int totalInds, float ssTotal);
		
		int calcInds(Dataset* set, float& ssTotal);
		
		int nnDepth, gramDepth;
		
};


#endif	/* _NNSOLUTION_H */
