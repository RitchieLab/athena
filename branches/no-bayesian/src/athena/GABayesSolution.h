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
/*
 * File:   GABayesSolution.h
 * Author: dudeksm
 *
 * Created on September 11, 2015, 3:38 PM
 */

#ifndef _GABAYESSOLUTION_H
#define	_GABAYESSOLUTION_H

#include "Solution.h"
#include "SolutionCalculator.h"

class GABayesSolution: public Solution{

public:

		/// Constructor
		GABayesSolution(){setName("Bayesian Network"); testScore=0;}

		virtual Solution* clone();

		/// returns the genotypes present in the solution
		virtual vector<int> getGenotypes(bool dummyEncoded=true){
			return genos;
		}

		/// returns covariates present in the solution
		virtual vector<int> getCovariates(){return contins;}

		/// outputs a more human-readable version of the network
		virtual void outputClean(std::ostream& os, data_manage::Dataholder& data,
			bool mapUsed,  bool ottDummy, bool continMapUsed);

		/// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
		virtual void adjustScoreOut(Dataset* trainSet, Dataset* testSet){}

		/// Adjusts output of the scores when needed for training set only
		virtual void adjustScoreOut(Dataset* trainSet){}

		/// Adjusts score passed and returns value
		virtual float adjustScoreOut(float score, int nIndsTested, float ssTotal){return 0.0;}

		/// Adjusts score passed and returns value
		virtual float adjustScoreOut(float score, int nIndsTested, float constant,
			std::string calcName){return 0.0;}

		/// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
		virtual void adjustScoreOut(Dataset* trainSet, Dataset* testSet, std::string calcName){}

		/// Adjusts output of the scores when needed
		virtual void adjustScoreOut(Dataset* trainSet, std::string calcName){}

		virtual void setCovariates(std::vector<int> c){contins=c;}
		virtual void setGenotypes(std::vector<int> g){genos=g;}

private:
		std::vector<int> contins, genos;
};


#endif	/* GABAYESSOLUTION */
