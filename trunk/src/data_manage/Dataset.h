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
#ifndef DATASET_H_
#define DATASET_H_

#include "Individual.h"
#include <iostream>

namespace data_manage
{

///
///	Dataset contains pointers to individual objects.  It can represent the entire
/// data set or can be used for a permuted set or a training, testing or crossvalidation
/// interval.
///
class Dataset
{

	friend std::ostream& operator<<(std::ostream& os, Dataset& d);

public:
	Dataset();
	Dataset(float coMissing, unsigned int genoMissing, bool caseControlStatus){
			missingCoValue = coMissing;
			missingGenotype = genoMissing;
			binaryStatusOnly = caseControlStatus;
			ssTotal = 0;
	}
	~Dataset();

	inline unsigned int numInds(){return inds.size();}

	inline void addInd(Individual* ind){inds.push_back(ind);}

	inline Individual* operator[](unsigned int index){return inds[index];}

	inline Individual* getInd(unsigned int index){return inds[index];}
	
	/// Returns number of genotypes in set
	inline unsigned int numGenos(){return inds[0]->numGenotypes();}
	
	/// Returns number of covariates in set
	inline unsigned int numCovariates(){return inds[0]->numCovariates();}

	void addInds(std::vector<Individual* >& new_inds);

	inline float getMissingCoValue(){return missingCoValue;}
	
	inline void setMissingCoValue(float miss){missingCoValue = miss;}
	
	inline int getMissingGenotype(){return missingGenotype;}
	
	inline void setMissingGenotype(int miss){missingGenotype = miss;}
	
	/// returns SStotal for this set (used in r-squared calculations)
	float getSSTotal(){return ssTotal;}
	
	/// Calculates SStotal
	void calcSSTotal();
	
	/// Sets constant value used in converting scores
	void setConstant(double c){constantValue = c;}
	
	/// Returns constant value 
	double getConstant(){return constantValue;}
	
	/// Returns whether set is case/control
	bool isCaseControl(){return binaryStatusOnly;}
	
	/// Sets case/control status
	void setCaseControlStatus(bool tf){binaryStatusOnly=tf;}
	
	Dataset operator+(Dataset& d);
	
	/// Gets number of levels for a continuous variable
	unsigned int getNumLevels(unsigned int varIndex){return continLevels[varIndex];}
	
	/// Sets number of levels for a continuous variable
	void setNumLevels(unsigned int varIndex, unsigned int nLevels){
		if(continLevels.empty())
			continLevels.resize(numCovariates());
		continLevels[varIndex]=nLevels;
	}
	
	/// Returns number of levels for phenotype
	unsigned int getNumStatusLevels(){return statusLevels;}
	
	/// Sets number of levels for phenotype
	void setNumStatusLevels(unsigned int nLevels){statusLevels=nLevels;}	
	
private:
	std::vector<Individual*> inds;
	std::vector<unsigned int> continLevels;
	float missingCoValue, ssTotal;
	int missingGenotype;
	bool binaryStatusOnly;
	double constantValue;
	unsigned int statusLevels;

};

}

#endif /*DATASET_H_*/
