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
#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include <vector>
#include <string>
#include <map>

namespace data_manage
{

/// Contains data and status for individuals in dataset
class Individual
{
public:
	Individual();

	~Individual();

	/// status of individual
	inline float status(){return indStatus;}

	/// returns genotype or environmental value
	inline float value(unsigned int index){
		if(index<numCovars) return covariates[index];
		else return float(genotype[index+numCovars]);
	}

	inline float getGenotype(unsigned int index){return genotype[index];}
	
	inline void setGenotype(unsigned int index, float value){genotype[index]=value;}

	inline void setStatus(float stat){indStatus = stat;}

	inline float getStatus(){return indStatus;}

	/// appends genotype to genotype lsit
	inline void addGenotype(float geno){genotype.push_back(geno);numLoci++;}
	
	/// appends environmental variable to list
	void addCovariate(float val){covariates.push_back(val);numCovars++;}

	/// returns covariate value
	float getCovariate(unsigned int index){return covariates[index];}

	/// sets covariate value
	void setCovariate(unsigned int index, float val){covariates[index] = val;}

	/// returns number of genotypes
	inline unsigned int numGenotypes(){return numLoci;}

	/// returns number of environmental variables
	inline unsigned int numCovariates(){return numCovars;}

	/// sets genotype vector to be equal to vector passed
	inline void setAllGenotypes(std::vector<char>& genos){genotype = genos; numLoci = genos.size();}

	/// returns ID value as string
	std::string getID(){return id;}
	
	/// sets ID for individual 
	void setID(std::string identification){id=identification;}
	
private:
	float indStatus;
	unsigned int numLoci, numCovars;
	std::string id;
	std::vector<float> covariates;
	std::vector<char> genotype;

};

}

#endif /*INDIVIDUAL_H_*/

