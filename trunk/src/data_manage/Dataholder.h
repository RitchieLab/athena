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
#ifndef DATAHOLDER_H_
#define DATAHOLDER_H_

#include "Individual.h"
#include "DataExcept.h"
#include <string>

using namespace std;

namespace data_manage
{

///
/// Holds data for analysis.  Actual analysis is run using the
/// Dataset class.  The Dataset class contains pointers back 
///
class Dataholder
{
public:
	Dataholder();
	~Dataholder();

	/// Adds an individual to the dataholder
	void addInd(Individual& ind);
	
	/// Adds an individual to the dataholder
	void addInd(Individual* ind);

	/// Returns pointer to indiated individual
	inline Individual* operator[] (unsigned int indIndex){return inds[indIndex];}

	/// Returns number of individuals in set
	inline unsigned int numInds(){return inds.size();}

	/// Add variable name
	inline void addCovarName(string varName, string scaleGroup="default"){
		covars.push_back(varName);
		covarsMap[varName]=covars.size()-1;
		if(covarsScaleGroup.find(scaleGroup) == covarsScaleGroup.end()){
			vector<int> newVector;
			covarsScaleGroup[scaleGroup]=newVector;
		}
		covarsScaleGroup[scaleGroup].push_back(covars.size()-1);
		covarsGroup[covars.size()-1]=covarsScaleGroup.find(scaleGroup);
	}
	
	/// returns group name based on covariate index
	inline string getCovarGroupName(unsigned int index){
		std::map<unsigned int, std::map<std::string, std::vector<int> >::iterator>::iterator iter = covarsGroup.find(index);
		return iter->second->first;
	}
	
	inline std::map<std::string, std::vector<int> >& getContinGroups(){return covarsScaleGroup;}

	/// Retrieve variable name
	inline string getCovarName(unsigned int index){return covars[index];}

	/// Retrieve index by name
	inline unsigned int getCovarIndex(string covarName){return covarsMap[covarName];}

	/// Number of covariates in set
	inline unsigned int numCovars(){return covars.size();}

	/// Add genotype name
	inline void addGenoName(string varName){genos.push_back(varName);genosMap[varName]=genos.size()-1;}

	/// Retrieve index for name
	inline unsigned int getGenoIndex(string genoName){return genosMap[genoName];}

	/// Retrieve variable name
	inline string getGenoName(unsigned int index){return genos[index];}

	/// Number of genotypes in set
	inline unsigned int numGenotypes(){return genos.size();}
	
	/// Sets maximum locus value in set
	inline void setMaxLocusValue(unsigned int maxLocusValue){maxLocus=maxLocusValue;}

	/// Returns maximum  locus value in set
	inline unsigned int getMaxLocusValue(){return maxLocus;}

	/// Sets indicator when any genotpyes are missing in set
	inline void anyMissingGenos(bool missing){anyMissing = missing;}

	/// Gets indicator true when mmissing genotypes in data
	inline bool anyMissingGenos(){return anyMissing;}

	/// Returns genotype for indicated individual at indicated locus
	inline float getGenotype(unsigned int currInd, unsigned int currLoc){return inds[currInd]->value(currLoc);}

	/// Returns number of genotypes in set
	inline unsigned int numGenos(){return inds[0]->numGenotypes();}

	/// Returns number of covariates in set
	inline unsigned int numCovariates(){return inds[0]->numCovariates();}

	/// Returns pointer to ind
	inline Individual* getInd(unsigned int index){return inds[index];}
	
	/// Returns whether data matches case/control set
	inline bool isCaseControl(){return binaryStatusOnly;}
	
	/// Sets whether data is case/control only
	inline void setCasecontrol(bool tf){binaryStatusOnly=tf;}

	/// Returns individual based on id
	Individual* getIndByID(string id){
		if(indsMap.find(id)==indsMap.end()) throw DataExcept("Unable to find individual with ID=" + id);
		return indsMap[id];
	}
	
	/// Gets value for continous variable missing
	float getMissingCoValue(){return missingCoValue;}
	
	/// Sets missing value for continuous variables
	void setMissingCoValue(float miss){missingCoValue = miss;}
	
	/// Gets missing value for genotypes
	int getMissingGenotype(){return missingGenotype;}
	
	/// Sets missing value for genotypes
	void setMissingGenotype(int miss){missingGenotype = miss;}
	
	/// Adds default snp names
	void addDefaultSnps();
	
	/// Adds default covariate names
	void addDefaultCovars();
	
	/// Indicates whether ott-dummy encoding used for genotypes
	void ottDummyEncoding(bool val){ottEncoded = val;}
	
	/// Returns status on ott_dummy_encoding
	bool ottDummyEncoding(){return ottEncoded;}
	
	/// Set the cut-off for splitting the training and testing when the 2 files are split
	void setTestSplit(int indnum){splitNum = indnum;}
	
	/// Returns the split number
	int getTestSplit(){return splitNum;}
	
	/// Check for variance in variables and create lists of ones that have zero variance
	void checkVariance();
	
	/// Returns list of excluded Genotypes
	vector<unsigned int> getExcludedGenotypes(){return excludedGenos;}
	
	/// Returns list of excluded Continuous variables
	vector<unsigned int> getExcludedContins(){return excludedContin;}
	
private:
	std::vector<Individual*> inds;
	std::vector<std::string> genos;
	std::vector<std::string> covars;
	std::vector<unsigned int> excludedGenos, excludedContin;
	std::map<std::string, unsigned int> genosMap, covarsMap; //key is geno name, value is index into genos array
	std::map<std::string, Individual*> indsMap;
	std::map<std::string, std::vector<int> > covarsScaleGroup;
	std::map<unsigned int, std::map<std::string, std::vector<int> >::iterator> covarsGroup;
	unsigned int maxLocus;
	bool anyMissing, ottEncoded, binaryStatusOnly;
	int missingGenotype, splitNum;
	float missingCoValue;

};

}

#endif /*DATAHOLDER_H_*/
