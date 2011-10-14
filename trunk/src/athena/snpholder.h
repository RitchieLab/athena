//
// C++ Implementation: snpholder.h
//
// Description: Represents an entity which contains one or more SNPs
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef _SNP_HOLDER_H
#define _SNP_HOLDER_H
#include "snpsnpmodel.h"
#include <map>
#include <vector>
namespace Biofilter {

/**
 * @brief Represents an object that aggregates 1 or more SNPs (by RS number)
 */
class SnpHolder {
public:
    SnpHolder(uint id);
    virtual ~SnpHolder();

	uint ID();

	virtual std::set<uint> SNPs();			///< Return a copy of the snp list
	/**
	 * @brief return list of rs IDs (lookup is pos->rsID)
	 */
	virtual std::set<uint> SNPs(std::map<uint, uint>& lookup);
	
	void AddSnp(uint id);					///< Add a snp to the holder
	
	/**
	 * Brief Generate models
	 */
	void GenerateModels(SnpModelCollection& models, SnpHolder& other, float implIndex);
	//void GenerateModels(SnpModelCollection& models, SnpHolder& other, std::map<uint, uint>& lookup, float implicationIndex);

	SnpSnpModel *RandomModel(SnpHolder& other, float implicationIndex);

	/**
	 * @brief Counts the number of RS Ids found int he set, snpList
	 */
	uint SnpCount(std::set<uint>& snpList);
	uint SnpCount();						///< Return snp count

	/**
	 * @brief Estimate snp/snp models when combined with another snpholder
	 */
	uint EstimateModelCount(SnpHolder* other);


	/**
	 * @Brief Returns the matching SNPs found in the snps and the local snp list
	 */
	virtual std::set<uint> GetSnpCoverage(std::set<uint>& snps);
protected:
	uint id;								///< Unique ID for item
    std::set<uint> snps;					///< SNPs
};

inline
SnpHolder::SnpHolder(uint id) : id(id) { }

inline
SnpHolder::~SnpHolder() { }


inline
std::set<uint> SnpHolder::SNPs() {
 	return snps;
}
inline
std::set<uint> SnpHolder::SNPs(std::map<uint, uint>& lookup) {
	std::set<uint> rsIDs;

	std::set<uint>::iterator itr = snps.begin();
	std::set<uint>::iterator end = snps.end();
	
	while (itr != end) {
		uint pos = *itr;
		if (lookup.find(pos) != lookup.end())
			rsIDs.insert(lookup[pos]);
	}
	return snps;
}

inline
void SnpHolder::AddSnp(uint id) {
	snps.insert(id);
}

inline
uint SnpHolder::ID() {
	return id;
}

inline
uint SnpHolder::SnpCount() {
	return snps.size();
}

inline
std::set<uint> SnpHolder::GetSnpCoverage(std::set<uint>& snpList) {
	std::set<uint> common;
	set_intersection(snps.begin(), snps.end(), snpList.begin(), snpList.end(), inserter(common, common.begin()));
	return common;
}

inline
uint SnpHolder::SnpCount(std::set<uint>& snpList) {
	return GetSnpCoverage(snpList).size();
}

inline
uint SnpHolder::EstimateModelCount(SnpHolder* other) {
	uint common = SnpCount(other->snps);
	return (snps.size() - common) * (other->snps.size() - common);
}

inline
SnpSnpModel *SnpHolder::RandomModel(SnpHolder& other, float implicationIndex) {
	std::set<uint> s = SNPs();
	std::vector<uint> localSnps;
	localSnps.insert(localSnps.end(), s.begin(), s.end());
	s = other.SNPs();
	std::vector<uint> otherSnps;
	otherSnps.insert(otherSnps.end(), s.begin(), s.end());
	uint lIdx = rand() % localSnps.size();
	uint rIdx = rand() % otherSnps.size();
	return new SnpSnpModel(localSnps[lIdx], otherSnps[rIdx], implicationIndex);
}

inline
void SnpHolder::GenerateModels(SnpModelCollection& models, SnpHolder& other, float implicationIndex) {
	std::set<uint> localSnps = SNPs();
	std::set<uint> otherSnps = other.SNPs();

	std::set<uint> left;
	set_difference(localSnps.begin(), localSnps.end(), otherSnps.begin(), otherSnps.end(), inserter(left, left.begin()));

	std::set<uint> right;
	set_difference(otherSnps.begin(), otherSnps.end(), localSnps.begin(), localSnps.end(), inserter(right, right.begin()));
	
	std::set<uint>::iterator localItr = left.begin();
	std::set<uint>::iterator localEnd = left.end();

	std::set<uint>::iterator otherItr = right.begin();
	std::set<uint>::iterator otherEnd = right.end();

	SnpSnpModel *model;
	while (localItr != localEnd) {
		uint snp1 = *localItr++;
		otherItr = right.begin();

		while (otherItr != otherEnd) {
			uint snp2 = *otherItr++;
			model = new SnpSnpModel(snp1, snp2, implicationIndex);
			if (models.find(model) == models.end())
				models.insert(model);
			else
				delete model;
		}
	}

}
}

#endif //_SNP_HOLDER_H

