/* 
 * File:   GeneGeneModel.h
 * Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
 * Copyright: See COPYING file that comes with this distribution
 * Created on October 19, 2009, 1:45 PM
 */
#include <string>
#include <sstream>
#include "region.h"
#include <map>
#include <iomanip>

#ifndef BIO_GENEGENEMODEL_H
#define	BIO_GENEGENEMODEL_H

namespace Biofilter {

typedef std::map<uint, Region*> RegionLookup;
/**
 * @brief Represents the gene/gene models 
 */
class GeneGeneModel : public BioModel {
public:
	GeneGeneModel();
    GeneGeneModel(Region* g1, Region* g2, uint pairID = (uint)-1);
    GeneGeneModel(const GeneGeneModel& orig);
    virtual ~GeneGeneModel();

	/* file buffer requirements */
	bool operator<(const GeneGeneModel& other) const;
	bool operator==(const GeneGeneModel& other) const;
	GeneGeneModel& operator=(const GeneGeneModel& other);
	void WriteBinary(std::ostream& file) const;
	void Write(std::ostream& os, bool isBinary) const;
	void MergeGroups(const GeneGeneModel& other);


	
	/**
	 * @brief Intiate a load
	 * @param file source of the archive
	 * @param isBinary True for binary formatted files
	 */
	bool Load(std::istream& file, bool isBinary);
	bool LoadBinary(std::istream& file);		///< Loads from binary archive

	std::set<uint> DdGroups() const;

	/**
	 * @brief Returns implication index
	 */
	float ImplicationIndex() const;

	/**
	 * @brief Generate snp-snp models associated with local genes
	 * @note This will insert models into the set if a match isn't already present
	 */
	void GenerateModels(SnpModelCollection& models);

	SnpSnpModel *RandomModel();
	
	/**
	 * @brief Purges contents of a SnpModelCollection (deallocating)
	 */
	void Purge(SnpModelCollection& models);

	uint EstimateModelCount();					///< Returns number of models estimated to be generated

	static RegionLookup geneLookup;				///< Gene lookups
	static float DuplicateDD_Weight;			///< Weight applied to duplicate disease dependent groups
	void WriteSummary(std::ostream& os);		///< Summarize a gene/gene model
	std::string ToString(std::set<uint> groups) const;
	void Reset();								///< Clears out the model
protected:
	Region *region1;
	Region *region2;
    std::set<uint> diGroups;					///< pairing information
};


struct ggmScoreGT {
	bool operator()(const GeneGeneModel*model1, const GeneGeneModel* model2) const {
		float ii1 = model1->ImplicationIndex();
		float ii2 = model2->ImplicationIndex();

		if (ii1 == ii2)
			return *model2 < *model1;
		else
			return ii1 > ii2;
	}
};


inline
void GeneGeneModel::Reset() {
	region1 = region2 = NULL;
	diGroups.clear();
}

inline
GeneGeneModel::GeneGeneModel() : BioModel(0.0), region1(NULL), region2(NULL) { }

inline
GeneGeneModel::GeneGeneModel(Region* g1, Region* g2, uint pairID) : BioModel(0.0) {
	uint id1 = g1->ID();
	uint id2 = g2->ID();

	if (pairID != (uint)-1)
		diGroups.insert(pairID);

	if (id1 > id2) {
		region1 = g2;
		region2 = g1;
	} else {
		region1 = g1;
		region2 = g2;
	}
}

inline
GeneGeneModel::GeneGeneModel(const GeneGeneModel& orig) {
	implicationIndex = orig.implicationIndex;
	region1 = orig.region1;
	region2 = orig.region2;
	diGroups = orig.diGroups;
}

inline
GeneGeneModel::~GeneGeneModel() {	}



/* file buffer requirements */
inline
bool GeneGeneModel::operator<(const GeneGeneModel& other) const {
	if (region1 == other.region1)
		return region2 < other.region2;
	else
		return region1 < other.region1;
}

inline
bool GeneGeneModel::operator==(const GeneGeneModel& other) const {
	return region1 == other.region1 && region2 == other.region2;
}


inline
GeneGeneModel& GeneGeneModel::operator=(const GeneGeneModel& other) {
	implicationIndex = other.implicationIndex;
	region1 = other.region1;
	region2 = other.region2;
	diGroups = other.diGroups;

	return *this;
}

inline
uint GeneGeneModel::EstimateModelCount() {
	return region1->EstimateModelCount(region2);
}

inline
void GeneGeneModel::WriteSummary(std::ostream& os) {
	std::set<uint> dd = region1->DdGroups();	//diGroups;

	std::set<uint> otherDD = region2->DdGroups();
	dd.insert(otherDD.begin(), otherDD.end());
	os<<std::setw(35)<<region1->RegionName()<<std::setw(8)<<region1->SnpCount()
		<<std::setw(35)<<region2->RegionName()<<std::setw(8)<<region2->SnpCount()
		<<std::setw(10)<<ImplicationIndex()<<std::setw(10)<<EstimateModelCount()<<"\t"
		<<ToString(diGroups)<<","<<ToString(dd)<<"\n";

}

inline
std::string GeneGeneModel::ToString(std::set<uint> groups) const {
	std::stringstream ss;
	std::set<uint>::iterator itr = groups.begin();
	std::set<uint>::iterator end = groups.end();

	uint count = 0;
	while (itr != end) {
		if (count++ > 0)
			ss<<"|";
		ss<<*itr++;
	}
	return ss.str();
}
inline
void GeneGeneModel::Write(std::ostream& os, bool isBinary=false)const {
	if (isBinary)
		WriteBinary(os);
	else {
		uint id1=region1->ID(), id2=region2->ID();
		os<<id1<<"\t"<<id2<<"\t"<<ImplicationIndex()<<"\t"<<ToString(diGroups)<<"\n";
	}
}
inline
void GeneGeneModel::WriteBinary(std::ostream& file) const {
	uint id1=region1->ID(), id2=region2->ID();
	uint count = diGroups.size();
	float score = ImplicationIndex();
	file.write((char*)&id1, sizeof(uint));
	file.write((char*)&id2, sizeof(uint));
	file.write((char*)&score, sizeof(float));
	file.write((char*)&count, sizeof(uint));
	std::set<uint>::iterator itr = diGroups.begin();
	std::set<uint>::iterator end = diGroups.end();

	while (itr != end) {
		uint group = *itr++;
		file.write((char*)&group, sizeof(uint));

	}

}

inline
bool GeneGeneModel::Load(std::istream& file, bool isBinary) {
	Reset();
	if (isBinary)
		return LoadBinary(file);
	else {
		uint id1=0, id2=0;
		float score;
		file>>id1>>id2>>score;
		region1 = region2 = NULL;
		if (id1 > 0 && id2 > 0) {
			region1 = geneLookup[id1];
			region2 = geneLookup[id2];

			std::string groups, group;
			file>>groups;
			std::stringstream groupList(groups);
			while (getline(groupList, group,  '|'))
				diGroups.insert(atoi(group.c_str()));
		}
	}
	return (region1 && region2);
}


inline
bool GeneGeneModel::LoadBinary(std::istream& file) {
	uint id1=0, id2=0, count;
	float score = 0.0;
	bool isGood = file.good() && !file.eof();

	if (isGood) {
		file.read((char*)&id1, sizeof(uint));
		file.read((char*)&id2, sizeof(uint));

		if (id1 > 0 && id2 > 0) {
			region1 = geneLookup[id1];
			region2 = geneLookup[id2];

			file.read((char*)&score, sizeof(float));
			file.read((char*)&count, sizeof(uint));

			diGroups.clear();
			for (uint i=0; i<count; i++) {
				uint group = 0;
				file.read((char*)&group, sizeof(uint));
				diGroups.insert(group);
			}

		}
	}
	return isGood && id1 > 0 && id2 > 0;
}

inline
void GeneGeneModel::MergeGroups(const GeneGeneModel& other) {
	diGroups.insert(other.diGroups.begin(), other.diGroups.end());
}

inline
std::set<uint> GeneGeneModel::DdGroups() const {
	std::set<uint> ddSingular;
	std::set<uint> dd = region1->DdGroups();
	ddSingular.insert(dd.begin(), dd.end());
	dd = region2->DdGroups();
	ddSingular.insert(dd.begin(), dd.end());
	return ddSingular;
}
inline
float GeneGeneModel::ImplicationIndex() const {
	std::set<uint> ddSingular = DdGroups();
	uint totalDD =region1->DdGroups().size() + region2->DdGroups().size();
	return (float)(diGroups.size() + ddSingular.size()) + ((float)(totalDD - ddSingular.size()) * DuplicateDD_Weight);
}

inline
SnpSnpModel *GeneGeneModel::RandomModel() {
	return region1->RandomModel(*region2, ImplicationIndex());
}

inline
void GeneGeneModel::GenerateModels(SnpModelCollection& models) {
	region1->GenerateModels(models, *region2, ImplicationIndex());
}

inline
void GeneGeneModel::Purge(SnpModelCollection& models) {
	SnpModelCollection::iterator itr = models.begin();
	SnpModelCollection::iterator end = models.end();

	while (itr != end)
		delete *itr++;
	models.clear();
}


}
#endif	/* BIO_GENEGENEMODEL_H */

