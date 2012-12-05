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
//
// C++ Implementation: snpsnpmodel.h
//
// Description: Defines functionality associated with Snp/Snp models
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __SNP_SNP_MODEL_H
#define __SNP_SNP_MODEL_H
#include "biomodel.h"
#include <set>
#include <iomanip>
#include <assert.h>
#ifndef uint 
typedef unsigned int uint;
#endif
namespace Biofilter {



class SnpSnpModel : public BioModel {
public:
	SnpSnpModel(uint snp1, uint snp2, float implicationIndex);
	virtual ~SnpSnpModel();

	/* file buffer requirements */
	bool operator<(const SnpSnpModel& other) const;
	bool EvalLT(const SnpSnpModel* other) const;
	bool operator==(const SnpSnpModel& other) const;
	SnpSnpModel& operator=(const SnpSnpModel& other);
	void WriteBinary(std::ostream& file) const;
	void Write(std::ostream& os, bool useBinary = true) const;
	bool LoadBinary(std::istream& file);

	std::set<uint> snps;
};

class SnpSnpModel;
struct LtSnpModelPointer {
	bool operator()(const SnpSnpModel *model1, const SnpSnpModel *model2) const {
		return *model1 < *model2;
	}
};
typedef std::set<SnpSnpModel*, LtSnpModelPointer> SnpModelCollection;

inline
SnpSnpModel::SnpSnpModel(uint snp1, uint snp2, float implicationIndex) : BioModel(implicationIndex) {
	snps.insert(snp1);
	snps.insert(snp2);
}

inline
SnpSnpModel::~SnpSnpModel() { }

inline
bool SnpSnpModel::operator<(const SnpSnpModel& other) const {
	return snps < other.snps;
}

inline
bool SnpSnpModel::EvalLT(const SnpSnpModel* other) const {
	return snps < other->snps;
}

inline
bool SnpSnpModel::operator==(const SnpSnpModel& other) const {
	return snps == other.snps;
}

inline
SnpSnpModel& SnpSnpModel::operator=(const SnpSnpModel& other) {
	implicationIndex = other.implicationIndex;
	snps = other.snps;
	return *this;
}


inline
void SnpSnpModel::Write(std::ostream& os, bool useBinary) const {
	if (useBinary) {
		WriteBinary(os);
		return;
	}
	std::set<uint>::iterator itr = snps.begin();
	std::set<uint>::iterator end = snps.end();
	//uint count = 0;
	while (itr != end) {
		os<<std::setw(11)<<std::right<<*itr++;
	}
	os<<std::setw(7)<<std::right<<std::setprecision(1)<<ImplicationIndex()<<"\n";
}

inline
void SnpSnpModel::WriteBinary(std::ostream& file) const {
	//uint modelSize = snps.size();
	file.write((char*)&implicationIndex, sizeof(float));           //# implication index

	std::set<uint>::iterator itr = snps.begin();
	std::set<uint>::iterator end = snps.end();
	
	while (itr != end) {
		uint locus = *itr++;
		file.write((char*)&locus, sizeof(uint));
	}
}

inline
bool SnpSnpModel::LoadBinary(std::istream& file) {
	uint modelSize = 2;
	file.read((char*)&implicationIndex, sizeof(float));

	for (size_t i=0; i<modelSize; i++) {
		uint locus = 0;
		file.read((char*)&locus, sizeof(uint));
		snps.insert(locus);
	}
	return modelSize > 0;
}
}

#endif //__SNP_SNP_MODEL_H
