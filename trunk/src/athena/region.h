//
// C++ Implementation: region.h
//
// Description: Represents a basic region (gene) object
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __REGION_H
#define __REGION_H

#include <fstream>
#include "snpholder.h"
#include <iostream>
namespace Biofilter {
/**
* @brief encapsulates a bit more information specific to a region(gene)
*/
class Region : public SnpHolder {
public:
	Region(const char *name, uint id = 0);
	virtual ~Region();


	bool IsPresent(uint snp);			///< Returns true if the SNP is a member of the gene

	/** 
	 * @brief Load from file 
	 * @param source stream being read
	 * @param useBinary True for binary format
	 */
	bool Load(std::istream& source, bool useBinary = false);
	/** 
	 * @brief write to file 
	 * @param source stream being written to
	 * @param useBinary True for binary format
	 */
	void Write(std::ostream& source, bool useBinary = false);

	void InsertDI(uint);				///< Add a disease independent group
	void InsertDD(uint);				///< Add a disease dependent group

	std::set<uint> DiGroups();			///< Return copy of the di groups
	std::set<uint> DdGroups();			///< return a copy of the dd groups

	virtual std::string RegionName();	///< Returns the region name
	void RegionName(const char *name);	///< Sets the region name
protected:
	void WriteBinary(std::set<uint>& items, std::ostream& s);
	void Write(std::set<uint>& items, std::ostream& s);
	void LoadBinary(std::set<uint>& items, std::istream& s);
	std::string name;					///< Common Name for the gene (for reporting purposes)
	std::set<uint> diGroups;			///< Identify which metagroups we are members of
	std::set<uint> ddGroups;			///< Identify which disease dependent groups we are members of
	/** Loads into items, values from ascii stream*/
	void Load(std::set<uint>& items, std::istream& s);

};

inline
Region::Region(const char *name, uint id) : SnpHolder(id), name(name) { }

inline
Region::~Region() { }

inline
void Region::InsertDI(uint group) {
	diGroups.insert(group);
}
inline
void Region::InsertDD(uint group) {
	ddGroups.insert(group);
}

inline
std::set<uint> Region::DiGroups() {
	return diGroups;
}

inline
std::set<uint> Region::DdGroups() {
	return ddGroups;
}

inline
std::string Region::RegionName() {
	return name;
}

inline
void Region::RegionName(const char *name) {
	this->name = name;
}
inline
bool Region::IsPresent(uint snp) {
	return snps.find(snp) != snps.end();
}

inline
void Region::WriteBinary(std::set<uint>& items, std::ostream& s) {
	uint count = items.size();
	uint item;
	s.write((char*)&count, sizeof(uint));
	std::set<uint>::iterator itr = items.begin();
	std::set<uint>::iterator end = items.end();
	while (itr != end) {
		item = *itr++;
		s.write((char*)&item, sizeof(uint));
	}
}
inline
void Region::Write(std::set<uint>& items, std::ostream& s) {
	uint count = items.size();
	uint item;

	std::set<uint>::iterator itr = items.begin();
	std::set<uint>::iterator end = items.end();
	count = 0;
	while (itr != end) {
		item = *itr++;
		if (count++ > 0)
			s<<"|";
		s<<item;
	}
}

inline
void Region::Load(std::set<uint>& items, std::istream& s) {
	std::string list, item;
	s>>list;
	std::stringstream listItems(list);
	while (getline(listItems, item, '|'))
		items.insert(atoi(item.c_str()));
	
}

inline
void Region::LoadBinary(std::set<uint>& items, std::istream& s) {
	uint count = 0, item = 0;
	s.read((char*)&count, sizeof(uint));
	for (uint i=0; i<count; i++) {
		item = 0;
		s.read((char*)item, sizeof(uint));
		items.insert(item);
	}
}

inline
bool Region::Load(std::istream& source, bool useBinary) {
	bool isSuccessful = source.good() && !source.eof();

	if (isSuccessful) {
		if (useBinary) {
			source>>name;
			source.read((char*)&id, sizeof(uint));
			LoadBinary(snps, source);
			LoadBinary(diGroups, source);
			LoadBinary(ddGroups, source);
		}
		else {
			char line[1048576];
			source.getline(line, 1048576);
			std::stringstream ss(line);
			ss>>name;
			ss>>id;
			Load(snps, ss);
			Load(diGroups, ss);
			Load(ddGroups, ss);
		}
	}

	return isSuccessful;
}

inline
void Region::Write(std::ostream& source, bool useBinary) {
	std::set<uint> snps = SNPs();
	if (snps.size() > 0) {
		source<<RegionName();
		if (useBinary) {
			source.write((char*)&id, sizeof(uint));
			WriteBinary(snps, source);
			WriteBinary(diGroups, source);
			WriteBinary(ddGroups, source);
			assert(ddGroups.size() == 0);
		}
		else {
			source<<"\t"<<id<<"\t";
			Write(snps,source);
			source<<"\t";
			Write(diGroups, source);
			source<<"\t";
			Write(ddGroups, source);
			source<<"\n";
		}
	}
}

}

#endif //__REGION_H
