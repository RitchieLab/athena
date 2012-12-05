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
// C++ Implementation: genegenemodelreader.h
//
// Description: Serializes gene gene models from binary file
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <assert.h>
#include "genegenemodel.h"
#include <fstream>

namespace Biofilter {

class GeneGeneModelReader {
public:
	GeneGeneModelReader(const char *geneSource , const char *modelFile, bool binaryStorage = true);
	~GeneGeneModelReader();

	/**
	 * @Brief Returns the region lookup (geneID->Region*) used by the reader
	 */
	RegionLookup& GetGeneLookup();

	/**
	 * @Brief Archives snp-snp models
	 * @param filename The name of the gene-gene model
	 * @param modelCount The max number of models allowed
	 * @param minImplicationIndex The minimum Implication index to be considered
	 * @param binaryStorage true indicates binary format
	 */
	std::map<float, uint> ArchiveSnpModels(const char *filename, uint modelCount, uint minImplicationIndex, bool binaryStorage = true);

	/**
	 * @brief This is used to provide a way to iterate over gene-gene models
	 * @note In addition to standard functionality, users can ask the iterator to iterate over itself and generate a large set of models
	 */
	struct iterator {
		friend class GeneGeneModelReader;
		iterator(const iterator& other);
		~iterator();
		GeneGeneModel operator*();			///< Acquire a copy of the model using dereference
		GeneGeneModel* operator->();			///< Return a reference to the current model
		iterator& operator++();				///< Prefix (this is the one I've used...hint hint)
		iterator operator++(int c);			///< Postfix....this should work?

		/**
		 * @Brief Returns a collection of SnpSnp models.
		 * @param models This is the collection which will be cleared and new models added to
		 * @param maxModelCount The collection will be at most maxModelCount, with the exception of times when the first gene/gene model produces more than maxModelCount
		 * @param minImplicationIndex models present in the snp collection will have an implication no lower than minImplicationIndex
		 **/
		uint GetModels(SnpModelCollection& models,uint maxModelCount, uint minImplicationIndex);

		/**
		 * @Brief Used for stl like iteration (itr == end)
		 */
		bool operator==(const iterator& other);
		bool operator!=(const iterator& other);
		iterator& operator=(const iterator& other);
	protected:
		std::streampos GetPos() const;

		bool LoadModel();
		/**
		 * @brief only GeneGenemodelReaders can create these
		 */
		iterator(const char *filename, bool binaryStorage, bool isEnd);
		bool binaryStorage;					///< True indicates binary formatting
		bool isEnd;							///< Used to recognize when we are at the end of the file
		std::ifstream source;					///< The actual data file
		GeneGeneModel curModel;					///< Current model
		std::streampos pos;
	};

	iterator begin();						///< STL style iterator acquisition
	iterator end();							///< STL style iterator stuff
protected:
	/**
	 * @brief We need a gene lookup which the user must provide (biofilter produces this when it writes  a gene-gene archive
	 */
	void LoadGeneLookup(const char *filename, bool binaryStorage);

	bool binaryStorage;						///< True is binary format

        std::string modelSource;					///< The file containing models
	//std::ifstream *source;					///< File we are reading from

};

GeneGeneModelReader::iterator& GeneGeneModelReader::iterator::operator=(const GeneGeneModelReader::iterator& other) {
  source.clear();
	binaryStorage	= other.binaryStorage;
	isEnd			= other.isEnd;
	std::streampos  pos = other.GetPos();
	source.seekg(pos);
	curModel		= other.curModel;
	return *this;
}

GeneGeneModelReader::iterator GeneGeneModelReader::begin() {
	return iterator(modelSource.c_str(), binaryStorage, false);
}

GeneGeneModelReader::iterator GeneGeneModelReader::end() {
	return iterator(modelSource.c_str(), binaryStorage, true);
}

inline
bool GeneGeneModelReader::iterator::LoadModel() {
	if (!isEnd) {
		isEnd = !curModel.Load(source, binaryStorage);
		pos = source.tellg();
	}
	return !isEnd;
}
inline
uint GeneGeneModelReader::iterator::GetModels(SnpModelCollection& models, uint maxModelCount, uint minImplIndex) {
	models.clear();
	uint count = 0;
	if (isEnd)
		return 0;
	do {
		curModel.GenerateModels(models);
		++(*this);
		count++;
	} while (!isEnd && models.size() < maxModelCount && curModel.ImplicationIndex() >= minImplIndex);
	return count;
}

inline
std::streampos GeneGeneModelReader::iterator::GetPos() const {
	return pos;
}
inline
GeneGeneModelReader::iterator::iterator(const GeneGeneModelReader::iterator& other): binaryStorage(other.binaryStorage), isEnd(other.isEnd), curModel(other.curModel) {
    source.seekg(other.GetPos());
}

inline
GeneGeneModelReader::iterator::iterator(const char *filename, bool binaryStorage, bool isEnd) : binaryStorage(binaryStorage), isEnd(isEnd) {

    source.open(filename, std::ios::binary|std::ios::in);

    if (!source.good()) {
	std::cerr<<"Well, the file wasn't opened!\n";
    }
    if (binaryStorage) {
	uint d=0;
	source.read((char*)&d, sizeof(uint));
	if (d != 0) {
	    std::cerr<<"Attempting to open text file in binary mode!";
	    exit(1);
	}
	uint modelCount;
	source.read((char*)&modelCount, sizeof(uint));
    }
    else {
	uint modelCount = 0;
	source >> modelCount;
    }

    isEnd = !source.good() || source.eof();

    if (!isEnd)
	LoadModel();
	
}

inline
GeneGeneModelReader::iterator::~iterator() { }

inline
GeneGeneModel GeneGeneModelReader::iterator::operator*() {
	assert(!isEnd);
	return curModel;
}

inline
GeneGeneModel *GeneGeneModelReader::iterator::operator->() {
	return &curModel;
}

inline
GeneGeneModelReader::iterator& GeneGeneModelReader::iterator::operator++() {
	LoadModel();
	return *this;
}

inline
GeneGeneModelReader::iterator GeneGeneModelReader::iterator::operator++(int c) {
	GeneGeneModelReader::iterator itr(*this);
	LoadModel();
	return itr;
}

inline
bool GeneGeneModelReader::iterator::operator!=(const iterator& other) {
    return !(*this == other);
}

inline
bool GeneGeneModelReader::iterator::operator==(const iterator& other) {
	if (isEnd == true || other.isEnd == true)
		return isEnd == other.isEnd;
  return pos == other.pos;
//	return source == other.source;
}

inline
GeneGeneModelReader::GeneGeneModelReader(const char *geneSource ,const char *modelFile, bool binaryStorage)
	: binaryStorage(binaryStorage), modelSource(modelFile) {
	LoadGeneLookup(geneSource, false);
}

inline
GeneGeneModelReader::~GeneGeneModelReader() {
}

inline
RegionLookup& GeneGeneModelReader::GetGeneLookup() {
	return GeneGeneModel::geneLookup;
}

inline
void GeneGeneModelReader::LoadGeneLookup(const char *filename, bool binaryStorage) {
	std::ifstream file(filename);
	if (binaryStorage) {
		file.close();
		file.open(filename, std::ios::binary);
	}
	GeneGeneModel::geneLookup.clear();
	while (file.good() && !file.eof()) {
		Region *region = new Region("", 0);
		if (region->Load(file, binaryStorage)) {
			GeneGeneModel::geneLookup[region->ID()]=region;
		}
		else
			delete region;
	}
}


inline
std::map<float, uint> GeneGeneModelReader::ArchiveSnpModels(const char *filename, uint modelCount, uint minImplicationIndex, bool binaryStorage) {
	GeneGeneModelReader::iterator genegenemodels = begin();
	SnpModelCollection models;
	genegenemodels.GetModels(models, modelCount, minImplicationIndex);
	SnpModelCollection::iterator itr = models.begin();
	SnpModelCollection::iterator end = models.end();
	modelCount = models.size();

	std::map<float, uint> scores;

	std::ofstream file(filename);
	if (binaryStorage) {
		file.close();
		file.open(filename, std::ios::binary);
		uint empty = 0;
		file.write((char*)&empty, sizeof(uint));
		file.write((char*)&modelCount, sizeof(uint));
	}
	else 
		file<<modelCount<<"\n";
	
	while (itr != end) {
		SnpSnpModel *model = *itr;
		float ii = model->ImplicationIndex();
		if (scores.find(ii) == scores.end())
			scores[ii] = 1;
		else
			scores[ii]++;
		model->Write(file, binaryStorage);
		itr++;
	}
	return scores;
}

}
