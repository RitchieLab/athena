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
//BioFilterBinReader.cpp

#include "BioFilterBinReader.h"
#include "AthenaExcept.h"
#include <sstream>
#include <iostream>

using namespace std;

///
/// Constructor
///
BioFilterBinReader::BioFilterBinReader(){
	initialize();
}

///
/// Initializes object
///
void BioFilterBinReader::initialize(){
	maximumReads = 500;
}


///
/// Fills vector with models from file
/// @param models Vector will contain models in order read from the file
/// @param filename File to read from
/// @param maxRead Maximum number of models to read
/// @return Number of models actually read
///
int BioFilterBinReader::getModels(std::vector<BioModel>& models, string filename, 
	unsigned int maxRead){
	
	maximumReads = maxRead;
	if(!reader.is_open()){
		reader.open(filename.c_str(), ios::binary);
		if(!reader.is_open()){
			throw AthenaExcept("Unable to open bio filter file " + filename);
		}
	}
	
	models.clear();
	
	unsigned int modelSize, groupCount, ddCount, ddRedundant, geneCount;
	unsigned int locus, group, ddGroup, gene;
	unsigned int impScore;
	
	while(!reader.eof() && (models.size() < maximumReads) && reader.good()){

		BioModel mod;
	
		// read in all data -- but only using ID of snps and implication index score
		reader.read((char*)&modelSize, sizeof(unsigned int));
		reader.read((char*)&groupCount, sizeof(unsigned int));
		reader.read((char*)&ddCount, sizeof(unsigned int));
		reader.read((char*)&ddRedundant, sizeof(unsigned int));
		reader.read((char*)&geneCount, sizeof(unsigned int));
		reader.read((char*)&impScore, sizeof(float));
		mod.implicationIndex = impScore;
		
		// get model id numbers ( rs numbers with 'rs' removed)
		// put 'rs' back on to ID before setting ID string in model
		for(size_t i=0; i<modelSize; i++){
			reader.read((char*)&locus, sizeof(unsigned int));
			stringstream ss;
			ss << locus;
			mod.idString.push_back("rs" + ss.str());
		}

		// get group IDs -- not currently used in ATHENA
		vector<unsigned int> groups(groupCount,0);
		for(size_t i=0; i<groupCount; i++){
			reader.read((char*)&group, sizeof(unsigned int));
			groups[i] = group;
		}
		
		vector<unsigned int> ddGroups(ddCount, 0);
		for(size_t i=0; i<ddCount; i++){
			reader.read((char*)&ddGroup, sizeof(unsigned int));
			ddGroups[i] = group;
		}
		
		vector<unsigned int> genes(geneCount, 0);
		for(size_t i=0; i<geneCount; i++){
			reader.read((char*)&gene, sizeof(unsigned int));
			genes[i] = gene;
		}
		models.push_back(mod);
	}
	if(models.size() < maximumReads || reader.eof() || !reader.good())
		reader.close();  
	return models.size();
}

