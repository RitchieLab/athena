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
//BioFilterReader.cpp

#include "BioFilterReader.h"
#include "AthenaExcept.h"
#include <sstream>
#include <iostream>

using namespace std;

///
/// Constructor
///
BioFilterReader::BioFilterReader(){
	initialize();
}

///
/// Initializes object
///
void BioFilterReader::initialize(){
	maximumReads = 500;
}



///
/// Fills vector with models from file
/// @param models Vector will contain models in order read from the file
/// @param filename File to read from
/// @param max_read Maximum number of models to read
/// @return Number of models actually read
///
int BioFilterReader::getModels(std::vector<BioModel>& models, string filename, unsigned int maxRead){
	
	maximumReads = maxRead;

	if(!reader.is_open()){
		reader.open(filename.c_str(), ios::in);
		if(!reader.is_open())
			throw AthenaExcept("Unable to open bio filter file " + filename);
	}
	
	models.clear();
	
	// assume 2 locus models for now
	string id, line;
	while(!reader.eof() && (models.size() < maximumReads)){
		getline(reader, line);
		if(line.find_first_of("0123456789") == string::npos)
			continue;
		stringstream ss(line);
		BioModel mod;
		ss >> id;
		mod.idString.push_back(id);
		ss >> id;
		mod.idString.push_back(id);
		ss >> mod.implicationIndex;
		models.push_back(mod);
	}

	// when all models done reading close the stream
	if(models.size() < maximumReads)
		reader.close();

	return models.size();
}
