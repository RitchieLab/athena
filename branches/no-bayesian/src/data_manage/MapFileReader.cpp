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
#include "MapFileReader.h"
#include <sstream>
#include <fstream>

#include <cstdlib>

namespace data_manage
{

///
/// Parses map file and stores loci names in the Dataholder
/// @param mapFile
/// @param dataholder Dataholder
///
void MapFileReader::parseMapFile(string mapFile, Dataholder* dataholder){
	std::ifstream mapStream(mapFile.c_str(), ios::in);

	if(!mapStream.is_open()){
		throw DataExcept("ERROR: Unable to open " + mapFile + "\n");
	}

	string line;

	string snpID;
//	unsigned int chrom;
//	unsigned int pos;
string chrom,pos;

	while(!mapStream.eof()){
		getline(mapStream, line);

		if(line.find_first_of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") == string::npos){
			continue;
		}

		stringstream ss(line);
		ss >> chrom >> snpID >> pos;
		dataholder->addGenoName(snpID);
	}

	mapStream.close();
}


}

