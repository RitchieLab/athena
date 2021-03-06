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
//BioFilterReader.h

#ifndef __BIOFILTERREADER_H__
#define __BIOFILTERREADER_H__

#include <vector>
#include <string>
#include <fstream>
#include "BioReader.h"


/// Reads models from file 
class BioFilterReader: public BioReader{

	public:
	
		BioFilterReader();
		
		/// Fills vector with models from file
		int getModels(std::vector<BioModel>& models, std::string filename, unsigned int maxRead);
	
	private:
	
		void initialize();
	
		std::string filename;
		unsigned int maximumReads;
		
		std::ifstream reader;
		
};

#endif


