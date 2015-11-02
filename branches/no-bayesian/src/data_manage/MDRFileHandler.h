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
#ifndef MDRFILEHANDLER_H_
#define MDRFILEHANDLER_H_

#include "FileHandler.h"

namespace data_manage
{

class MDRFileHandler : public data_manage::FileHandler
{
public:

	MDRFileHandler();

	virtual ~MDRFileHandler();

	/// provides interface for filling dataholder object with genotypes
	void parseFile(std::string filename, Dataholder* holder, int missingValue,
		float statusMissingValue, bool containsID=false);

	/// parses testing and training files
	void parseFile(std::string trainFile, std::string testFile, Dataholder* holder,
	  int missingValue, float statusMissingValue, bool containsID=false);

	/// provides interface for writing out dataholder object
	void writeFile(std::string filename, Dataholder* holder);
	
private:
	int dummyID;

};

}

#endif /*MDRFILEHANDLER_H_*/
