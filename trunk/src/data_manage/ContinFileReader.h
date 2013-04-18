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
#ifndef CONTINFILEREADER_H_
#define CONTINFILEREADER_H_

#include "Dataholder.h"
#include "DataExcept.h"

namespace data_manage
{

///
/// Reads continuous variables from file and stores in Dataholder.  The
/// order of individuals must match that of the original data file.
///
class ContinFileReader
{
public:
	ContinFileReader();
	~ContinFileReader();

	void readContinFile(std::string filename, Dataholder* holder, float missingValue,
					bool containsID=false);
					
	void readContinFile(std::string trainFile, std::string testFile, Dataholder* holder, 
		float missingValue, bool containsID=false);  

private:
	int dummyID;
};

}

#endif /*CONTINFILEREADER_H_*/
