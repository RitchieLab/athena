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

#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_

#include "Dataholder.h"
#include "DataExcept.h"

namespace data_manage
{

/// Provides interface for multiple data files
class FileHandler
{
public:
	/// virtual destructor
	virtual ~FileHandler();
	
	/// provides interface for filling dataholder object with genotypes
	virtual void parse_file(std::string filename, Dataholder* holder, int missingValue,
           float statusMissingValue, bool contains_id=false)=0; 
	
	/// provides interface for filling dataholder object with genotypes from a testing
	/// and a training file
	virtual void parse_file(std::string train_file, std::string test_file, Dataholder* holder,
	  int missingValue, float statusMissingValue, bool contains_id=false)=0;
	
	/// provides interface for writing out dataholder object
	virtual void write_file(std::string filename, Dataholder* holder)=0;
	
};

}

#endif /*FILEHANDLER_H_*/
