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
#ifndef MAPFILEREADER_H_
#define MAPFILEREADER_H_

#include "Dataholder.h"
#include "DataExcept.h"

namespace data_manage
{

class MapFileReader
{
public:

	/// parses the map file and stores the names in the Dataholder
	void parseMapFile(std::string mapFile, Dataholder* holder);

};

}

#endif /*MAPFILEREADER_H_*/
