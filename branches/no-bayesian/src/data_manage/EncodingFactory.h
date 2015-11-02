//EncodingFactory.h

/*
Copyright Marylyn Ritchie 2013

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
#ifndef ENCODINGFACTORY_H_
#define ENCODINGFACTORY_H_

#include "DummyConvert.h"
#include "DataExcept.h"

namespace data_manage{

class EncodingFactory{

	public:

		static DummyConvert* createEncoder(string encodeType);
		
	private:
		
		enum EncodeType{
			NoMatch,
			NoEncode,
			Custom,
			Add_quad,
			Additive
		};

		static map<string, EncodeType> encodeMap;
		
		/// Sets the map
		static void setEncodeMap();
		
};


}

#endif /*ENCODINGFACTORY_H_*/
