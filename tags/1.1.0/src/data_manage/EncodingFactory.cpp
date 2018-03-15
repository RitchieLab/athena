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
#include "EncodingFactory.h"
#include "StephenDummyConvert.h"
#include "OttDummyConvert.h"
#include "CustomDummyConvert.h"
#include "Stringmanip.h"

namespace data_manage{

map<string, EncodingFactory::EncodeType> EncodingFactory::encodeMap;

///
/// Function creates a data scaler based on the name passed
/// @param scaleType name of scaling type
/// @return pointer to new ScaleData
/// @throws DataExcept when no match
///
DummyConvert* EncodingFactory::createEncoder(string encodeType){
	
	if(encodeMap.empty()){
		setEncodeMap();
	}
	
	DummyConvert* newEncoder;
	
	vector<string> tokens = Stringmanip::split(encodeType, ':');
	vector<int> conv;
	encodeType = tokens[0];
	
	switch(encodeMap[encodeType]){
		case NoMatch:
			throw DataExcept("No data encoder matching " + encodeType);
			break;
		case NoEncode:
			newEncoder = new DummyConvert;
			break;
		case Add_quad:
			newEncoder = new OttDummyConvert;
			break;
		case Additive:
			newEncoder = new StephenDummyConvert;
			break;
		case Custom:
			for(vector<string>::iterator iter=tokens.begin()+1; iter != tokens.end(); iter++){
				conv.push_back(Stringmanip::stringToNumber<int>(*iter));
			}
			newEncoder = new CustomDummyConvert(conv);
			break;

		default:
			throw DataExcept("No data encoder matching " + encodeType);
	}
	
	return newEncoder;
	
}


///
/// Establishes the map for use in creating solutions
/// @return 
///
void EncodingFactory::setEncodeMap(){
	encodeMap["ADD_QUAD"]=Add_quad;
	encodeMap["OTT"]=Add_quad;
	encodeMap["TRUE"]=Add_quad;
	encodeMap["STEPHEN"]=Additive;
	encodeMap["ADDITIVE"]=Additive;
	encodeMap["CUSTOM"]=Custom;
	encodeMap["NONE"]=NoEncode;
	encodeMap["FALSE"]=NoEncode;

}

}
