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

#include "OptGrammarFactory.h"

const std::string& OptGrammarFactory::registerCalc(const std::string& key, createOptFunc* ptr){
	creation_map[key] = ptr;
	return key;
}

OptGrammar* OptGrammarFactory::create(const std::string& key){
	std::map<std::string, createOptFunc*>::const_iterator it=creation_map.find(key);
	if(it != creation_map.end()){
		return (*it).second();
	}else{
		throw AthenaExcept(key + " is not a valid optimization type");
	}
}



