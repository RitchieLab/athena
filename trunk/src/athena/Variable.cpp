/*
Copyright Marylyn Ritchie 2015

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
#include "Variable.h"

Variable::Variable(int idx){
	index=idx;
}

GenoVariable::GenoVariable(int index):Variable(index){
	gType = true;
	// one level for missing
	nLevels=4;
}

///
/// Return name of variable
/// @param holder Dataholder
/// @returns name
///
string GenoVariable::getName(data_manage::Dataholder* holder){
	return holder->getGenoName(index);
}


///
/// Returns value for individual
/// @param ind Individual
///
float GenoVariable::getValue(data_manage::Individual* ind){
	return ind->getGenotype(index);
}

///
/// Returns missing for variable
/// @returns missing value;
///
// float GenoVariable::getMissingVal(){
// 	return 3.0;
// }

ConVariable::ConVariable(int index):Variable(index){
	gType=false;
}

///
/// Returns value for individual for continuous variable
/// @param ind Individual
///
float ConVariable::getValue(data_manage::Individual* ind){
	return ind->getCovariate(index);
}

///
/// Return name of variable
/// @param holder Dataholder
/// @returns name
///
string ConVariable::getName(data_manage::Dataholder* holder){
	return holder->getCovarName(index);
}


///
/// Returns missing for variable
/// @returns missing value;
///
// float ConVariable::getMissingVal(){
// 	return holder->get
// }
