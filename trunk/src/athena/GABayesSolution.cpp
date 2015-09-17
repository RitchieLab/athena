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
#include "GABayesSolution.h"
#include "CalculatorFactory.h"
#include <Stringmanip.h>
#include <sstream>

#include <iostream>

///
/// Clones current solution
///
Solution* GABayesSolution::clone(){
		GABayesSolution* newGASol = new GABayesSolution;
		newGASol->setGenotypes(genos);
		newGASol->setCovariates(contins);

		Solution* newSol = newGASol;
		Solution* thiscopy = this;

		newSol->copy(thiscopy);
		return newSol;
}


///
/// Cleans up output to remove extra rule information. Changes Concat operator
/// into corresponding numbers for output
/// @param os ostream to write to
/// @param data
/// @param mapUsed
///
void GABayesSolution::outputClean(std::ostream& os, data_manage::Dataholder& data,
			bool mapUsed, bool ottDummy, bool continMapUsed){

	for(size_t i=0; i<symbols.size(); i++){
		os << symbols[i];
	}
}


