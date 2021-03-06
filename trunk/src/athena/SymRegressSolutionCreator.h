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
// SymRegressSolutionCreator.h

#ifndef _SYMBREGRESSSOLUTIONCREATOR_H
#define	_SYMBREGRESSSOLUTIONCREATOR_H

#include "NNSolutionCreator.h"

///
/// Inherits from the NNSolutionCreator class.  
/// Adds functionality allowing which assigns worst fitness to any network
/// that doesn't contain all the variables in the grammar.
///

class SymRegressSolutionCreator: public NNSolutionCreator{

	public:
		/// creates solution from vector of strings
		virtual void establishSolution(vector<string>& symbols, Dataset* set);
 
		 /// creates solution from vector of strings
		virtual void establishSolution(vector<string>& symbols);

		/// writes a dot compatible text file representing the network
		virtual void graphical_output(ostream& os, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy){}
		
		/// writes output as an equation	
		virtual void equationOutput(ostream& os, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed);	
		
		virtual void establishSolutionEquation(std::vector<std::string>& symbols);
		
	protected:
	
		std::string alterLabel(data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, string label, bool continMapUsed, bool equationOut);
	
		 /// evaluates single individual and returns value for that individual
		float evaluateInd(Individual* ind);
		
		vector<std::string> equationStack;
		
};

#endif
