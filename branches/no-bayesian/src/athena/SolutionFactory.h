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
/*
 * File:   SolutionFactory.h
 * Author: dudeksm
 *
 * Created on November 13, 2008, 4:24 PM
 */

#ifndef _SOLUTIONFACTORY_H
#define	_SOLUTIONFACTORY_H

#include <map>
#include "SolutionCreator.h"
#include "AthenaExcept.h"

///
/// Creates and returns solution based on name passed
///
class SolutionFactory{

public:

		static SolutionCreator* createSolution(std::string solution_name);
		static SolutionCreator* createSolution(std::string solution_name, std::vector<std::string>& vars);

private:

		enum SolutionType{
				/// For no match
				MissingSolutionType,
				/// For neural network
				NNSolutionType,
				/// For neural network requiring all variables to be fit
				NNSolutionAllType,
				/// For neural network restricting the networks to only one occurrence of each variable
				NNSolutionOnceType,
				/// For Symbolic Regression
				SymRegressSolutionType
		};

		static std::map<std::string, SolutionType> solutionMap;

		/// Sets map for use in Solution creation
		static void setSolutionMap();

};

#endif	/* _SOLUTIONFACTORY_H */

