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
#include "Solution.h"
#include "Terminals.h"

///    
/// outputs as a text file
///
void Solution::outputSolution(std::ostream& os){
		std::vector<std::string>::iterator iter;
		for(iter = symbols.begin(); iter != symbols.end(); ++iter){
				os << *iter << " ";
		}
		os << std::endl;
}


/// outputs a more human-readable version of the network
void Solution::outputClean(std::ostream& os,  data_manage::Dataholder& data, bool mapUsed,
	bool ottDummy, bool continMapUsed){
	outputSolution(os);
}


///
/// outputs as a graphviz compatible file
///
void Solution::outputGraph(std::ostream& os){
	 os << "Not implemented for " << solutionName << std::endl;
}


///
/// copies contents of 2 Solutions
///
void Solution::copy(Solution* other){
		symbols = other->symbols;
		solutionName = other->solutionName;
		solFitness = other->solFitness;
		testScore = other->testScore;
		complexity = other->complexity;
}


///
/// Returns a clone of this solution.  It is dynamically allocated
/// and memory management for it is the responsibility of the 
/// caller
///
Solution* Solution::clone(){
		Solution* newSol = new Solution;
		newSol->copy(this);
		return newSol;
}

