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
void Solution::output_solution(std::ostream& os){
    // simply output each symbol stored in the solution
    std::vector<std::string>::iterator iter;
    for(iter = symbols.begin(); iter != symbols.end(); ++iter){
        os << *iter << " ";
    }
    os << std::endl;
}


/// outputs a more human-readable version of the network
void Solution::output_clean(std::ostream& os,  data_manage::Dataholder& data, bool map_used,
  bool ott_dummy, bool continmap_used){
  output_solution(os);
}

///
/// outputs as a graphviz compatible file
///
void Solution::output_graph(std::ostream& os){
   os << "Not implemented for " << solution_name << std::endl;
}


///
/// copies contents of 2 Solutions
///
void Solution::copy(Solution* other){
    
//    cout << "copying this=" << this << " other= " << other << endl;
    symbols = other->symbols;
    solution_name = other->solution_name;
    sol_fitness = other->sol_fitness;
    test_score = other->test_score;   
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

