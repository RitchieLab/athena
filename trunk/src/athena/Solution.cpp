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

