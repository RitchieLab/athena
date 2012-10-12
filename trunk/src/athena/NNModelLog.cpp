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
//NNModelLog.cpp

#include "NNModelLog.h"
#include <iomanip>

///
/// Destructor
///
NNModelLog::~NNModelLog(){
    if(log_stream.is_open()){
        close_log();
    }
}

void NNModelLog::header(ostream& os, std::string fitness_name){
  os << "GEN\tRANK\t" << fitness_name  << 
    " FITNESS\tGRAM_DEPTH\tNN_DEPTH\tNUM_G\tNUM_C";
  if(detailed)
    os << "\tMODEL";
  os << "\n";
}

void NNModelLog::open_log(std::string filename, std::string fitness_name){
  log_stream.open(filename.c_str(), ios::out);
  if(!log_stream.is_open()){
    throw AthenaExcept(filename + " unable to open for writing results");
  }
  header(log_stream, fitness_name);
}

///
/// Closes output stream
///
void NNModelLog::close_log(){
    log_stream.close();
}
    
///
/// Add solution to log file
///
void NNModelLog::write_solution(NNSolution & solution, int generation, int rank){
    log_stream << generation << "\t";
    log_stream << rank << "\t";
    log_stream << solution.fitness() << "\t";
    log_stream << solution.get_gram_depth() << "\t";
    log_stream << solution.get_nn_depth() << "\t";
    log_stream << solution.get_genotypes().size() << "\t";
    log_stream << solution.get_covariates().size() << "\t";
    if(detailed)
        solution.output_solution(log_stream);
    else
        log_stream << "\n";
}
