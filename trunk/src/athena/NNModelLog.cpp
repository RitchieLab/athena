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

void NNModelLog::header(ostream& os){
  os << "GEN\tRANK\tFITNESS\tGRAM_DEPTH\tNN_DEPTH\tNUM_G\tNUM_C";
  if(detailed)
    os << "\tMODEL";
  os << "\n";
}

void NNModelLog::open_log(std::string filename){
  log_stream.open(filename.c_str(), ios::out);
  if(!log_stream.is_open()){
    throw AthenaExcept(filename + " unable to open for writing results");
  }
  header(log_stream);
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
