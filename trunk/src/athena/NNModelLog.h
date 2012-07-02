// NNModelLog.h

#ifndef _NNMODELLOG_H
#define	_NNMODELLOG_H

#include "NNSolution.h"
#include <fstream>

///
/// Outputs each model in the population and lists the details kept in
/// in the solutions about the model.  The output resembles the .best files
/// but with more information.  Some of the information will be grammar depth,
/// network depth, fitness and the terminal representation of the network
///

class NNModelLog{

public:
    NNModelLog(){};
    
    ~NNModelLog();

    void open_log(std::string filename);
    
    void close_log();
    
    void write_solution(NNSolution & solution, int generation, int rank);
    
private:
    
    void header(ostream& os);
    
    std::ofstream log_stream;

};

#endif