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
    NNModelLog(){detailed=false;}
    
    ~NNModelLog();

    void open_log(std::string filename, std::string fitness_name);
    
    void close_log();
    
    void write_solution(NNSolution & solution, int generation, int rank);
    
    void set_detailed(bool include){detailed=include;}
    
    bool get_detailed(){return detailed;}
    
private:
    
    bool detailed;
    
    void header(ostream& os, std::string fitness_name);
    
    std::ofstream log_stream;

};

#endif