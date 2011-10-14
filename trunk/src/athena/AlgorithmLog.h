//AlgorithmLog.h

#ifndef _ALGORITHMLOG_H
#define	_ALGORITHMLOG_H

#include <iostream>

#ifdef PARALLEL
  #include "mpi.h"
#endif

///
/// Stores data for display on algorithm processing for investigating results.
///
class AlgorithmLog{

  public:

    AlgorithmLog(int num_snps){}

    virtual ~AlgorithmLog(){}

    /// Outputs 
    virtual void output_log(std::ostream& os)=0;
    
    #ifdef PARALLEL
      virtual void sendLog()=0; // for slaves
      virtual void receiveLogs(int nprocs)=0; // for master
    #endif

};


#endif
