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
