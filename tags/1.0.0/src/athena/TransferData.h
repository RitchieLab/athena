/* 
 * File:   TransferData.h
 * Author: dudeksm
 *
 * Created on January 21, 2008, 2:59 PM
 */

#ifndef _TRANSFERDATA_H
#define	_TRANSFERDATA_H

#ifdef PARALLEL

#include "mpi.h"
#include <Dataholder.h>


///
/// Transfers data from master to slaves using MPI
///
class TransferData{
    
public:

  /// sends all data to slaves 
  void sendData(data_manage::Dataholder& holder);
  
  /// receives all data from slaves
  void receiveData(data_manage::Dataholder& holder);
     
private:
    
    
};

#endif /* if PARALLEL defined */

#endif	/* _CONFIG_H */

