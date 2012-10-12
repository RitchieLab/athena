/* 
 * File:   TransferData.h
 * Author: dudeksm
 *
 * Created on January 21, 2008, 2:59 PM
 */
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

