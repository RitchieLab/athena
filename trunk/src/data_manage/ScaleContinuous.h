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
#ifndef SCALECONTINUOUS_H_
#define SCALECONTINUOUS_H_

#include "ScaleData.h"

namespace data_manage
{

class ScaleContinuous : public data_manage::ScaleData
{
public:
  ScaleContinuous();
  ~ScaleContinuous();

  void adjust_contin(Dataholder* holder, unsigned int var_index);
  
  /// Scales status by dividing all 
  void adjust_status(Dataholder* holder);

  /// Returns original value for a scaled status value
  float get_original_status(float status){return (status*statmax)-statadjust;}
  
  /// Writes scaling information to ostream
  string output_scale_info();
  
private:
  
  float statmin, statmax, statadjust;
  

};

}

#endif /*SCALECONTINUOUS_H_*/
