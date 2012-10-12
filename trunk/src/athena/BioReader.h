//BioReader.h
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

#ifndef __BIOREADER_H__
#define __BIOREADER_H__

#include <vector>
#include <string>
#include <fstream>

/// Contains information related to Biofilter models provided
struct BioModel{
  std::vector<std::string> idstring;
  float implication_index;
};


///
/// Abstract base class for bio filter readers
///
class BioReader{

  public:
    virtual ~BioReader(){}

    /// Fills vector with models from file
    virtual int GetModels(std::vector<BioModel>& models, std::string filename, unsigned int max_read)=0;

};

#endif
