/*
Copyright Marylyn Ritchie 2014

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
/*
 * File:   SumFileReader.h
 * Author: dudeksm
 *
 * Created on Wed Jun  4 09:22:49 CDT 2014
 */

#ifndef _SUMFILEREADER_H
#define	_SUMFILEREADER_H

#include "Population.h"

class SumFileReader{

  public:
    void readSumFile(std::string filename);

    vector<Solution*>& getModelPopulation(){return modelPop;}

  private:
    vector<Solution*> modelPop;

};


#endif
