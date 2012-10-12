/* 
 * File:   OutputManager.h
 * Author: dudeksm
 *
 * Created on December 29, 2008, 3:37 PM
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

#ifndef _OUTPUTMANAGER_H
#define	_OUTPUTMANAGER_H

#include "Population.h"
#include "AthenaExcept.h"
#include <vector>
#include <fstream>
#include "Algorithm.h"


///
/// Ouputs best models from the analysis as a summary
/// and as individual files showing the results.
///

class OutputManager{
    
public:
    
    /// sets basename for output
    void setBasename(std::string base){basename = base;}
    
    /// outputs summary of best models
    void outputSummary(std::vector<Population>& pops, data_manage::Dataholder& data, 
      bool mapfile_used=false, bool dummy_encoded=true, bool continmap_used=false,
      std::string fitness_name=" ");
    
    /// outputs a file for each best model
    void outputBestModels(std::vector<Population>& pops, int nmodels, std::string scaleInfo,
      data_manage::Dataholder& data, bool map_used, bool ott_dummy, bool continmap_used);
    
    /// returns a stream for writing
    std::ostream& getStream(std::string filename);
    
    /// closes the provided stream
    void closeStream(){if(log_stream.is_open()){ log_stream.close();}}
    
    /// output graphical representation as defined in algorithm
    void outputGraphic(Algorithm* alg,  std::vector<Population>& pops,std::string basename, int numModels,
      data_manage::Dataholder& data, bool map_used, bool ott_dummy, bool continmap_used);
    
private:
    
    std::string basename;
    std::ofstream log_stream;
    
};

#endif	/* _OUTPUTMANAGER_H */

