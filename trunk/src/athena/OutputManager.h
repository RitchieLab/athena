/* 
 * File:   OutputManager.h
 * Author: dudeksm
 *
 * Created on December 29, 2008, 3:37 PM
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
      bool mapfile_used=false, bool dummy_encoded=true, bool continmap_used=true);
    
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

