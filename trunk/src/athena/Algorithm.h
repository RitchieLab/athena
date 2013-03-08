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
/* 
 * File:   Algorithm.h
 * Author: dudeksm
 *
 * Created on November 10, 2008, 1:18 PM
 */

#ifndef _ALGORITHM_H
#define	_ALGORITHM_H

#include <Dataset.h>
#include <Stringmanip.h>
#include "Config.h"
#include "Population.h"
#include "AthenaExcept.h"
#include "AlgorithmLog.h"
#include <Dataholder.h>

using namespace data_manage;

///
/// Abstract base class for algorithms
///

class Algorithm{
    
public:
    Algorithm();
    
    virtual ~Algorithm();
    
    /// Set the parameters for the algorithm
    virtual void set_params(AlgorithmParams& alg_params, int numExchanges, int numGenos, int numContin)=0;
    
    /// Set current Dataset for running algorithm
    virtual void set_dataset(Dataset* new_set){
        set = new_set;
    }
    
    void set_dummy_encoding(bool ott_dummy_encoded){dummy_encoded = ott_dummy_encoded;}
    
    /// Returns population 
    Population get_pop(){return pop;}
    
    /// Runs alorithm
    virtual void run()=0;
    
    /// Runs algorithm for a step which is whatever duration set by user, returns true if completed
    virtual int step()=0;
    
    /// Sets testing values for best solutions
    virtual void test_solution(Dataset* test_set, int nproc)=0;
    
    /// Outputs individual evaluations to stream
    virtual void output_ind_evals(Dataset* set, ostream& os, int model)=0;
    
    /// Initializes the algorithm for each new dataset to test
    virtual void initialize()=0;
    
    /// Sets random seed
    virtual void setrand(unsigned int seed)=0;
    
    /// Returns population of results
    Population& getPopulation(){return pop;}
    
    /// Writes output to stream
    virtual void writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy, bool continmap_used)=0;
    
    /// Returns extension for graphical representation
    virtual std::string getGraphicalFileExt()=0;
    
    /// Saves the log
    virtual void saveLog()=0;
    
    /// Starts log
    virtual void startLog(int num_snps)=0;

    /// Writes the log
    virtual void writeLog()=0;
    
    /// Clears the logs
    virtual void clearLogs()=0;
    
    /// Finish log and pulls together model information
    virtual void finishLog(std::string basename, int cv)=0;
    
    /// Sets log type 
    virtual void setLogType(LogType ltype){logTypeSelected = ltype;}
    
    /// Prepares log files
    virtual void prepareLog(std::string basename, int cv)=0;
    
    /// Close log files
    virtual void CloseLog()=0;
    
    /// Retrieves the models from BioFilter and stores the information in the algorithm
    virtual void getBioModels(std::string filename, std::string biofiletype, 
      data_manage::Dataholder* holder)=0;
    
    /// Retrieves the models from BioFilter archive and stores in algorithm
    virtual void getBioModelsArchive(string genegeneFile, string archiveFile, 
      data_manage::Dataholder* holder)=0;
    
    virtual void tempoOutputName(string outname)=0;
    
    virtual std::string get_fitness_name(){return fitness_name;}
    
    virtual void set_fitness_name(std::string fname){fitness_name = fname;}
    
    #ifdef PARALLEL
      virtual void setRank(int rank){myRank = rank;}
      int getRank(){return myRank;}
      void setTotalNodes(int total){totalNodes = total;}
      int getTotalNodes(){return totalNodes;}
    #endif
    
protected:
    bool dummy_encoded;
    data_manage::Dataset* set;
    LogType logTypeSelected;
    Population pop;
    std::string fitness_name;
    
    vector<AlgorithmLog*> logs;
    int myRank;
    
//    #ifdef PARALLEL
      int totalNodes;
//    #endif
    
};


#endif	/* _ALGORITHM_H */

