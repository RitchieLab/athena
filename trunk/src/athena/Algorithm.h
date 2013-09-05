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
		virtual void setParams(AlgorithmParams& algParams, int numExchanges, int numGenos, 
			int numContin)=0;
		
		/// Set current Dataset for running algorithm
		virtual void setDataset(Dataset* newSet){
				set = newSet;
		}
		
		void setDummyEncoding(bool ottDummyEncoded){dummyEncoded = ottDummyEncoded;}
		
		/// Returns population 
		Population getPop(){return pop;}
		
		/// Runs alorithm
		virtual void run()=0;
		
		/// Runs algorithm for a step which is whatever duration set by user, returns true if completed
		virtual int step()=0;
		
		/// Sets testing values for best solutions
		virtual void testSolution(Dataset* testSet, int nproc)=0;
		
		/// Outputs individual evaluations to stream
		virtual void outputIndEvals(Dataset* set, ostream& os, int model)=0;
		
		/// Initializes the algorithm for each new dataset to test
		virtual void initialize()=0;
		
		/// Sets random seed
		virtual void setRand(unsigned int seed)=0;
		
		/// Returns population of results
		Population& getPopulation(){return pop;}
		
		/// Writes output to stream
		virtual void writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed)=0;
		
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
		virtual void closeLog()=0;
		
		/// Retrieves the models from BioFilter and stores the information in the algorithm
		virtual void getBioModels(std::string fileName, std::string bioFileType, 
			data_manage::Dataholder* holder)=0;
		
		/// Retrieves the models from BioFilter archive and stores in algorithm
		virtual void getBioModelsArchive(string geneGeneFile, string archiveFile, 
			data_manage::Dataholder* holder)=0;
		
		virtual std::string getFitnessName(){return fitnessName;}
		
		/// Get column names for output
		virtual vector<std::string> getAdditionalOutputNames(){
			vector<std::string> temp; 
			return temp;
		}
		
		/// Get final values for best model for reporting (such as AUC)
		virtual	vector<std::string> getAdditionalFinalOutput(Dataset* set){
			vector<std::string> temp; 
			return temp;
		}
		
		virtual void setFitnessName(std::string fname){fitnessName = fname;}
		
		virtual void setConfigDefaults(Config& configuration,AlgorithmParams& algParam)=0;
		
		/// Select best model from models passed and return
		virtual void selectBestModel(std::vector<Solution*>& solutions, data_manage::Dataholder * holder,
			Dataset* set, Config& configuration)=0;
		
		#ifdef PARALLEL
			virtual void setRank(int rank){myRank = rank;}
			int getRank(){return myRank;}
			void setTotalNodes(int total){totalNodes = total;}
			int getTotalNodes(){return totalNodes;}
		#endif
		
protected:
		bool dummyEncoded;
		data_manage::Dataset* set;
		LogType logTypeSelected;
		Population pop;
		std::string fitnessName;
		
		vector<AlgorithmLog*> logs;
		int myRank;
		
		int totalNodes;
};


#endif	/* _ALGORITHM_H */

