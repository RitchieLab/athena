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
 * File:   GEBayes.h
 * Author: dudeksm
 *
 * Created on March 5, 2014, 11:03 AM
 */
 
 
#ifndef _GEBAYES_H
#define	_GEBAYES_H

#include "Algorithm.h"
#include "AlgorithmFactory.h"
#include "BioFilterModelCollection.h"
#include "Config.h"
#include "BayesSolution.h"
#include "NNLog.h"
#include "BayesModelLog.h"

#ifdef PARALLEL
	#include "GenomeTransfer.h"
#endif

class GEBayes:public Algorithm{

public:

	GEBayes();
		
	~GEBayes();

	/// Set the parameters for the algorithm
	virtual void setParams(AlgorithmParams& algParams, int numExchanges, int numGenos, 
		int numContin, vector<unsigned int>& excludedGenos, 
		vector<unsigned int>& excludedContins);
		
	/// Run algorithm 
	void run();
		
	/// Runs a step of the algorithm
	virtual int step();
		
	/// Sets random seed
	void setRand(unsigned int seed);
		
	void setDataset(Dataset* newSet);
	
	/// Initializes the algorithm for each new dataset to test
	void initialize();
	
	/// Sets configuration to any defaults specific to this algorithm
	void setConfigDefaults(Config& configuration,AlgorithmParams& algParam);
	
	/// Sets testing values for best solutions
	void testSolution(Dataset* testSet, int nproc);
	
	/// Select best model from models passed and return
	virtual void selectBestModel(std::vector<Solution*>& solutions, data_manage::Dataholder * holder,
		Dataset* set, Config& configuration){}
	
	/// Outputs individual evaluations to stream
	virtual void outputIndEvals(Dataset* set, ostream& os, int model){}
	
	/// Writes output to stream
	virtual void writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
		bool mapUsed, bool ottDummy, bool continMapUsed);
			
	virtual void writeEquation(ostream& os, Solution* sol, data_manage::Dataholder* holder,
		bool mapUsed, bool ottDummy, bool continMapUsed);

	virtual void produceGraphic(std::string inputGraphic, std::string outputGraphic, 
		std::string imgWriter){}

	/// Returns extension for graphical representation
	virtual std::string getGraphicalFileExt();

	/// Saves the log
	virtual void saveLog(){};

	/// Starts log
	virtual void startLog(int num_snps);

	/// Writes log
	void writeLog();

	/// Clears the logs
	virtual void clearLogs(){}

	/// Finish log and pulls together model information
	virtual void finishLog(std::string basename, int cv);
		
	/// Prepares log files
	void prepareLog(std::string basename, int cv);
		
	/// Close log files
	virtual void closeLog();
	
	virtual vector<Solution*> runValidation(std::string sumFile);
	
	virtual vector<std::string> getAdditionalOutputNames();
	
	/// Return formatted output for display for final models;
	virtual	void getAdditionalFinalOutput(Dataset* set);	

	#ifdef PARALLEL
		virtual void setRank(int rank);
	#endif
	
	
protected:
	
	void resetCrossover();
	
	/// Sets GA for run
	virtual void setGAParams(vector<unsigned int>& excludedGenos, 
			vector<unsigned int>& excludedContins);
		
	/// Sets default values for parameters
	void initializeParams(); 
	
	void fillPopulation();
	
	void fillLog();
	
	BayesSolution* convertGenome(GAGenome& ind);
	
	enum GEBayesParams{
		noMatchParam,
		parentRange,
		childRange
	};
	
	NNLog* geLog;
	BayesModelLog* modelLog;
	std::string mainLogFilename, fitnessLogFilename, snpnameLogFilename;
	
	std::map<std::string, GEBayesParams> paramMap;	
	unsigned int restrictStepsDone, minNumParents, maxNumParents, minNumChildren, 
		maxNumChildren;

	#ifdef PARALLEL
		GenomeTransfer popMigrator;
	#endif

};


#endif /* _GEBAYES_H */
