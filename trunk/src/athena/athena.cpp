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
 * File:   athena.cpp
 * Author: dudeksm
 *
 * Created on November 10, 2008, 2:10 PM
 */

//#include <stdlib.h>

#include "ConfigFileReader.h"
#include <Dataholder.h>
#include "MDRFileHandler.h"
#include "ContinFileReader.h"
#include "MapFileReader.h"
#include "ContinMapFileReader.h"
#include "CrossValidator.h"
#include "OttDummyConvert.h"
#include "AlgorithmFactory.h"
#include "OutputManager.h"
#include <iostream>
#include <sstream>
#include "OutputSet.h"
#include <ScaledDataFactory.h>
#include <StephenDummyConvert.h>
#include <time.h>

void exitApp(AthenaExcept& he, int myRank);
void exitApp(DataExcept& de, int myRank);
void adjustSeed(Config& config, int orig_seed, int cv, int nproc, int myRank);
std::string timeDiff(double dif);

int main(int argc, char** argv) {

	int nproc = 1; // only one processor when not running in parallel

 bool mapFileUsed = false, continMapUsed = false;
 int myRank = 0;
	 
#ifdef PARALLEL
 

	// set up MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
#endif /* end PARALLEL code block */
	 
		string versionDate = "4/7/13";
		string execName = "ATHENA";
		string version = "1.0.1";
		 time_t start,end;
		 
		if(argc < 2){
				AthenaExcept he("\n\tATHENA\n\tv" + version + "\n\t" + versionDate +  
						"\n\n\tUsage: ATHENA <config>\n\n");
				exitApp(he, myRank);
		}
		else{
#ifdef PARALLEL
	if(myRank==0){
#endif
				time (&start);
				cout << endl << "\t" << execName << ":\t" << versionDate << endl << endl;
#ifdef PARALLEL
				}
#endif
		}
		
		string configFile = argv[1];
		ConfigFileReader configRead;
		Config config;
		ScaleData* scaler = NULL;

		// read config file
		try{
			config = configRead.readConfig(configFile);
		}catch(AthenaExcept he){
				exitApp(he, myRank);
		}
		
		// fill dataholder with data
		data_manage::Dataholder data;

		try{
			// read in genotype data
				data_manage::MDRFileHandler mdrReader;
				if(config.getTrainFile().size() ==0)
					mdrReader.parseFile(config.getDataSetName(), &data, config.getMissingValue(),
							 config.getStatusMissingValue(),config.getIDinData());
				else
					mdrReader.parseFile(config.getTrainFile(), config.getTestFile(), 
						&data, config.getMissingValue(), config.getStatusMissingValue(), 
						config.getIDinData());
				
				// read in continuous data if any
				if(config.getContinTrainFile().size() > 0){
					data_manage::ContinFileReader continReader;
					continReader.readContinFile(config.getContinTrainFile(), config.getContinTestFile(),
						&data, config.getContinMiss(), config.getIDinData());          
				}
				else if(config.getContinFileName().size() > 0){
						data_manage::ContinFileReader continReader;
						continReader.readContinFile(config.getContinFileName(), &data,
										config.getContinMiss(), config.getIDinData());
				}
				
				// if present read map file
				if(config.getMapName().size() > 0){
						mapFileUsed = true;
						data_manage::MapFileReader mapReader;
						mapReader.parseMapFile(config.getMapName(), &data);
				}
				else{
						mapFileUsed = false;
						data.addDefaultSnps();
				}
				
				if(config.getContinMapName().size() > 0){
						continMapUsed = true;
						data_manage::ContinMapFileReader continMapReader;
						continMapReader.parseMapFile(config.getContinMapName(), &data);
				}
				else{
						continMapUsed = false;
						data.addDefaultCovars();
				}
			 
				// convert data to ott dummy representation if needed
				if(config.getOttEncoded()){
					data_manage::OttDummyConvert ott;
					ott.convertGenotypes(&data);
				}
				else if(config.getEncodeType() == Config::StephenDummy){
					data_manage::StephenDummyConvert st;
					st.convertGenotypes(&data);
				}

				// alter continuous variables and status value
				scaler = data_manage::ScaledDataFactory::createScaler(config.getStatusAdjust());
				scaler->adjustStatus(&data);
				for(unsigned int c=0; c < data.numCovariates(); c++){
					scaler->adjustContin(&data, c);
				}

		}catch(AthenaExcept he){
				exitApp(he, myRank);
		}
 
		// set random seed  before splitting
		srand(config.getRandSeed());

		// construct crossvalidation sets to use in running algorithm
		CrossValidator cvMaker;
		CVSet cvSet;

	if(config.getSplitFile()==""){
	    cvSet = cvMaker.splitData(config.getNumCV(), &data);
#ifdef PARALLEL
		if(myRank==0)
#endif
	    cvMaker.saveSplits(config.getOutputName() + ".cvsplit");
	}
	else{
	    try{
			cvSet = cvMaker.loadSplits(config.getSplitFile(), &data);
		}catch(DataExcept de){
				exitApp(de, myRank);
		}
	}
	
		// run crossvalidations and store the populations
		int numCV = cvSet.numIntervals();
			 
		// create algorithm
		vector<AlgorithmParams> algParams= config.getAlgorithmParams();
		Algorithm* alg = AlgorithmFactory::createAlgorithm(algParams[0].name);
#ifdef PARALLEL
		alg->setRank(myRank);
		alg->setTotalNodes(nproc);
#endif
		alg->setRand(config.getRandSeed());
		alg->setDummyEncoding(config.getOttEncoded());
		alg->setLogType(config.getLogType());
		
		alg->setParams(algParams[0], config.getNumExchanges(), 
			data.numGenos(), data.numCovariates());

		/// store results in a population vector
		vector<Population> pops;
		string cvFileName = "cv";

		OutputManager writer;
		writer.setBasename(config.getOutputName());
		int currCV=config.getStartCV()-1;
	if(currCV==0)
		writer.setFiles(mapFileUsed, alg->getFitnessName());
	  int originalSeed = config.getRandSeed();
		for(; currCV < numCV; currCV++){
			adjustSeed(config, originalSeed,currCV, nproc, myRank);
			alg->setRand(config.getRandSeed());
#ifdef PARALLEL
	if(myRank==0){
#endif
			cout << "Beginning Cross-validation " << currCV + 1 << "...";
			cout.flush();
#ifdef PARALLEL
	}
#endif

			alg->setDataset(&(cvSet.getInterval(currCV).getTraining()));

#ifdef PARALLEL
		if(myRank==0){  // only have master output cv when desired
#endif
			// output cv interval when requested
			if(config.getCVOutput()){
				OutputSet oSet;
				oSet.outputCV(cvFileName, cvSet.getInterval(currCV).getTraining(), currCV+1);
			}
#ifdef PARALLEL
	} /* end check for master writing CV splits */
#endif
			alg->startLog(data.numGenos());

		if(config.getBioFilterFile().size() > 0){
			alg->getBioModels(config.getBioFilterFile(), config.getBioFileType(), &data);
		}
		else if(config.getBioGeneFile().size() > 0){
			alg->getBioModelsArchive(config.getBioGeneFile(), config.getBioArchiveFile(), &data);
		}
			// prepare logs -- these are written as job progresses
			alg->prepareLog(config.getOutputName(), currCV+1);
			try{
					alg->initialize();
		}catch(AthenaExcept& ae){
				exitApp(ae, myRank);
		}
			
			for(int step=0; step < config.getNumExchanges(); step++){
				 if(alg->step()){
						break; // can complete early
				 }
			}

		alg->closeLog();
			if(numCV > 1)
				alg->testSolution(&(cvSet.getInterval(currCV).getTesting()), nproc);
			
			// check population values
	  Population pop = alg->getPopulation();

		int currProc = 0;
#ifdef PARALLEL
	if(myRank==0){
		int lastProc = 1;
		if(config.outputAllNodesBest())
			lastProc = nproc;
		for(currProc=0; currProc < lastProc; currProc++){
#endif
			if(config.getIndOutput()){
				stringstream ss;
				ss << config.getOutputName() << "." << currCV+1 << "." << currProc+1 << ".ind_results.txt";
				ostream & os = writer.getStream(ss.str());
				os << "Ind ID\tPredicted\tOriginal\n";
				alg->outputIndEvals(&(cvSet.getInterval(currCV).getTraining()), os, currProc);
				if(numCV > 1)
					alg->outputIndEvals(&(cvSet.getInterval(currCV).getTesting()), os, currProc);
				writer.closeStream();
			}
#ifdef PARALLEL
	}
#endif 
		cout << " Completed" << endl;
	alg->finishLog(config.getOutputName(),currCV+1);
		int nModels=1;
		
 #ifdef PARALLEL
		if(myRank==0){

			if(config.outputAllNodesBest())
				nModels = nproc;
#endif
	  // update output when needed
		if(pop.getConvertScores()){
			if(numCV > 1)
					pop.convertScores(&(cvSet.getInterval(currCV).getTraining()), 
						&(cvSet.getInterval(currCV).getTesting()));
			else
				pop.convertScores(&(cvSet.getInterval(0).getTraining()));
		}
		writer.outputSummary(pop, currCV, data, mapFileUsed, config.getOttEncoded(), continMapUsed,
				alg->getFitnessName());    
		switch(config.getSummaryOnly()){
			case Config::False:
				writer.outputGraphic(alg, pop, currCV, config.getOutputName(), nModels, data, 
					mapFileUsed, config.getOttEncoded(), continMapUsed);
			case Config::Best:
				writer.outputBestModels(pop, nModels, currCV,scaler->outputScaleInfo(), data, 
					mapFileUsed, config.getOttEncoded(), continMapUsed); 
			default:
				;
		}
		
#ifdef PARALLEL
} /* end of output */
#endif
		}
		
#ifdef PARALLEL
		}   /* ends master processing of output */
		MPI_Finalize();
 #endif
		delete scaler;

#ifdef PARALLEL
		if(myRank==0){
#endif
		time(&end);
		double dif = difftime (end,start);
		cout << "\n\tAnalysis took " << timeDiff(dif) << endl << endl;
#ifdef PARALLEL
		}
#endif

		return (EXIT_SUCCESS);
}

///
/// Returns string formatted to indicate passage of time
/// @param dif time in seconds
///
std::string timeDiff(double dif){
		double elapsed;
		string period;
		if(dif < 60){
				elapsed = dif;
				period = " seconds";
		}
		else if(dif < 3600){
				elapsed = dif/60;
				period = " minutes";
		}
		else{
				elapsed = dif/3600;
				period = " hours";
		}
		return Stringmanip::numberToString(elapsed) + period;
}

///
/// Outputs message in exception and exits program
/// @param he AthenaExcept
///
void exitApp(AthenaExcept& he, int myRank){
		if(myRank==0)
			cout << he.what() << endl << endl;;
#ifdef PARALLEL
	MPI_Finalize();
#endif
	exit(EXIT_FAILURE);    
} 



///
/// Outputs message in exception and exits program
/// @param he AthenaExcept
///
void exitApp(DataExcept& de, int myRank){
		if(myRank==0)
			cout << de.what() << endl << endl;;
#ifdef PARALLEL
	MPI_Finalize();
#endif
	exit(EXIT_FAILURE);    
} 


///
/// Adjust seed based on current CV and number of processors in run
/// @param config Configuration
/// @param cv Current cross-validation
/// @param nproc Total number of processors
/// @param myRank Rank of this process
///
void adjustSeed(Config& config, int origSeed, int cv, int nproc, int myRank){
	int newSeed = origSeed + nproc * cv + myRank;
	config.setRandSeed(newSeed);
}
