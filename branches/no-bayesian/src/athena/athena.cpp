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
#include "AlgorithmFactory.h"
#include "OutputManager.h"
#include <iostream>
#include <sstream>
#include "OutputSet.h"
#include <ScaledDataFactory.h>
#include <EncodingFactory.h>
#include <time.h>
#include <unistd.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

void exitApp(AthenaExcept& he, int myRank);
void exitApp(DataExcept& de, int myRank);
void adjustSeed(Config& config, int orig_seed, int cv, int nproc, int myRank);
std::string timeDiff(double dif);
void reportExcluded(vector<unsigned int> genotypes, vector<unsigned int> contins);

int main(int argc, char** argv) {

	int nproc = 1; // only one processor when not running in parallel

 bool mapFileUsed = false, continMapUsed = false;
 int myRank = 0;

#ifdef HAVE_CXX_MPI
	// set up MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif /* end HAVE_CXX_MPI code block */

		string versionDate = "11/12/2015";
		string execName = "ATHENA";
		string version = "1.1.0";
		 time_t start,end;

		if(argc < 2){
				AthenaExcept he("\n\tATHENA\n\tv" + version + "\n\t" + versionDate +
						"\n\n\tUsage: ATHENA <config>\n\n");
				exitApp(he, myRank);
		}
		else{
#ifdef HAVE_CXX_MPI
	if(myRank==0){
#endif
				time (&start);
				cout << endl << "\t" << execName << ":\t" << version << "\t(" <<  versionDate << ")"<< endl << endl;
#ifdef HAVE_CXX_MPI
				}
#endif
		}

		string configFile = argv[1];
		ConfigFileReader configRead;
		Config config;
		ScaleData* scaler = NULL, *continScaler=NULL;
		stringstream performanceStream;

		// read config file
		try{
			config = configRead.readConfig(configFile);
		}catch(AthenaExcept he){
				exitApp(he, myRank);
		}

		vector<AlgorithmParams> algParams= config.getAlgorithmParams();
		Algorithm* alg = AlgorithmFactory::createAlgorithm(algParams[0].name);
		alg->setConfigDefaults(config,algParams[0]);
		// fill dataholder with data
		data_manage::Dataholder data;
		CrossValidator cvMaker;
		CVSet cvSet;

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

		// set random seed  before splitting
		srand(config.getRandSeed());

		// construct crossvalidation sets to use in running algorithm

	if(config.getSplitFile()==""){
			if(config.getValidationSumFile().empty())
		    cvSet = cvMaker.splitData(config.getNumCV(), &data);
		  else
		  	cvSet = cvMaker.splitData(1, &data);
#ifdef HAVE_CXX_MPI
		if(myRank==0){
#endif
			if(config.getValidationSumFile().empty())
	    	cvMaker.saveSplits(config.getOutputName() + ".cvsplit");
#ifdef HAVE_CXX_MPI
		}
#endif
	}
	else{
			cvSet = cvMaker.loadSplits(config.getSplitFile(), &data);
	}

	// check variance of input variables
				data.checkVariance(cvSet);
#ifdef HAVE_CXX_MPI
		if(myRank==0){
#endif
				if(!data.getExcludedGenotypes().empty() || !data.getExcludedContins().empty()){
					reportExcluded(data.getExcludedGenotypes(), data.getExcludedContins());
				}
#ifdef HAVE_CXX_MPI
}
#endif

				// convert data if needed
				try{
					data_manage::DummyConvert * encoder = data_manage::EncodingFactory::createEncoder(config.getEncodeType());
					encoder->convertGenotypes(&data);
					config.setOttEncoded(encoder->getMultipleVars());
					delete encoder;
				}
				catch(DataExcept& de){
					exitApp(de, myRank);
				}

		}catch(AthenaExcept& ae){
				exitApp(ae, myRank);
		}
		catch(DataExcept& de){
			exitApp(de, myRank);
		}

					// alter continuous variables and status value
				scaler = data_manage::ScaledDataFactory::createScaler(config.getStatusAdjust());
				scaler->adjustStatus(&data);
				continScaler = data_manage::ScaledDataFactory::createScaler(config.getContinAdjust());
				continScaler->adjustContin(&data);

		// run crossvalidations and store the populations
		int numCV = cvSet.numIntervals();

		// create algorithm
#ifdef HAVE_CXX_MPI
		alg->setRank(myRank);
		alg->setTotalNodes(nproc);
#endif
		alg->setRand(config.getRandSeed());
		alg->setDummyEncoding(config.getOttEncoded());
		alg->setLogType(config.getLogType());

		try{
			vector<unsigned int> exGenos=data.getExcludedGenotypes();
			vector<unsigned int> exContins=data.getExcludedContins();
			alg->setParams(algParams[0], config.getNumExchanges(),
				data.numGenos(), data.numCovariates(), exGenos, exContins);
		}
		catch(AthenaExcept ex){
			exitApp(ex, myRank);
		}

		/// store results in a population vector
 		vector<Population> pops;
		vector<Solution*> bestSolutions;
		string cvFileName = "cv";

		OutputManager writer;
		writer.setBasename(config.getOutputName());
				// for validation of existing models
		if(!config.getValidationSumFile().empty()){
#ifdef HAVE_CXX_MPI
			if(myRank==0){
#endif
			try{
				alg->setDataset(&(cvSet.getInterval(0).getTraining()));
			}catch(AthenaExcept& ae){
				exitApp(ae, myRank);
			}
			vector<Solution*> models = alg->runValidation(config.getValidationSumFile());
			if(alg->getPopulation().getConvertScores()){
				for(size_t i=0; i<models.size(); i++){
					models[i]->adjustScoreOut(&(cvSet.getInterval(0).getTraining()), alg->getFitnessName());
				}
			}
			writer.writeValidation(alg->getFitnessName(), alg->getAdditionalOutputNames(),
				models, data, mapFileUsed, config.getOttEncoded(), continMapUsed, alg);
			if(config.getIndOutput()){
				stringstream tempss;
				vector<std::stringstream*> indss;
				for(size_t i=0; i<models.size(); i++){
					std::stringstream * newss = new std::stringstream;
					indss.push_back(newss);
				}
				alg->validationIndOutput(indss, models);
				writer.validationIndOutput(indss, config.getOutputName());
				for(size_t i=0; i<indss.size(); i++){
					delete indss[i];
				}
			}
#ifdef HAVE_CXX_MPI
			}
#endif
		}
		else{

		int currCV=config.getStartCV()-1;
	if(currCV==0 and config.getSummaryOnly()!=Config::Suppress){
		writer.setFiles(mapFileUsed, alg->getFitnessName(), alg->getAdditionalOutputNames());
	}

	  int originalSeed = config.getRandSeed();
		for(; currCV < numCV; currCV++){
			adjustSeed(config, originalSeed,currCV, nproc, myRank);
			alg->setRand(config.getRandSeed());
#ifdef HAVE_CXX_MPI
	if(myRank==0){
#endif
			cout << "Beginning Cross-validation " << currCV + 1 << "...";
			cout.flush();
#ifdef HAVE_CXX_MPI
	}
#endif

			try{
				alg->setDataset(&(cvSet.getInterval(currCV).getTraining()));
			}catch(AthenaExcept& ae){
				exitApp(ae, myRank);
			}
#ifdef HAVE_CXX_MPI
		if(myRank==0){  // only have master output cv when desired
#endif
			// output cv interval when requested
			if(config.getCVOutput()){
				OutputSet oSet;
				oSet.outputCV(config.getOutputName() + ".train.cv", cvSet.getInterval(currCV).getTraining(), currCV+1);
				if(numCV > 1)
				oSet.outputCV(config.getOutputName() + ".test.cv", cvSet.getInterval(currCV).getTesting(), currCV+1);
			}
#ifdef HAVE_CXX_MPI
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

		alg->getAdditionalFinalOutput(&(cvSet.getInterval(currCV).getTraining()));

		if(numCV > 1){
			alg->testSolution(&(cvSet.getInterval(currCV).getTesting()), nproc);
			alg->getAdditionalFinalOutput(&(cvSet.getInterval(currCV).getTesting()),
				&(cvSet.getInterval(currCV).getTraining()), &data, mapFileUsed, config.getOttEncoded(),
				continMapUsed);
		}
			// check population values
	  pops.push_back(alg->getPopulation());
	  Population& pop = pops.back();
	  bestSolutions.push_back(pop[0]->clone());
		// update output when needed
		if(myRank ==0 or config.getSummaryOnly()==Config::All){
			if(pop.getConvertScores()){
				if(numCV > 1)
						pop.convertScores(&(cvSet.getInterval(currCV).getTraining()),
							&(cvSet.getInterval(currCV).getTesting()), alg->getFitnessName());
				else
					pop.convertScores(&(cvSet.getInterval(0).getTraining()), alg->getFitnessName());
			}
			if(config.getSummaryOnly()==Config::All){
				writer.outputAllModels(alg, pop, myRank, currCV,scaler->outputScaleInfo(), data,
					mapFileUsed, config.getOttEncoded(), continMapUsed, numCV > 1);
			}
		}

		int currProc = 0;
#ifdef HAVE_CXX_MPI
	if(myRank==0){
		int lastProc = 1;
		if(config.outputAllNodesBest())
			lastProc = nproc;
		for(currProc=0; currProc < lastProc; currProc++){
#endif
			if(config.getIndOutput()){
				performanceStream << "CV " << Stringmanip::numberToString(currCV+1) << endl;
				performanceStream << "Training" << endl;
				alg->outputIndEvals(&(cvSet.getInterval(currCV).getTraining()), performanceStream, currProc);
				if(numCV > 1){
					performanceStream << "Testing" << endl;
					alg->outputIndEvals(&(cvSet.getInterval(currCV).getTesting()), performanceStream, currProc);
				}

			}
#ifdef HAVE_CXX_MPI
	}
#endif
		cout << " Completed" << endl;
	alg->finishLog(config.getOutputName(),currCV+1);
		int nModels=1;



 #ifdef HAVE_CXX_MPI
		if(myRank==0){

			if(config.outputAllNodesBest())
				nModels = nproc;
#endif

		if(config.getSummaryOnly()!=Config::Suppress){
			writer.outputSummary(pop, currCV, data, alg, mapFileUsed, config.getOttEncoded(), continMapUsed,
				 alg->getFitnessName());
			switch(config.getSummaryOnly()){
				case Config::All:
				case Config::False:
					writer.outputGraphic(alg, pop, currCV, config.getOutputName(), nModels, data,
						mapFileUsed, config.getOttEncoded(), continMapUsed, config.getImgWriter());
				case Config::Best:
					writer.outputBestModels(pop, nModels, currCV,scaler->outputScaleInfo(), data,
						mapFileUsed, config.getOttEncoded(), continMapUsed);
				default:
				;
			}
		}

#ifdef HAVE_CXX_MPI
} /* end of output */
#endif
#ifdef HAVE_CXX_MPI
		}   /* ends master processing of output */
#endif
		}

	// add equation output to summary
#ifdef HAVE_CXX_MPI
if(myRank==0){
#endif
	if(config.getSummaryOnly()!=Config::Suppress){
		writer.outputEquations(alg, bestSolutions, data, mapFileUsed, config.getOttEncoded(),
			continMapUsed);
	}
#ifdef HAVE_CXX_MPI
}
#endif

	// if chosen run best model selection from list of best models
	try{
	if(config.selectBestModel()){
		Dataset selectSet;
		if(numCV > 1)
			selectSet = cvSet.getInterval(0).getTraining() + cvSet.getInterval(0).getTesting();
		else
			selectSet = cvSet.getInterval(0).getTraining();
		alg->selectBestModel(bestSolutions, &data, &selectSet, config);
#ifdef HAVE_CXX_MPI
	if(myRank==0){
#endif

		Population bestPop = alg->getPopulation();
		alg->getAdditionalFinalOutput(&selectSet);
		writer.outputBest(bestPop[0],data,mapFileUsed,config.getOttEncoded(),continMapUsed,
			alg->getFitnessName());
			if(config.getIndOutput()){
				performanceStream << "CV Best" << endl;
				alg->outputIndEvals(&selectSet, performanceStream, 0);
			}
#ifdef HAVE_CXX_MPI
	}
#endif
	}
	}
	catch(AthenaExcept ae){
#ifdef HAVE_CXX_MPI
	if(myRank==0)
#endif
		cout << ae.what() << endl;
	}

#ifdef HAVE_CXX_MPI
		if(myRank==0)
#endif
	if(config.getIndOutput())
			writer.outputInds(performanceStream, config.getOutputName(), alg->getFitnessName());
	} // end for standard run (no input validation file)
#ifdef HAVE_CXX_MPI
		MPI_Finalize();
#endif
		delete scaler;
		delete continScaler;

#ifdef HAVE_CXX_MPI
		if(myRank==0){
#endif
		if(config.getSummaryOnly()==Config::All){
			for(int cv=0; cv < numCV; cv++){
				writer.combineAllModels(nproc, cv, alg);
			}
		}
		time(&end);
		double dif = difftime (end,start);
		cout << "\n\tAnalysis took " << timeDiff(dif) << endl << endl;
#ifdef HAVE_CXX_MPI
		}
#endif

	for(size_t i=0; i<bestSolutions.size(); i++){
		delete bestSolutions[i];
	}

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
			cout << "\n" << he.what() << endl << endl;;
#ifdef HAVE_CXX_MPI
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
			cout << "\nERROR: " << de.what() << endl << endl;;
#ifdef HAVE_CXX_MPI
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

///
/// Output any excluded variables
/// @param genotypes
/// @param contins
///
void reportExcluded(vector<unsigned int> genotypes, vector<unsigned int> contins){
	vector<unsigned int>::iterator iter;
	cout << "\nExcluded variables:";
	set<string> excluded;
	for(iter=genotypes.begin(); iter != genotypes.end(); ++iter){
		excluded.insert("G" + Stringmanip::numberToString(*iter+1));
	}
	for(iter=contins.begin(); iter != contins.end(); ++iter){
		excluded.insert("C" + Stringmanip::numberToString(*iter+1));
	}

	for(set<string>::iterator iter=excluded.begin(); iter != excluded.end();
		++iter){
		cout << " " << *iter;
	}
	cout << "\n" << endl;
}

