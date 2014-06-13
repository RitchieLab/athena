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
#include <ga/ga.h>
#include "GENNGrammarAdjuster.h"
#include "InitGEgenome.h"
#include "BioFilterModelCollection.h"

#ifdef PARALLEL
#define MAX_GENOME_SIZE 100000
#endif

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
			int numContin);
		
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
			
		virtual void writeEquation(ostream& os, Solution* sol, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed)=0;
			
		virtual void produceGraphic(std::string inputGraphic, std::string outputGraphic, 
			std::string imgWriter)=0;
		
		/// Returns extension for graphical representation
		virtual std::string getGraphicalFileExt()=0;
		
		/// Saves the log
		virtual void saveLog()=0;
		
		/// Starts log
		virtual void startLog(int num_snps){}

		/// Writes the log
		virtual void writeLog(){}
		
		/// Clears the logs
		virtual void clearLogs(){}
		
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
			data_manage::Dataholder* holder);
		
		/// Retrieves the models from BioFilter archive and stores in algorithm
		virtual void getBioModelsArchive(string geneGeneFile, string archiveFile, 
			data_manage::Dataholder* holder);
		
		virtual std::string getFitnessName(){return fitnessName;}
		
		/// Get column names for output
		virtual vector<std::string> getAdditionalOutputNames(){
			vector<std::string> temp; 
			return temp;
		}
		
		/// Get final values for best model for reporting (such as AUC)
		virtual	void getAdditionalFinalOutput(Dataset* set){}
		
		virtual void setFitnessName(std::string fname){fitnessName = fname;}
		
		virtual void setConfigDefaults(Config& configuration,AlgorithmParams& algParam)=0;
		
		/// Select best model from models passed and return
		virtual void selectBestModel(std::vector<Solution*>& solutions, data_manage::Dataholder * holder,
			Dataset* set, Config& configuration)=0;
			
		/// Runs models from summary file against validation set
		virtual vector<Solution*> runValidation(std::string sumFile)=0;			
		
		#ifdef PARALLEL
			virtual void setRank(int rank){myRank = rank;}
			int getRank(){return myRank;}
			void setTotalNodes(int total){totalNodes = total;}
			int getTotalNodes(){return totalNodes;}
		#endif
		
protected:

		/// Sets default values for parameters
		virtual void initializeParams();
		
		/// Sets GA parameters
		virtual void setGAParams();
		
		virtual void freeMemory();

		void setInitParams();

		void expandVariables();
		
		void setMapperPrefs(AthenaGrammarSI& athenaMapper);
		
		void setBioModels(BioFilterModelCollection& collection, data_manage::Dataholder* holder);

		enum AlgParams{
				minSizeParam,
				maxSizeParam,
				tailRatioParam,
				growRateParam,
				maxDepthParam,
				tailSizeParam,
				sensibleInitParam,
				popSizeParam,
				probCrossParam,
				probMutParam,
				gramFileParam,
				stepSizeParam,
				calcType,
				useEffectiveXO,
				useAllSnps,
				useAllCovariates,
				bioModelSelection,
// 				requireAll,
// 				requireAllOnce,
// 				bioInitFract,
// 				restrictVarGens,
// 				bioModelSelection,
				blockCrossGens,
// 				resetVarsAtMigration,
// 				bpfreq,
// 				bpstart,
				gaSelection,
			#ifdef ATHENA_BLOAT_CONTROL
				doubleTournF,
				doubleTournD,
				doubleTournFitFirst,
				prunePlantFract,
			#endif
				fitGoal,
				bestCVThresh,
				bestCorrThresh,
// 				constantSpan
		};
		
		enum GASelectionType{
			NoMatchSelector,
			#ifdef ATHENA_BLOAT_CONTROL
			DoubleTournamentSelection,
			#endif
			RouletteWheelSelection,
			ParetoFrontSelection,
			ParetoRankSelection
		};

		enum BioSelectionType{
			NoMatchSelection,
			rouletteSelect,
			orderedSelect
		};

		std::map<std::string, AlgParams> algParamMap;
		std::map<std::string, GASelectionType> gaSelectorMap;
		std::map<std::string, BioSelectionType> bioModelSelectionMap;
		GASelectionType gaSelector;

		// BioFilter parameters
		BioSelectionType biofilterSelectorType;
		
		bool effectiveXO, useAllVars, useAllCovars;
		unsigned int wrapEvents, randSeed;
		unsigned int minSize, maxSize;
		float tailRatio, growRate;
		unsigned int maxDepth, tailSize;
		bool sensibleInit, maxBest;
		float pruneAndPlantFract;
		unsigned int popSize, numGenerations, stepSize, ngensBlockCross;
		
		// Genetic algorithm
		GASimpleGA* ga;
		
		std::string grammarFile, calculatorName;
		double probCross, probMut, initBioFract, fitnessGoal, bestCorrThreshold;
		int numGenotypes, numContinuous, bestCVThreshold;		

		bool dummyEncoded;
		data_manage::Dataset* set;
		LogType logTypeSelected;
		Population pop;
		std::string fitnessName;
		
		vector<AlgorithmLog*> logs;
		int myRank;
		AthenaGrammarSI mapper;
		GENNGrammarAdjuster adjuster;
		
		int totalNodes;
};


#endif	/* _ALGORITHM_H */

