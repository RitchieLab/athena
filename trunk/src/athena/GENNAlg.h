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
 * File:   GENNAlg.h
 * Author: dudeksm
 *
 * Created on November 10, 2008, 2:03 PM
 */

#ifndef _GENNALG_H
#define	_GENNALG_H

#include "Algorithm.h"
#include "InitGEgenome.h"
#include "AlgorithmFactory.h"
#include "NNLog.h"
#include "NNModelLog.h"
#include "GENNGrammarAdjuster.h"
#include "BioFilterModelCollection.h"
#include "Config.h"

#ifdef PARALLEL
#define MAX_GENOME_SIZE 100000
#endif


class GENNAlg:public Algorithm{
		
public:

		GENNAlg();
		
		~GENNAlg();
		
		/// Set the parameters for the algorithm
		virtual void setParams(AlgorithmParams& algParams, int numExchanges, int numGenos, int numContin);
		
		/// Run algorithm 
		void run();
		
		/// Runs a step of the algorithm
		virtual int step();
		
		/// Sets random seed
		void setRand(unsigned int seed);
		
		void setDataset(Dataset* newSet);
		
		/// Sets testing values for best solutions
		void testSolution(Dataset* testSet, int nproc);
		
		/// Outputs individual evaluations to stream
		void outputIndEvals(Dataset* set, ostream& os, int model);
		
		/// Initializes the algorithm for each new dataset to test
		void initialize();
		
		/// Writes output to stream
		void writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed);\
			
		/// Writes model as equation to stream
		void writeEquation(ostream& os, Solution* sol, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed);
		
		/// Returns extension for graphical representation
		std::string getGraphicalFileExt();
		
		/// Saves the log
		void saveLog();
		
		/// Starts log
		void startLog(int num_snps);
	 
		/// Writes log
		void writeLog();
		
		/// Clears log
		void clearLogs();
		
		/// Close log
		void closeLog();
		
		/// Finish log and pulls together model information
		void finishLog(std::string basename, int cv);
		
		/// Prepares log files
		void prepareLog(std::string basename, int cv);
		
		/// Returns covariates and snps in best network 
		vector<string> getBestVariables();
		
		/// Retrieves the models from BioFilter and stores the information in the algorithm
		void getBioModels(std::string filename, std::string bioFileType, data_manage::Dataholder* holder);
		
		/// Retrieves the models from the BioFilter archive files and stores information in algorithm
		void getBioModelsArchive(string genegeneFile, string archiveFile, data_manage::Dataholder* holder);
		
		/// Return additional output names used for final models
		virtual vector<std::string> getAdditionalOutputNames();
		
		/// Return formatted output for display for final models;
		virtual	vector<std::string> getAdditionalFinalOutput(Dataset* set);
		
		/// Return single best model
		virtual void selectBestModel(std::vector<Solution*>& solutions, data_manage::Dataholder * holder,
			Dataset* set, Config& config);
			
		virtual void setConfigDefaults(Config& configuration, AlgorithmParams& algParam);
		
		#ifdef PARALLEL
			struct genomeMPI{
				float genomeParams[7];
				int codons[MAX_GENOME_SIZE];
			};
		
			typedef struct genomeMPI structMPI;
			
			virtual void setRank(int rank);
			void updateWithMigration(float* stats, int* codons, int totalNodes, int myRank, int max_length=0);
			void updateWithMigration(structMPI * genomes, int totalNodes, int myRank);
			void sendAndReceive(int totalNodes, int myRank);
			void sendAndReceiveStruct(int totalNodes, int myRank);
			int nodesCompleted(int complete);
			void exchangeBestVariables(int totalNodes, int myRank, vector<int>& genotypes,
				vector<int>& contins);
		#endif
	
		
protected:
		
		/// Sets GA for run
		void setGAParams();
		
		/// Sets default values for parameters
		void initializeParams();   
		
		void fillPopulation();
		
		void fillLog();
		
		void outputGenome(GAGenome& g);
		
		void convertNetworks(AthenaGrammarSI& currentMapper, AthenaGrammarSI& newMapper);
		
		void setMapperPrefs(AthenaGrammarSI& athenaMapper);
		
		void resetCrossover();
		
		void setRestrictedGrammar(bool clearVariables);
		
		void setRestrictedGrammar(bool clearVariables, vector<int>& genos,
		  vector<int>& contins);
		
		void runBackPropagation();
		
		void setBioModels(BioFilterModelCollection& collection, data_manage::Dataholder* holder);
		
		NNSolution* convertGenome(GAGenome& ind);
		
		int reachedGoal();
		
		enum GENNParams{
				noMatchParam,
				minSizeParam,
				maxSizeParam,
				tailRatioParam,
				growRateParam,
				maxDepthParam,
				tailSizeParam,
				sensibleInitParam,
				popSizeParam,
				numGenParam,
				probCrossParam,
				probMutParam,
				gramFileParam,
				stepSizeParam,
				calcType,
				useEffectiveXO,
				useAllSnps,
				useAllCovariates,
				requireAll,
				requireAllOnce,
				bioInitFract,
				restrictVarGens,
				bioModelSelection,
				blockCrossGens,
				resetVarsAtMigration,
				bpfreq,
				bpstart,
				gaSelection,
				doubleTournF,
				doubleTournD,
				doubleTournFitFirst,
				prunePlantFract,
				fitGoal,
				bestCVThresh,
				bestCorrThresh
		};
		
		
		enum BioSelectionType{
			NoMatchSelection,
			rouletteSelect,
			orderedSelect
		};
		
		enum GASelectionType{
			NoMatchSelector,
			#ifdef ATHENA_BLOAT_CONTROL
			DoubleTournamentSelection,
			#endif
			RouletteWheelSelection
		};
	
		
		void freeMemory();
		void outputAlgInds();
		
		bool parametersSet;
		std::string oName;
		
		std::map<std::string, GENNParams> paramMap;
		std::map<std::string, BioSelectionType> bioModelSelectionMap;
		std::map<std::string, GASelectionType> gaSelectorMap;
		
		void makeRestrictedGrammar(bool clearVariables, AthenaGrammarSI& useMapper);
		
		// BioFilter parameters
		BioSelectionType biofilterSelectorType;
		
		// GA Selection parameters
		GASelectionType gaSelector;
		bool fitFirst;
		int doubleTourneyF;
		float doubleTourneyD;
		
		//Random initialization parameters
		unsigned int minSize, maxSize;
		
		//Sensible initialization parameters
		float tailRatio, growRate;
		unsigned int maxDepth, tailSize;
		bool sensibleInit;
		
		//Prune and plant parameters
		float pruneAndPlantFract;
		
		//General algorithm parameters
		bool effectiveXO, useAllVars, useAllCovars, requireAllVars, requireAllVarsOnce, maxBest,
			resetRestrictedAtMigration;
		unsigned int wrapEvents, randSeed;
		std::string grammarFile, calculatorName, mainLogFilename, fitnessLogFilename, 
				snpnameLogFilename;
		unsigned int popSize, numGenerations, stepSize, ngensVarRestrict, restrictStepsDone,
			ngensBlockCross;
		double probCross, probMut, initBioFract, fitnessGoal, bestCorrThreshold;
		int numGenotypes, numContinuous, bpFirstGen, bpFreqGen, bpNextOpt, bestCVThreshold;
		
		NNLog* geLog;
		NNModelLog* modelLog;
		
		// GE parameters
		AthenaGrammarSI mapper, restrictMapper;
		
		// Genetic algorithm
		GASimpleGA* ga;
		GENNGrammarAdjuster adjuster;
		
		#ifdef PARALLEL
			int genomeInfo, genomeArray;
		#endif
		
};


#endif	/* _GENNALG_H */

