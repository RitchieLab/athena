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
#include "AlgorithmFactory.h"
#include "NNLog.h"
#include "NNModelLog.h"

#include "Config.h"
#ifdef HAVE_CXX_MPI
	#include "GenomeTransfer.h"
#endif

class GENNAlg:public Algorithm{

public:

		GENNAlg();

		~GENNAlg();

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

		/// Sets testing values for best solutions
		void testSolution(Dataset* testSet, int nproc);

		/// Outputs individual evaluations to stream
		void outputIndEvals(Dataset* set, ostream& os, int model);

		/// Initializes the algorithm for each new dataset to test
		void initialize();

		/// Writes output to stream
		void writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed);\

		/// Produce graphical file
		virtual void produceGraphic(std::string inputGraphic, std::string outputGraphic,
			std::string imgWriter);

		/// Writes model as equation to stream
		void writeEquation(ostream& os, Solution* sol, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed);

		/// Returns extension for graphical representation
		virtual std::string getGraphicalFileExt();

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

		/// Return additional output names used for final models
		virtual vector<std::string>  getAdditionalOutputNames();

		/// Return formatted output for display for final models;
		virtual	void getAdditionalFinalOutput(Dataset* set);

		/// Return single best model
		virtual void selectBestModel(std::vector<Solution*>& solutions, data_manage::Dataholder * holder,
			Dataset* set, Config& config);

		virtual void setConfigDefaults(Config& configuration, AlgorithmParams& algParam);

		virtual vector<Solution*> runValidation(std::string sumFile);

		virtual void validationIndOutput(vector<std::stringstream*>& indss, vector<Solution*>& models);

 		#ifdef HAVE_CXX_MPI
 			virtual void setRank(int rank);
 		#endif


protected:

		/// Sets GA for run
		virtual void setGAParams(vector<unsigned int>& excludedGenos,
			vector<unsigned int>& excludedContins);

		/// Sets default values for parameters
		void initializeParams();

		void fillPopulation();

		void fillLog();

		void outputGenome(GAGenome& g);

		void convertNetworks(AthenaGrammarSI& currentMapper, AthenaGrammarSI& newMapper);

		void resetCrossover();

		void setRestrictedGrammar(bool clearVariables);

		void setRestrictedGrammar(bool clearVariables, vector<int>& genos,
		  vector<int>& contins);

		void runBackPropagation();

		NNSolution* convertGenome(GAGenome& ind);

		int reachedGoal();

		enum GENNParams{
				noMatchParam,
				calcType,
				requireAll,
				requireAllOnce,
				bioInitFract,
				restrictVarGens,
				resetVarsAtMigration,
				bpfreq,
				bpstart,
				constantSpan
		};

		enum GASelectionType{
			NoMatchSelector,
			#ifdef ATHENA_BLOAT_CONTROL
			DoubleTournamentSelection,
			#endif
			RouletteWheelSelection
		};

		void outputAlgInds();

		std::string oName;

		std::map<std::string, GENNParams> paramMap;

		void makeRestrictedGrammar(bool clearVariables, AthenaGrammarSI& useMapper);

		//General algorithm parameters
		bool requireAllVars, requireAllVarsOnce, resetRestrictedAtMigration, simpleConstants;
		std::string mainLogFilename, fitnessLogFilename, snpnameLogFilename;
		unsigned int ngensVarRestrict, restrictStepsDone;
		double initBioFract;
		float minConstant, maxConstant, constantInterval;
		int numGenotypes, numContinuous, bpFirstGen, bpFreqGen, bpNextOpt, bestCVThreshold;

		NNLog* geLog;
		NNModelLog* modelLog;

		// GE parameters
		AthenaGrammarSI restrictMapper;

		#ifdef HAVE_CXX_MPI
			GenomeTransfer popMigrator;
		#endif

};


#endif	/* _GENNALG_H */

