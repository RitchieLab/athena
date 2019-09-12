/*
Copyright Marylyn Ritchie 2015

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

#ifndef _GEDISCRIMBAYES_H
#define	_GEDISCRIMBAYES_H

#include "GEBayes.h"

#ifdef HAVE_CXX_MPI
	#include "GenomeTransfer.h"
#endif

#define MAXMODELSTRING 2048

class GEDiscrimBayes:public Algorithm{

public:

	GEDiscrimBayes();

	~GEDiscrimBayes();

	/// Sets random seed
	void setRand(unsigned int seed);

	void setDummyEncoding(bool ottDummyEncoded);

	/// Set the parameters for the algorithm
	virtual void setParams(AlgorithmParams& algParams, int numExchanges, int numGenos,
		int numContin, vector<unsigned int>& excludedGenos,
		vector<unsigned int>& excludedContins);

	void setDataset(Dataset* newSet);

	/// Starts log
	virtual void startLog(int num_snps);

	/// Prepares log files
	void prepareLog(std::string basename, int cv);

	/// closes logs
	void closeLog();

	/// Runs a step of the algorithm
	virtual int step();

	/// Unused in this algorithm
	virtual	void getAdditionalFinalOutput(Dataset* set){}

	/// Selects best models and does discriminant analysis
	virtual void getAdditionalFinalOutput(Dataset* testing, Dataset* training,
		data_manage::Dataholder* holder, bool mapUsed, bool ottDummy, bool continMapUsed);

	/// Unused in this algorithm
	void testSolution(Dataset* testSet, int nproc){}

	/// For setting any defaults
	void setConfigDefaults(Configuration& configuration, AlgorithmParams& algParam);

	/// Retrieves the models from BioFilter and stores the information in the algorithm
	void getBioModels(std::string fileName, std::string bioFileType,
			data_manage::Dataholder* holder);

	/// Retrieves the models from BioFilter archive and stores in algorithm
	void getBioModelsArchive(string geneGeneFile, string archiveFile,
			data_manage::Dataholder* holder);

	/// Initializes the algorithm for each new dataset to test
	void initialize();

	/// Select best model from models passed and return
	virtual void selectBestModel(std::vector<Solution*>& solutions, data_manage::Dataholder * holder,
		Dataset* set, Configuration& configuration){}

	virtual vector<Solution*> runValidation(std::string sumFile);

		/// Saves the log
	virtual void saveLog(){};

	/// Finish log and pulls together model information
	virtual void finishLog(std::string basename, int cv);

	/// Returns extension for graphical representation
	virtual std::string getGraphicalFileExt();

	/// Run algorithm
	void run();

	/// Outputs individual evaluations to stream
	virtual void outputIndEvals(Dataset* set, ostream& os, int model){}

	/// Writes output to stream
	virtual void writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
		bool mapUsed, bool ottDummy, bool continMapUsed);

	virtual void writeEquation(ostream& os, Solution* sol, data_manage::Dataholder* holder,
		bool mapUsed, bool ottDummy, bool continMapUsed);

	virtual void produceGraphic(std::string inputGraphic, std::string outputGraphic,
		std::string imgWriter){}

	#ifdef HAVE_CXX_MPI
		virtual void setRank(int rank){myRank = rank; caseAlg->setRank(rank);
			controlAlg->setRank(rank);}
		virtual int getRank(){return myRank;}
		virtual void setTotalNodes(int total){
			totalNodes = total;
			caseAlg->setTotalNodes(total);
			controlAlg->setTotalNodes(total);
			}
		virtual int getTotalNodes(){return totalNodes;}
	#endif


protected:

	/// Sets GA for run
	virtual void setGAParams(vector<unsigned int>& excludedGenos,
			vector<unsigned int>& excludedContins);

	struct ConditTable{
		bool genoNode;
		int nodeIndex;
		vector<unsigned int> parentIndexes;
		vector<bool> parentGenos;
		vector<vector<double> > probs;
		vector<int> cumulativeLevels;
	};

	struct modScores{
		int count;
		float score;
		vector<ConditTable> tables;
		std::set<string> varsWithParents;
	};

	struct IndResults{
		string indID;
		vector<long double> scores;
		int phenotype;
		double predicted;
	};

	struct mScores{
		string mString;
		modScores* mPtr;
	};

	static bool sortByScore(const mScores& lhs, const mScores& rhs){
		return lhs.mPtr->score > rhs.mPtr->score;
	}

	/// Counts unique models and stores in modelTotals set
	void totalModels(Algorithm* alg,  std::map<string,modScores>& modelTotals, data_manage::Dataholder* holder,
		bool mapUsed, bool ottDummy, bool continMapUsed, bool caseMods);

  /// Re-write equation string in sorted order so identical equations are matched
  std::string rewriteEquation(string equation);

  /// Sorts parent string from bayesian equation
  string sortParentStr(string parentString);

	// calculate probability tables for each variable when no parent node
	void calcProbTables(Dataset* dset, map<string, vector<double> >& caseOrphanProbs,
		data_manage::Dataholder* holder);

	/// split equation into tokens
	vector<string> tokenizeEquation(string equation);

	void setConditionalTables(Dataset* dset, map<string, modScores>& topModels,
		Dataholder* holder, double missValue);

	void writeUniqueFiles(ostream& outstream, map<float, vector<string> >& sortedModels,
		map<string, modScores>& modelHolder);

	int configParentData(vector<int>& parentValues, vector<unsigned int> &parents,
		vector<bool>& isGeno, Dataset* dSet, vector<int>& cumulativeLevels);

	void setIndModScores(Dataset* dset, map<string,modScores>& models,
		vector<IndResults>& scores, map<std::string, vector<double> >& orphanProbs);

	/// returns AUC
	double setPredictedScores(vector<IndResults>& indScores, map<string,modScores>& caseModels,
		map<string,modScores>& controlModels, double caseRatio);

	/// Sets default values for parameters
	void initializeParams();

	string getLabel(string modelString, Dataholder* holder,
		bool mapUsed, bool continMapUsed);

	void writeIndScores(string filename, vector<IndResults>& scores);

	GEBayes* caseAlg, *controlAlg;
	Dataset * caseDataset, *controlDataset;
	int topModelsUsed, currCV;
	string outputName;

	enum GEDiscrimBayesParams{
		noMatchParam,
		modelsToUse
	};

	std::map<std::string, GEDiscrimBayesParams> paramMap;

	#ifdef HAVE_CXX_MPI
		GenomeTransfer popMigrator;
		struct uniqueModelMPI{
			char modelEquation[MAXMODELSTRING];
			int count;
			float score;
		};

			void gatherModelInformation(map<string, modScores>& models);
	#endif

};


#endif /* _GEDISCRIMBAYES_H */
