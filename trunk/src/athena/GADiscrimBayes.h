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


#ifndef _GADISCRIMBAYES_H
#define	_GADISCRIMBAYES_H

#include "Algorithm.h"
#include "Variable.h"
#define MAXMODELTRANSFER 16384

class GADiscrimBayes:public Algorithm{

public:

	GADiscrimBayes();

	~GADiscrimBayes();

	/// Set the parameters for the algorithm
	virtual void setParams(AlgorithmParams& algParams, int numExchanges, int numGenos,
		int numContin, vector<unsigned int>& excludedGenos,
		vector<unsigned int>& excludedContins);

	/// add current dataset
	void setDataset(Dataset* newSet);

	/// Initialize GA algorithm
	void initialize();

	/// Runs a step of the algorithm
	virtual int step();

	/// Run algorithm
	void run();

	/// Unused in this algorithm
	void testSolution(Dataset* testSet, int nproc){}

	/// Outputs individual evaluations to stream
	void outputIndEvals(Dataset* set, ostream& os, int model){}

		/// Sets random seed
	void setRand(unsigned int seed);

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

	/// Finish log and pulls together model information
	virtual void finishLog(std::string basename, int cv);

	/// Prepares log files
	void prepareLog(std::string basename, int cv){}

	/// closes logs
	void closeLog(){}

	/// For setting any defaults
	void setConfigDefaults(Config& configuration, AlgorithmParams& algParam);

	/// Select best model from models passed and return
	virtual void selectBestModel(std::vector<Solution*>& solutions, data_manage::Dataholder * holder,
		Dataset* set, Config& configuration){}

	virtual vector<Solution*> runValidation(std::string sumFile);

	void getAdditionalFinalOutput(Dataset* testing, Dataset* training,
		data_manage::Dataholder* holder, bool mapUsed, bool ottDummy, bool continMapUsed);

	void finalFromFile(Dataset* testing, Dataset* training,
		data_manage::Dataholder* holder, bool mapUsed, bool ottDummy, bool continMapUsed);

protected:
	/// Sets GA for run
	virtual void setGAParams(vector<unsigned int>& excludedGenos,
			vector<unsigned int>& excludedContins);

	void addGenotypes(vector<unsigned int>& excludedGenos,
			int nVars);

	void addContin(vector<unsigned int>& excludedContin,
			int nVars);

	void initializeParams();

	enum GABayesParams{
		noMatchParam,
		initConnectProb,
		modelsToUse,
		maximumChildren,
		maximumParents,
		limitMethod,
		caseAllFileName,
		controlAllFileName
	};

	struct ConditionalTable{
		bool genoNode;
		int nodeIndex;
		vector<unsigned int> parentIndexes;
		vector<bool> parentGenos;
		vector<vector<double> > probs;
		vector<int> cumulativeLevels;
	};

	struct ModelScores{
		int count;
		float score;
		vector<ConditionalTable> tables;
		std::set<int> varsWithParents;
		vector<vector<int> > originalConns;
	};

	struct connComparison{
		vector<vector<int> > originalConns, newConns;
		float newScore;
	};

	struct IndivResults{
		string indID;
		vector<long double> scores;
		int phenotype;
		double predicted;
	};

	struct mdScores{
		string mString;
		ModelScores* mPtr;
	};

	static bool sortByScore(const mdScores& lhs, const mdScores& rhs){
		return lhs.mPtr->score > rhs.mPtr->score;
	}

	void writeIndScores(string filename, vector<IndivResults>& scores);

	void writeUniqueFiles(ostream& outstream, map<float,  vector<vector<vector<int> >  > >& sortedModels,
		map<vector<vector<int> >, ModelScores>& modelHolder, Dataholder* holder,
		bool genoMapUsed, bool continMapUsed);

	void writeDotFiles(multimap<float, connComparison>& sortedModels, Dataholder* holder,
		bool genoMapUsed, bool continMapUsed, string endName);

	void writeReducedFile(multimap<float, connComparison>& sortedModels, Dataholder* holder,
		bool genoMapUsed, bool continMapUsed, string fileName);

	string constructBayesStr(vector<vector<int> > network,
		Dataholder* holder,bool genoMapUsed, bool continMapUsed);

	void setConditionalTables(Dataset* dset, map<vector<vector<int> > , ModelScores>& topModels,
		Dataholder* holder, double missValue);

	int configParentData(vector<int>& parentValues, vector<unsigned int> &parents,
		Dataset* dSet, vector<int>& cumulativeLevels);

	void	calcProbTables(Dataset* dset, vector<vector<double> >& orphanProbs,
		data_manage::Dataholder* holder);

	vector<vector<int> > constructEquation(GA2DArrayGenome<int>& genome);

	double setPredictedScores(vector<IndivResults>& indScores, map<vector<vector<int> >,ModelScores>& caseModels,
		map<vector<vector<int> >,ModelScores>& controlModels, double caseRatio);

	void configGA(GASimpleGA* ga);

	void totalModels(GASimpleGA* alg,  map<vector<vector<int> >, ModelScores>& topModels,
		data_manage::Dataholder* holder,bool caseMods,bool genoMapUsed, bool continMapUsed);

	void setIndModScores(Dataset* dset, map<vector<vector<int> >,ModelScores>& models,
		vector<IndivResults>& indScores, vector<vector<double> >& orphanProbs);

	void pruneModels(map<vector<vector<int> >, ModelScores>& caseModels,
		map<vector<vector<int> >, ModelScores>& controlModels, data_manage::Dataholder* holder,
		bool genoMapUsed, bool continMapUsed);

void writeGenoNet(vector<vector<int> >& eq);
	void readAllFile(string allFileName,map<vector<vector<int> >, ModelScores>& models,
		map<string,int>& nameToIndex,data_manage::Dataholder* holder,
		bool caseMods, bool genoMapUsed, bool continMapUsed);

	void selectTopModels(map<vector<vector<int> >, ModelScores>& modelHolder,
		map<vector<vector<int> >, ModelScores>& topModels, data_manage::Dataholder* holder,
		bool caseMods, bool genoMapUsed, bool continMapUsed);

	void runDiscriminantAnalysis(map<vector<vector<int> >, ModelScores>& caseModels,
		map<vector<vector<int> >, ModelScores>& controlModels, Dataset* testing, Dataset* training,
		data_manage::Dataholder* holder, bool mapUsed, bool continMapUsed);

	#ifdef HAVE_CXX_MPI
		void sendAndReceiveGenomes(int totalNodes, int myRank, GASimpleGA* ga);
		void updateWithMigration(unsigned int* newGenes, float * recvScores,
			int genomeSize,	int totalNodes, int myRank, GASimpleGA* ga);
		struct modelMPI{
			short modelSeq[MAXMODELTRANSFER];
			int count;
			float score;
		};
		void constructModelVec(vector<vector<int> >& model, short* modelSeq, int seqSize);
		bool constructTransferSeq(vector<vector<int> > model, short* seq,	int seqSize);
		void gatherModelInformation(map<vector<vector<int> >, ModelScores>& models);
	#endif

	std::map<std::string, GABayesParams> paramMap;

	GASimpleGA * caseGA, *controlGA;
	Dataset * caseDataset, *controlDataset;
	vector<Variable*> varList;
	float initProbConn;
	size_t totalVars;
	string outputName, limitMethodType, caseAllFile, controlAllFile;
	int topModelsUsed, currCV, maxChildren, maxParents;

};

#endif