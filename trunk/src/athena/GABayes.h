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


#ifndef _GABAYES_H
#define	_GABAYES_H

#include "Algorithm.h"
#include "Variable.h"
#include "GABayesSolution.h"
#include "Athena2DArrayGenome.h"

class GABayes:public Algorithm{

public:

	GABayes();

	~GABayes();

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
	void testSolution(Dataset* testSet, int nproc);

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
	void setConfigDefaults(Configuration& configuration, AlgorithmParams& algParam);

	/// Select best model from models passed and return
	virtual void selectBestModel(std::vector<Solution*>& solutions, data_manage::Dataholder * holder,
		Dataset* set, Configuration& configuration){}

	virtual vector<Solution*> runValidation(std::string sumFile);

	/// Return formatted output for display for final models;
	virtual	void getAdditionalFinalOutput(Dataset* set);

	void getAdditionalFinalOutput(Dataset* testing, Dataset* training,
		data_manage::Dataholder* holder, bool mapUsed, bool ottDummy, bool continMapUsed);

	vector<string> getAdditionalOutputNames();

protected:
	/// Sets GA for run
	virtual void setGAParams(vector<unsigned int>& excludedGenos,
			vector<unsigned int>& excludedContins);

	void addGenotypes(vector<unsigned int>& excludedGenos,
			int nVars);

	void addContin(vector<unsigned int>& excludedContin,
			int nVars);

	void fillPopulation();

	void constructSymbols(Athena2DArrayGenome<int>& genome,
	 	Solution* sol);

	GABayesSolution* convertGenome(GAGenome& ind);

	void initializeParams();

	enum GABayesParams{
		noMatchParam,
		initConnectProb,
		maximumChildren,
		maximumParents,
		limitMethod
	};

// 	enum NodeLimitMethods{
// 		limitRandom,
// 		limitMI
// 	};

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

	string constructBayesStr(vector<vector<int> > network,
		Dataholder* holder,bool genoMapUsed, bool continMapUsed);

	vector<vector<int> > constructEquation(Athena2DArrayGenome<int>& genome);

	void configGA(GASimpleGA* ga);

	void totalModels(GASimpleGA* alg,  map<vector<vector<int> >, ModelScores>& topModels,
		data_manage::Dataholder* holder,bool caseMods,bool genoMapUsed, bool continMapUsed);

std::string timeDiff(double dif);

void writeGenoNet(vector<vector<int> >& eq);

	#ifdef HAVE_CXX_MPI
		void sendAndReceiveGenomes(int totalNodes, int myRank, GASimpleGA* ga);
		void updateWithMigration(int * newGenes, float * recvScores,
			int genomeSize,	int totalNodes, int myRank, GASimpleGA* ga);
	#endif

	std::map<std::string, GABayesParams> paramMap;
// 	std::map<std::string, NodeLimitMethods> limitMethodMap;
	vector<Variable*> varList;
	float initProbConn;
	size_t totalVars;
	string outputName, limitMethodType;
	int currCV, maxChildren, maxParents;
	bool mapfileUsed, continMapfileUsed;
// 	NodeLimitMethods limitMethodType;
};

#endif
