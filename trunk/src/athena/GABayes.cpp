#include "GABayes.h"
#include "GAFunct.h"
#include <set>
#include "ModelLogParser.h"
#include "SumFileReader.h"
#include <algorithm>
#include <unistd.h>
using namespace std;

///
/// Constructor
///
GABayes::GABayes(){
	initializeParams();
}

///
/// Destructor
///
GABayes::~GABayes(){
		freeMemory();
		for(size_t i=0; i<varList.size(); i++){
			delete varList[i];
		}
}


///
/// Initializes parameters for basic run of algorithm.  These
/// can be modified using the set_params function
///
void GABayes::initializeParams(){

	 myRank = 0;
	 modelSortCol = 2;

	 popSize = 100;
	 numGenerations = 100;
	 probCross = 0.9;
	 probMut = 0.01;
	 fitnessGoal = 1.0;
	 bestCorrThreshold = 0.8;
	 bestCVThreshold = 2;

	 stepSize = 100;
	 randSeed = 7;

	 initProbConn = 0.5;

	 calculatorName = "K2";
	 currCV=1;

	 // establish map for parameters
	 paramMap["INITCONNECTP"] = initConnectProb;

	 gaSelectorMap["ROULETTE"] = RouletteWheelSelection;
	 gaSelectorMap["PARETO"] = ParetoFrontSelection;
	 gaSelectorMap["PARETORANK"] = ParetoRankSelection;

	 useAllVars = false;
	 numGenotypes = 0;
	 numContinuous = 0;

	 ga = NULL;
	 maxBest = true;
	 gaSelector = RouletteWheelSelection;
// 	 BayesSolution* sol = (BayesSolution*)GAFunct::getBlankSolution();
// 	 pop.insert(sol);
}


///
/// Sets parameters within the algorithm for the
/// analysis.
/// @param algParam AlgorithmParams
/// @param numExchanges total number of times the best model will be passed to other algorithms
/// when more than one algorithm
/// @param numGenos number of Genotypes in set
/// @param numContin number of Continuous variables in set
///
void GABayes::setParams(AlgorithmParams& algParam, int numExchanges, int numGenos,
	int numContin, vector<unsigned int>& excludedGenos, vector<unsigned int>& excludedContins){


		Algorithm::setParams(algParam, numExchanges, numGenos, numContin, excludedGenos,
			excludedContins);
		map<string, string>::iterator mapIter;
		vector<string> tokens;

		for(mapIter = algParam.params.begin(); mapIter != algParam.params.end();
			mapIter++){
				switch(paramMap[mapIter->first]){
						case noMatchParam:
							if(paramMap.find(mapIter->first) == paramMap.end())
								throw AthenaExcept("No match for parameter " + mapIter->first +
												" in Algorithm GADescrimBayes");
								break;
// 							case modelsToUse:
// 								topModelsUsed = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
// 								break;
// 						case parentRange:
// 							// split on ':' for min:max
// 							tokens = Stringmanip::split(mapIter->second, ':');
// 							minNumParents = Stringmanip::stringToNumber<unsigned int>(tokens[0]);
// 							maxNumParents = Stringmanip::stringToNumber<unsigned int>(tokens[1]);
// 							if(minNumParents < 1)
// 								throw AthenaExcept("PARENTRANGE minimum must be 1 or greater");
// 							break;
// 						case childRange:
// 							tokens = Stringmanip::split(mapIter->second, ':');
// 							minNumChildren = Stringmanip::stringToNumber<unsigned int>(tokens[0]);
// 							maxNumChildren = Stringmanip::stringToNumber<unsigned int>(tokens[1]);
// 							if(minNumChildren < 1)
// 								throw AthenaExcept("CHILDRANGE minimum must be 1 or greater");
// 							break;
						case initConnectProb:
							initProbConn = Stringmanip::stringToNumber<float>(mapIter->second);
							break;
// 						case bastart:
// 							balAccStart = Stringmanip::stringToNumber<int>(mapIter->second);
// 							break;
						default:
							if(paramMap.find(mapIter->first) == paramMap.end())
								throw AthenaExcept("No match for parameter " + mapIter->first +
												" in Algorithm GADescrimBayes");
				}
		}

		// first optimization of backpropagation
// 		baNextOpt = balAccStart;
		setGAParams(excludedGenos, excludedContins);
}


///
/// Sets parameters for use with GAlib
/// @param excludedGenos Vector of excluded SNP variables
/// @param excludedContins List of excluded continuous variables
/// @throws AthenaExcept on error
///
void GABayes::setGAParams(vector<unsigned int>& excludedGenos,
			vector<unsigned int>& excludedContins){
		GARandomSeed(randSeed);
		srand(randSeed);
		// first free ga memory if already run
		freeMemory();
		varList.clear();

		addGenotypes(excludedGenos, numGenotypes);
		addContin(excludedContins, numContinuous);
		GAFunct::setSolutionType(calculatorName);

	 if(calculatorName.find("K2") != string::npos){
		 pop.setConvertScores(false);
	 }

}


///
/// Set list of included variables
/// @param excludedGenos Vector of excluded variables
/// @param nVars Total number of variables
///
void GABayes::addGenotypes(vector<unsigned int>& excludedGenos,
			int nVars){

			std::set<int> excluded;
			for(int i=0; i<excludedGenos.size(); i++){
				excluded.insert(excludedGenos[i]);
			}

			for(int i=0; i<nVars; i++){
				if(excluded.find(i) == excluded.end()){
					varList.push_back(new GenoVariable(i));
				}
			}
}

///
/// Set list of included variables
/// @param excludedGenos Vector of excluded variables
/// @param nVars Total number of variables
///
void GABayes::addContin(vector<unsigned int>& excludedContin,
			int nVars){

			std::set<int> excluded;
			for(int i=0; i<excludedContin.size(); i++){
				excluded.insert(excludedContin[i]);
			}

			for(int i=0; i<nVars; i++){
				if(excluded.find(i) == excluded.end()){
					varList.push_back(new ConVariable(i));
// 					varList.back()->setNumLevels(controlDataset->getNumLevels(i));
				}
			}
}



///
/// Split Dataset into case/control sets
/// @param newSet Dataset to split
///
void GABayes::setDataset(Dataset* newSet){
	set = newSet;
	for(size_t i=0; i<varList.size(); i++){
		if(!varList[i]->isGeno()){
			varList[i]->setNumLevels(set->getNumLevels(varList[i]->getIndex()));
		}
	}
	GAFunct::setDataset(set, varList);
}

///
/// Establishes test data set values for the population
/// @param test_set Dataset containing individuals for the test set
///
void GABayes::testSolution(Dataset* testSet, int nproc){

	 // make sure test set is categorical
// 	 ScaleCategorical scaler;
// 	 scaler.adjustContin(testSet);
// 	 scaler.adjustStatus(testSet);

	// use first dataset
	set = testSet;
	GAFunct::setDataset(set,varList);
  int n = ga->population().size();

  for(int i=0; i<n; i++){
		GA2DBinaryStringGenome bestGenome = (GA2DBinaryStringGenome&)ga->population().best(i);
		float testScore = GAFunct::GACaseObjective(bestGenome);
// 		((GA2DBinaryStringGenome&)ga->population().best(i)).setTestValue(testScore);
		pop[i]->testVal(testScore);
	}
}


///
/// Set up GA
///
void GABayes::configGA(GASimpleGA* ga){
		GASelectionScheme* selector;
		switch(gaSelector){
			case NoMatchSelector:
			case DoubleTournamentSelection:
			case ParetoFrontSelection:
			case ParetoRankSelection:
			case RouletteWheelSelection:
				selector = new GARouletteWheelSelector;
				break;
		};

		ga->selector(*selector);
		// safe to delete selector as algorithm clones it
		delete selector;
		ga->scaling(GANoScaling());
		// individuals in population
		ga->populationSize(popSize);
		ga->pMutation(probMut);
		ga->pCrossover(probCross);
		ga->nGenerations(numGenerations);
		// always maximize these scores
			ga->maximize();
			maxBest = true;
		ga->initialize();
}


///
/// Initializes the algorithm
///
void GABayes::initialize(){
		// free ga if already established
		freeMemory();

		// total number of variables
		totalVars = varList.size();
		GA2DBinaryStringGenome genome(totalVars,totalVars, GAFunct::GACaseObjective);
		// evaluator set as part of constructor
		GAFunct::setInitConP(initProbConn);
		genome.initializer(GAFunct::initCase);
		ga = new GASimpleGA(genome);

		configGA(ga);
}

///
/// Establishes the algorithm for the run based on the parameters
/// set in this algorithm.
///
void GABayes::run(){
		while(!ga->done()){
				ga->step();
		}
}

///
/// Runs an indicated interval for the algorithm. This interval is set
/// in the stepSize variable.
///
int GABayes::step(){
		int completed =0;
		// run each separately
	for(unsigned int i=0; i < stepSize; i++){
		ga->step();
	}

#ifdef HAVE_CXX_MPI
	// perform migration
	sendAndReceiveGenomes(totalNodes,myRank,ga);
#endif

		fillPopulation();
		return completed;
}

///
/// Fills the population after a step to allow for exchange of models with
/// other algorithms.
///
void GABayes::fillPopulation(){
		unsigned int numInds = ga->population().size();
		// clear solution population
		pop.clear();


// GAGenome& ind=ga->population().individual(0);
// cout << "best score=" << ind.score() << endl;
// GA2DBinaryStringGenome& genome = (GA2DBinaryStringGenome&) ind;
// for(int i=0; i<genome.height(); i++){
// 	for(int j=0; j<genome.width(); j++){
// 	cout << " " << genome.gene(i,j);
// 	}
// 	cout << endl;
// }
// cout << "----------------" << endl;

		for(unsigned int currInd = 0; currInd < numInds; currInd++){
				pop.insert(convertGenome(ga->population().individual(currInd)));
		}
}

///
/// Converts an individual from population into a Solution and returns it
/// @param GAGenome ind
/// @return BayesSolution*
///
GABayesSolution* GABayes::convertGenome(GAGenome& ind){
	GA2DBinaryStringGenome& genome = (GA2DBinaryStringGenome&) ind;
	GABayesSolution* sol = (GABayesSolution*)GAFunct::getBlankSolution();
// 	vector<vector<int> > = constructEquation(genome);
	constructSymbols(genome, sol);
	sol->fitness(genome.score());
// 	sol->setSymbols(symbols);
// 	sol->fitness(genome.score());
// 	sol->testVal(genome.getTestValue());
// 	sol->setGramDepth(genome.getGramDepth());
// 	sol->setNNDepth(genome.getDepth());
// 	sol->setComplexity(genome.getComplexity());
	return sol;
}


/// write out equation
void GABayes::writeGenoNet(vector<vector<int> >& eq){
	for(size_t i=0; i<eq.size(); i++){
		cout << "[G" << i+1;
		if(!eq[i].empty()){
			cout << "|G" << (eq[i][0]+1);
		}
		for(size_t j=1; j<eq[i].size(); j++){
			cout << ":G" << (eq[i][j]+1);
		}
		cout << "]";
	}
	cout << "\n";
}

///
/// Sets random seed
/// @param seed
///
void GABayes::setRand(unsigned int seed){
	GARandomSeed(seed);
	srand(seed);
}


/// Writes graphical representation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GABayes::writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
	bool mapUsed, bool ottDummy, bool continMapUsed){
// 	os << "implement writeGraphical" << endl;
	os << "digraph G{\n";
	os << "\tgraph [ dpi = 300 ];\n";
	os << "\tsize=\"8.5,11.0\";\n";
	os << "\tdir=\"none\";\n";
	os << "\trankdir=\"LR\";\n";
	os << "\torientation=\"portrait\";\n";

	std::set<std::string> added;

	vector<string> symbols=sol->getSymbols();
	for(vector<string>::iterator symbIter=symbols.begin(); symbIter != symbols.end();
		++symbIter){

		// remove first and last ('[' and ']')
		string symb=(*symbIter).substr(1,(*symbIter).length()-2);
		// split on "|"
		vector<string> pcs=data_manage::Stringmanip::split(symb, '|');
		// check for parents
		if(pcs.size() > 1){
			vector<string> parents=data_manage::Stringmanip::split(pcs[1],':');
			for(size_t i=0; i<parents.size(); i++){
				if(added.find(parents[i])==added.end()){
					os << "\t" << parents[i] << " [shape=\"box\" style=\"filled\" label=\"" <<
						parents[i] << "\"];\n";
					added.insert(parents[i]);
				}
				os << "\t" << parents[i] <<  "->" << pcs[0] << ";\n";
			}
		}
		if(added.find(pcs[0])==added.end()){
			os << "\t" << pcs[0] << " [shape=\"box\" style=\"filled\" label=\"" <<
				pcs[0] << "\"];\n";
			added.insert(pcs[0]);
		}
	}
	os << "}" << endl;
}


///
/// Writes model represented as an equation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GABayes::writeEquation(ostream& os, Solution* sol, data_manage::Dataholder* holder,
	bool mapUsed, bool ottDummy, bool continMapUsed){
	sol->outputClean(os,*holder,mapUsed,ottDummy,continMapUsed);
}

///
/// Returns graphical file extension appropriate to solution
/// @return extension
///
std::string GABayes::getGraphicalFileExt(){
	return ".dot";
// 	return GEObjective::getGraphicalExt();
}

///
/// Finishes logs and compiles them into a single file
/// in cases where there are multiple processors
/// @param basename
/// @param cv
///
void GABayes::finishLog(string basename, int cv){
		if(myRank==0){
				if(logTypeSelected==LogNone){
					return;
				}
				vector<string> filenames;
				for(int rank=0; rank < totalNodes; rank++){
						string modelLogName = basename + ".rank." + Stringmanip::numberToString<int>(rank) + ".cv." +
								Stringmanip::numberToString<int>(cv) + ".models.log";
						filenames.push_back(modelLogName);
				}
				string outName = basename + ".cv." +
								Stringmanip::numberToString<int>(cv) + ".models.log";
				ModelLogParser parser;
				if(logTypeSelected == LogVariables)
					parser.compileVariableFiles(filenames, outName);
				else if(logTypeSelected != LogOverview)
					parser.compileFiles(filenames, outName, GAFunct::getWorstScore());
		}
}


///
/// Run dataset against models contained in previously generated summary file
/// @param sumFile ATHENA summary file
///
vector<Solution*> GABayes::runValidation(std::string sumFile){
	SumFileReader reader;
	reader.readSumFile(sumFile);
	vector<Solution*> models = reader.getModelPopulation();
// 	for(size_t i=0; i<models.size(); i++){
// 		GEObjective::calcFitness(models[i]);
// 		models[i]->setAdditionalOutput(GEObjective::calcAdditionalFinalOutput(models[i]));
// 		models[i]->testVal(0);
// 	}
	return models;
}


///
/// Sets values in main configuration to defaults needed by
/// GADiscrimBayes Algorithm
/// @param configuration Config
///
void GABayes::setConfigDefaults(Config& configuration, AlgorithmParams& algParam){
	// algorithm will write its own output
	configuration.setLogType("NONE");
// 	configuration.setSummaryOnly("SUPPRESS");
	configuration.setContinAdjust("MAKECATEGORIAL");
	configuration.setStatusAdjust("MAKECATEGORIAL");
	outputName = configuration.getOutputName();
	mapfileUsed = (configuration.getMapName().length() > 1);
	continMapfileUsed = (configuration.getContinMapName().length() > 1);
}

///
/// Return additional output names (if any such as AUC)
///
vector<string> GABayes::getAdditionalOutputNames(){
  return GAFunct::getAdditionalOutputNames();
}

///
/// Performs discriminant analysis
///
void GABayes::getAdditionalFinalOutput(Dataset* testing, Dataset* training,
	data_manage::Dataholder* holder, bool mapUsed, bool ottDummy, bool continMapUsed){

	GABayes::setDataset(testing);
	unsigned int numInds = ga->population().size();
  for(unsigned int currInd = 0; currInd < numInds; currInd++){
		pop[currInd]->setAdditionalOutput(GAFunct::getAdditionalFinalOutput(ga->population().individual(currInd)));
	}

}

///
/// Calculate and return additional output (whether model is better than unconnected
///  bayesian network) for best model
///
void GABayes::getAdditionalFinalOutput(Dataset* set){
	GABayes::setDataset(set);

	unsigned int numInds = ga->population().size();
  for(unsigned int currInd = 0; currInd < numInds; currInd++){
		pop[currInd]->setAdditionalOutput(GAFunct::getAdditionalFinalOutput(ga->population().individual(currInd)));
	}
}


///
/// Construct equation string from genome
/// @genome GA2DBinaryStringGenome
///
vector<vector<int> > GABayes::constructEquation(GA2DBinaryStringGenome& genome){
	vector<int> empty;
	vector<vector<int> > conns(varList.size(), empty);
	for(int i=0; i<genome.height(); i++){
		for(int j=0; j<genome.width(); j++){
			if(genome.gene(i,j)==1){
				conns[j].push_back(i);
			}
		}
	}
	return conns;
}

///
/// Constructs bnnlearn compatible symbols only for nodes thare
/// have either parents or children (or both)
/// @param genome
/// @param solution
/// @returns vector of symbols
void GABayes::constructSymbols(GA2DBinaryStringGenome& genome,
	Solution* sol){
  vector<string> symbols;
  std::set<int> genos, contins;
  int height=genome.height();
  int width=genome.width();

	string genoPrefix ="";
	string continPrefix ="";
  if(!mapfileUsed){
  	genoPrefix="G";
  }
  if(!continMapfileUsed){
  	continPrefix="C";
  }

  Dataholder* holder=set->getHolder();
  for(int ch=0; ch<width; ch++){
		vector<int> parents;
		bool include=false;
		for(int p=0; p<height; p++){
			if(genome.gene(p,ch)){
				parents.push_back(p);
			}
		}
		// if parents empty, check to see if this node is a parent itself
		if(parents.empty()){
			for(int i=0; i<width; i++){
				if(genome.gene(ch,i)){
					include=true;
					break;
				}
			}
		}
		else{
			include=true;
		}
		if(include){
			//construct symbol
			string prefix;

			if(varList[ch]->isGeno())
				prefix=genoPrefix;
			else
				prefix=continPrefix;
			string symb="[" + prefix + varList[ch]->getName(holder);
			if(!parents.empty()){
				if(varList[parents[0]]->isGeno())
					prefix=genoPrefix;
				else
					prefix=continPrefix;
				symb += "|" + prefix + varList[parents[0]]->getName(holder);
				for(size_t par=1; par<parents.size(); par++){
					if(varList[parents[par]]->isGeno())
						prefix=genoPrefix;
					else
						prefix=continPrefix;
					symb += ":" + varList[parents[par]]->getName(holder);
				}
			}
			symb += "]";
			symbols.push_back(symb);
			if(varList[ch]->isGeno())
				genos.insert(varList[ch]->getIndex()+1);
			else
				contins.insert(varList[ch]->getIndex()+1);
		}
  }

  if(symbols.empty()){
  	symbols.push_back(" No connections");
  }

	sol->setSymbols(symbols);
	vector<int> g,c;

	// set genotypes and continuous variable indexes
	std::copy(genos.begin(), genos.end(), std::back_inserter(g));
	std::copy(contins.begin(), contins.end(), std::back_inserter(c));
	sol->setGenotypes(g);
	sol->setCovariates(c);
	sol->setComplexity(g.size()+c.size());

}

///
/// Constructs bayes network representation compatible with bnlearn package in R
/// @param network 2-D vector showing parents for each node
/// @param holder Dataholder
/// @returns network string
///
string GABayes::constructBayesStr(vector<vector<int> > network,
	Dataholder* holder, bool genoMapUsed, bool continMapUsed){
	string netString;
	vector<string> prefixes(2, "");
	if(!genoMapUsed){
		prefixes[1]="G";
	}
	if(!continMapUsed){
		prefixes[0]="C";
	}
	for(size_t i=0; i<network.size(); i++){
		netString += "[";
		netString += prefixes[varList[i]->isGeno()] + varList[i]->getName(holder);
		if(!network[i].empty()){
			netString += "|" + prefixes[varList[network[i][0]]->isGeno()]  + varList[network[i][0]]->getName(holder);
			for(size_t j=1; j<network[i].size(); j++){
				netString += ":" + prefixes[varList[network[i][j]]->isGeno()] + varList[network[i][j]]->getName(holder);
			}
		}
		netString += "]";
	}
	return netString;
}


///
///  Writes list of all unique models and counts
///
void GABayes::writeUniqueFiles(ostream& outstream, map<float,  vector<vector<vector<int> >  > >& sortedModels,
	map<vector<vector<int> >, ModelScores>& modelHolder, Dataholder* holder, bool genoMapUsed,
		bool continMapUsed){

	outstream << "Network\tScore\tCount\n";
	for(map<float,  vector<vector<vector<int> > > >::reverse_iterator sortedIter=sortedModels.rbegin();
		sortedIter != sortedModels.rend(); ++sortedIter){
// outstream << "num strings=" << sortedIter->second.size() << "\n";
		for(vector<vector<vector<int> > >::iterator strIter = sortedIter->second.begin(); strIter != sortedIter->second.end();
			++strIter){
// 			outstream << *strIter << "\t" << sortedIter->first << "\t" << modelHolder[*strIter].count << "\n";
			outstream << constructBayesStr(*strIter, holder, genoMapUsed, continMapUsed) << "\t" << sortedIter->first
				<< "\t" << modelHolder[*strIter].count << "\n";
		}
	}
}

#ifdef HAVE_CXX_MPI

void GABayes::sendAndReceiveGenomes(int totalNodes, int myRank, GASimpleGA* ga){

	GA2DBinaryStringGenome& genome=(GA2DBinaryStringGenome &)ga->statistics().bestIndividual();

	int height = genome.height();
	int width = genome.width();
	int sendSize = height * width;

// if(myRank==1){
// cout << "totalnodes=" << totalNodes << endl;
// for(int i=0; i<genome.height(); i++){
// 	cout << " " << i;
// }
// cout << "\n-------- SENT ---------------\n";
// for(int i=0; i<genome.height(); i++){
// 	cout << "i=" << i << " ";
// 	for(int j=0; j<genome.width(); j++){
// 	cout <<  genome.gene(i,j) << " ";
// 	}
// 	cout << endl;
// }
// }

	// transfer entire binary genome
	unsigned char* send = new unsigned char[sendSize];
	// send new scores
	float sendScore = genome.score();
	float * recvScore = new float[totalNodes];
	MPI_Allgather(&sendScore, 1, MPI_FLOAT, recvScore, 1, MPI_FLOAT, MPI_COMM_WORLD);


// if(myRank==1){
// 	cout << myRank << " recvScore=>" << recvScore[0] << " " << recvScore[1] << endl;
// 	cout << myRank << " creating genome to send" << endl;
// }
	int index=0;
	// fill send buffer
	for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
			send[index]=genome.gene(i,j);
			index++;
		}
	}
	int recvSize = sendSize * totalNodes;
	unsigned char* recv = new unsigned char[recvSize];
	MPI_Allgather(send, sendSize, MPI_UNSIGNED_CHAR, recv, sendSize, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);

	updateWithMigration(recv, recvScore, sendSize, totalNodes, myRank, ga);

	delete [] recvScore;
	delete [] send;
	delete [] recv;
}


void GABayes::updateWithMigration(unsigned char* newGenes, float * recvScores,
	int genomeSize,	int totalNodes, int myRank, GASimpleGA* ga){
	GAPopulation pop(ga->population());
	int geneIndex=0;

	for(int node=0; node < totalNodes; node++){
		if(myRank==node){
			geneIndex+=genomeSize;
			continue;
		}
		GAGenome *tmpInd = ga->population().individual(0).clone();
		GA2DBinaryStringGenome& genome = (GA2DBinaryStringGenome&)*tmpInd;
		genome.score(recvScores[node]);
		for(int i=0; i<genome.height(); i++){
			for(int j=0; j<genome.width(); j++){
				genome.gene(i,j,newGenes[geneIndex]);
				geneIndex++;
			}
		}

// 		for(int i=0; i<genomeSize; i++){
// 			genome.gene(i, mpiGenomes[node].codons[i]);
// 		}
		pop.add(genome);


// if(myRank==0){
// for(int i=0; i<genome.height(); i++){
// 	cout << " " << i;
// }
// cout << "\n-------- RECEIVED ---------------\n";
// for(int i=0; i<genome.height(); i++){
// 	cout << "i=" << i << " ";
// 	for(int j=0; j<genome.width(); j++){
// 	cout <<  genome.gene(i,j) << " ";
// 	}
// 	cout << endl;
// }
// }

		delete tmpInd;
	}

	for(int i=0; i < totalNodes-1; i++)
		pop.destroy();

	ga->population(pop);
}

#endif

