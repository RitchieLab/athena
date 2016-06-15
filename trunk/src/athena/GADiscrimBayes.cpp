#include "GADiscrimBayes.h"
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
GADiscrimBayes::GADiscrimBayes(){
	initializeParams();
}

///
/// Destructor
///
GADiscrimBayes::~GADiscrimBayes(){
// 		if(geLog != NULL){
// 				delete geLog;
// 				geLog=NULL;
// 		}
		freeMemory();
// 	for(size_t i=0; i<categoryGAs.size(); i++){
// 		delete categoryGAs[i];
// 	}
// 	for(size_t i=0; i<categoryDatasest.size(); i++){
// 		delete categoryDatasest[i];
// 	}
// 	if(caseGA != NULL)
// 		delete caseGA;
// 	if(controlGA != NULL)
// 		delete controlGA;
// 	if(caseDataset != NULL)
// 		delete caseDataset;
// 	if(controlDataset != NULL)
// 		delete controlDataset;
}


void GADiscrimBayes::freeMemory(){
	for(size_t i=0; i<categoryGAs.size(); i++){
		delete categoryGAs[i];
	}
	categoryGAs.clear();
	for(size_t i=0; i<categoryDatasets.size(); i++){
		delete categoryDatasets[i];
	}
	categoryDatasets.clear();
}

///
/// Initializes parameters for basic run of algorithm.  These
/// can be modified using the set_params function
///
void GADiscrimBayes::initializeParams(){

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
	 topModelsUsed = 5;

		maxChildren = 5;
		maxParents = 3;
		limitMethodType = "BESTMI"; // "RANDOM" is also valid

	 // establish map for parameters
	 paramMap["INITCONNECTP"] = initConnectProb;
	 paramMap["NUMTOPMODELS"] = modelsToUse;
	 paramMap["MAXPARENTS"] = maximumParents;
	 paramMap["MAXCHILDREN"] = maximumChildren;
	 paramMap["LIMITMETHOD"] =limitMethod;
// 	 paramMap["CASEALLFILE"] = caseAllFileName;
// 	 paramMap["CONTROLALLFILE"] = controlAllFileName;
	 paramMap["CATEGORYALLFILES"] = categoryAllFileNames;

	 gaSelectorMap["ROULETTE"] = RouletteWheelSelection;
	 gaSelectorMap["PARETO"] = ParetoFrontSelection;
	 gaSelectorMap["PARETORANK"] = ParetoRankSelection;

	 useAllVars = false;
	 outputDotFiles = true;
	 categorizedPhenos = false;
	 numGenotypes = 0;
	 numContinuous = 0;
	 categoryAllFiles="";

	 ga = NULL;
	 maxBest = true;
	 gaSelector = RouletteWheelSelection;
// 	 caseDataset = controlDataset = NULL;
// 	 caseGA=controlGA=NULL;
	 GABayesSolution* sol = (GABayesSolution*)GAFunct::getBlankSolution();
	 pop.insert(sol);

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
void GADiscrimBayes::setParams(AlgorithmParams& algParam, int numExchanges, int numGenos,
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
							case modelsToUse:
								topModelsUsed = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
// cout << "topModelsUsed=" << topModelsUsed << endl;
								break;
							case maximumChildren:
								maxChildren = Stringmanip::stringToNumber<int>(mapIter->second);
								break;
							case maximumParents:
								maxParents = Stringmanip::stringToNumber<int>(mapIter->second);
								break;
							case limitMethod:
								limitMethodType = mapIter->second;
								break;
							case	categoryAllFileNames:
								categoryAllFiles = mapIter->second;
								break;
// 							case	controlAllFileName:
// 								controlAllFile = mapIter->second;
// 								break;
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
void GADiscrimBayes::setGAParams(vector<unsigned int>& excludedGenos,
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
		 pop.setConvertScores(true);
	 }

}


///
/// Set list of included variables
/// @param excludedGenos Vector of excluded variables
/// @param nVars Total number of variables
///
void GADiscrimBayes::addGenotypes(vector<unsigned int>& excludedGenos,
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
void GADiscrimBayes::addContin(vector<unsigned int>& excludedContin,
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
void GADiscrimBayes::setDataset(Dataset* newSet){
	set = newSet;

	if(!categorizedPhenos){
		set->getHolder()->reNumberStatus();
		categorizedPhenos = true;
	}

// 	if(caseDataset != NULL)
// 		delete caseDataset;
// 	if(controlDataset != NULL)
// 		delete controlDataset;

// 	for(size_t i=0; i<categoryDatasets.size(); i++){
// 		delete categoryDatasets[i];
// 	}
// 	categoryDatasets.clear();
	freeMemory();
// 	vector<Dataset*> splitSets = newSet->splitCaseControl();
	categoryDatasets = newSet->splitCategories();
// 	controlDataset=splitSets[0];
// 	caseDataset = splitSets[1];
	for(size_t i=0; i<varList.size(); i++){
		if(!varList[i]->isGeno()){
			varList[i]->setNumLevels(newSet->getNumLevels(varList[i]->getIndex()));
		}
	}
	GAFunct::setDatasets(categoryDatasets, varList, categoryAllFiles.size() < 1);
}

///
/// Set up GA
///
void GADiscrimBayes::configGA(GASimpleGA* ga){
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
		if(categoryAllFiles.size() < 1)
			ga->initialize();
}


///
/// Initializes the algorithm
///
void GADiscrimBayes::initialize(){
	// total number of variables
	totalVars = varList.size();
	GAFunct::setInitConP(initProbConn);
	for(size_t i=0; i<categoryDatasets.size(); i++){
		Athena2DArrayGenome<int> genome(maxParents,totalVars, GAFunct::GAObjective);
		genome.category(i);
		genome.initializer(GAFunct::init);
		genome.mutator(GAFunct::mutate);
		categoryGAs.push_back(new GASimpleGA(genome));
		configGA(categoryGAs.back());
	}
}

///
/// Establishes the algorithm for the run based on the parameters
/// set in this algorithm.
///
void GADiscrimBayes::run(){
		while(!ga->done()){
				ga->step();
		}
}

///
/// Runs an indicated interval for the algorithm. This interval is set
/// in the stepSize variable.
///
int GADiscrimBayes::step(){
	int completed =0;
// 	if(controlAllFile.length() > 0 && caseAllFile.length() > 0){
	if(categoryAllFiles.length() > 0){
		return 0;
	}

	// run each separately
	for(unsigned int i=0; i < stepSize; i++){
// 		caseGA->step();
// GA2DBinaryStringGenome genome = (GA2DBinaryStringGenome&)caseGA->population().best(0);
// vector<vector<int> > eq = constructEquation(genome);
// writeGenoNet(eq);
// cout << "step=" << i << " case best score=" << genome.score() << "\n";
// 		controlGA->step();
		for(size_t i=0; i<categoryGAs.size(); i++){
			categoryGAs[i]->step();
		}
	}

#ifdef HAVE_CXX_MPI
	// perform migration
	for(size_t i=0; i<categoryGAs.size(); i++){
		sendAndReceiveGenomes(totalNodes,myRank,categoryGAs[i]);
	}
// 	sendAndReceiveGenomes(totalNodes,myRank,caseGA);
// 	sendAndReceiveGenomes(totalNodes,myRank,controlGA);
#endif

	return completed;
}

/// write out equation
void GADiscrimBayes::writeGenoNet(vector<vector<int> >& eq){
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
void GADiscrimBayes::setRand(unsigned int seed){
	GARandomSeed(seed);
	srand(seed);
}


/// Writes graphical representation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GADiscrimBayes::writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
	bool mapUsed, bool ottDummy, bool continMapUsed){
	os << "implement writeGraphical" << endl;
// 	GEObjective::outputModel(os, sol, holder, mapUsed, ottDummy, continMapUsed);
}

///
/// Writes model represented as an equation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GADiscrimBayes::writeEquation(ostream& os, Solution* sol, data_manage::Dataholder* holder,
	bool mapUsed, bool ottDummy, bool continMapUsed){
	os << "implement writeEquation" << endl;
// 	GEObjective::outputEquation(os, sol, holder, mapUsed, ottDummy, continMapUsed);
}

///
/// Returns graphical file extension appropriate to solution
/// @return extension
///
std::string GADiscrimBayes::getGraphicalFileExt(){
	return ".dot";
// 	return GEObjective::getGraphicalExt();
}

///
/// Finishes logs and compiles them into a single file
/// in cases where there are multiple processors
/// @param basename
/// @param cv
///
void GADiscrimBayes::finishLog(string basename, int cv){
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
vector<Solution*> GADiscrimBayes::runValidation(std::string sumFile){
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
void GADiscrimBayes::setConfigDefaults(Config& configuration, AlgorithmParams& algParam){

	// algorithm will write its own output
	configuration.setLogType("NONE");
	if(configuration.getSummaryOnly()==Config::True){
		outputDotFiles=false;
	}
	configuration.setSummaryOnly("SUPPRESS");
	configuration.setContinAdjust("MAKECATEGORIALRECODE");
	configuration.setMultiCategory(true);
	configuration.setStatusAdjust("NONE");
	outputName = configuration.getOutputName();
	currCV = configuration.getStartCV();
}


///
/// Extract list of files from .all files of GABayes run
///
void GADiscrimBayes::finalFromFile(Dataset* testing, Dataset* training,
	data_manage::Dataholder* holder, bool mapUsed, bool ottDummy, bool continMapUsed){

	// create map using variable list
// 	vector<Variable*> varList;
	string genoPrefix, continPrefix;
	map<string,int> nameToIndex;
	if(mapUsed)
		genoPrefix="";
	else
		genoPrefix="G";

	if(continMapUsed)
		continPrefix="";
	else
		continPrefix="C";
	string varName;

	for(size_t i=0; i<varList.size(); i++){
		if(varList[i]->isGeno()){
			varName = genoPrefix + varList[i]->getName(holder);
		}
		else{
			varName = continPrefix + varList[i]->getName(holder);
		}
		nameToIndex[varName]=i;
	}

	// loop through files
	vector<string> fileInfo = Stringmanip::split(categoryAllFiles, ' ');
	map<int, string> filenameMap;
	// should be outcome category <space> filename
	for(size_t i=0; i<fileInfo.size(); i+=2){
		filenameMap[Stringmanip::stringToNumber<int>(fileInfo[i])]=fileInfo[i+1];
	}
	vector<map<vector<vector<int> >, ModelScores> > models;
	int category=0;
	// iterate through filename map and read files for appropriate category
	for(map<int,string>::iterator iter=filenameMap.begin(); iter != filenameMap.end(); ++iter){
		map<vector<vector<int> >, ModelScores> fileModels;
		readAllFile(iter->second, fileModels, nameToIndex, holder, category, mapUsed, continMapUsed);
		models.push_back(fileModels);
		category++;
	}

	runDiscriminantAnalysis(models, testing, training, holder, mapUsed,
		continMapUsed);
}


void GADiscrimBayes::readAllFile(string allFileName,map<vector<vector<int> >, ModelScores>& models,
	map<string,int>& nameToIndex, data_manage::Dataholder* holder,
		int category, bool genoMapUsed, bool continMapUsed){
	ifstream allIn(allFileName.c_str(), ios::in);
	if(!allIn.is_open()){
		throw AthenaExcept("ERROR: Unable to open " + allFileName + "\n");
	}
	string line;
	getline(allIn, line);
	getline(allIn, line);
	getline(allIn, line);
	string delim = "\t";
	size_t nVariables = nameToIndex.size();
	vector<int> empty;
	map<vector<vector<int> >, ModelScores> modelHolder;

	do{
		vector<string> columns = Stringmanip::split(line, delim[0]);
		vector<vector<int> >conns(nVariables, empty);
		string modelStr = columns[0];
// cout << "modelStr=" << columns[0] << endl;
		vector<string> nodes = Stringmanip::split(modelStr,']');
		for(vector<string>::iterator iter=nodes.begin(); iter != nodes.end(); ++iter){
			string symb=(*iter).substr(1,(*iter).length()-1);
			vector<string> pcs=data_manage::Stringmanip::split(symb, '|');
			// check for parents and add to connections
			if(pcs.size() > 1){
// cout << "pcs[0]=" << pcs[0] << endl;
				int childIdx = nameToIndex[pcs[0]];
// cout << "childIdx=" << childIdx << endl;
				vector<string> parents=data_manage::Stringmanip::split(pcs[1],':');
				for(vector<string>::iterator parIter=parents.begin(); parIter != parents.end();
					++parIter){
// cout << "parIter=" << *parIter << " index=" << nameToIndex[*parIter] << endl;
					conns[childIdx].push_back(nameToIndex[*parIter]);
				}
			}
		}

		if(modelHolder.find(conns) == modelHolder.end()){
			modelHolder[conns].count=1;
			modelHolder[conns].score=Stringmanip::stringToNumber<float>(columns[2]);
		}
		else{
			modelHolder[conns].count++;
		}
		getline(allIn,line);
	}while(!allIn.eof());

	allIn.close();
	selectTopModels(modelHolder, models, holder, category, genoMapUsed, continMapUsed);
}

///
/// Performs discriminant analysis
///
void GADiscrimBayes::runDiscriminantAnalysis(vector<map<vector<vector<int> >, ModelScores> >& models,
	Dataset* testing, Dataset* training, data_manage::Dataholder* holder, bool mapUsed,
	bool continMapUsed){

	pruneModels(models, holder, mapUsed, continMapUsed);

	int totalInds=0;
	vector<vector<double> > emptyVec;
	vector<vector<vector<double> > > orphanProbs(models.size(), emptyVec);
	for(size_t i=0; i<models.size(); i++){
		calcProbTables(categoryDatasets[i], orphanProbs[i], holder);
		totalInds += categoryDatasets[i]->numInds();
	}
	vector<double> missingValues(models.size(), 0.0);
	vector<double> categoryRatios;
	for(size_t i=0; i<models.size(); i++){
// cout << "num cat=" << i << " is " << categoryDatasets[i]->numInds() << endl;
		double ratio = categoryDatasets[i]->numInds() / double(totalInds);
// cout << "ratio=" << ratio << endl;
		int numCat = categoryDatasets[i]->numInds() + int(testing->numInds()*ratio+0.5);
		categoryRatios.push_back(ratio);
		missingValues[i] = double(1) / numCat;
		setConditionalTables(categoryDatasets[i], models[i], holder, missingValues[i]);
	}

	// for each individual calculate total probability score by multiplying all the probabilities
	// for all variables in set (using conditional probs when appropriate)
	// Use the probs for the model no matter the status of the individual.
	// Need to loop through the categories making each category the "case" and the others
	// the control.  When only 2 categories, just use the category=1 as the case

	IndivResults emptyResult;

	size_t caseIdx=0;
	if(models.size() == 2)
		caseIdx=1;
	bool headerNeeded=true;
	for(;caseIdx<models.size();caseIdx++){
// cout << "caseIdx=" << caseIdx << " original status=" << holder->getOriginalStatus(caseIdx) << endl;
		vector<IndivResults> trainingScores(training->numInds(),emptyResult),
			testingScores(testing->numInds(),emptyResult);

		setIndModScores(training, models[caseIdx], trainingScores, orphanProbs[caseIdx], caseIdx);
		setIndModScores(testing, models[caseIdx], testingScores, orphanProbs[caseIdx], caseIdx);

		for(size_t conIdx=0; conIdx < models.size(); conIdx++){
			if(caseIdx==conIdx)
				continue;
// cout << "conIndex " << conIdx << endl;
			setIndModScores(training, models[conIdx], trainingScores, orphanProbs[conIdx], caseIdx);
			setIndModScores(testing, models[conIdx], testingScores, orphanProbs[conIdx], caseIdx);
// 			categoryDatasets[conIdx]->numInds();
		}
		double trainingAUC=setPredictedScores(trainingScores, models, caseIdx,categoryRatios[caseIdx]);
		double testingAUC=setPredictedScores(testingScores, models, caseIdx,categoryRatios[caseIdx]);
// cout << "trainingAUC=" << trainingAUC << " testingAUC=" << testingAUC << endl;
#ifdef HAVE_CXX_MPI
if(myRank == 0){
#endif
		if(caseIdx==0 || models.size()==2){
			outputProbFiles(trainingScores, training, testingScores, testing, models);
		}

		string trainingFile = outputName + ".cv." + Stringmanip::numberToString(currCV) +
			".case." + Stringmanip::numberToString(holder->getOriginalStatus(caseIdx)) + ".train.inds.txt";
		writeIndScores(trainingFile, trainingScores);
		string testingFile = outputName + ".cv." + Stringmanip::numberToString(currCV) +
				".case." + Stringmanip::numberToString(holder->getOriginalStatus(caseIdx)) + ".test.inds.txt";
		writeIndScores(testingFile, testingScores);

		// add this CV result to .sum file
		string sumFile = outputName + ".GADBN.sum";
		ofstream sumstream;
		if(currCV > 1 || !headerNeeded){
			sumstream.open(sumFile.c_str(), ios::app);
		}
		else{
			sumstream.open(sumFile.c_str(), ios::out);
			if(headerNeeded){
				sumstream << "CV\tCase Category\tCase Models\tCount";

				for(size_t i=0; i<models.size()-1; i++){
					sumstream << "\tCategory\tControl Models\tCount";
				}
				sumstream << "\tTrain AUC\tTest AUC\n";
				headerNeeded=false;
			}
		}

		mdScores tmpScore;
		vector<mdScores> eVec;
		vector<vector<mdScores> > categoryMods(models.size(), eVec);
		// sort models by score
		for(size_t i=0; i<models.size(); i++){
			for(map<vector<vector<int> >, ModelScores>::iterator mIter=models[i].begin();
				mIter != models[i].end(); ++mIter){
					tmpScore.mString = constructBayesStr(mIter->first,holder, mapUsed, continMapUsed, true);
					tmpScore.mPtr = &(mIter->second);
					categoryMods[i].push_back(tmpScore);
			}
			sort(categoryMods[i].begin(), categoryMods[i].end(), sortByScore);
		}

		int modIndex=0;
		sumstream << currCV << "\t" << holder->getOriginalStatus(caseIdx)<< "\t"<< categoryMods[caseIdx][modIndex].mString << "\t"
			<< categoryMods[caseIdx][modIndex].mPtr->count;
		for(size_t i=0; i<models.size(); i++){
			if(i==caseIdx)
				continue;
			sumstream << "\t" <<  holder->getOriginalStatus(i) << "\t" << categoryMods[i][modIndex].mString << "\t"
				<< categoryMods[i][modIndex].mPtr->count;
		}
		sumstream << "\t" << trainingAUC << "\t" << testingAUC << "\n";

		for(modIndex=1; modIndex<topModelsUsed; modIndex++){
			if(modIndex < categoryMods[caseIdx].size()){
				sumstream << "\t" << holder->getOriginalStatus(caseIdx) << "\t"<< categoryMods[caseIdx][modIndex].mString << "\t"
					<< categoryMods[caseIdx][modIndex].mPtr->count;
			}
			else{
				sumstream <<  "\t\t\t";
			}
			for(size_t i=0; i<models.size(); i++){
				if(i==caseIdx)
					continue;
				if(modIndex < categoryMods[i].size()){
					sumstream << "\t" <<  holder->getOriginalStatus(i)  << "\t" << categoryMods[i][modIndex].mString << "\t"
						<< categoryMods[i][modIndex].mPtr->count;
				}
				else{
					sumstream << "\t\t\t";
				}
			}
			sumstream << "\t\t\n";
// 			if(modIndex < caseMods.size()){
// 				sumstream << "\t" << caseMods[modIndex].mString
// 				<<  "\t" << caseMods[modIndex].mPtr->count;
// 			}
// 			else{
// 				sumstream << "\t\t";
// 			}
// 			if(modIndex < controlMods.size()){
// 				sumstream << "\t" << controlMods[modIndex].mString
// 					<< "\t" << controlMods[modIndex].mPtr->count << "\t\t\n";
// 			}
// 			else{
// 				sumstream << "\t\t\t\t\n";
// 			}
		}

		sumstream.close();
#ifdef HAVE_CXX_MPI
}
#endif
	}
	currCV++;
		// sort models by score
// 		map<vector<vector<int> >, ModelScores>::iterator caseIter=caseModels.begin();
// 		map<vector<vector<int> >, ModelScores>::iterator conIter=controlModels.begin();
// 		mdScores tmpScore;
// 		vector<mdScores> caseMods, controlMods;
//
// 		for(int i=0; i<topModelsUsed; i++){
// 			if(caseIter != caseModels.end()){
// 				tmpScore.mString = constructBayesStr(caseIter->first,holder, mapUsed, continMapUsed);
// 				tmpScore.mPtr = &(caseIter->second);
// 				caseMods.push_back(tmpScore);
// 				++caseIter;
// 			}
// 			if(conIter != controlModels.end()){
// 				tmpScore.mString = constructBayesStr(conIter->first, holder, mapUsed, continMapUsed);
// 				tmpScore.mPtr = &(conIter->second);
// 				controlMods.push_back(tmpScore);
// 				++conIter;
// 			}
// 		}
//
// 		sort(caseMods.begin(), caseMods.end(), sortByScore);
// 		sort(controlMods.begin(), controlMods.end(), sortByScore);

//
// 		sumstream << currCV << "\t" << caseMods[modIndex].mString << "\t"
// 			<< caseMods[modIndex].mPtr->count << "\t"
// 			<< controlMods[modIndex].mString << "\t"
// 			<< controlMods[modIndex].mPtr->count << "\t"
// 			<< trainingAUC << "\t" << testingAUC << "\n";

}


///
/// Outputs probability tables for top models used in analysis
///
void GADiscrimBayes::outputProbFiles(vector<IndivResults>& trainingScores, data_manage::Dataset* training,
	vector<IndivResults>&  testingScores, data_manage::Dataset* testing,
	vector<map<vector<vector<int> >, ModelScores> >& models){
	int currModelIdx=0;
	for(size_t category=0; category<models.size(); category++){
		size_t nModels = models[category].size();
		string filename = outputName + ".cv." + Stringmanip::numberToString(currCV) +
			".cat." + Stringmanip::numberToString(training->getHolder()->getOriginalStatus(category)) + ".train.probtable";
// cout << "category=" << category << " probtable cat=" << training->getHolder()->getOriginalStatus(category) <<endl;
		ofstream probFile(filename.c_str(), ios::out);
		probFile << "ind ID";
		for(size_t i=0; i<nModels; i++){
			probFile << "\tV" << i+1;
		}
		probFile << "\n";
		// for each individual in the sets output individual ID and then the
		// probability for each model
		for(size_t i=0; i<training->numInds(); i++){
			probFile << (*training)[i]->getID();
// cout << "i=" << i << endl;
			for(size_t j=currModelIdx; j<(nModels+currModelIdx); j++){
				probFile << "\t" << trainingScores[i].scores[j];
			}
			probFile << "\n";
		}
		probFile.close();

		filename = outputName + ".cv." + Stringmanip::numberToString(currCV) +
			".cat." + Stringmanip::numberToString(testing->getHolder()->getOriginalStatus(category)) + ".test.probtable";
		probFile.open(filename.c_str(), ios::out);
		probFile << "ind ID";
		for(size_t i=0; i<nModels; i++){
			probFile << "\tV" << i+1;
		}
		probFile << "\n";
		for(size_t i=0; i<testing->numInds(); i++){
			probFile << (*testing)[i]->getID();
			for(size_t j=currModelIdx; j<(nModels+currModelIdx); j++){
				probFile << "\t" << testingScores[i].scores[j];
			}
			probFile << "\n";
		}
		probFile.close();
		currModelIdx += nModels;
	}
}


///
/// Performs discriminant analysis
///
void GADiscrimBayes::getAdditionalFinalOutput(Dataset* testing, Dataset* training,
	data_manage::Dataholder* holder, bool mapUsed, bool ottDummy, bool continMapUsed){

// 	if(controlAllFile.length() > 0 && caseAllFile.length() > 0){
	if(categoryAllFiles.length() > 0){
// 		if(currCV==1)
		finalFromFile(testing, training, holder, mapUsed, ottDummy, continMapUsed);
		return;
	}

	map<vector<vector<int> >, ModelScores > empty;
	vector<map<vector<vector<int> >, ModelScores > > models(categoryGAs.size(), empty);
	for(size_t i=0; i<categoryGAs.size(); i++){
		totalModels(categoryGAs[i], models[i], holder, i, mapUsed, continMapUsed);
	}

	runDiscriminantAnalysis(models, testing, training, holder, mapUsed, continMapUsed);
}


void GADiscrimBayes::pruneModels(vector<map<vector<vector<int> >, ModelScores> >& models,
	data_manage::Dataholder* holder, bool genoMapUsed, bool continMapUsed){

	string endName, outName;
// cout << "start pruneModels" << endl;
	for(size_t i=0; i<models.size(); i++){
// cout << "i=" << i << endl;
		map<vector<vector<int> >, ModelScores> temp;
		multimap<float, connComparison> sortedModels;
		connComparison comparison;
			for(map<vector<vector<int> >, ModelScores>::iterator iter=models[i].begin(); iter != models[i].end();
				++iter){
				vector<vector<int> > conns = iter->first;
				ModelScores modScore = iter->second;
// cout << "call prune" << endl;
				modScore.score = GAFunct::prune(conns, i);
				comparison.originalConns = iter->first;
				comparison.newConns = conns;
				comparison.newScore=modScore.score;
				if(temp.find(conns) != temp.end()){
					temp[conns].count += iter->second.count;
				}
				else{
					temp[conns]=modScore;
				}
				// sorted by old score
				sortedModels.insert(std::pair<float, connComparison>(iter->second.score, comparison));
			}
			models[i] = temp;
			#ifdef HAVE_CXX_MPI
				if(myRank == 0){
			#endif
// cout << "myRank=" << myRank  << "sorted models" << endl;
				outName = outputName + ".cv." + Stringmanip::numberToString(currCV) + ".cat." +
				Stringmanip::numberToString(holder->getOriginalStatus(i)) + ".reduced";
				endName = ".cat." +  Stringmanip::numberToString(holder->getOriginalStatus(i)) + ".dot";
				if(outputDotFiles)
					writeDotFiles(sortedModels, holder, genoMapUsed, continMapUsed, endName);
				writeReducedFile(sortedModels, holder, genoMapUsed, continMapUsed, outName);
#ifdef HAVE_CXX_MPI
}
#endif
	}

}


// void GADiscrimBayes::pruneModels(map<vector<vector<int> >, ModelScores>& caseModels,
// 		map<vector<vector<int> >, ModelScores>& controlModels, data_manage::Dataholder* holder,
// 		bool genoMapUsed, bool continMapUsed){
// 	map<vector<vector<int> >, ModelScores> temp;
// 	multimap<float, connComparison> sortedModels;
// 	connComparison comparison;
// 	for(map<vector<vector<int> >, ModelScores>::iterator iter=caseModels.begin(); iter != caseModels.end();
// 		++iter){
// 		vector<vector<int> > conns = iter->first;
// 		ModelScores modScore = iter->second;
// // cout << "original score=" << modScore.score << endl;
// 		modScore.score = GAFunct::pruneCase(conns);
// 		comparison.originalConns = iter->first;
// 		comparison.newConns = conns;
// 		comparison.newScore=modScore.score;
// 		if(temp.find(conns) != temp.end()){
// 			temp[conns].count += iter->second.count;
// 		}
// 		else{
// 			temp[conns]=modScore;
// 		}
// 		// sorted by old score
// 		sortedModels.insert(std::pair<float, connComparison>(iter->second.score, comparison));
// 	}
// 	caseModels = temp;
// 	string endName, outName;
//
// #ifdef HAVE_CXX_MPI
// if(myRank == 0){
// #endif
// // cout << "myRank=" << myRank  << "sorted models" << endl;
// 	outName = outputName + ".cv." + Stringmanip::numberToString(currCV) + ".case.reduced";
// 	endName = ".case.dot";
// 	writeDotFiles(sortedModels, holder, genoMapUsed, continMapUsed, endName);
// 	writeReducedFile(sortedModels, holder, genoMapUsed, continMapUsed, outName);
// #ifdef HAVE_CXX_MPI
// }
// #endif
//
// 	temp.clear();
// 	sortedModels.clear();
// 	for(map<vector<vector<int> >, ModelScores>::iterator iter=controlModels.begin(); iter != controlModels.end();
// 		++iter){
// 		vector<vector<int> > conns = iter->first;
// 		ModelScores modScore = iter->second;
// // cout << "control original score=" << modScore.score << endl;
// 		modScore.score = GAFunct::prune(conns);
// 		modScore.originalConns = iter->first;
// 		comparison.originalConns = iter->first;
// 		comparison.newConns = conns;
// 		comparison.newScore=modScore.score;
// 		if(temp.find(conns) != temp.end()){
// 			temp[conns].count += iter->second.count;
// 		}
// 		else{
// 			temp[conns]=modScore;
// 		}
// 		sortedModels.insert(std::pair<float, connComparison>(iter->second.score, comparison));
// 	}
// 	controlModels = temp;
//
// #ifdef HAVE_CXX_MPI
// if(myRank == 0){
// #endif
// 	endName = ".control.dot";
// 	writeDotFiles(sortedModels, holder, genoMapUsed, continMapUsed, endName);
// 	outName = outputName + ".cv." + Stringmanip::numberToString(currCV) + ".control.reduced";
// 	writeReducedFile(sortedModels, holder, genoMapUsed, continMapUsed, outName);
// #ifdef HAVE_CXX_MPI
// }
// #endif
//
// }

///
/// Create a .reduced file showing how every model was reduced
///
void GADiscrimBayes::writeReducedFile(multimap<float, connComparison>& sortedModels, Dataholder* holder,
	bool genoMapUsed, bool continMapUsed, string fileName){

	ofstream os(fileName.c_str(), ios::out);
	os << "Original\tReduced\tOrig Score\tReduced Score\n";
	for(multimap<float, connComparison>::iterator iter = sortedModels.begin(); iter != sortedModels.end();
		++iter){
		string original = constructBayesStr(iter->second.originalConns,  holder, genoMapUsed, continMapUsed);
		string newmodel = constructBayesStr(iter->second.newConns, holder, genoMapUsed, continMapUsed);
		os << original << "\t" << newmodel << "\t" << iter->first << "\t" << iter->second.newScore << "\n";

	}
	os.close();

}

///
/// Create a .dot file for every pruned model in the set
///
void GADiscrimBayes::writeDotFiles(multimap<float, connComparison>& sortedModels, Dataholder* holder,
	bool genoMapUsed, bool continMapUsed, string endName){

	multimap<float, connComparison>::iterator iter = sortedModels.begin();
	int modNum=1;
	while(iter != sortedModels.end()){
		string outName = outputName + ".cv." + Stringmanip::numberToString(currCV) + "." +
			Stringmanip::numberToString(modNum) + endName;
		string network = constructBayesStr(iter->second.newConns, holder, genoMapUsed,
			continMapUsed, true);
		// tokenize
		vector<string> nodes = Stringmanip::split(network, ']');
		ofstream os;
		os.open(outName.c_str(), ios::out);
		os << "digraph G{\n";
		os << "\tgraph [ dpi = 300 ];\n";
		os << "\tsize=\"8.5,11.0\";\n";
		os << "\tdir=\"none\";\n";
		os << "\trankdir=\"LR\";\n";
		os << "\torientation=\"portrait\";\n";

		std::set<std::string> include;
		std::set<std::string> added;
		for(vector<string>::iterator symbIter=nodes.begin(); symbIter != nodes.end();
			++symbIter){
			// remove '['
			string symb=(*symbIter).substr(1,(*symbIter).length()-1);
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
		os.close();
		++iter;
		modNum++;
	}
}

///
/// Set predicted phenotype for each individual in set
/// @return AUC
///
double GADiscrimBayes::setPredictedScores(vector<IndivResults>& indScores,
	vector<map<vector<vector<int> > ,ModelScores> >& models,
	int caseIdx, double caseRatio){

	long double caseScore, conScore;
	int sIndex, totalCount;
	stat::TestResult tempResult;
	std::vector<stat::TestResult> results;
	for(vector<IndivResults>::iterator indIter=indScores.begin(); indIter != indScores.end();
		++indIter){
		caseScore = conScore = 0.0;
		sIndex=totalCount=0;
		for(map<vector<vector<int> >,ModelScores>::iterator caseModIter=models[caseIdx].begin(); caseModIter != models[caseIdx].end();
			++caseModIter){
			caseScore += indIter->scores[sIndex++] * caseModIter->second.count;
			totalCount += caseModIter->second.count;
		}
		caseScore /= (long double)totalCount;
// cout << "indIter->scores.size()=" << indIter->scores.size() << endl;
// cout <<  "caseScore=" << caseScore << endl;
		totalCount=0;
		for(size_t conIdx=0; conIdx < models.size(); conIdx++){
			if(conIdx==caseIdx)
				continue;
			for(map<vector<vector<int> >,ModelScores>::iterator conModIter=models[conIdx].begin(); conModIter != models[conIdx].end();
				++conModIter){
				conScore += indIter->scores[sIndex++] * conModIter->second.count;
				totalCount += conModIter->second.count;
			}
		}
		conScore /= (long double)totalCount;
// cout << "conScore=" << conScore << endl;
// cout << "caseRatio=" << caseRatio << endl;
		indIter->predicted = caseRatio * caseScore / (caseRatio * caseScore +
			(1-caseRatio) * conScore);
		tempResult.score = indIter->predicted;
		tempResult.status = indIter->phenotype;
// cout << "tempResult.score=" << tempResult.score << endl;
// cout << "tempResult.status=" << tempResult.status << endl;
		results.push_back(tempResult);
	}

	//calculate and return AUC
	return stat::AUCCalc::calculateAUC(results);
}

///
/// Sets the individual scores for every individual in set for the models passed
///
void GADiscrimBayes::setIndModScores(Dataset* dset, map<vector<vector<int> >,ModelScores>& models,
	vector<IndivResults>& indScores, vector<vector<double> >& orphanProbs, int caseValue){

	data_manage::Individual * ind;
	for(size_t i=0; i<dset->numInds(); i++){
		ind=(*dset)[i];
		indScores[i].indID = (*dset)[i]->getID();
		indScores[i].phenotype = (*dset)[i]->getStatus();
		// set any that don't match case value to be a control
		if(indScores[i].phenotype == caseValue)
			indScores[i].phenotype = 1;
		else
			indScores[i].phenotype = 0;
		// loop through each model and assign a value to each individual for each model
		for(map<vector<vector<int> >, ModelScores>::iterator modIter=models.begin(); modIter != models.end();
			++modIter){
			indScores[i].scores.push_back(1.0);
			for(vector<ConditionalTable>::iterator tableIter=modIter->second.tables.begin();
				tableIter != modIter->second.tables.end(); ++tableIter){
					// have to get the probability and include it
					int value = 0;

					for(size_t j=0; j<tableIter->parentIndexes.size(); j++){
						value += varList[tableIter->parentIndexes[j]]->getValue(ind) *
							tableIter->cumulativeLevels[j];
					}
					indScores[i].scores.back() *= tableIter->probs[varList[tableIter->nodeIndex]->getValue(ind)][value];
			}

			// add in all the scores for the other variables in the set
			for(size_t vIndex=0; vIndex < varList.size(); vIndex++){
				if(modIter->second.varsWithParents.find(vIndex) == modIter->second.varsWithParents.end()){
					indScores[i].scores.back() *= orphanProbs[vIndex][int(varList[vIndex]->getValue(ind))];
				}
			}
		}
	}

}

void GADiscrimBayes::selectTopModels(map<vector<vector<int> >, ModelScores>& modelHolder,
	map<vector<vector<int> >, ModelScores>& topModels, data_manage::Dataholder* holder,
	int category, bool genoMapUsed, bool continMapUsed){
	// insert into map with score as key
	map<float, vector<vector<vector<int> >  >  > sortedModels;
	for(map<vector<vector<int> >, ModelScores>::iterator iter=modelHolder.begin(); iter != modelHolder.end();
		++iter){
		sortedModels[iter->second.score].push_back(iter->first);
	}

#ifdef HAVE_CXX_MPI
if(myRank == 0){
#endif
	string outName;
	outName = outputName + ".cv." + Stringmanip::numberToString(currCV) + ".cat." +
		Stringmanip::numberToString(holder->getOriginalStatus(category)) + ".uniq";
	ofstream os;
	os.open(outName.c_str(), ios::out);
	writeUniqueFiles(os,sortedModels,modelHolder, holder, genoMapUsed, continMapUsed);
	os.close();
#ifdef HAVE_CXX_MPI
}
#endif

	int i=0;
	map<float,  vector<vector<vector<int> >  > >::reverse_iterator sortedIter=sortedModels.rbegin();
	while(i<topModelsUsed){
		for(vector<vector<vector<int> >  >::iterator modIter=sortedIter->second.begin();
			modIter != sortedIter->second.end(); ++modIter){
			if(i++ < topModelsUsed){
				topModels[*modIter] = modelHolder[*modIter];
			}
			else{
				break;
			}
		}
		++sortedIter;
	}
}

void GADiscrimBayes::totalModels(GASimpleGA* ga,  map<vector<vector<int> >, ModelScores>& topModels,
	data_manage::Dataholder* holder,int category, bool genoMapUsed, bool continMapUsed){
	topModels.clear();

	int nGenomes = ga->population().size();
	map<vector<vector<int> >, ModelScores> modelHolder;
	vector<vector<int> > mNodes;
	for(int i=0; i<nGenomes; i++){
		Athena2DArrayGenome<int> genome = (Athena2DArrayGenome<int>&)ga->population().best(i);
		// construct the model string
		mNodes = constructEquation(genome);
		if(modelHolder.find(mNodes)==modelHolder.end()){
			modelHolder[mNodes].count=1;
			modelHolder[mNodes].score=genome.score();
		}
		else{
			modelHolder[mNodes].count++;
		}
	}

	// gather all models
	#ifdef HAVE_CXX_MPI
		gatherModelInformation(modelHolder);
	#endif
	selectTopModels(modelHolder, topModels, holder, category, genoMapUsed, continMapUsed);

}


///
/// Construct equation string from genome
/// @genome GA2DBinaryStringGenome
///
vector<vector<int> > GADiscrimBayes::constructEquation(Athena2DArrayGenome<int>& genome){
	vector<int> empty;
	vector<vector<int> > conns(varList.size(), empty);
// 	for(int i=0; i<genome.height(); i++){
// 		for(int j=0; j<genome.width(); j++){
// 			if(genome.gene(i,j)==1){
// 				conns[j].push_back(i);
// 			}
// 		}
// 	}
// 	return conns;

	for(int y=0; y<genome.height(); y++){
		for(int x=0; x<genome.width(); x++){
			if(genome.gene(x,y) != -1){
				conns[y].push_back(genome.gene(x,y));
			}
		}
	}
	return conns;

}


///
/// calculates probability tables for every variable in variable list
///
void	GADiscrimBayes::calcProbTables(Dataset* dset, vector<vector<double> >& orphanProbs,
	data_manage::Dataholder* holder){

	orphanProbs.clear();
	// iterate through varlist and set probabilities
	for(size_t i=0; i<varList.size(); i++){
		int nLevels = varList[i]->getNumLevels();
		vector<int> totals(nLevels,0);
		vector<double> table(nLevels,0.0);
		for(unsigned int ind=0; ind < dset->numInds(); ind++){
			totals[int(varList[i]->getValue(dset->getInd(ind)))]++;
		}
		double allTotal = 0.0;
		// exclude highest level (it is for individuals who are missing data)
		for(int l=0; l<nLevels-1; l++){
			allTotal += totals[l];
		}
		for(int i=0; i<nLevels-1; i++){
			table[i] = totals[i] / allTotal;
		}
		// set missing probability to 1.0
		// when set to 1.0 it means that the missing data will not affect overall
		// probability score of the individual
		table[nLevels-1]=1.0;
		orphanProbs.push_back(table);
	}
}


///
/// Determines which values correspond to one or more
/// missing data indicators in the variables
/// @param parents vector of Variable pointers
/// @param missingVals set which will contain the values for missing
/// @param nLevels number of levels for each
///
void	GADiscrimBayes::setMissingIndexes(vector<unsigned int> &parents, std::set<int>& missingVals,
	vector<int>& nLevels, vector<int>& cumulativeLevels){
	// 3 sets with max of 4 on each
	vector<int> max(parents.size(),0);
	for(size_t i=0; i<max.size(); i++){
		max[i]=nLevels[i]-1;
	}

	// track the current index on each level
	// set to max
	vector<int> index(max.size(),0);
	for(size_t i=0; i<index.size(); i++){
		index[i]=max[i];
	}
	// current index to change (last one -- or inner most one)
	int currIndex=index.size()-1;
	bool endLoop=false;
	int value;

	do{
		value=0;
		for(size_t i=0; i<index.size(); i++){
// 			cout << index[i] << " ";
			if(index[i]==max[i]){
// 				cout << " <-- " << index[i] <<" is max";
				for(size_t j=0; j<index.size(); j++){
					value += index[j] * cumulativeLevels[j];
				}
				missingVals.insert(value);
// 			cout << "insert " << value << endl;
				break;
			}
		}
// 		cout << endl;
		index[currIndex]--;
		while(index[currIndex] < 0){
			index[currIndex]=max[currIndex];
			if(currIndex==0)
				endLoop=true;
			else{
				currIndex--;
				index[currIndex]--;
			}
		}
		currIndex=index.size()-1;
	}while(!endLoop);
}

///
/// Create parent data combination
/// @param parentValues
/// @param parents contains indexes to parent genotypes
/// @returns number of different levels(factors) in the parent combined values
///
int GADiscrimBayes::configParentData(vector<int>& parentValues, vector<unsigned int> &parents,
	Dataset* dSet, vector<int>& cumulativeLevels, std::set<int>& missingDataVals){
	// assume three levels (to hold SNP data)
// 	int constLevels = 3; // this has to be changed also
// 	int snpLevels = 3;

	// set number of levels for each parent
	vector<int> nLevels(parents.size(), 0);
	int nl = 1;
	for(size_t i=0; i<parents.size(); i++){
		nLevels[i]=varList[parents[i]]->getNumLevels();
		nl *= nLevels[i];
	}

	cumulativeLevels.assign(parents.size(), 1);
	for(unsigned int i=1; i<cumulativeLevels.size(); i++){
		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels[i-1];
	}

	// create set of missing values
	missingDataVals.clear();
	setMissingIndexes(parents, missingDataVals, nLevels, cumulativeLevels);

// 	missingDataVals.clear();
// 	for(unsigned int i=0; i<parents.size(); i++){
// 		// set each to missing and store all combinations with others set to all values
// 		int val=varList[parents[i]]->getMissingVal() * cumulativeLevels[i];
// 		for(unsigned int j=0; i<parents.size(); j++){
// 			if(i==j)
// 				continue;
// 			for(unsigned int k=0; k<varList[parents[j]]->getNumLevels(); k++){
// 				val += k * cumulativeLevels[j];
// 			}
// 		}
// 		missingDataVals.insert(val);
// 	}

// 	deque<float> args;
	unsigned int nParents = parents.size();
	Individual* ind;

	for(unsigned int i=0; i < dSet->numInds(); i++){
 		ind = (*dSet)[i];
// 		IndividualTerm::setInd(ind);

		int value = 0;
		for(unsigned int j=0; j < nParents; j++){
			value += varList[parents[j]]->getValue(ind) * cumulativeLevels[j];
		}
		parentValues.push_back(value);
	}
	return nl;
}


///
/// Set conditional probability tables for models
///
void GADiscrimBayes::setConditionalTables(Dataset* dset, map<vector<vector<int> > , ModelScores>& topModels,
	Dataholder* holder, double missValue){

	Individual* ind;
	for(map<vector<vector<int> >, ModelScores>::iterator modIter=topModels.begin(); modIter != topModels.end();
		++modIter){
// 		int tableIndex=0;
		modIter->second.tables.clear();

// vector<vector<int> > v=modIter->firscondt;
// cout << constructBayesStr(v,holder,false,false) << endl;
// cout << "score=" << modIter->second.score << endl;
		ConditionalTable emptyTable;
		for(size_t nodeIdx=0; nodeIdx < modIter->first.size(); nodeIdx++){
			// if no parents then skip
			if(modIter->first.at(nodeIdx).empty())
				continue;
// cout << "nodeIdx=" << nodeIdx << endl;
			Variable * childVar = varList[nodeIdx];
			modIter->second.tables.push_back(emptyTable);
			modIter->second.tables.back().nodeIndex = nodeIdx;
			modIter->second.varsWithParents.insert(nodeIdx);
			vector<int> parentVec = modIter->first.at(nodeIdx);
// 			for(vector<int>::iterator parentIter=modIter->first.at(nodeIdx).begin();
// 				parentIter!=modIter->first.at(nodeIdx).end(); ++parentIter){
			for(vector<int>::iterator parentIter=parentVec.begin();
				parentIter!=parentVec.end(); ++parentIter){
				// get each parent
				modIter->second.tables.back().parentIndexes.push_back(*parentIter);
			}


// Individual* ind2;
// cout << "\n";
// for(unsigned int i=0; i<dset->numInds(); i++){
// ind2=(*dset)[i];
// cout << childVar->getValue(ind2);
// for(unsigned int j=0; j<modIter->second.tables.back().parentIndexes.size(); j++){
// cout << "\t" << varList[modIter->second.tables.back().parentIndexes[j]]->getValue(ind2);
// }
// cout << "\n";
// }

			vector<int> parentValues;
			std::set<int> missingDataVals;
			int parentLevels = configParentData(parentValues, modIter->second.tables.back().parentIndexes,
					dset, modIter->second.tables.back().cumulativeLevels, missingDataVals);
			int nodeLevels = childVar->getNumLevels();
			vector<int> inner(parentLevels,0);
			vector<vector<int> > totals(nodeLevels, inner);
			vector<int> nodeTotals(nodeLevels, 0);
					// cycle through individuals and total
// 			int val;
			for(unsigned int i=0; i<dset->numInds(); i++){
				ind=(*dset)[i];
				totals[int(childVar->getValue(ind))][parentValues[i]]++;
// 				nodeTotals[int(childVar->getValue(ind))]++;
			}
			for(unsigned int lvl=0; lvl<nodeLevels-1; lvl++){
				for(unsigned int parLvl=0; parLvl < parentLevels; parLvl++){
					if(missingDataVals.find(parLvl) == missingDataVals.end()){
						nodeTotals[lvl] += totals[lvl][parLvl];
					}
				}
			}

			vector<double> innerProbs(parentLevels, 0.0);
			modIter->second.tables.back().probs.assign(nodeLevels, innerProbs);
			for(int i=0; i<nodeLevels; i++){
				for(int j=0; j<parentLevels; j++){
				if(missingDataVals.find(j) != missingDataVals.end() || i == childVar->getMissingVal())
					modIter->second.tables.back().probs[i][j] = 1.0;
				else if(totals[i][j] > 0)
					modIter->second.tables.back().probs[i][j] = totals[i][j]/double(nodeTotals[i]);
				else
					modIter->second.tables.back().probs[i][j] = missValue;
				}
			}
// cout << "table" << endl;
// for(size_t i=0; i<modIter->second.tables.back().probs.size(); i++){
// 	for(size_t j=0; j<modIter->second.tables.back().probs[i].size(); j++){
// 		cout << "i=" << i << " j=" << j << " p=" << modIter->second.tables.back().probs[i][j] << endl;
// 	}
// }
// exit(1);
		}
// exit(1);
	}
}


///
/// Constructs bayes network representation compatible with bnlearn package in R
/// @param network 2-D vector showing parents for each node
/// @param holder Dataholder
/// @param excludeIndependent Default is false.  Exclude nodes that are independent.
/// @returns network string
///
string GADiscrimBayes::constructBayesStr(vector<vector<int> > network,
	Dataholder* holder, bool genoMapUsed, bool continMapUsed, bool excludeIndependent){
	string netString;
	vector<string> prefixes(2, "");
	if(!genoMapUsed){
		prefixes[1]="G";
	}
	if(!continMapUsed){
		prefixes[0]="C";
	}

	std::set<int> excludedNodes;
	if(excludeIndependent){
		std::set<int> isParent;
		for(size_t i=0; i<network.size(); i++){
			for(size_t j=0; j<network[i].size(); j++){
				isParent.insert(network[i][j]);
			}
		}

		for(size_t i=0; i<network.size(); i++){
			if(network[i].empty() && isParent.find(i) == isParent.end()){
				excludedNodes.insert(i);
			}
		}
	}


	for(size_t i=0; i<network.size(); i++){
		if(excludedNodes.find(i) == excludedNodes.end()){
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
	}
	return netString;
}


///
///  Writes list of all unique models and counts
///
void GADiscrimBayes::writeUniqueFiles(ostream& outstream, map<float,  vector<vector<vector<int> >  > >& sortedModels,
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

///
/// Output individual results with predicted scores
///
void GADiscrimBayes::writeIndScores(string filename, vector<IndivResults>& scores){
	ofstream out;
	out.open(filename.c_str(), ios::out);
	out << "ID\tPHENO\tPREDICTED\n";
	for(vector<IndivResults>::iterator iter=scores.begin(); iter != scores.end();
		++iter){
		out << iter->indID << "\t" << iter->phenotype << "\t" << iter->predicted << "\n";
	}
	out.close();
}

#ifdef HAVE_CXX_MPI

void GADiscrimBayes::sendAndReceiveGenomes(int totalNodes, int myRank, GASimpleGA* ga){

	Athena2DArrayGenome<int>& genome=(Athena2DArrayGenome<int> &)ga->statistics().bestIndividual();

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
	unsigned int* send = new unsigned int[sendSize];
	// send new scores
	float * sendInfo = new float[2];
// 	float sendScore = genome.score();
	sendInfo[0] = genome.score();
	sendInfo[1] = genome.category();
	float * recvScore = new float[totalNodes*2];
	MPI_Allgather(sendInfo, 2, MPI_FLOAT, recvScore, 2, MPI_FLOAT, MPI_COMM_WORLD);


// if(myRank==1){
// 	cout << myRank << " recvScore=>" << recvScore[0] << " " << recvScore[1] << endl;
// 	cout << myRank << " creating genome to send" << endl;
// }
	int index=0;
	// fill send buffer
	for(int y=0; y<height; y++){
		for(int x=0; x<width; x++){
			send[index]=genome.gene(x,y);
			index++;
		}
	}
	int recvSize = sendSize * totalNodes;
	unsigned int* recv = new unsigned int[recvSize];
	MPI_Allgather(send, sendSize, MPI_INT, recv, sendSize, MPI_INT, MPI_COMM_WORLD);

	updateWithMigration(recv, recvScore, sendSize, totalNodes, myRank, ga);

	delete [] sendInfo;
	delete [] recvScore;
	delete [] send;
	delete [] recv;
}


void GADiscrimBayes::updateWithMigration(unsigned int* newGenes, float * recvScores,
	int genomeSize,	int totalNodes, int myRank, GASimpleGA* ga){
	GAPopulation pop(ga->population());
	int geneIndex=0;

	for(int node=0; node < totalNodes; node++){
		if(myRank==node){
			geneIndex+=genomeSize;
			continue;
		}
		GAGenome *tmpInd = ga->population().individual(0).clone();
		Athena2DArrayGenome<int>& genome = (Athena2DArrayGenome<int>&)*tmpInd;
		genome.score(recvScores[node*2]);
		genome.category(recvScores[node*2+1]);
		for(int y=0; y<genome.height(); y++){
			for(int x=0; x<genome.width(); x++){
				genome.gene(x,y,newGenes[geneIndex]);
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


///
/// converts 2D vector into 1D array.  Transfers only information for
/// those indexes where there are parents in the network.
/// @returns false when constructed sequence too long
///
bool GADiscrimBayes::constructTransferSeq(vector<vector<int> > model, short* seq,
	int seqSize){

	int seqIndex=0;
	size_t i=0;
	for(; i<model.size(); i++){
		if(model[i].empty())
			continue;
		seq[seqIndex++]=i;
		if(seqIndex == seqSize)
			return false;
		for(size_t j=0; j<model[i].size(); j++){
			seq[seqIndex++]=model[i][j];
			if(seqIndex == seqSize)
				return false;
		}
		seq[seqIndex++]=-1;
		if(seqIndex == seqSize)
			return false;
	}
// if(myRank==1){
// cout << "SENT" << endl;
// for(int x=0; x<seqIndex; x++){
// 	cout << seq[x] << " ";
// }
// cout << endl;
// for(size_t i=0; i<model.size(); i++){
// 	cout << i << "=>";
// 	for(size_t j=0; j<model[i].size(); j++){
// 			cout << model[i][j] << ",";
// 	}
// 	cout << "|";
// }
// cout << "-=======" << endl;
// }
	for(;seqIndex<seqSize;seqIndex++){
		seq[seqIndex]=-1;
	}



	return true;
}


void GADiscrimBayes::constructModelVec(vector<vector<int> >& model, short* modelSeq,
	int seqSize){

	vector<int> empty;
	model.assign(varList.size(), empty);

	for(size_t i=0; i<seqSize; i++){
// 		if(modelSeq[seqSize] == -1){
// 			varIndex++;
// 		}
		int currVar=modelSeq[i];
		if(currVar != -1){
			i++;
			while(modelSeq[i] != -1 and i<seqSize){
				model[currVar].push_back(modelSeq[i]);
				i++;
			}
		}
	}


// cout << "RECEIVED" << endl;
// for(int x=0; x<seqSize; x++){
// 	cout << modelSeq[x] << " ";
// }
// cout << endl;
// for(size_t i=0; i<model.size(); i++){
// 	cout << i << "=>";
// 	for(size_t j=0; j<model[i].size(); j++){
// 			cout << model[i][j] << ",";
// 	}
// 	cout << "|";
// }
// cout << "-=======" << endl;


}


void GADiscrimBayes::gatherModelInformation(map<vector<vector<int> >, ModelScores>& models){
	modelMPI empty;
	modelMPI* modelSend = new modelMPI[popSize];

	int i=0;
	for(map<vector<vector<int> >, ModelScores>::iterator iter=models.begin(); iter != models.end();
		++iter){
		modelSend[i].score=iter->second.score;
		modelSend[i].count=iter->second.count;
		if(!constructTransferSeq(iter->first,modelSend[i].modelSeq, MAXMODELTRANSFER)){
			modelSend[i].score=GAFunct::getWorstScore();
		}
		i++;
	}

	modelMPI* modelRecv=NULL;
	int recvSize = popSize * totalNodes;
	if(myRank==0){
		modelRecv = new modelMPI[recvSize];
	}

	MPI_Gather(modelSend, sizeof(empty)*popSize, MPI_BYTE, modelRecv, sizeof(empty)*popSize, MPI_BYTE, 0, MPI_COMM_WORLD);

	// merge all
	if(myRank==0){
		models.clear();
		for(int modIndex=0; modIndex < recvSize; modIndex++){
			vector<vector<int> > modVec;
			ModelScores modScore;
			constructModelVec(modVec, modelRecv[modIndex].modelSeq,MAXMODELTRANSFER);
			if(modelRecv[modIndex].count > 0){
				if(models.find(modVec) == models.end()){
					models[modVec].count = modelRecv[modIndex].count;
					models[modVec].score = modelRecv[modIndex].score;
				}
				else{
					models[modVec].count += modelRecv[modIndex].count;
				}
			}
// cout << "rank=" << myRank << " finished merging models " << endl;
		}
	}

	delete [] modelSend;
	if(myRank==0){
		delete [] modelRecv;
	}
}

#endif

