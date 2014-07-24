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
//Algorithm.cpp

#include "Algorithm.h"
#include "GEObjective.h"
// #include <set>


Algorithm::Algorithm(){
		fitnessName=" ";
		initializeParams();
}

Algorithm::~Algorithm(){
		for(unsigned int i=0; i<logs.size(); i++){
			delete logs[i];
		}
}


void Algorithm::setParams(AlgorithmParams& algParam, int numExchanges, int numGenos, 
	int numContin, vector<unsigned int>& excludedGenos, vector<unsigned int>& excludedContins){
		map<string, string>::iterator mapIter;
		vector<string> tokens;
		
		numGenotypes = numGenos;
		numContinuous = numContin;
		
		for(mapIter = algParam.params.begin(); mapIter != algParam.params.end(); 
			mapIter++){     
				switch(algParamMap[mapIter->first]){
// 						case noMatchParam:
// 								throw AthenaExcept("No match for parameter " + mapIter->first +
// 												" in Algorithm");
// 								break;
						case minSizeParam:
								minSize = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case maxSizeParam:
								maxSize = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								#ifdef PARALLEL
									if(maxSize > MAX_GENOME_SIZE)
										throw AthenaExcept(mapIter->first + " can be no more than " + 
											Stringmanip::numberToString<int>(MAX_GENOME_SIZE));
								#endif
								break;
						case tailRatioParam:
								tailRatio = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case growRateParam:
								growRate = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case maxDepthParam:
								maxDepth = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case tailSizeParam:
								tailSize = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case sensibleInitParam:
								sensibleInit = Stringmanip::check_true_false(mapIter->second);
								break;
						case popSizeParam:
								popSize = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case probCrossParam:
								probCross = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case bioModelSelection:
								if(bioModelSelectionMap.find(Stringmanip::to_upper(mapIter->second)) == 
										bioModelSelectionMap.end())
									throw AthenaExcept("No match for bio model selection type " + mapIter->second);
								else
									biofilterSelectorType = bioModelSelectionMap[Stringmanip::to_upper(mapIter->second)];
								break;
						case probMutParam:
								probMut = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case gramFileParam:
								grammarFile = mapIter->second;
								break;
						case stepSizeParam:
								stepSize = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case calcType:
								calculatorName = Stringmanip::to_upper(mapIter->second);
								break;        
						case useEffectiveXO:
								effectiveXO = Stringmanip::check_true_false(mapIter->second);
								break;
						case useAllSnps:
								useAllVars = Stringmanip::check_true_false(mapIter->second);
								break;
						case useAllCovariates:
								useAllCovars = Stringmanip::check_true_false(mapIter->second);
								break;
						case fitGoal:
								fitnessGoal = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case gaSelection:
								if(gaSelectorMap.find(Stringmanip::to_upper(mapIter->second)) ==
									gaSelectorMap.end())
									throw AthenaExcept("No match for GA selection type " + mapIter->second);
								else
									gaSelector = gaSelectorMap[Stringmanip::to_upper(mapIter->second)];
									break;
#ifdef ATHENA_BLOAT_CONTROL
						case doubleTournF:
								doubleTourneyF = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case doubleTournD:
								doubleTourneyD = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case doubleTournFitFirst:
								fitFirst = Stringmanip::check_true_false(mapIter->second);
								break;
#endif
						case blockCrossGens:
								ngensBlockCross = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case bestCVThresh:
						    bestCVThreshold = Stringmanip::stringToNumber<int>(mapIter->second);
						    break;
						case bestCorrThresh:
						    bestCorrThreshold = Stringmanip::stringToNumber<double>(mapIter->second);
						    break;
#ifdef ATHENA_BLOAT_CONTROL
						case prunePlantFract:
								break;
#endif
// 						default:
// 								throw AthenaExcept("No match for parameter " + mapIter->first +
// 												" in Algorithm GENN");               
				}
		}
		
		if(stepSize > numGenerations){
				stepSize = numGenerations;
		}
		numGenerations = stepSize * numExchanges;
		
		setFitnessName(calculatorName);
}


///
/// Initializes parameters for basic run of algorithm.  These
/// can be modified using the set_params function
///
void Algorithm::initializeParams(){
	 
	 myRank = 0;
	 totalNodes = 1;
	 minSize = 50; // Minimum size for Random Initialization
	 maxSize = 200; // Maximum size for Random Initialization

	 tailRatio = 0.0;
	 growRate = 0.5;

	 sensibleInit = false;
	 maxDepth = 10;
	 tailSize = 0;

	 grammarFile = "";
	 wrapEvents = 0;

	 popSize = 100;
	 numGenerations = 100;
	 probCross = 0.9;
	 probMut = 0.01;
	 fitnessGoal = 1.0;
	 bestCorrThreshold = 0.8;
	 bestCVThreshold = 2;
	 
	 stepSize = 100;

	 effectiveXO = false;
	 randSeed = 7;
 
	 calculatorName = "BALANCEDACC";
	 
	 // establish map for parameters
	 algParamMap["MINSIZE"] = minSizeParam;
	 algParamMap["MAXSIZE"] = maxSizeParam;
	 algParamMap["TAILRATIO"] = tailRatioParam;
	 algParamMap["GROWRATE"] = growRateParam;
	 algParamMap["MAXDEPTH"] = maxDepthParam;
	 algParamMap["TAILSIZE"] = tailSizeParam;
	 algParamMap["SENSIBLEINIT"] = sensibleInitParam;
	 algParamMap["POPSIZE"] = popSizeParam;
	 algParamMap["PROBCROSS"] = probCrossParam;
	 algParamMap["PROBMUT"] = probMutParam;
	 algParamMap["GRAMMARFILE"] = gramFileParam;
	 algParamMap["GENSPERSTEP"] = stepSizeParam;
	 algParamMap["CALCTYPE"] = calcType;
	 algParamMap["EFFECTIVEXO"] = useEffectiveXO;
	 algParamMap["INCLUDEALLSNPS"] = useAllSnps;
	 algParamMap["BLOCKCROSSGENS"] = blockCrossGens;
	 algParamMap["GASELECTION"] = gaSelection;
	 algParamMap["FITNESSGOAL"] = fitGoal;
	 algParamMap["BESTCVTHRESH"] = bestCVThresh;
	 algParamMap["BESTCORRTHRESH"] = bestCorrThresh;
	 algParamMap["BIOMODELSELECTION"] = bioModelSelection;
#ifdef ATHENA_BLOAT_CONTROL
	 algParamMap["DOUBTOURNF"] = doubleTournF;
	 algParamMap["DOUBTOURND"] = doubleTournD;
	 algParamMap["DOUBTOURNFITFIRST"] = doubleTournFitFirst;
	 algParamMap["PRUNEPLANT"] = prunePlantFract;
#endif
	 
#ifdef ATHENA_BLOAT_CONTROL
	 gaSelectorMap["DOUBLE"] = DoubleTournamentSelection;
#endif
	 gaSelectorMap["ROULETTE"] = RouletteWheelSelection;
	 gaSelectorMap["PARETO"] = ParetoFrontSelection;
	 gaSelectorMap["PARETORANK"] = ParetoRankSelection;
	 
// 	 biofilterSelectorType = orderedSelect;
#ifdef ATHENA_BLOAT_CONTROL
	 pruneAndPlantFract = 0.0;
#endif

	 useAllVars = false;
	 numGenotypes = 0;
	 numContinuous = 0;
	 // default is to add any additional variables to migration instead of replacing
	 ngensBlockCross = 0;
	 
	 ga = NULL;
	 maxBest = true;
	 
	 gaSelector = RouletteWheelSelection;
#ifdef ATHENA_BLOAT_CONTROL
	 fitFirst = false;
	 doubleTourneyF = 7;
	 doubleTourneyD = 1.4;
#endif
	 
}


///
/// Frees memory
///
void Algorithm::freeMemory(){
		if(ga != NULL){
				delete ga;
				ga = NULL;
		}
}


void Algorithm::setInitParams(){
	 InitGEgenome::setMinSize(minSize);
	 InitGEgenome::setMaxSize(maxSize);
	 InitGEgenome::setMapper(&mapper);
}


void Algorithm::expandVariables(){
	// adjust grammar when ott dummy encoding has been used
	// or for shorthand ways of specifying variables
	adjuster.readGrammarFile(grammarFile);
	if(useAllVars){
		adjuster.includeAllVars(numGenotypes, numContinuous);
	}
	else{
		adjuster.expandVariables();
		if(dummyEncoded)
			adjuster.doubleGenotypeGrammar();
		}
}

///
/// Exclude variables from grammar so they will not appear in any model
///
void Algorithm::excludeVariables(vector<unsigned int>& excludedGenos, 
			vector<unsigned int>& excludedContins){
	
	std::set<string> varSet;
	string name;
	// construct a map of all excluded variables
	for(vector<unsigned int>::iterator iter=excludedContins.begin();
		 iter!=excludedContins.end(); ++iter){
		name = "C" + Stringmanip::numberToString(*iter+1);
		varSet.insert(name);
	}
	
	for(vector<unsigned int>::iterator iter=excludedGenos.begin();
		 iter!=excludedGenos.end(); ++iter){
		if(dummyEncoded){
			unsigned int num = (*iter*2)+1;
			name = "G" + Stringmanip::numberToString(num);
			varSet.insert(name);
			name = "G" + Stringmanip::numberToString(num+1);
			varSet.insert(name);
		}
		else{
			name = "G" + Stringmanip::numberToString(*iter+1);
			varSet.insert(name);
		}
	}
	
	adjuster.excludeVariables(varSet);
	
}


/// Sets parameters for use with GAlib
/// @param alg_params AlgorithmParams
/// @throws AthenaExcept on error
///
void Algorithm::setGAParams(vector<unsigned int>& excludedGenos, 
	vector<unsigned int>& excludedContins){   
	GARandomSeed(randSeed);
	srand(randSeed);   
	// first free ga memory if already run
	freeMemory();
		
	//Initialize the GE mapper
	//Set maximum number of wrapping events per mapping
	mapper.setMaxWraps(wrapEvents);
	expandVariables();
	adjuster.setMapper(mapper);

	setMapperPrefs(mapper);

	 // add mapper to objective function
	 GEObjective::setMapper(&mapper);
	 if(logTypeSelected == LogNone)
		 GEObjective::addLogging(false);
	 else
		 GEObjective::addLogging(true);

	setInitParams();
}


void Algorithm::setMapperPrefs(AthenaGrammarSI& athenaMapper){
		// Mapper settings.
	  athenaMapper.setGenotypeMaxCodonValue(INT_MAX);
	  athenaMapper.setPopSize(popSize);
	  athenaMapper.setGrow(growRate);
	  athenaMapper.setMaxDepth(maxDepth);
	  
	  if(tailSize){
		  athenaMapper.setTailSize(tailSize);
	  }
	  else{
		  athenaMapper.setTailRatio(tailRatio);
	  }
	  athenaMapper.setRestrictRule("<v>");
	  
	  // set up reverse grammar that can be used with resetting genome based
	  // on optimization (such as backpropagation)
		athenaMapper.constructReverseGrammar();
}



///
/// Retrieves the models from BioFilter and stores the information in the algorithm
/// @param filename File with biofilter models
/// @param bioFileType Type of biological filter file (BINARY or TEXT)
/// @param holder Dataholder that contains map file info to convert model names
/// into indexes in the 
///
void Algorithm::getBioModels(std::string filename, std::string bioFileType, data_manage::Dataholder* holder){
	BioFilterModelCollection collection(filename, 100000, bioFileType);
	setBioModels(collection, holder);
}



///
/// Fills biomodel collection from archive files 
/// @param genegeneFile genegene filename
/// @param archiveFile arcchive filename
/// @param holder Dataholder that contains map file info to convert model names
/// into indexes
///
void Algorithm::getBioModelsArchive(string genegeneFile, string archiveFile, data_manage::Dataholder* holder){
	BioFilterModelCollection collection(genegeneFile, archiveFile, 10000);
	setBioModels(collection, holder);
}



///
/// Sets the bio models based on collection passed
///
void Algorithm::setBioModels(BioFilterModelCollection& collection, data_manage::Dataholder* holder){

	vector<BioModel> bioModels;
	int numModelsNeeded = (unsigned int)(initBioFract * popSize);

	if(biofilterSelectorType == rouletteSelect){  
		for(int i=0; i<numModelsNeeded; i++){
			bioModels.push_back(collection.getRandomModel());
		}
	}
	else{
		// set starting point for this algorithm -- in parallel mode this will give each 
		// node a different set
		collection.setStartModel(myRank * numModelsNeeded);
		for(int i=0; i<numModelsNeeded; i++){
			bioModels.push_back(collection.getNextModel());
		}
	}

	mapper.clearModels();

	// from models find the variables to use and store in mapper for use in initializing dataset
	for(vector<BioModel>::iterator iter=bioModels.begin(); iter != bioModels.end(); iter++){
		vector<int> indexes;
		int index;
		for(unsigned int i=0; i<iter->idString.size(); i++){
			index = holder->getGenoIndex(iter->idString[i]);
			indexes.push_back(index);
		}
		mapper.addModel(indexes);
	}
	  mapper.setModelCodons(dummyEncoded);
}

