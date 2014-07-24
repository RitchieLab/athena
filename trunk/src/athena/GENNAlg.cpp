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

#include "GEObjective.h"
#include "GENNAlg.h"
#include "Terminals.h"
#include "GENNGrammarAdjuster.h"
#include "ModelLogParser.h"
#include "BestModelSelector.h"
#include "GAParetoSelector.h"
#include "GAParetoRankSelector.h"
#include "SumFileReader.h"
#include <ga/ga.h>
#include <ctime>
#include <set>
#include <algorithm>

///
/// Constructor
///
GENNAlg::GENNAlg(){
	initializeParams();
}


///
/// Destructor
///
GENNAlg::~GENNAlg(){
		if(geLog != NULL){
				delete geLog;
				geLog=NULL;
		}
		freeMemory();
}


///
/// Frees memory
///
// void GENNAlg::freeMemory(){
// 		if(ga != NULL){
// 				delete ga;
// 				ga = NULL;
// 		}
// }

///
/// Set current Dataset for running algorithm
/// @param new_set Dataset
/// 
void GENNAlg::setDataset(Dataset* newSet){
	 set = newSet;
	 set->calcSSTotal();
	 GEObjective::setDataset(set);
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
void GENNAlg::setParams(AlgorithmParams& algParam, int numExchanges, int numGenos, 
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
												" in Algorithm GENN");
								break;
						case requireAll:
								requireAllVars = Stringmanip::check_true_false(mapIter->second);
								break;
						case requireAllOnce:
								requireAllVarsOnce = Stringmanip::check_true_false(mapIter->second);
								break;
						case bioInitFract:
								initBioFract = Stringmanip::stringToNumber<double>(mapIter->second);     
								break;
						case restrictVarGens:
								ngensVarRestrict = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case bioModelSelection:
								if(bioModelSelectionMap.find(Stringmanip::to_upper(mapIter->second)) == 
										bioModelSelectionMap.end())
									throw AthenaExcept("No match for bio model selection type " + mapIter->second);
								else
									biofilterSelectorType = bioModelSelectionMap[Stringmanip::to_upper(mapIter->second)];
								break;
						case blockCrossGens:
								ngensBlockCross = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case resetVarsAtMigration:
								resetRestrictedAtMigration = Stringmanip::check_true_false(mapIter->second);
								break;
						case bpstart:
								bpFirstGen = Stringmanip::stringToNumber<int>(mapIter->second);
								break;
						case bpfreq:    
								bpFreqGen = Stringmanip::stringToNumber<int>(mapIter->second);
								break;
						case bestCVThresh:
						    bestCVThreshold = Stringmanip::stringToNumber<int>(mapIter->second);
						    break;
						case bestCorrThresh:
						    bestCorrThreshold = Stringmanip::stringToNumber<double>(mapIter->second);
						    break;
						case constantSpan:
								tokens = Stringmanip::split(mapIter->second, ':');
								minConstant = Stringmanip::stringToNumber<float>(tokens[0]);
								maxConstant = Stringmanip::stringToNumber<float>(tokens[1]);
								constantInterval = Stringmanip::stringToNumber<float>(tokens[2]);
								simpleConstants=true;
								break;
#ifdef ATHENA_BLOAT_CONTROL
						case prunePlantFract:
#endif
						default:
						  if(paramMap.find(mapIter->first) == paramMap.end())
								throw AthenaExcept("No match for parameter " + mapIter->first +
												" in Algorithm GENN");               
				}
		}
		
		// first optimization of backpropagation 
		bpNextOpt = bpFirstGen;
			 
		setGAParams(excludedGenos, excludedContins);
		
}


///
/// Sets values in main configuration to defaults needed by 
/// GENN Algorithm
/// @param configuration Config
///
void GENNAlg::setConfigDefaults(Config& configuration, AlgorithmParams& algParam){
  if(algParam.params["CALCTYPE"].compare("RSQUARED")==0){
    configuration.setStatusAdjust("MINMAX");
  }
}


///
/// Initializes parameters for basic run of algorithm.  These
/// can be modified using the set_params function
///
void GENNAlg::initializeParams(){
	 
	 initBioFract = 0.0;
 
	 calculatorName = "BALANCEDACC";
	 
	 // establish map for parameters
	 paramMap["REQUIREALLVARS"] = requireAll;
	 paramMap["REQUIREALLONCE"] = requireAllOnce;
	 paramMap["BIOFILTERFRACT"] = bioInitFract;
	 paramMap["NUMGENSRESTRICTVARS"] = restrictVarGens;

	 paramMap["RESETVARSMIGRATION"] = resetVarsAtMigration;
	 paramMap["BACKPROPFREQ"] = bpfreq;
	 paramMap["BACKPROPSTART"] = bpstart;
	 paramMap["CONSTANTS"]=constantSpan;
	 
	 bioModelSelectionMap["ROULETTE"] = rouletteSelect;
	 bioModelSelectionMap["ORDERED"] = orderedSelect;
	 
	 biofilterSelectorType = orderedSelect;

	 useAllVars = false;
	 requireAllVars = false;
	 requireAllVarsOnce = false;
	 ngensVarRestrict = 0;
	 // default is to add any additional variables to migration instead of replacing
	 resetRestrictedAtMigration = false;
	 ngensBlockCross = 0;
	 
	 geLog = NULL;
	 
	 minConstant = -1.0;
	 maxConstant = 1.0;
	 constantInterval = 0.01;
	 simpleConstants = false;
	 
	 // these control when first backpropagation occurs and how often thereafter
	 bpFirstGen = -1;
	 bpFreqGen = 0;
	 bpNextOpt = -1;
	 
	 logTypeSelected = LogNone;
}



///
/// Sets random seed 
/// @param seed 
///
void GENNAlg::setRand(unsigned int seed){
	GARandomSeed(seed);
	srand(seed);
}


///
/// Sets parameters for use with GAlib
/// @param alg_params AlgorithmParams
/// @throws AthenaExcept on error
///
void GENNAlg::setGAParams(vector<unsigned int>& excludedGenos, 
			vector<unsigned int>& excludedContins){   
    
		GARandomSeed(randSeed);
		srand(randSeed);   
		// first free ga memory if already run
		freeMemory();
		
		//Initialize the GE mapper
	 //Set maximum number of wrapping events per mapping
	 mapper.setMaxWraps(wrapEvents);
	 expandVariables();
		 
	  // set constants
		if(simpleConstants)
			adjuster.setConstants(minConstant, maxConstant, constantInterval); 
		 
		adjuster.setMapper(mapper);

		setMapperPrefs(mapper);
	  
	  if(!requireAllVars && !requireAllVarsOnce)
			GEObjective::setSolutionType("NN", calculatorName);
		else if(requireAllVars){
			vector<string> varStrings = adjuster.getVariables();
			GEObjective::setSolutionType("NNALL", calculatorName, varStrings);
		}
		else if(requireAllVarsOnce){
			vector<string> varStrings = adjuster.getVariables();
			GEObjective::setSolutionType("NNONCE", calculatorName, varStrings);
		}

	if(simpleConstants)
		GEObjective::addConstants(adjuster.getIncludedConstants());

	 if(calculatorName.compare("RSQUARED")==0 || calculatorName.compare("MEANABSOLUTE")==0){
//    if(calculatorName.compare("RSQUARED")==0){
		 pop.setConvertScores(true);
	 }

	 // add mapper to objective function
	 GEObjective::setMapper(&mapper);
	 if(logTypeSelected == LogNone)
		 GEObjective::addLogging(false);
	 else
		 GEObjective::addLogging(true);
		 
    setInitParams();
}



///
/// sets mapper preferences
/// @param athenaMapper AthenaGrammarSI
///
// void GENNAlg::setMapperPrefs(AthenaGrammarSI& athenaMapper){
// 		// Mapper settings.
// 	  athenaMapper.setGenotypeMaxCodonValue(INT_MAX);
// 	  athenaMapper.setPopSize(popSize);
// 	  athenaMapper.setGrow(growRate);
// 	  athenaMapper.setMaxDepth(maxDepth);
// 	  
// 	  if(tailSize){
// 		  athenaMapper.setTailSize(tailSize);
// 	  }
// 	  else{
// 		  athenaMapper.setTailRatio(tailRatio);
// 	  }
// 	  athenaMapper.setRestrictRule("<v>");
// 	  
// 	  // set up reverse grammar that can be used with resetting genome based
// 	  // on optimization (such as backpropagation)
// 		athenaMapper.constructReverseGrammar();
// }



///
/// Compare bet individual against goal fitness.  If met or exceeded, then
/// return 1.  If not met, return 0.
/// @return 1 if goal reached, and 0 otherwise
///
int GENNAlg::reachedGoal(){
	GE1DArrayGenome& bestGenome = (GE1DArrayGenome&) ga->population().individual(0);
	
	bool goalDone = 0;
	double fitness;
	// when running r-squared need to convert the mean-squared error into
	// r-squared result
	if(pop.getConvertScores()){
		NNSolution* bestSolution = convertGenome(ga->population().individual(0));
		bestSolution->adjustScoreOut(set, calculatorName);
		fitness = double(bestSolution->fitness());
		delete bestSolution;
	}
	else{
		fitness = bestGenome.fitness();
	}
	
	if(fitness >= fitnessGoal){
		goalDone=1;
	}
	return goalDone;
}


///
/// Establishes the algorithm for the run based on the parameters
/// set in this algorithm.
///
void GENNAlg::run(){
		while(!ga->done()){
				ga->step();
		}
}


///
/// Runs an indicated interval for the algorithm. This interval is set
/// in the stepSize variable.  
///
int GENNAlg::step(){
		int completed =0;
		GE1DArrayGenome::myRank = myRank;
		for(unsigned int i=0; i < stepSize; i++){
				if(!ga->done()){
						if(logTypeSelected!=LogNone){
								geLog->addGeneration();
						}
				 ga->step();
				 restrictStepsDone++;
					if(ngensVarRestrict && restrictStepsDone == ngensVarRestrict){
					 // have to convert the networks back to original mapping
					 GEObjective::setMapper(&mapper);
					 GE1DArrayGenome::setMapper(&mapper);
					 convertNetworks(restrictMapper, mapper);
				 }
				 
				 // check for need to change crossovers
				 if(ngensBlockCross && restrictStepsDone == ngensBlockCross){
						resetCrossover();
				 }
 
				 // check to see if the back propagation should be done at this generation
				 if(int(restrictStepsDone) == bpNextOpt && bpFirstGen >= 0){
					runBackPropagation();
					bpNextOpt += bpFreqGen;
				 }
				 fillLog();
				 // stop analysis when best model reaches or exceeds the fitness goal
				 completed = reachedGoal();
					#ifdef PARALLEL
							completed = popMigrator.nodesCompleted(completed);
					#endif
					if(completed){
						break;
					}
				}
		}

		#ifdef PARALLEL
			// if restricted variables has been used and are still in effect
			// need to convert all networks back to original grammar and exchange
			// then construct new grammar restricted to only variables in the 
			// new population and continue to process
			if(ngensVarRestrict && restrictStepsDone < ngensVarRestrict){
				// have to convert the networks back to original mapping
				 GEObjective::setMapper(&mapper);
				 GE1DArrayGenome::setMapper(&mapper);
				 convertNetworks(restrictMapper, mapper);
			} 

			popMigrator.sendAndReceiveStruct(totalNodes, myRank, ga);

			// when running in parallel transfer around populations
			if(ngensVarRestrict && restrictStepsDone < ngensVarRestrict){
				// after transfer construct new grammar and set for use
				setRestrictedGrammar(resetRestrictedAtMigration);
			}
		#endif
		
		// only need to fill population at this point not at end of each generation
		fillPopulation();
		return completed;
}


///
/// Run dataset against models contained in previously generated summary file
/// @param sumFile ATHENA summary file
///
vector<Solution*> GENNAlg::runValidation(std::string sumFile){
	SumFileReader reader;
	reader.readSumFile(sumFile);
	vector<Solution*> models = reader.getModelPopulation();
	for(size_t i=0; i<models.size(); i++){
		GEObjective::calcFitness(models[i]);
		models[i]->setAdditionalOutput(GEObjective::calcAdditionalFinalOutput(models[i]));
		models[i]->testVal(0);
	}
	return models;
}

///
/// Writes individual output information to stringstreams
/// @param indss
/// @param models
///
void GENNAlg::validationIndOutput(vector<std::stringstream*>& indss, vector<Solution*>& models){
  for(size_t i=0; i<indss.size(); i++){
    GEObjective::calcFitnessOut(models[i], *(indss[i]));
  }
}


///
/// Fills the algorithm log with data from the population
///
void GENNAlg::fillLog(){
		
		// don't fill log if not being kept
		if(logTypeSelected==LogNone){
			return;
		}


		// add to detailed log files
		fillPopulation();


		unsigned int numInds = ga->population().size();	
		float worstScore = GEObjective::getWorstScore();
		
		for(unsigned int currInd =0; currInd < numInds; currInd++){

			GE1DArrayGenome& genome = (GE1DArrayGenome&)(ga->population().individual(currInd));
			if(worstScore != genome.score()){
				geLog->addNetwork();
				geLog->addNNSize(genome.getEffectiveSize());
				geLog->addNNDepth(genome.getDepth());
				if(!pop.getConvertScores()){
					geLog->addFitness(genome.score(), genome.getGenos(), genome.getCovars());
				}
				else{  // add converted score to log file
					geLog->addFitness(pop[0]->adjustScoreOut(genome.score(), genome.getNumIndsEvaluated(),
						genome.getSSTotal(), getFitnessName()), genome.getGenos(), genome.getCovars());
				}
				geLog->addNumGenos(genome.getNumGenes());
				geLog->addNumCovars(genome.getNumCovars());
				geLog->addEpochs(genome.getNumEpochsTrained());
				// zero out epochs`
				genome.setNumEpochsTrained(0);
			}
		}
		geLog->completeGen();

		#ifdef PARALLEL
			geLog->sendReceiveLogs(totalNodes, myRank);
		#endif
		
		// output log information -- appended to existing log files
		if(myRank==0){
				writeLog();
		}

		if(logTypeSelected==LogVariables){
			NNSolution * solution;
			int nSolutions = pop.numSolutions();
			modelLog->addGeneration(geLog->getCurrentGen());
			for(int i=0; i<nSolutions; i++){
				solution = (NNSolution*)pop[i];
				modelLog->writeVariables(*solution, dummyEncoded);
			}
			return;
		}

		if(logTypeSelected != LogOverview){
			NNSolution * solution;
			int nSolutions = pop.numSolutions();
			for(int i=0; i<nSolutions; i++){
				solution = (NNSolution*)pop[i];
				modelLog->writeSolution(*solution, geLog->getCurrentGen(), i+1);
			}
		}
	
}


///
/// Gets log and saves in vector
///
void GENNAlg::saveLog(){

	// need to get all slaves logs when running parallel
	#ifdef PARALLEL
		if(myRank==0)
			geLog->receiveLogs(totalNodes);
		else
			geLog->sendLog();
	#endif

	logs.push_back(geLog);
}


///
/// Return additional output names (if any such as AUC)
/// 
vector<string> GENNAlg::getAdditionalOutputNames(){
  return GEObjective::getAdditionalOutputNames();
}


///
/// Calculate and return additional output (like AUC) for best model
///
void GENNAlg::getAdditionalFinalOutput(Dataset* set){
	GEObjective::setDataset(set);
	unsigned int numInds = ga->population().size();
  for(unsigned int currInd = 0; currInd < numInds; currInd++){
		pop[currInd]->setAdditionalOutput(GEObjective::getAdditionalFinalOutput(ga->population().individual(currInd)));
	}
	
}


///
/// Close the log files
///
void GENNAlg::closeLog(){
	if(logTypeSelected != LogNone){
		modelLog->closeLog();
	}
}


///
/// Starts fresh log 
///
void GENNAlg::startLog(int numSnps){
	if(geLog != NULL){
		delete geLog;
		if(modelLog != NULL)
	    delete modelLog;
		geLog=NULL;
		modelLog=NULL;
	}
	geLog = new NNLog(numSnps);
	geLog->setMaxBest(GEObjective::logMaxBest());
}



///
/// Outputs current population of genome 
/// 
void GENNAlg::outputAlgInds(){

	 unsigned int numInds = ga->population().size();
	 
		for(unsigned int currInd = 0; currInd < numInds; currInd++){
			GE1DArrayGenome& genome = (GE1DArrayGenome&)(ga->population().individual(currInd));
			mapper.setGenotype(genome);  
			Phenotype const *phenotype=mapper.getPhenotype();
			unsigned int phenoSize=(*phenotype).size();
			vector<string> symbols(phenoSize, "");
		 
			for(unsigned int i=0; i<phenoSize; ++i){
				 symbols[i] = *((*phenotype)[i]);
			}
			cout << "IND in ga: " << currInd+1 << " score: " << genome.score() << endl;
		}
}



///
/// Converts an individual from population into a Solution and returns it
/// @param GAGenome ind
/// @return NNSolution*
///
NNSolution* GENNAlg::convertGenome(GAGenome& ind){
	
	GE1DArrayGenome& genome = (GE1DArrayGenome&) ind;
	NNSolution* sol = (NNSolution*)GEObjective::getBlankSolution();
	mapper.setGenotype(genome);
	Phenotype const *phenotype=mapper.getPhenotype();
	unsigned int phenoSize=(*phenotype).size();
	vector<string> symbols(phenoSize, "");			
	for(unsigned int i=0; i<phenoSize; ++i){
		symbols[i] = *((*phenotype)[i]);
	}

	sol->setSymbols(symbols);
	sol->fitness(genome.score());
	sol->testVal(genome.getTestValue());
	sol->setGramDepth(genome.getGramDepth());
	sol->setNNDepth(genome.getDepth());
	sol->setComplexity(genome.getComplexity());
	return sol;
}



///
/// Fills the population after a step to allow for exchange of models with 
/// other algorithms.
/// 
void GENNAlg::fillPopulation(){
		unsigned int numInds = ga->population().size();
		// clear solution population
		pop.clear();
		
		for(unsigned int currInd = 0; currInd < numInds; currInd++){
				pop.insert(convertGenome(ga->population().individual(currInd)));
		}
}



///
/// Establishes test data set values for the population
/// @param test_set Dataset containing individuals for the test set
///
void GENNAlg::testSolution(Dataset* testSet, int nproc){
	 
	 // use first dataset
	 set = testSet;
	 GEObjective::setDataset(set);
  int n = ga->population().size();
 
  for(int i=0; i<n; i++){
		GE1DArrayGenome bestGenome = (GE1DArrayGenome&)ga->population().best(i);
		float testScore = GEObjective::GEObjectiveFunc(bestGenome);
		((GE1DArrayGenome&)ga->population().best(i)).setTestValue(testScore);
		pop[i]->testVal(testScore);
	}
}



///
/// Evaluates and outputs dataset to indicated output stream.
/// @param set Dataset
/// @param os ostream to write to
/// @param model index of model to output
/// 
void GENNAlg::outputIndEvals(Dataset* set, ostream& os, int model){
	GEObjective::setDataset(set);
	GE1DArrayGenome bestGenome = (GE1DArrayGenome&)ga->population().best(model);
	GEObjective::GEObjectiveFuncOut(bestGenome, os);
}



///
/// Returns the snps and covariates present in the network
///
vector<string> GENNAlg::getBestVariables(){
	vector<string> symbols;
	return symbols;
}


///
/// Sets all genomes in a population to either effective crossover or one point
///
void GENNAlg::resetCrossover(){
		if(effectiveXO){
			ga->crossover(GE1DArrayGenome::effCrossover);
		}
		else{
			ga->crossover(GE1DArrayGenome::OnePointCrossover);
		}
}



///
/// Initializes the algorithm
///
void GENNAlg::initialize(){
		restrictStepsDone=0;
		
		// free ga if already established
		freeMemory();
		
		// establish genome type and assign intialization and objective functions
		GE1DArrayGenome genome(50);

		// need to add objective function
		genome.evaluator(GEObjective::GEObjectiveFunc);
		// set function for establishing basic values without full evaluation
		genome.establishinator(GEObjective::GEObjectiveInit);
		
		if(sensibleInit){
				genome.initializer(InitGEgenome::initFuncSI);
		}
		else{
				genome.initializer(InitGEgenome::initFuncRandom);
		}
		
		// assign crossover type
		if(ngensBlockCross > 0)
			genome.crossover(GE1DArrayGenome::blockCrossover);
		else if(!effectiveXO){
				genome.crossover(GE1DArrayGenome::OnePointCrossover);
		}
		else{
				genome.crossover(GE1DArrayGenome::effCrossover);
		}
		
		// use point mutator for GEListGenome
		genome.mutator(GE1DArrayGenome::codonMutator);
		
		// controls allowable sizes
		genome.resizeBehaviour(1, maxSize);
		GEObjective::setMaxGenomeSize(maxSize);
		
		// Set up the algorithm
		ga = new GASimpleGA(genome);
		
#ifdef ATHENA_BLOAT_CONTROL
		ga->setPrunePlant(pruneAndPlantFract);
		ga->pruneplant(GE1DArrayGenome::pruneAndPlant);
#endif

		// Set selector type
		GASelectionScheme* selector;
		switch(gaSelector){
			case NoMatchSelector:
			case RouletteWheelSelection:
				selector = new GARouletteWheelSelector;
				break;
			case ParetoFrontSelection:
 			  selector = new GAParetoSelector;
			  break;
			case ParetoRankSelection:
 			  selector = new GAParetoRankSelector;
			  break;			  
#ifdef ATHENA_BLOAT_CONTROL
			case DoubleTournamentSelection:
				selector = new GADoubleTournamentSelector(doubleTourneyD, doubleTourneyF, fitFirst);
				break;
#endif
		};
		
		ga->selector(*selector);
		// safe to delete selector as algorithm clones it
		delete selector; 
		// scaling
		ga->scaling(GANoScaling());
		// individuals in population
		ga->populationSize(popSize);
		ga->pMutation(probMut);        
		ga->pCrossover(probCross);
		ga->nGenerations(numGenerations);

		if(GEObjective::maxBest()){
			ga->maximize();
			maxBest = true;
		}
		else{
			ga->minimize();
			pop.sortAscending();
			maxBest = false;
			geLog->setMaxBest(GEObjective::logMaxBest());
		}
	 
		mapper.resetGrammarModels();
		
		mapper.setVariableCodonMap();
		GE1DArrayGenome::setMapper(&mapper);
		ga->initialize();

		// when restricted variables are requested 
		// fill new mapper with grammar with variables from original mapper
		if(ngensVarRestrict > 0){
			setRestrictedGrammar(true);
		}
		
		// first optimization of backpropagation 
		bpNextOpt = bpFirstGen;

		// run optimization after initialization when indicated
		if(bpNextOpt == 0){
			runBackPropagation();
			bpNextOpt += bpFreqGen;
		}
		fillPopulation();
		fillLog();
}



///
/// Runs backpropagation on current population of genomes
///
void GENNAlg::runBackPropagation(){
	
	unsigned int numInds = ga->population().size();
	
	for(unsigned int currInd = 0; currInd < numInds; currInd++){
		GEObjective::optimizeSolution(ga->population().individual(currInd));
	}
	ga->evaluatePop();
}


///
/// Establish restricted grammar based on current population
/// @param clearVariables when true only variables in current populations will
/// be included in the new restricted grammar.  When false it will keep existing
/// variables and add new ones from population.
///
void GENNAlg::setRestrictedGrammar(bool clearVariables){
	if(clearVariables)
		adjuster.clearVariables();
	fillPopulation();
	vector<string> variables;
	for(Solution* sol=pop.best(); sol != NULL; sol=pop.GetNext()){
		adjuster.addVariables(sol->getSymbols());
	}
	
	makeRestrictedGrammar(clearVariables, restrictMapper);
	
	convertNetworks(mapper, restrictMapper);	
	
// 	adjuster.editOnlyVarIncluded();
// 	adjuster.setMapper(restrictMapper, myRank);
// 	setMapperPrefs(restrictMapper);
// 	restrictMapper.setVariableCodonMap();
// 	GEObjective::setMapper(&restrictMapper);
// 	GE1DArrayGenome::setMapper(&restrictMapper);
// 	convertNetworks(mapper, restrictMapper);
}


///
/// Establish restricted grammar based on current population
/// @param clearVariables when true only variables in current populations will
/// be included in the new restricted grammar.  When false it will keep existing
/// variables and add new ones from population.
///
void GENNAlg::makeRestrictedGrammar(bool clearVariables,
  AthenaGrammarSI& useMapper){

// 	adjuster.addVariables(variables);
	adjuster.editOnlyVarIncluded();
	adjuster.setMapper(useMapper, myRank);
	setMapperPrefs(useMapper);
	useMapper.setVariableCodonMap();
	GEObjective::setMapper(&useMapper);
	GE1DArrayGenome::setMapper(&useMapper);
// 	convertNetworks(mapper, restrictMapper);	 
}

///
///
///
void GENNAlg::setRestrictedGrammar(bool clearVariables, vector<int>& genos,
  vector<int>& contins){
  if(clearVariables)
		adjuster.clearVariables();
  vector<string> variables;
  vector<int>::iterator iter;
  for(iter=genos.begin(); iter!=genos.end(); iter++){
    variables.push_back("G" + Stringmanip::numberToString(*iter));
  }
  for(iter=contins.begin(); iter!=contins.end(); iter++){
    variables.push_back("C" + Stringmanip::numberToString(*iter));
  }

  adjuster.addVariables(variables);
  makeRestrictedGrammar(clearVariables,mapper);
  InitGEgenome::setMapper(&mapper);
// 	fillPopulation();
// 	for(Solution* sol=pop.best(); sol != NULL; sol=pop.GetNext()){
// 		adjuster.addVariables(sol->getSymbols());
// 	}
// 	adjuster.editOnlyVarIncluded();
// 	adjuster.setMapper(restrictMapper, myRank);
// 	setMapperPrefs(restrictMapper);
// 	restrictMapper.setVariableCodonMap();
// 	GEObjective::setMapper(&restrictMapper);
// 	GE1DArrayGenome::setMapper(&restrictMapper);  
  
}


///
/// Converts networks from one mapping to other
/// @param currentMapper
/// @param newMapper
///
void GENNAlg::convertNetworks(AthenaGrammarSI& currentMapper, AthenaGrammarSI& newMapper){
		unsigned int numInds = ga->population().size();
		for(unsigned int currInd = 0; currInd < numInds; currInd++){
			GE1DArrayGenome& genome = (GE1DArrayGenome&)(ga->population().individual(currInd));            
			currentMapper.convertGenomeVariables(newMapper, genome);
			// reset the genome codons to the new values
			unsigned int genomeSize = currentMapper.getGenotype()->size();
			genome.resize(genomeSize);
			int i=0;

 	    // Now copy genotype onto genome
			Genotype::const_iterator genIt=(currentMapper.getGenotype())->begin();
			while(genIt!=(currentMapper.getGenotype())->end()){
	      genome.gene(i, *genIt);
			  genIt++;
	      i++;
			}
		}
		ga->statistics().setBestIndividual(ga->population());
}



void GENNAlg::outputGenome(GAGenome& g){
	GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);
	cout << "Score: " << genome.score() << " ";
	
	int length = genome.length();
	for(int i=0; i<length; i++){
		cout << (unsigned int)(genome.gene(i)) << " ";
	}
	cout << endl;
}


///
/// Finishes logs and compiles them into a single file
/// in cases where there are multiple processors
/// @param basename
/// @param cv
///
void GENNAlg::finishLog(string basename, int cv){
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
					parser.compileFiles(filenames, outName, GEObjective::getWorstScore());					
		}
}



///
/// Prepares for logging
/// @param outname
/// @param cv
///
void GENNAlg::prepareLog(string basename, int cv){
//     int totalPopSize = popSize;
		
		mainLogFilename = basename + ".cv." + Stringmanip::numberToString<int>(cv) + ".log";
		fitnessLogFilename =  basename + ".cv." + Stringmanip::numberToString<int>(cv) +
				".fitness.log";    
		snpnameLogFilename = basename + ".cv." + Stringmanip::numberToString<int>(cv) + 
				".snpsize.log";
				
		modelLog = new NNModelLog;
		
		string modelLogName = basename + ".rank." + Stringmanip::numberToString<int>(myRank) + ".cv." + 
				Stringmanip::numberToString<int>(cv) + ".models.log";
	
	if(logTypeSelected != LogNone && logTypeSelected != LogOverview){
		modelLog->openLog(modelLogName, GEObjective::calculatorName(), logTypeSelected==LogVariables);
		if(logTypeSelected == LogDetailed)
				modelLog->setDetailed(true);
	}
		
	#ifdef PARALLEL
		if(myRank==0){
	#endif
	
	ofstream outFile;
	switch(logTypeSelected){
		case LogDetailed:   
		case LogVariables:
		case LogOverview: 
		case LogSummary:
			outFile.open(mainLogFilename.c_str(), ios::out);
				if(!outFile.is_open()){
						throw AthenaExcept(mainLogFilename + " unable to open for writing log file");
				}  
			geLog->outputMainHeaders(outFile); 
			outFile.close();
		case LogNone:
		break;
	}
	
	
	#ifdef PARALLEL
		}
	#endif
}


///
/// Process and fill population with best models
///
void GENNAlg::selectBestModel(std::vector<Solution*>& solutions, data_manage::Dataholder * holder,
			Dataset* set, Config& config){

  setDataset(set);
// get list of variables
// edit grammar file
// run entire dataset with best variables
// return single best solution
  BestModelSelector bestSelector;
  bestSelector.setCorrThreshold(bestCorrThreshold);
  bestSelector.setCVThreshold(bestCVThreshold);
  bestSelector.selectBestVariables(solutions, holder);
  ofstream corrfile((config.getOutputName() + ".correlation.txt").c_str());
  corrfile << bestSelector;
  corrfile.close();
  
  vector<int> bestGenos, bestContin;
  if(myRank==0){
  	bestGenos = bestSelector.getIncludedGenos();
	  bestContin = bestSelector.getIncludedContins();
	}

#ifdef PARALLEL
//   exchangeBestVariables(totalNodes, myRank, bestGenos, bestContin);
  popMigrator.exchangeBestVariables(totalNodes, myRank, bestGenos, bestContin);
#endif
	
	if(bestGenos.empty() && bestContin.empty()){
	  throw  AthenaExcept("\nNo genotypes or continuous variables meet criteria for best model selection");
	}
	
	for_each (bestGenos.begin(), bestGenos.end(), data_manage::Utility::increment);
	for_each (bestContin.begin(), bestContin.end(), data_manage::Utility::increment);
	
	// edit grammar
  setRestrictedGrammar(true, bestGenos, bestContin);

  // set up run with this restricted grammar for same number of generations
  // as in the original run
  startLog(set->numGenos());
  prepareLog(config.getOutputName(), 0);
	initialize();
	// don't restrict variables as that has already been done
	ngensVarRestrict=0;
	for(int i=0; i < config.getNumExchanges(); i++){
	  if(step()){
		  break; // can complete early
			}
	}
  closeLog();

}


///
/// Writes current log to output
/// @param outname
///
void GENNAlg::writeLog(){

		ofstream outFile;
		outFile.open(mainLogFilename.c_str(), ios::app);
		if(!outFile.is_open()){
				throw AthenaExcept(mainLogFilename + " unable to open for writing log file");
		 }
		 geLog->outputLog(outFile);
		 outFile.close();
}



///
/// Clears logs from the algorithm
///
void GENNAlg::clearLogs(){
	for(unsigned int i=0; i<logs.size(); i++){
		delete logs[i];
	}
	logs.clear();
}


///
/// Writes graphical representation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GENNAlg::writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
	bool mapUsed, bool ottDummy, bool continMapUsed){
	GEObjective::outputModel(os, sol, holder, mapUsed, ottDummy, continMapUsed);
}


///
/// Produce graphical representation from writer specified
///
void GENNAlg::produceGraphic(std::string inputGraphic, std::string outputGraphic, 
			std::string imgWriter){
	string command = imgWriter + " -Tpng " + inputGraphic + " -o " + outputGraphic + ".png";
	system(command.c_str());
}


///
/// Writes model represented as an equation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GENNAlg::writeEquation(ostream& os, Solution* sol, data_manage::Dataholder* holder,
	bool mapUsed, bool ottDummy, bool continMapUsed){
	GEObjective::outputEquation(os, sol, holder, mapUsed, ottDummy, continMapUsed);
}

///
/// Returns graphical file extension appropriate to solution
/// @return extension
///
std::string GENNAlg::getGraphicalFileExt(){
	return GEObjective::getGraphicalExt();
}



///
/// Retrieves the models from BioFilter and stores the information in the algorithm
/// @param filename File with biofilter models
/// @param bioFileType Type of biological filter file (BINARY or TEXT)
/// @param holder Dataholder that contains map file info to convert model names
/// into indexes in the 
///
// void GENNAlg::getBioModels(std::string filename, std::string bioFileType, data_manage::Dataholder* holder){
// 	BioFilterModelCollection collection(filename, 100000, bioFileType);
// 	setBioModels(collection, holder);
// }



///
/// Fills biomodel collection from archive files 
/// @param genegeneFile genegene filename
/// @param archiveFile arcchive filename
/// @param holder Dataholder that contains map file info to convert model names
/// into indexes
///
// void GENNAlg::getBioModelsArchive(string genegeneFile, string archiveFile, data_manage::Dataholder* holder){
// 	BioFilterModelCollection collection(genegeneFile, archiveFile, 10000);
// 	setBioModels(collection, holder);
// }



///
/// Sets the bio models based on collection passed
///
// void GENNAlg::setBioModels(BioFilterModelCollection& collection, data_manage::Dataholder* holder){
// 
// 	vector<BioModel> bioModels;
// 	int numModelsNeeded = (unsigned int)(initBioFract * popSize);
// 
// 	if(biofilterSelectorType == rouletteSelect){  
// 		for(int i=0; i<numModelsNeeded; i++){
// 			bioModels.push_back(collection.getRandomModel());
// 		}
// 	}
// 	else{
// 		// set starting point for this algorithm -- in parallel mode this will give each 
// 		// node a different set
// 		collection.setStartModel(myRank * numModelsNeeded);
// 		for(int i=0; i<numModelsNeeded; i++){
// 			bioModels.push_back(collection.getNextModel());
// 		}
// 	}
// 
// 	mapper.clearModels();
// 
// 	// from models find the variables to use and store in mapper for use in initializing dataset
// 	for(vector<BioModel>::iterator iter=bioModels.begin(); iter != bioModels.end(); iter++){
// 		vector<int> indexes;
// 		int index;
// 		for(unsigned int i=0; i<iter->idString.size(); i++){
// 			index = holder->getGenoIndex(iter->idString[i]);
// 			indexes.push_back(index);
// 		}
// 		mapper.addModel(indexes);
// 	}
// 	  mapper.setModelCodons(dummyEncoded);
// }



#ifdef PARALLEL

void GENNAlg::setRank(int rank){
  popMigrator.setRank(rank);
 			Algorithm::setRank(rank); 
// 				InitGEgenome::setrank(rank);
// 				GEObjective::setrank(rank);
}


///
/// Returns sum of all the fitness goal checks on all nodes.  A 0 means no node has reached
/// the fitness.  A value > 0 means that many populations found a good enough fit.
/// @param rank
/// @param complete
///
// int GENNAlg::nodesCompleted(int complete){
// 	int sum_completed=0;
// 	MPI_Allreduce(&complete, &sum_completed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
// 	return sum_completed;
// }

///
/// Alternate MPI messaging taking advantage of MPI_Allgather functionality using struct
/// @param totalNodes
/// @param myRank
///
// void GENNAlg::sendAndReceiveStruct(int totalNodes, int myRank){
// 
// 	structMPI * send = new structMPI;
// 	
// 	GE1DArrayGenome& genome = (GE1DArrayGenome&)ga->statistics().bestIndividual();
// 	// package genome Info
// 	send->genomeParams[0] = genome.length();
// 	send->genomeParams[1] = genome.score();;
// 	send->genomeParams[2] = genome.getEffectiveSize();
// 	send->genomeParams[3] = genome.getNumGenes();
// 	send->genomeParams[4] = genome.getNumCovars();
// 	send->genomeParams[5] = genome.getNumIndsEvaluated();
// 	send->genomeParams[6] = genome.getSSTotal();
// 	send->genomeParams[7] = genome.getNumNodes();
// 	
// 	// package codons
// 	for(int i=0; i<genome.length(); i++){
// 		send->codons[i] = genome.gene(i);
// 	}
// 	
// 	// prepare receiving array
// 	structMPI * recv = new structMPI[totalNodes];
// 	
// 	MPI_Allgather (send, sizeof(*send), MPI_BYTE, recv, sizeof(*send), MPI_BYTE, MPI_COMM_WORLD);
// 	
// 	updateWithMigration(recv, totalNodes, myRank);
// 	
// 	delete send;
// 	delete [] recv;
// 	
// }

///
/// Broadcast best genotypes and continuous variables to all other nodes
///
// void GENNAlg::exchangeBestVariables(int totalNodes, int myRank, vector<int>& genotypes,
//   vector<int>& contins){
// 
//   int nVars = 1000;
//   int * variables = new int[nVars];
// 
//   vector<int>::iterator iter;
//   if(myRank==0){
//     int currVar=0;
//     variables[currVar++]=int(genotypes.size());
//     for(iter=genotypes.begin(); iter!=genotypes.end(); iter++){
//       variables[currVar++]=*iter;
//     }
//     variables[currVar++]=int(contins.size());
//     for(iter=contins.begin(); iter!=contins.end(); iter++){
//       variables[currVar++]=*iter;
//     }    
//   }
//   MPI_Bcast(variables, nVars, MPI_INT, 0, MPI_COMM_WORLD);
//   
//   genotypes.clear();
//   contins.clear();
//   
//   int currVar=0;
//   int n = variables[currVar++];
//   for(int i=0; i<n; i++){
//     genotypes.push_back(variables[currVar++]);
//   }
//   n=variables[currVar++];
//   for(int i=0; i<n; i++){
//     contins.push_back(variables[currVar++]);
//   }
// }



///
/// Incorporates migration into population
/// @param genomes 
/// @param totalNodes
/// @param myRank
///
// void GENNAlg::updateWithMigration(structMPI* mpiGenomes, int totalNodes, int myRank){
// 	GAPopulation pop(ga->population());
// 	
// 	for(int node=0; node < totalNodes; node++){
// 		if(myRank==node){
// 			continue;
// 		}
// 		
// 		GAGenome *tmpInd = ga->population().individual(0).clone();
// 		GE1DArrayGenome& genome = (GE1DArrayGenome&)*tmpInd;
// 		int len = mpiGenomes[node].genomeParams[0];
// 		genome.length(len);
// 		genome.setEffectiveSize(mpiGenomes[node].genomeParams[2]);
// 		genome.setNumGenes(mpiGenomes[node].genomeParams[3]);
// 		genome.setNumCovars(mpiGenomes[node].genomeParams[4]);
// 		genome.setNumIndsEvaluated(mpiGenomes[node].genomeParams[5]);
// 		genome.setSSTotal(mpiGenomes[node].genomeParams[6]);
// 		genome.setNumNodes(mpiGenomes[node].genomeParams[7]);
// 		for(int i=0; i<len; i++){
// 			genome.gene(i, mpiGenomes[node].codons[i]);
// 		}
// 		genome.score(mpiGenomes[node].genomeParams[1]);
// 
// 		pop.add(genome);
// 		
// 		delete tmpInd;
// 	}
// 
// 	// remove worst individuals from population
// 	for(int i=0; i < totalNodes-1; i++)
// 		pop.destroy();
// 		
// 	ga->population(pop);
// 
// }


#endif
