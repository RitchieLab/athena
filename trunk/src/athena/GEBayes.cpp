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
#include "GEBayes.h"
#include "Terminals.h"
// #include "GENNGrammarAdjuster.h"
#include "ModelLogParser.h"
// #include "BestModelSelector.h"
#include "GAParetoSelector.h"
#include "GAParetoRankSelector.h"
#include "SumFileReader.h"
#include <ga/ga.h>
#include <ctime>
#include <set>
#include <algorithm>
#include <ScaleCategorical.h>

///
/// Constructor
///
GEBayes::GEBayes(){
	initializeParams();
}


///
/// Destructor
///
GEBayes::~GEBayes(){
		if(geLog != NULL){
				delete geLog;
				geLog=NULL;
		}
		freeMemory();
}

///
/// Initializes parameters for basic run of algorithm.  These
/// can be modified using the set_params function
///
void GEBayes::initializeParams(){
 
	calculatorName = "K2";
	logTypeSelected = LogNone;
	geLog = NULL;
	  
	minNumParents = 1;
	minNumChildren = 1;
	maxNumParents = 3;
	maxNumChildren = 2;
	balAccStart=-1;
	balAccFreq=0;

	paramMap["CHILDRANGE"] = childRange;
	paramMap["PARENTRANGE"] = parentRange;
	paramMap["BAFREQ"] = bafreq;
	paramMap["BABEGIN"] = bastart;
	
}


///
/// Sets all genomes in a population to either effective crossover or one point
///
void GEBayes::resetCrossover(){
		if(effectiveXO){
			ga->crossover(GE1DArrayGenome::effCrossover);
		}
		else{
			ga->crossover(GE1DArrayGenome::OnePointCrossover);
		}
}


///
/// Sets values in main configuration to defaults needed by 
/// GEBAyes Algorithm
/// @param configuration Config
///
void GEBayes::setConfigDefaults(Config& configuration, AlgorithmParams& algParam){
//   if(algParam.params["CALCTYPE"].compare("RSQUARED")==0){
//     configuration.setStatusAdjust("MINMAX");
//   }
	configuration.setContinAdjust("MAKECATEGORIAL");
	configuration.setStatusAdjust("MAKECATEGORIAL");
}


///
/// Sets random seed 
/// @param seed 
///
void GEBayes::setRand(unsigned int seed){
	GARandomSeed(seed);
	srand(seed);
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
void GEBayes::setParams(AlgorithmParams& algParam, int numExchanges, int numGenos, 
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
												" in Algorithm GEBayes");
								break;
						case parentRange:
							// split on ':' for min:max
							tokens = Stringmanip::split(mapIter->second, ':');
							minNumParents = Stringmanip::stringToNumber<unsigned int>(tokens[0]);
							maxNumParents = Stringmanip::stringToNumber<unsigned int>(tokens[1]);
							if(minNumParents < 1)
								throw AthenaExcept("PARENTRANGE minimum must be 1 or greater");
							break;
						case childRange:
							tokens = Stringmanip::split(mapIter->second, ':');
							minNumChildren = Stringmanip::stringToNumber<unsigned int>(tokens[0]);
							maxNumChildren = Stringmanip::stringToNumber<unsigned int>(tokens[1]);
							if(minNumChildren < 1)
								throw AthenaExcept("CHILDRANGE minimum must be 1 or greater");
							break;
						case bafreq:
							balAccFreq = Stringmanip::stringToNumber<int>(mapIter->second);
							break;
						case bastart:
							balAccStart = Stringmanip::stringToNumber<int>(mapIter->second);
							break;
#ifdef ATHENA_BLOAT_CONTROL
						case prunePlantFract:
#endif
						default:
							if(paramMap.find(mapIter->first) == paramMap.end())
								throw AthenaExcept("No match for parameter " + mapIter->first +
												" in Algorithm GEBayes");               
				}
		}
		
		// first optimization of backpropagation 
// 		bpNextOpt = bpFirstGen;
		baNextOpt = balAccStart;
			 
		setGAParams(excludedGenos, excludedContins);
		
}


///
/// Set current Dataset for running algorithm
/// @param new_set Dataset
/// 
void GEBayes::setDataset(Dataset* newSet){
	 set = newSet;
// 	 ScaleCategorical scaler;
// 	 scaler.adjustContin(set);
// 	 scaler.adjustStatus(set);
	 GEObjective::setDataset(set);
}


///
/// Initializes the algorithm
///
void GEBayes::initialize(){
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
// 			geLog->setMaxBest(GEObjective::logMaxBest());
		}
		mapper.resetGrammarModels();
		
		mapper.setVariableCodonMap();
		GE1DArrayGenome::setMapper(&mapper);

		ga->initialize();
		baNextOpt = balAccStart;
		// run optimization after initialization when indicated
		if(baNextOpt == 0){
			runBalancedAccuracyOptimization();
			baNextOpt += balAccFreq;
		}
			
		fillPopulation();
		fillLog();
}


///
/// Sets parameters for use with GAlib
/// @param alg_params AlgorithmParams
/// @throws AthenaExcept on error
///
void GEBayes::setGAParams(vector<unsigned int>& excludedGenos, 
			vector<unsigned int>& excludedContins){   
    
		GARandomSeed(randSeed);
		srand(randSeed);   
		// first free ga memory if already run
		freeMemory();
		
		//Initialize the GE mapper
	 //Set maximum number of wrapping events per mapping
	 mapper.setMaxWraps(wrapEvents);
	 expandVariables();
	 // ADJUSTING GRAMMAR -- can be added back for the new discriminant network grammar
// 	 adjuster.setBayesianSize(minNumParents, maxNumParents, minNumChildren, 
// 			maxNumChildren);
		// remove any excluded SNPs and/or continuous variables
		if(!excludedGenos.empty() || !excludedContins.empty())
			excludeVariables(excludedGenos, excludedContins);
			
		adjuster.setMapper(mapper);
		setMapperPrefs(mapper);
			
		GEObjective::setSolutionType("BAYES", calculatorName);

	 if(calculatorName.find("K2") != string::npos){
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
/// Fills the population after a step to allow for exchange of models with 
/// other algorithms.
/// 
void GEBayes::fillPopulation(){
		unsigned int numInds = ga->population().size();
		// clear solution population
		pop.clear();		
		for(unsigned int currInd = 0; currInd < numInds; currInd++){
				pop.insert(convertGenome(ga->population().individual(currInd)));
		}
}


///
/// Converts an individual from population into a Solution and returns it
/// @param GAGenome ind
/// @return BayesSolution*
///
BayesSolution* GEBayes::convertGenome(GAGenome& ind){
	GE1DArrayGenome& genome = (GE1DArrayGenome&) ind;
	BayesSolution* sol = (BayesSolution*)GEObjective::getBlankSolution();
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
// 	sol->setComplexity(genome.getNumNodes());
	return sol;
}


///
/// Return additional output names (if any such as AUC)
/// 
vector<string> GEBayes::getAdditionalOutputNames(){
  return GEObjective::getAdditionalOutputNames();
}

///
/// Calculate and return additional output (whether model is better than unconnected 
///  bayesian network) for best model
///
void GEBayes::getAdditionalFinalOutput(Dataset* set){
	GEObjective::setDataset(set);
	GEObjective::setRefDataset(set);
// cout << "getAdditionalFinalOutput" << endl;

	unsigned int numInds = ga->population().size();
  for(unsigned int currInd = 0; currInd < numInds; currInd++){
		pop[currInd]->setAdditionalOutput(GEObjective::getAdditionalFinalOutput(ga->population().individual(currInd)));
	}

	// if used the balanced accuracy optimization use that as the sorting for final population 
// 	if(int(restrictStepsDone) == (baNextOpt - balAccFreq)){
// 		setBABest();
// 	}
	
}


///
/// Calculate and return additional output (whether model is better than unconnected 
///  bayesian network) for best model
///
void GEBayes::getAdditionalFinalOutput(Dataset* testing, Dataset* training,
		data_manage::Dataholder* holder, bool mapUsed, bool ottDummy, bool continMapUsed){

	GEObjective::setDataset(testing);
	GEObjective::setRefDataset(training);
	unsigned int numInds = ga->population().size();
;
  for(unsigned int currInd = 0; currInd < numInds; currInd++){
		pop[currInd]->setAdditionalOutput(GEObjective::getAdditionalFinalOutput(ga->population().individual(currInd)));
// cout << "addoutput: ";
// for(unsigned int i=0; i<pop[currInd]->getAdditionalOutput().size(); i++){
// cout << pop[currInd]->getAdditionalOutput()[i] << " ";
// }
// cout << endl;
	}

	// if used the balanced accuracy optimization use that as the sorting for final population 
	if(int(restrictStepsDone) == (baNextOpt - balAccFreq)){
		setBABest();
	}
}


///
/// Set final best model based on the balanced accuracy if that feature was used
///
void GEBayes::setBABest(){
	modelSortCol = 5;
	float value;
	vector<Solution*> sols;
	vector<string> values;
	vector<float> scores;
	
	Solution* sol = pop.GetFirst();
	while(sol != NULL){
		sols.push_back(sol);
		values = sol->getAdditionalOutput();
		if(!values.empty()){
			value = Stringmanip::stringToNumber<float>(values[1]);
// if(values.size() > 2)
// cout << "training=" << values[1] << " testing=" << values[3] << endl;
// else
// cout << "training=" << values[1] << endl;
			scores.push_back(value);
		}
		else{
			scores.push_back(0.0);
		}
		sol = pop.GetNext();
	}
	unsigned int numInds = ga->population().size();
	pop.clearTree();
	for(unsigned int currInd=0; currInd < numInds; currInd++){
		pop.insert(sols[currInd], scores[currInd]);
	}
// cout << "end setBABest"<<endl;
}


///
/// Run dataset against models contained in previously generated summary file
/// @param sumFile ATHENA summary file
///
vector<Solution*> GEBayes::runValidation(std::string sumFile){
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
/// Fills the algorithm log with data from the population
///
void GEBayes::fillLog(){
		
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
				// no need to convert scores for Bayesian network
					geLog->addFitness(genome.score(), genome.getGenos(), genome.getCovars());
				}
				else{  // add converted score to log file
					geLog->addFitness(pop[0]->adjustScoreOut(genome.score(), genome.getNumIndsEvaluated(),
						set->getConstant(), getFitnessName()), genome.getGenos(), genome.getCovars());
				}
				geLog->addNumGenos(genome.getNumGenes());
				geLog->addNumCovars(genome.getNumCovars());
				geLog->addEpochs(genome.getNumEpochsTrained());
				// zero out epochs`
				genome.setNumEpochsTrained(0);
			}
		}
		geLog->completeGen();

		#ifdef HAVE_CXX_MPI
			geLog->sendReceiveLogs(totalNodes, myRank);
		#endif
		
		// output log information -- appended to existing log files
		if(myRank==0){
				writeLog();
		}

		if(logTypeSelected==LogVariables){
			BayesSolution * solution;
			int nSolutions = pop.numSolutions();
			modelLog->addGeneration(geLog->getCurrentGen());
			for(int i=0; i<nSolutions; i++){
				solution = (BayesSolution*)pop[i];
				modelLog->writeVariables(*solution, dummyEncoded);
			}
			return;
		}

		if(logTypeSelected != LogOverview){
			BayesSolution * solution;
			int nSolutions = pop.numSolutions();
			for(int i=0; i<nSolutions; i++){
				solution = (BayesSolution*)pop[i];
				modelLog->writeSolution(*solution, geLog->getCurrentGen(), i+1);
			}
		}
	
}


///
/// Close the log files
///
void GEBayes::closeLog(){
	if(logTypeSelected != LogNone){
		modelLog->closeLog();
	}
}


///
/// Starts fresh log 
///
void GEBayes::startLog(int numSnps){
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
/// Writes current log to output
/// @param outname
///
void GEBayes::writeLog(){

		ofstream outFile;
		outFile.open(mainLogFilename.c_str(), ios::app);
		if(!outFile.is_open()){
				throw AthenaExcept(mainLogFilename + " unable to open for writing log file");
		 }
		 geLog->outputLog(outFile);
		 outFile.close();
}


///
/// Prepares for logging
/// @param outname
/// @param cv
///
void GEBayes::prepareLog(string basename, int cv){
//     int totalPopSize = popSize;
		
		mainLogFilename = basename + ".cv." + Stringmanip::numberToString<int>(cv) + ".log";
		fitnessLogFilename =  basename + ".cv." + Stringmanip::numberToString<int>(cv) +
				".fitness.log";    
		snpnameLogFilename = basename + ".cv." + Stringmanip::numberToString<int>(cv) + 
				".snpsize.log";
				
		modelLog = new BayesModelLog;
		
		string modelLogName = basename + ".rank." + Stringmanip::numberToString<int>(myRank) + ".cv." + 
				Stringmanip::numberToString<int>(cv) + ".models.log";
	
	if(logTypeSelected != LogNone && logTypeSelected != LogOverview){
		modelLog->openLog(modelLogName, GEObjective::calculatorName(), logTypeSelected==LogVariables);
		if(logTypeSelected == LogDetailed)
				modelLog->setDetailed(true);
	}
		
	#ifdef HAVE_CXX_MPI
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
	
	
	#ifdef HAVE_CXX_MPI
		}
	#endif
}


///
/// Finishes logs and compiles them into a single file
/// in cases where there are multiple processors
/// @param basename
/// @param cv
///
void GEBayes::finishLog(string basename, int cv){
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


/// Writes graphical representation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GEBayes::writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
	bool mapUsed, bool ottDummy, bool continMapUsed){
	GEObjective::outputModel(os, sol, holder, mapUsed, ottDummy, continMapUsed);
}

///
/// Writes model represented as an equation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GEBayes::writeEquation(ostream& os, Solution* sol, data_manage::Dataholder* holder,
	bool mapUsed, bool ottDummy, bool continMapUsed){
	GEObjective::outputEquation(os, sol, holder, mapUsed, ottDummy, continMapUsed);
}


///
/// Establishes the algorithm for the run based on the parameters
/// set in this algorithm.
///
void GEBayes::run(){
		while(!ga->done()){
				ga->step();
		}
}

///
/// Runs an indicated interval for the algorithm. This interval is set
/// in the stepSize variable.  
///
int GEBayes::step(){
		int completed =0;
		GE1DArrayGenome::myRank = myRank;
// cout << "myRank=" << myRank << " in GEBayes step()" << endl;
		for(unsigned int i=0; i < stepSize; i++){
// cout << "myRank=" << myRank << " step i=" << i;
				if(!ga->done()){
					if(logTypeSelected!=LogNone){
						geLog->addGeneration();
					}
					ga->step();
// cout << " completed" << endl;
				 restrictStepsDone++;
				 // check for need to change crossovers
				 if(ngensBlockCross && restrictStepsDone == ngensBlockCross){
						resetCrossover();
				 }
// cout << "restrictStepsDone=" << restrictStepsDone << endl;
				 // check to see if the back propagation should be done at this generation
				 if(int(restrictStepsDone) == baNextOpt && balAccStart >= 0){
// cout << "restrictStepsDone=" << restrictStepsDone << " running optimization\n";
						runBalancedAccuracyOptimization();
						baNextOpt += balAccFreq;
				 }				 
				 
				 fillLog();
				}	
		}
// cout << "myRank=" << myRank << " after steps done in GEBayes step()" << endl;
// if(int(restrictStepsDone) == bpNextOpt && bpFirstGen >= 0){
// runBalancedAccuracyOptimization();
// bpNextOpt += bpFreqGen;
// }

		#ifdef HAVE_CXX_MPI
// cout << "myRank=" <<  myRank << " call sendAndReceive" << endl;
			popMigrator.sendAndReceiveStruct(totalNodes, myRank,ga);

			// when running in parallel transfer around populations
// 			if(ngensVarRestrict && restrictStepsDone < ngensVarRestrict){
// 				// after transfer construct new grammar and set for use
// 				setRestrictedGrammar(resetRestrictedAtMigration);
// 			}
		#endif
		
// cout << "myRank=" <<  myRank << " after sendAndReceive" << endl;
		// only need to fill population at this point not at end of each generation
		fillPopulation();
// cout << "myRank=" <<  myRank << " after filling population" << endl;
		return completed;
}


///
/// Runs backpropagation on current population of genomes
///
void GEBayes::runBalancedAccuracyOptimization(){
	
	unsigned int numInds = ga->population().size();

	vector<string> values;
	float accuracy, maxAcc=0.0;
	int bestAccInd = 0;
	
// 	GEObjective::setDataset(set);
	GEObjective::setRefDataset(set);
	for(unsigned int currInd = 0; currInd < numInds; currInd++){
// cout << "currInd=" << currInd << " total=" << numInds << endl;
		GEObjective::optimizeSolution(ga->population().individual(currInd));
// cout << "finished optimize" << endl;
		// get additional results -- second will be the training balanced accuracy
		values = GEObjective::getAdditionalFinalOutput(ga->population().individual(currInd));
		if(!values.empty()){
			accuracy = Stringmanip::stringToNumber<float>(values[1]);
			if(accuracy > maxAcc){
				bestAccInd = currInd;
				maxAcc = accuracy;
			}
		}
	}
	ga->evaluatePop();
// cout << "bestAcc=" << maxAcc << " best genome index=" << bestAccInd << endl;	
	GE1DArrayGenome bestOptGenome = (GE1DArrayGenome &)ga->population().individual(bestAccInd);
	
	/// get best optimized individual and create a population to hold it
	GAPopulation tmpPop(bestOptGenome);
	
	/// set bestIndividual in stats to match it so it will be passed to nodes
	ga->statistics().setBestIndividual(tmpPop, gaFalse);
}


///
/// Establishes test data set values for the population
/// @param test_set Dataset containing individuals for the test set
///
void GEBayes::testSolution(Dataset* testSet, int nproc){
	 
	 // make sure test set is categorical
// 	 ScaleCategorical scaler;
// 	 scaler.adjustContin(testSet);
// 	 scaler.adjustStatus(testSet);
	 
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
/// Returns graphical file extension appropriate to solution
/// @return extension
///
std::string GEBayes::getGraphicalFileExt(){
	return GEObjective::getGraphicalExt();
}


#ifdef HAVE_CXX_MPI
void GEBayes::setRank(int rank){
  popMigrator.setRank(rank);
                        Algorithm::setRank(rank);
}
#endif
