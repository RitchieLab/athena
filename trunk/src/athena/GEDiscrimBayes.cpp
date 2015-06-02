#include "GEDiscrimBayes.h"
#include "GEObjective.h"
#include "SumFileReader.h"
#include "ModelLogParser.h"
#include <algorithm>

///
/// Constructor
///
GEDiscrimBayes::GEDiscrimBayes(){
	initializeParams();
}

///
/// Destructor
///
GEDiscrimBayes::~GEDiscrimBayes(){
// 		if(geLog != NULL){
// 				delete geLog;
// 				geLog=NULL;
// 		}
// 		freeMemory();
		delete caseAlg;
		delete controlAlg;
	if(caseDataset != NULL)
		delete caseDataset;
	if(controlDataset != NULL)
		delete controlDataset;		
}


///
/// Initializes parameters for basic run of algorithm.  These
/// can be modified using the set_params function
///
void GEDiscrimBayes::initializeParams(){
 
 	caseAlg = new GEBayes;
 	controlAlg = new GEBayes;
 	topModelsUsed = 5;
  
// 	calculatorName = "K2";
// 	logTypeSelected = LogNone;
// 	geLog = NULL;
	ga = NULL;
	caseDataset = controlDataset = NULL;
	currCV=1;
	
	paramMap["NUMTOPMODELS"] = modelsToUse;
// 	  
// 	minNumParents = 1;
// 	minNumChildren = 1;
// 	maxNumParents = 3;
// 	maxNumChildren = 2;
// 	balAccStart=-1;
// 	balAccFreq=0;
// 
// 	paramMap["CHILDRANGE"] = childRange;
// 	paramMap["PARENTRANGE"] = parentRange;
// 	paramMap["BAFREQ"] = bafreq;
// 	paramMap["BABEGIN"] = bastart;
	
}


///
/// Sets random seed 
/// @param seed 
///
void GEDiscrimBayes::setRand(unsigned int seed){
	GARandomSeed(seed);
	srand(seed);
}

///
/// Sets values in main configuration to defaults needed by 
/// GEDiscrimBayes Algorithm
/// @param configuration Config
///
void GEDiscrimBayes::setConfigDefaults(Config& configuration, AlgorithmParams& algParam){
	// algorithm will write its own output
	configuration.setLogType("NONE");
	configuration.setSummaryOnly("SUPPRESS");
	outputName = configuration.getOutputName();
}

///
/// Sets encoding value for genotypes for all algorithms
/// @param ottDummyEncoded
///
void GEDiscrimBayes::setDummyEncoding(bool ottDummyEncoded){
	dummyEncoded = ottDummyEncoded;
	caseAlg->setDummyEncoding(ottDummyEncoded);
	controlAlg->setDummyEncoding(ottDummyEncoded);
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
void GEDiscrimBayes::setParams(AlgorithmParams& algParam, int numExchanges, int numGenos, 
	int numContin, vector<unsigned int>& excludedGenos, vector<unsigned int>& excludedContins){
// cout << "setting params" << endl;	
		caseAlg->setParams(algParam, numExchanges, numGenos, numContin, excludedGenos, excludedContins);
		controlAlg->setParams(algParam, numExchanges, numGenos, numContin, excludedGenos, excludedContins);
		
		Algorithm::setParams(algParam, numExchanges, numGenos, numContin, excludedGenos,
			excludedContins);
		
// cout << "call params here" << endl;
		map<string, string>::iterator mapIter;
		vector<string> tokens;
		
		for(mapIter = algParam.params.begin(); mapIter != algParam.params.end(); 
			mapIter++){     
				switch(paramMap[mapIter->first]){
					case noMatchParam:
						break;
					case modelsToUse:
						topModelsUsed = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
// cout << "topModelsUsed=" << topModelsUsed << endl;
					break;
				}
			}
		/*
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
		*/
}


///
/// Sets parameters for use with GAlib
/// @param alg_params AlgorithmParams
/// @throws AthenaExcept on error
///
void GEDiscrimBayes::setGAParams(vector<unsigned int>& excludedGenos, 
			vector<unsigned int>& excludedContins){
    /*
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
		 
		// can probably drop this
    setInitParams();
    */
}


///
/// Split Dataset into case/control sets and pass to appropriate
/// GEBAyes algorithm
/// @param newSet Dataset to split
///
void GEDiscrimBayes::setDataset(Dataset* newSet){
// cout << "start setDataset" << endl;
	set = newSet;
	if(caseDataset != NULL)
		delete caseDataset;
	if(controlDataset != NULL)
		delete controlDataset;
		
// cout << z[0] << " " << z[1] << endl;		
// cout << "first => " << (*newSet)[0]->getStatus() << "\n";		
	vector<Dataset*> splitSets = newSet->splitCaseControl();
// cout << "first immediate=> " << (*newSet)[0]->getStatus() << "\n";	
	controlAlg->setDataset(splitSets[0]);
// cout << "first after set=> " << (*newSet)[0]->getStatus() << "\n";
	caseAlg->setDataset(splitSets[1]);
// cout << "first after second set=> " << (*newSet)[0]->getStatus() << "\n";
	caseDataset = splitSets[1];
	controlDataset = splitSets[0];
// cout << "caseDataset has " << caseDataset->numInds() << endl;
// cout << "controlDataset has " << controlDataset->numInds() << endl;
	
// cout << "first after => " << (*newSet)[0]->getStatus() << "\n";
	
}

///
/// Starts fresh log 
///
void GEDiscrimBayes::startLog(int numSnps){
	caseAlg->startLog(numSnps);
	controlAlg->startLog(numSnps);
}

///
/// Retrieves the models from BioFilter and stores the information in the algorithm
/// @param filename File with biofilter models
/// @param bioFileType Type of biological filter file (BINARY or TEXT)
/// @param holder Dataholder that contains map file info to convert model names
/// into indexes in the 
///
void GEDiscrimBayes::getBioModels(std::string filename, std::string bioFileType, data_manage::Dataholder* holder){
	caseAlg->getBioModels(filename, bioFileType, holder);
	controlAlg->getBioModels(filename, bioFileType, holder);
}



///
/// Fills biomodel collection from archive files 
/// @param genegeneFile genegene filename
/// @param archiveFile arcchive filename
/// @param holder Dataholder that contains map file info to convert model names
/// into indexes
///
void GEDiscrimBayes::getBioModelsArchive(string genegeneFile, string archiveFile, data_manage::Dataholder* holder){
	caseAlg->getBioModelsArchive(genegeneFile, archiveFile, holder);
	controlAlg->getBioModelsArchive(genegeneFile, archiveFile, holder);
}


///
/// Prepares for logging
/// @param outname
/// @param cv
///
void GEDiscrimBayes::prepareLog(string basename, int cv){
	caseAlg->prepareLog(basename, cv);
	controlAlg->prepareLog(basename, cv);
}

/// Initializes the algorithm for each new dataset to test
void GEDiscrimBayes::initialize(){
	caseAlg->initialize();
	controlAlg->initialize();
	
	// initialize structures here for next CV
	// to do after
	
}


///
/// Runs an indicated interval for the algorithm. This interval is set
/// in the stepSize variable.  
///
int GEDiscrimBayes::step(){
		int completed =0;
		// run each separately
// cout << "myRank=" << myRank << "run case step" << endl;
		GEObjective::setDataset(caseDataset);
		caseAlg->step();
// cout << "myRank=" << myRank << "run control step" << endl;		
		GEObjective::setDataset(controlDataset); 
		controlAlg->step();
// cout << "myRank=" << myRank << "finished case and control steps" << endl;
		return completed;
}


///
/// Starts fresh log 
///
void GEDiscrimBayes::closeLog(){
	caseAlg->closeLog();
	controlAlg->closeLog();
}

///
/// Performs discriminant analysis
///
void GEDiscrimBayes::getAdditionalFinalOutput(Dataset* testing, Dataset* training,
	data_manage::Dataholder* holder, bool mapUsed, bool ottDummy, bool continMapUsed){

// 	cout << "Fill getAdditionalFinalOutput here" << endl;

	// store in set -- use model equation as the key (sorted?)
	// get each model convert to model equation 
	// store in set and increment counter
	map<string, modScores> caseModels;
	map<string, modScores> controlModels;
	totalModels(caseAlg, caseModels, holder, mapUsed, ottDummy, continMapUsed, true);
// 	cout << "finished compiling case model totals" << endl;
	totalModels(controlAlg, controlModels, holder, mapUsed, ottDummy, continMapUsed, false);
// 	cout << "finished compiling control model totals" << endl;
	
	map<unsigned int, vector<double> > caseOrphanProbs, controlOrphanProbs;
	// calculate probability tables for each variable in case training set
	// for situation where there are no parents for that variable in the model
	calcProbTables(caseDataset, caseOrphanProbs, holder);
	// calculate probability tables for each variable in control training set
	// for situation where there are no parents in model
	calcProbTables(controlDataset, controlOrphanProbs, holder);
	
	double case2conRatio = caseDataset->numInds() / double(controlDataset->numInds()+caseDataset->numInds());
	int numCases = caseDataset->numInds() + int(testing->numInds()*case2conRatio+0.5);
	int numControls = controlDataset->numInds() + int(testing->numInds()*case2conRatio+0.5);
	double caseMissCell = double(1) / numCases;
	double conMissCell = double(1) / numControls;
		
	// loop through each model (calculate conditional probability for cases/controls)
	// calculate conditional probabilities for case models
	setConditionalTables(caseDataset, caseModels, holder, caseMissCell);
	// calculate conditional probabilities for control models
	setConditionalTables(controlDataset, controlModels, holder, conMissCell);
	
	// for each individual calculate total probability score by multiplying all the probabilities
	// for all variables in set (using conditional probs where appropriate)
	// when model is a case use the case probs and when model is a control use
	// the control probabilities no matter the status of the individual
	IndResults emptyResult;
	vector<IndResults> trainingScores(training->numInds(),emptyResult), 
		testingScores(testing->numInds(),emptyResult);
	// each individual will end up with a score for each top model
	// store results for testing and training sets
// 	cout << "training set used" << endl;
	setIndModScores(training, caseModels, trainingScores, caseOrphanProbs);
	setIndModScores(training, controlModels, trainingScores, controlOrphanProbs);
// 	cout << "testing set used" << endl;
	setIndModScores(testing, caseModels, testingScores, caseOrphanProbs);
	setIndModScores(testing, controlModels, testingScores, controlOrphanProbs);
	
	double trainingAUC = setPredictedScores(trainingScores, caseModels, controlModels, case2conRatio);
	double testingAUC = setPredictedScores(testingScores, caseModels, controlModels, case2conRatio);
	 
#ifdef HAVE_CXX_MPI
if(myRank == 0){
#endif	 
	string trainingFile = outputName + ".cv." + Stringmanip::numberToString(currCV) + ".train.inds.txt";
	writeIndScores(trainingFile, trainingScores);
	string testingFile = outputName + ".cv." + Stringmanip::numberToString(currCV) + ".test.inds.txt";
	writeIndScores(testingFile, testingScores);
	
	// add this CV result to .sum file
	string sumFile = outputName + ".GEBN.sum";
	ofstream sumstream;
	if(currCV > 1){
		sumstream.open(sumFile.c_str(), ios::app);
	}
	else{
		sumstream.open(sumFile.c_str(), ios::out);
		sumstream << "CV\tCase Models\tControl Models\tTrain AUC\tTest AUC\n";
	}
	
	map<string, modScores>::iterator caseIter=caseModels.begin();
	map<string, modScores>::iterator conIter=controlModels.begin();
	sumstream << currCV << "\t" << caseIter->first << "\t" << conIter->first << "\t" <<
		trainingAUC << "\t" << testingAUC << "\n";
	for(int i=1; i<topModelsUsed; i++){
		++caseIter;++conIter;
		sumstream << "\t" << caseIter->first << "\t" << conIter->first << "\t\t\n";
	}
	sumstream.close();
#ifdef HAVE_CXX_MPI
}
#endif
	currCV++;
	
	// below is just a placeholder
	pop=caseAlg->getPopulation();
}

///
/// Set predicted phenotype for each individual in set
/// @return AUC
///
double GEDiscrimBayes::setPredictedScores(vector<IndResults>& indScores, map<string,modScores>& caseModels,
	map<string,modScores>& controlModels, double caseRatio){
	
	double caseScore, conScore;
	int sIndex, totalCount;
	stat::TestResult tempResult;
	std::vector<stat::TestResult> results;
	for(vector<IndResults>::iterator indIter=indScores.begin(); indIter != indScores.end();
		++indIter){
		caseScore = conScore = 0.0;
		sIndex=totalCount=0;
		for(map<string,modScores>::iterator caseModIter=caseModels.begin(); caseModIter != caseModels.end();
			++caseModIter){
			caseScore += indIter->scores[sIndex++] * caseModIter->second.count;
			totalCount += caseModIter->second.count;
		}
		caseScore /= double(totalCount);
		
		totalCount=0;
		for(map<string,modScores>::iterator conModIter=controlModels.begin(); conModIter != controlModels.end();
			++conModIter){
			conScore += indIter->scores[sIndex++] * conModIter->second.count;
			totalCount += conModIter->second.count;
		}
		conScore /= double(totalCount);
		indIter->predicted = caseRatio * caseScore / (caseRatio * caseScore + 
			(1-caseRatio) * conScore);
		tempResult.score = indIter->predicted;
		tempResult.status = indIter->phenotype;
		results.push_back(tempResult);
	} 
	
	//calculate and return AUC
	return stat::AUCCalc::calculateAUC(results);
}


///
/// Sets the individual scores for every individual in set for the models passed
///
void GEDiscrimBayes::setIndModScores(Dataset* dset, map<string,modScores>& models,
	vector<IndResults>& indScores, map<unsigned int, vector<double> >& orphanProbs){
	unsigned int totalGenos = dset->numGenos();
// vector<int> t(2,0);
// cout << "start setIndModScores set size=" << dset->numInds() << endl;
	for(size_t i=0; i<dset->numInds(); i++){
		indScores[i].indID = (*dset)[i]->getID();
		indScores[i].phenotype = (*dset)[i]->getStatus();
// t[indScores[i].phenotype]++;
//cout << "i=" << i << " ID=" << indScores[i].indID  << " " << indScores[i].phenotype << endl;
		// loop through each model and assign a value to each individual for each model
		for(map<string, modScores>::iterator modIter=models.begin(); modIter != models.end();
			++modIter){
// cout << "model is " << modIter->first << endl;
			indScores[i].scores.push_back(1.0);
			for(vector<ConditTable>::iterator tableIter=modIter->second.tables.begin();
				tableIter != modIter->second.tables.end(); ++tableIter){
					// have to get the probability and include it
					int value = 0;
					for(size_t j=0; j<tableIter->parentIndexes.size(); j++){
// cout << "j=" << j << " genotype=" << (*dset)[i]->getGenotype(tableIter->parentIndexes[j]) 
// 	<< " SNP=" << tableIter->parentIndexes[j] << endl;
						value += (*dset)[i]->getGenotype(tableIter->parentIndexes[j]) *
							tableIter->cumulativeLevels[j]; 
					}
// cout << "value=" << value << endl;
// cout << "nodeIndex=" << tableIter->nodeIndex << endl;
// 					indScores[i].scores.back() *= tableIter->probs[indScores[i].phenotype][index];
// 					indScores[i].scores.back() *= tableIter->probs[tableIter->nodeIndex][value];
					indScores[i].scores.back() *= tableIter->probs[(*dset)[i]->getGenotype(tableIter->nodeIndex)][value];
			}
			// add in all the scores for the other variables in the set
			for(unsigned int geno=0; geno < totalGenos; geno++){
				if(modIter->second.varsWithParents.find(geno) == modIter->second.varsWithParents.end()){
					indScores[i].scores.back() *= orphanProbs[geno][(*dset)[i]->getGenotype(geno)];
				}
			}
		}
	}	
// cout << "cases=" << t[1] << " controls=" << t[0] << "\n";
}

///
/// Output individual results with predicted scores
///
void GEDiscrimBayes::writeIndScores(string filename, vector<IndResults>& scores){
	ofstream out;
	out.open(filename.c_str(), ios::out);
	out << "ID\tPHENO\tPREDICTED\n";
	for(vector<IndResults>::iterator iter=scores.begin(); iter != scores.end();
		++iter){
		out << iter->indID << "\t" << iter->phenotype << "\t" << iter->predicted << "\n";
	}	
	out.close();
}

///
/// Set conditional probability tables for models
///
void GEDiscrimBayes::setConditionalTables(Dataset* dset, map<string, modScores>& topModels,
	Dataholder* holder, double missValue){
	for(map<string, modScores>::iterator modIter=topModels.begin(); modIter != topModels.end();
		++modIter){
		int tableIndex=0;
		modIter->second.tables.clear();
		ConditTable emptyTable;

// for(size_t i=0; i<dset->numGenos(); i++){
// 	cout << "i=" << i << " genoname=" << holder->getGenoName(i) << endl;
// }
// exit(1);
		
		// need to split model string and then use only those that have parents
		vector<string> nodes = tokenizeEquation(modIter->first);
		for(vector<string>::iterator nodeIter=nodes.begin(); nodeIter != nodes.end();
			++nodeIter){
			size_t cPos = (*nodeIter).find("|");

			if(cPos != string::npos){
// cout << "node => " << *nodeIter << endl;
				modIter->second.tables.push_back(emptyTable);
				string nodeID = (*nodeIter).substr(1, cPos-1);
				cPos += 1;
				string substr = nodeIter->substr(cPos, (*nodeIter).length()-cPos-1);
				vector<string> parentIDs = data_manage::Stringmanip::split(substr, ':');
// 				calcParentCondTable(dset, *nodeIter);
// 				unsigned int nodeIndex = holder->getGenoIndex(nodeID.substr(1,nodeID.length()-cPos+1));
				modIter->second.tables.back().nodeIndex = holder->getGenoIndex(nodeID.substr(1,nodeID.length()-cPos+1));
				modIter->second.varsWithParents.insert(modIter->second.tables.back().nodeIndex);
// cout << nodeID.substr(1,nodeID.length()-cPos+1) <<  " nodeIndex=" << modIter->second.tables.back().nodeIndex << endl;
// 				vector<unsigned int> parentIndexes;
// 				modIter->second.tables[tableIndex].parentIndexes.clear();
				// split it here and get indexes for parent genotypes
				for(vector<string>::iterator pIter=parentIDs.begin(); pIter != parentIDs.end();
					++pIter){
// cout << *pIter << " = " << holder->getGenoIndex((*pIter).substr(1,(*pIter).length()-1)) << endl;
					modIter->second.tables.back().parentIndexes.push_back(holder->getGenoIndex((*pIter).substr(1,(*pIter).length()-1)));
				}
				
				vector<int> parentValues;
				int parentLevels = configParentData(parentValues, modIter->second.tables.back().parentIndexes, 
					dset, modIter->second.tables.back().cumulativeLevels);
				
				// for the conditional table for each identical combination of the parents
				// calculate the occurrence of the child 
				// will be 3 different for each parent combo
				// scale the child to 1.0 for the 3 occurrences
				int nodeLevels=3;
				vector<int> inner(parentLevels,0);
				vector<vector<int> > totals(nodeLevels, inner);
// 				vector<int> parentTotals(parentLevels, 0);
				vector<int> nodeTotals(nodeLevels, 0);
				
				// cycle through individuals and total 
				for(unsigned int i=0; i<dset->numInds(); i++){
					totals[((*dset)[i])->getGenotype(modIter->second.tables.back().nodeIndex)][parentValues[i]]++;
					nodeTotals[((*dset)[i])->getGenotype(modIter->second.tables.back().nodeIndex)]++;
				}
				
				vector<double> innerProbs(parentLevels, 0.0);
// 				vector<vector<double> > conditProbs(nodeLevels, innerProbs);
				modIter->second.tables.back().probs.assign(nodeLevels, innerProbs);
				// create conditional probability by dividing each cell total
				// by the number in that row
				for(int i=0; i<nodeLevels; i++){
					for(int j=0; j<parentLevels; j++){
						if(totals[i][j] > 0)
							modIter->second.tables.back().probs[i][j] = totals[i][j]/double(nodeTotals[i]);
						else
							modIter->second.tables.back().probs[i][j] = missValue;
// cout << "i=" << i << " j=" << j << " total=" << totals[i][j] << " nodeTotals=" << 
// nodeTotals[i] << " conditProbs=" << modIter->second.tables.back().probs[i][j]  << endl;
					}
				}
				
				// so the structure needs 1)a table for each node with parents
				// 2) the table will have the conditProbs
				// 3) also need the list of parent indexes
				// 4) need node index
				// 5) need way to convert parents to value (need cumulativeLevels to do that)
// 				struct conditTable{
// 					int nodeIndex;
// 					vector<int> parentIndexes;
// 					vector<vector<double> > probs;
// 					vector<int> cumulativeLevels;
// 				};
				
				
// 				modIter->second.condTable[]=calcCondTable(dset, *nodeIter);
// How to store the conditional table?  a tree structure?  top row is the child and then a parent
// at each level below?  -- Probably better and faster to do something like the 
// k2calcwithparent in BayesSolutionCreator
			}
		}		
	}	
}


///
/// calculates probability tables for every variable
/// currently assumes using only SNPs
///
void	GEDiscrimBayes::calcProbTables(Dataset* dset, map<unsigned int, vector<double> >& orphanProbs,
	data_manage::Dataholder* holder){
	
	orphanProbs.clear();
	for(unsigned int geno=0; geno < dset->numGenos(); geno++){
		vector<int> totals(3,0);
		vector<double> table(3,0.0);
		for(unsigned int ind=0; ind < dset->numInds(); ind++){
			totals[int(dset->getInd(ind)->getGenotype(geno))]++;
		}
		double allTotal = totals[0]+totals[1]+totals[2];
		for(unsigned int i=0; i<3; i++){
			table[i] = totals[i] / allTotal;
		}
		orphanProbs[geno]=table;
	}
}


///
/// Calculates conditional probability table for node that has parents
///
// void GEDiscrimBayes::calcParentCondTable(){
// 	
// }

///
/// Create parent data combination 
/// @param parentValues
/// @param parents contains indexes to parent genotypes
/// @returns number of different levels(factors) in the parent combined values
///
int GEDiscrimBayes::configParentData(vector<int>& parentValues, vector<unsigned int> &parents,
	Dataset* dSet, vector<int>& cumulativeLevels){
	// assume three levels (to hold SNP data)
	int constLevels = 3; // this has to be changed also
	
	// set number of levels for each parent
	vector<int> nLevels(parents.size(), 0);
	int nl = 1;
	for(size_t i=0; i<parents.size(); i++){
// 		nLevels[i] = parents[i]->getNumLevels(currentSet);
		nLevels[i] = constLevels;
		nl *= nLevels[i];
	}
	
	cumulativeLevels.assign(parents.size(), 1);
	for(unsigned int i=1; i<cumulativeLevels.size(); i++){
// 		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels;
		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels[i-1];
	}
// 	int nl = cumulativeLevels.back() * nLevels;
	deque<float> args;
	unsigned int nParents = parents.size();
	Individual* ind;
	
	for(unsigned int i=0; i < dSet->numInds(); i++){
// 		ind = (*currentSet)[i];
// 		IndividualTerm::setInd(ind);

		int value = 0;	
		for(unsigned int j=0; j < nParents; j++){
// 			value += parents[j]->evaluate(args) * cumulativeLevels[j];
			value += (*dSet)[i]->getGenotype(parents[j]) * cumulativeLevels[j];
		}
		parentValues.push_back(value);
	}
	return nl;
}


void GEDiscrimBayes::totalModels(Algorithm* alg,  map<string, modScores>& topModels, data_manage::Dataholder* holder,
		bool mapUsed, bool ottDummy, bool continMapUsed, bool caseMods){
// cout << "myRank=" << myRank << endl;		
	topModels.clear();	
	map<string, modScores> modelHolder;
	Population algPop = alg->getPopulation();
// cout << "best score = " << algPop[0]->fitness() << endl;
// cout << "number of models=" << algPop.numSolutions() << endl;
	for(size_t i=0; i<algPop.numSolutions(); i++){
		stringstream ss;
		try{
			GEObjective::outputEquation(ss, algPop[i], holder, mapUsed, ottDummy, continMapUsed);
// cout << ss.str() << " => ";
			string sortedStr = rewriteEquation(ss.str());
// cout << sortedStr << "\n";
			if(modelHolder.find(sortedStr) == modelHolder.end()){
				modelHolder[sortedStr].count = 1;
				modelHolder[sortedStr].score = algPop[i]->fitness();
			}
			else{
				modelHolder[sortedStr].count++;
			}
		}
		catch(AthenaExcept & ae){
			// skip incomplete models
// 			cout << "incomplete model -- not counted" << endl;
		}
	}
	
	// gather all models 
	#ifdef HAVE_CXX_MPI
// cout << "myRank=" << myRank << " about to call gatherModelInformation for caseModels=" << caseMods << endl;
		gatherModelInformation(modelHolder);
	#endif
	
	// insert into map with score as key
	map<float, vector<string> > sortedModels;
	for(map<string, modScores>::iterator iter=modelHolder.begin(); iter != modelHolder.end();
		++iter){
		sortedModels[iter->second.score].push_back(iter->first);
	}

#ifdef HAVE_CXX_MPI
if(myRank == 0){
#endif	
// cout << "myRank=" << myRank  << "sorted models" << endl;
	string outName;
	if(caseMods){
		outName = outputName + ".cv." + Stringmanip::numberToString(currCV) + ".case.uniq";
	}
	else{
		outName = outputName + ".cv." + Stringmanip::numberToString(currCV) + ".control.uniq";
	}
	ofstream os;
	os.open(outName.c_str(), ios::out);
	writeUniqueFiles(os,sortedModels,modelHolder);
	os.close();
#ifdef HAVE_CXX_MPI
}
#endif
	
// now iterate through models and save only the top N results
// cout << "Model\tscore\tcount\n";
// for(map<float, vector<string> >::reverse_iterator iter=sortedModels.rbegin(); iter!= sortedModels.rend();
// 	++iter){
// 	for(vector<string>::iterator viter=(iter->second).begin(); viter!=(iter->second).end();
// 		++viter){
// 		cout << *viter << "\t" << iter->first << "\t" << modelHolder[*viter].count <<"\n";
// 	}
// }

// 	int i=0;
// 	map<float, vector<string> >::iterator sortedIter=sortedModels.begin();
// 	while(i<topModelsUsed){
// 		for(vector<string>::iterator modIter=sortedIter->second.begin(); 
// 			modIter != sortedIter->second.end(); ++modIter){
// 			if(i++ < topModelsUsed){
// 				topModels[*modIter] = modelHolder[*modIter];			
// 			}
// 			else{
// 				break;
// 			}
// 		}
// 		++sortedIter;
// 	}
	
	int i=0;
	map<float, vector<string> >::reverse_iterator sortedIter=sortedModels.rbegin();
	while(i<topModelsUsed){
		for(vector<string>::iterator modIter=sortedIter->second.begin(); 
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
	
// cout << "myRank=" << myRank << " number top models=" << topModels.size() << endl;

}


///
/// Split model equation into tokens
///
vector<string> GEDiscrimBayes::tokenizeEquation(string equation){
	vector<string> tokens;
	// split string into tokens
	string tok;
	for(size_t i=0; i<equation.size(); i++){
		if(equation[i] == ']'){
			tok += equation[i];
			// end last token
// cout << "push back tok=" << tok << endl;
			if(tok.find(":") != string::npos){
				tok = sortParentStr(tok);
			}
			tokens.push_back(tok);
			tok="";
		}
		else{
			tok += equation[i];
		}
	}
	return tokens;
}


///
///  Writes list of all unique models and counts
/// 
void GEDiscrimBayes::writeUniqueFiles(ostream& outstream, map<float, vector<string> >& sortedModels,
	map<string, modScores>& modelHolder){

	outstream << "Network\tScore\tCount\n";
	for(map<float, vector<string> >::reverse_iterator sortedIter=sortedModels.rbegin();
		sortedIter != sortedModels.rend(); ++sortedIter){
// outstream << "num strings=" << sortedIter->second.size() << "\n";
		for(vector<string>::iterator strIter = sortedIter->second.begin(); strIter != sortedIter->second.end();
			++strIter){
			outstream << *strIter << "\t" << sortedIter->first << "\t" << modelHolder[*strIter].count << "\n";
		}
	}
}


///
/// Re-write equation string in sorted order so identical equations are matched
///
string GEDiscrimBayes::rewriteEquation(string equation){
	vector<string> tokens = tokenizeEquation(equation);
// split string into tokens
// 	string tok;
// 	for(size_t i=0; i<equation.size(); i++){
// 		if(equation[i] == ']'){
// 			tok += equation[i];
// 			// end last token
// // cout << "push back tok=" << tok << endl;
// 			if(tok.find(":") != string::npos){
// 				tok = sortParentStr(tok);
// 			}
// 			tokens.push_back(tok);
// 			tok="";
// 		}
// 		else{
// 			tok += equation[i];
// 		}
// 	}
	
	string sortedStr;
	//sort tokens
	sort(tokens.begin(), tokens.end());
	for(size_t i=0; i<tokens.size(); i++){
// 	for(vector<string>::iterator iter=tokens.begin(); iter != tokens.end(); ++iter){
		sortedStr+=tokens[i];
	}
// cout << "before " << equation << " => " << sortedStr << endl;
	return sortedStr;
}

///
/// sorts parent string in bayesian equation
/// @return string with parents sorted 
///
string GEDiscrimBayes::sortParentStr(string parentString){
// cout << "sorting parents " << parentString << endl;
	size_t startPos = parentString.find("|")+1;

	// skip the first '[' and last ']'
	string substr = parentString.substr(startPos, parentString.length()-startPos-1);
// cout << "substr=" << substr << endl;
	vector<string> tokens = data_manage::Stringmanip::split(substr, ':');
	sort(tokens.begin(), tokens.end());
	string sorted = parentString.substr(0,startPos);
	vector<string>::iterator iter=tokens.begin();
	sorted += *iter;
	++iter;
	for(; iter != tokens.end(); ++iter){
		sorted += ":" + *iter;
	}
	sorted += "]";
// cout << "after sorting " << sorted << endl;
  return sorted;
}

///
/// Run dataset against models contained in previously generated summary file
/// @param sumFile ATHENA summary file
///
vector<Solution*> GEDiscrimBayes::runValidation(std::string sumFile){
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
/// Finishes logs and compiles them into a single file
/// in cases where there are multiple processors
/// @param basename
/// @param cv
///
void GEDiscrimBayes::finishLog(string basename, int cv){
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
/// Returns graphical file extension appropriate to solution
/// @return extension
///
std::string GEDiscrimBayes::getGraphicalFileExt(){
	return GEObjective::getGraphicalExt();
}

///
/// Establishes the algorithm for the run based on the parameters
/// set in this algorithm.
///
void GEDiscrimBayes::run(){
		while(!ga->done()){
				ga->step();
		}
}

/// Writes graphical representation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GEDiscrimBayes::writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
	bool mapUsed, bool ottDummy, bool continMapUsed){
	GEObjective::outputModel(os, sol, holder, mapUsed, ottDummy, continMapUsed);
}


///
/// Writes model represented as an equation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GEDiscrimBayes::writeEquation(ostream& os, Solution* sol, data_manage::Dataholder* holder,
	bool mapUsed, bool ottDummy, bool continMapUsed){
	GEObjective::outputEquation(os, sol, holder, mapUsed, ottDummy, continMapUsed);
}

#ifdef HAVE_CXX_MPI

	void GEDiscrimBayes::gatherModelInformation(map<string, modScores>& models){
		// population size is maximum possible * number of nodes to hold them all
		//create array of models big enough to hold them all
// cout << "rank=" << myRank << " started gatherModelInformation" << endl;
		uniqueModelMPI empty;
		uniqueModelMPI* modelSend = new uniqueModelMPI[popSize];
// cout << "sizeof modelSend=" << sizeof(
		int i=0;
		for(map<string, modScores>::iterator iter=models.begin(); iter != models.end();
			++iter){
			modelSend[i].score=iter->second.score;
			modelSend[i].count=iter->second.count;
			strcpy(modelSend[i].modelEquation, iter->first.c_str());
// if(myRank==1){
// cout << "i=" << i << " eq=" << modelSend[i].modelEquation << " count=" <<
// modelSend[i].count << endl;
// }
			i++;
		}
// cout << "rank=" << myRank << " included " << i+1 << " models to send" << endl;
		while(i<popSize){
			modelSend[i].count=0;
			i++;
		}
// cout << "rank=" << myRank << " gatherModelInformation built send array" << endl;		
		uniqueModelMPI* modelRecv=NULL;
		int recvSize = popSize * totalNodes;
		if(myRank==0){
			modelRecv = new uniqueModelMPI[recvSize];
		}
// cout << "rank=" << myRank << " calling MPI_Gather " << endl;		
		MPI_Gather(modelSend, sizeof(empty)*popSize, MPI_BYTE, modelRecv, sizeof(empty)*popSize, MPI_BYTE, 0, MPI_COMM_WORLD);
// cout << "rank=" << myRank << " finished MPI_Gather " << endl;
		// merge all 
		if(myRank==0){
			models.clear();
// cout << "recvSize=" << recvSize << endl;
			for(int modIndex=0; modIndex < recvSize; modIndex++){
				if(modelRecv[modIndex].count > 0){
// cout << "modelEquation=" << modelRecv[modIndex].modelEquation << endl;
					if(models.find(modelRecv[modIndex].modelEquation) == models.end()){
						models[modelRecv[modIndex].modelEquation].count = modelRecv[modIndex].count;
						models[modelRecv[modIndex].modelEquation].score = modelRecv[modIndex].score;
					}
					else{
						models[modelRecv[modIndex].modelEquation].count += modelRecv[modIndex].count;
					}
				}
			}
// cout << "rank=" << myRank << " finished merging models " << endl;	
		}
		
		delete [] modelSend;
		if(myRank==0){
			delete [] modelRecv;
		}

	}

#endif
