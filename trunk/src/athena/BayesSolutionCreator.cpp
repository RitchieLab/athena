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
#include "BayesSolutionCreator.h"
#include "TerminalSymbol.h"
#include "Terminals.h"
#include <stack>
#include <cmath>
#include <Stringmanip.h>
#include "AthenaGrammarSI.h"

///
/// Constructor
///
BayesSolutionCreator::BayesSolutionCreator(){
		initialize();
}

///
/// Alternative constructor
///
BayesSolutionCreator::BayesSolutionCreator(vector<string>& symbols){
	initialize();
}


///
/// Initializes variables
///
void BayesSolutionCreator::initialize(){
	terminalsSet = false;
	calculator = NULL;
	graphicExtension = ".dot";
	
	leftOptBound = '(';
	rightOptBound = ')';
	numNodes=0;
}


///
/// Creates neural network for evaluation
/// sets up a postfix stack for evaluation of the neural network
/// @param pheno Phenotype from libGE that can be turned into a neural network
/// @param set Dataset
///
void BayesSolutionCreator::establishSolution(vector<string>& symbols, Dataset* set){
		if(!terminalsSet)
			setVariables(set);
			
		// sets values for this set for all variables
		// the map will be key of pointer and score values as the value of map
		// parentScoresSet will reset to false every the dataset changes		
		if(currentSet != set){
			currentSet = set;
			setParentScores(set);
		}
		
		establishSolution(symbols);
}

///
/// Creates neural network for output as equation
/// sets up a postfix stack for evaluation of the neural network
/// @param pheno Phenotype from libGE that can be turned into a neural network
/// @param set Dataset
///
void BayesSolutionCreator::establishSolutionEquation(vector<string>& symbols){
	establishSolution(symbols);
}


///
/// Creates Bayesian network for evaluation
/// @param symbols Phenotype from libGE that can be turned into a Bayesian network
///
void BayesSolutionCreator::establishSolution(vector<string>& symbols){
	
	int varCount=0;
	for(unsigned int i=0; i<symbols.size(); i++){
		if(symbols[i][0] == 'G'){
			varCount++;
		}
	}
	if(varCount > currentSet->numGenos() ){
		throw AthenaExcept("not enough variables");
	}

	// make list of children and parents of pheno (can't be both) -- so switch any children
	// to a second variable
	// also can't have duplicate parents or children -- switch those to new variable
	// need to propagate changes back to the original genome -- tricky part
  // change vector of symbols	 
	 
	set<TerminalSymbol*> included, parentIncluded;

	 stack<TerminalInfo> children;
	 
	 network.clearNodes();
	 
	 // pop from holder -- should be ')'
	 // grab the variables within -- all will go onto the children stack until reach '('
	 // if find a '->' then grab variables until ')' 
	 // point any of these to the child in question (top of children stack)
	 // after terminating '(' for children should have a list of children on stack

	int currPheno=int(symbols.size()-1);
	TerminalSymbol* currTerm = termHolder.getTerm(symbols[currPheno]);
	children.push(TerminalInfo(currTerm, currPheno));
	--currPheno;
	
	set<TerminalSymbol*> phenoChildren;
	// key is original phenotype position
	changedVariables.clear();
	
	while(currPheno >= 0){
		currTerm = termHolder.getTerm(symbols[currPheno]);	
		if(currTerm->getTermType() == TerminalSymbol::Genotype or 
			currTerm->getTermType() == TerminalSymbol::Covariate or
			currTerm->getTermType() == TerminalSymbol::Phenotype){
			children.push(TerminalInfo(currTerm, currPheno));
		}
		else if(currTerm == termHolder.connector()){
			// get parent(s) and apply to child/children
			// first will be ')'
			// grab variables until reach '('
			vector<TerminalInfo> parents;
			--currPheno;--currPheno;
			currTerm = termHolder.getTerm(symbols[currPheno]);
			while(currTerm != termHolder.leftParen()){
				parents.push_back(TerminalInfo(currTerm, currPheno));
				--currPheno;
				currTerm = termHolder.getTerm(symbols[currPheno]);
			}
			
			// here -- need to replace any duplicates with another variable for the parents
			if(parents.size() > 1){
				parentIncluded.clear();
				parentIncluded.insert(parents.back().term);
				for(int i=int(parents.size())-2; i>= 0; i--){
					if(parentIncluded.find(parents[i].term) != parentIncluded.end()){
						// need to get another parent to include
						TerminalSymbol* newTerm;
						int newCodon; 
						do{
 						 string newSymb = mapper->getNewVariable(newCodon);
 						 newTerm = termHolder.getTerm(newSymb);
						}while(parentIncluded.find(newTerm) != parentIncluded.end());
						parents[i].term=newTerm;
						parents[i].newValue = newCodon;
						changedVariables[parents[i].phenoIndex] = parents[i];
					}
					// whether changed or not insert here
					parentIncluded.insert(parents[i].term);
				}
			}
		
			// get children
			TerminalInfo child = children.top();
			children.pop();
			vector<TerminalInfo> tmpchildren;
			if(child.term == termHolder.leftParen()){
				child = children.top();
				children.pop();
				while(child.term != termHolder.rightParen()){
					tmpchildren.push_back(child);
					child = children.top();
					children.pop();
				}
			}
			else{
				tmpchildren.push_back(child);
			}	
			// replace any tmpchildren duplicates with another variable
			if(tmpchildren.size() > 1){
				included.clear();
				included.insert(tmpchildren.back().term);
				for(int i=int(tmpchildren.size()-2); i>=0; i--){
					if(included.find(tmpchildren[i].term) != included.end()){
						TerminalSymbol* newTerm; 
							int newCodon;
						do{
							string newSymb = mapper->getNewVariable(newCodon);
							newTerm = termHolder.getTerm(newSymb);
						}while(included.find(newTerm) != included.end());
 						tmpchildren[i].term=newTerm;
						tmpchildren[i].newValue=newCodon;
						changedVariables[tmpchildren[i].phenoIndex] = tmpchildren[i];
					}
 					// whether changed or not insert here
					included.insert(tmpchildren[i].term);			
				}	
			}

	
			// check that parents and children are not the same
			for(vector<TerminalInfo>::iterator parIter=parents.begin(); parIter != parents.end();
				++parIter){
				for(vector<TerminalInfo>::iterator childIter=tmpchildren.begin(); childIter != tmpchildren.end();
					++childIter){
					// pick new variable when a variable is pointing to itself
					if(parIter->term == childIter->term){
						TerminalSymbol* newTerm;
						int newCodon; 
						do{
 						 string newSymb = mapper->getNewVariable(newCodon);
 						 newTerm = termHolder.getTerm(newSymb);
						}while(parentIncluded.find(newTerm) != parentIncluded.end() || newTerm == childIter->term);
						parIter->term=newTerm;
						parIter->newValue = newCodon;
						changedVariables[parIter->phenoIndex] = *parIter;
					}
				}
			}
	
			// set parents to point to children
			// CODE HERE
			for(vector<TerminalInfo>::iterator parIter=parents.begin(); parIter != parents.end();
				++parIter){
				for(vector<TerminalInfo>::iterator childIter=tmpchildren.begin(); childIter != tmpchildren.end();
					++childIter){
						if(parIter->term->getTermType()==TerminalSymbol::Phenotype){
							phenoChildren.insert(childIter->term);
						}
						else if(childIter->term->getTermType()==TerminalSymbol::Phenotype){
							// can't be both a child and a parent of the phenotype
							if(phenoChildren.find(parIter->term) != phenoChildren.end()){
							// select new parent here to replace one already used
								TerminalSymbol* newTerm;
								int newCodon;
								do{
									string newSymb = mapper->getNewVariable(newCodon);
									newTerm = termHolder.getTerm(newSymb);
								}while(phenoChildren.find(newTerm) != phenoChildren.end());
								phenoChildren.insert(parIter->term);
								parIter->term=newTerm;
								parIter->newValue = newCodon;
								changedVariables[parIter->phenoIndex]=TerminalInfo(parIter->term,
									parIter->phenoIndex, parIter->newValue);
							}
						}
						else{ // non-pheno parent and child
							
						}						
						network.addConnection(parIter->term, childIter->term);
				}
			}
		

			// if only one child place back on stack
			if(tmpchildren.size()==1)
				children.push(tmpchildren[0]);
			// if next symbol is a connector -- push parent on stack
			if(currPheno-1 >= 0 && (termHolder.getTerm(symbols[currPheno-1]) == termHolder.connector())){
				children.push(parents[0]);
			}
			
		}
		else if(currTerm == termHolder.leftParen()){
			children.push(TerminalInfo(currTerm, currPheno));
		}
		--currPheno;
	}
 
		
	// iterate through stack and store covars and genotypes for use in checking
	// for missing data
	genos.clear();
	covars.clear();

	for(GraphNodeIter iter=network.begin(); iter != network.end(); ++iter){
		if((*iter)->term->getTermType() == TerminalSymbol::Genotype){
			genos[(*iter)->term]=1;
		}
		else if((*iter)->term->getTermType() == TerminalSymbol::Genotype){
			covars[(*iter)->term]=1;
		}
	}

}

///
/// placeholder for optimization method for Bayes Networks
/// @param symbols vector of strings the can be converted into network
/// @param set Dataset for optimizing solution
/// @returns indicator of optimization
///
int BayesSolutionCreator::optimizeSolution(std::vector<std::string>& symbols, Dataset* set){
	return 0;
}

///
/// Returns vector with indexes of SNPs in the latest solution
/// @return vector of indexes
///
vector<int> BayesSolutionCreator::getGeneIndexes(){
	vector<int> indexes;
	map<TerminalSymbol*, int>::iterator iter;
	GenotypeTerm* g;
	int i;
	for(iter = genos.begin(); iter != genos.end(); iter++){
		g = dynamic_cast<GenotypeTerm*>(iter->first);
		for(i=0; i<iter->second; i++){
			indexes.push_back(g->getIndex());
		}
	}
	return indexes;
}

///
/// Returns vector with indexes of Covariates in the latest solution
/// @return vector of indexes
///
inline vector<int> BayesSolutionCreator::getCovarIndexes(){
	vector<int> indexes;
	map<TerminalSymbol*, int>::iterator iter;
	ContinVariable* c;
	int i;
	for(iter = covars.begin(); iter != covars.end(); iter++){
		c = dynamic_cast<ContinVariable*>(iter->first);
		for(i=0; i<iter->second; i++){
			indexes.push_back(c->getIndex());
		}
	}  
	return indexes;
}


///
/// Creates variable terminals in terminal holder.
/// @param
///
void BayesSolutionCreator::setVariables(Dataset* set){
		unsigned int numCovars = (*set)[0]->numCovariates();
		unsigned int numGenos = (*set)[0]->numGenotypes();
		termHolder.createTerminals(numGenos, numCovars);
		terminalsSet = true;
}


///
/// Evaluates single individual and returns value
/// @param ind Individual to evaluate
/// returns value
///
float BayesSolutionCreator::evaluateInd(Individual* ind){

		termHolder.setInd(ind);
		
		return 0.0; // temporary placeholder
}


///
/// Checks to see if the individual has any missing data
/// @param ind
/// @param set
/// @return true if the individual should be used
///
bool BayesSolutionCreator::useInd(Individual* ind, Dataset* set){
		deque<float> args;
		
		map<TerminalSymbol*, int>::iterator iter;
 
		for(iter=covars.begin(); iter != covars.end(); iter++){
			if(iter->first->evaluate(args) == set->getMissingCoValue()){
				return false;
			}
		}
		
		for(iter=genos.begin(); iter != genos.end(); iter++){
			if(iter->first->evaluate(args) == set->getMissingGenotype()){
				return false;
			}
		}

		return true;
}



///
/// Evaluates the bayesian network and returns the K2 score
/// of the bayesian network.
/// @param set Dataset to use
///
float BayesSolutionCreator::evaluate(Dataset* set){
	 
	 nIndsEvaluated = set->numInds();
	 double score = 0.0;
	 vector<IndividualTerm*> parents;
	 calculator->reset();
	 double nParams,cumulativeBaseScore=0.0;
	 
	 
	 // now when a node has no parents no need to do anything as it will not
	 // change the score of the network.  It is the same as in the reference
	 // network where all nodes are included.
	 
	 for(GraphNodeIter iter=network.begin(); iter != network.end();
	 	++iter){
	 		// values for nodes without parents are cached to speed evaluation
	 		if((*iter)->parents.empty()){
				// set second argument to zero as the penalty is already applied
	 			calculator->addIndScore(parentScore[(*iter)->term],0);
	 			cumulativeBaseScore += parentScore[(*iter)->term];
	 		}
	 		else{
	 			parents.clear();
	 			for(GraphNodeIter parentIter=(*iter)->parents.begin();
	 				parentIter != (*iter)->parents.end(); ++parentIter){
	 					parents.push_back(static_cast<IndividualTerm*>((*parentIter)->term));
	 			}
	 			score = k2calcWithParent(static_cast<IndividualTerm*>((*iter)->term), parents, nParams);
	 			calculator->addIndScore(score, nParams);
	 			cumulativeBaseScore += parentScore[(*iter)->term];
	 		}
	 }

		// score will be the difference between the cumulativeBaseScore
		score = -(cumulativeBaseScore + calculator->getScore());	
		if(score < 0.0){
			score = 0.0;
		}
	
		return score;
}



///
/// For Networks, the detailed logging involves any 
/// additional processing needed.
///
void BayesSolutionCreator::detailedLogging(){
// 	ExpressionTree extree;
// 	extree.convertPostFix(postFixStack);
// 	nnDepth = extree.getMaxDepth();
}

///
/// Return detailed information 
///
unsigned int BayesSolutionCreator::getDetailedLog(){
// 		return nnDepth;
	return 0;
}

///
/// writes a dot compatible text file representing the network
/// @param os ostream to write to
///
void BayesSolutionCreator::graphicalOutput(ostream& os, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed){
	
	// network will be part of this tree
	network.outputDot(os,holder,mapUsed,ottDummy, continMapUsed);
}

void BayesSolutionCreator::equationOutput(ostream& os, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed){
			
		// traverse the network and output it
		// find the phenotype
		GraphNodeIter phenoIter = network.begin();
		for(; phenoIter != network.end(); ++phenoIter){
			if((*phenoIter)->term == termHolder.phenotype()){
				break;
			}
		}
		
		GraphNodeIter phenoParentsEnd = (*phenoIter)->parents.end();
		bool first = true;
		string starter = "(", finish="";
		// first output all the phenotype parents
		for(GraphNodeIter iter=(*phenoIter)->parents.begin(); iter != phenoParentsEnd; 
			++iter){
			if(first){
				os << starter;
				first=false;
			}
			os << getLabel(iter, holder) << " ";
			finish = ")->";
		}
		os << finish << "pheno";
		if((*phenoIter)->children.size() > 0){
			os << "->(";
			for(GraphNodeIter childIter=(*phenoIter)->children.begin(); 
				childIter != (*phenoIter)->children.end(); ++childIter){			
				first = true;
				finish ="";
				if((*childIter)->parents.size() > 1){
					finish = ")->";
					for(GraphNodeIter childParentIter=(*childIter)->parents.begin(); 
						childParentIter != (*childIter)->parents.end(); ++childParentIter){	
						if((*childParentIter)->term == (*phenoIter)->term)
							continue;
						if(first){
							os << starter;
							first=false;
						}
						os << getLabel(childParentIter, holder) << " ";
					}
					os << finish;
				}
				os << getLabel(childIter, holder) << " ";
			}
			os << ")";
		}
		
}

///
/// @param node GraphNodeIter
/// @returns variable index
///
string BayesSolutionCreator::getLabel(GraphNodeIter& node, Dataholder* holder){ 
	string name = (*node)->term->getName();
	stringstream ss(name.substr(1,name.length()-1));	
	int num;
	ss >> num;
	if(holder != NULL)
		return holder->getGenoName(num-1);
	else
		return ss.str();
}

/// 
/// Performs evaluation on the dataset.  Each individual has their result passed
/// to the output stream provided.
/// @param set Dataset 
/// @param os ostream 
/// @param results TestResult vector for holding result values
///
float BayesSolutionCreator::evaluateWithOutput(Dataset* set, ostream& os){

		// to check for missing data use the list of variables gotten from the
		// creation of the stack and then check each ind to  make sure
		// there is a value in each
		Individual * ind;

		calculator->reset();
		
		// when missing skip that ind
		for(unsigned int i=0; i < set->numInds(); i++){
				ind = (*set)[i];
							 
				ContinVariable::setInd(ind);
				GenotypeTerm::setInd(ind);       
				if(!useInd(ind, set)){
						os << ind->getID() << "\tMissing data\t" << ind->getStatus() << endl;
						continue;
				}
				float score = evaluateInd(ind);
				os << ind->getID() << "\t" << score << "\t" << ind->getStatus() << endl;
				calculator->addIndScore(score, ind->getStatus());
		}
 
		return calculator->getScore();
}


///
/// stores evaluation results as TestResult
/// @param set Dataset
/// @param results TestResult vector to contain values
///
void BayesSolutionCreator::evaluateForOutput(Dataset* set){
		Individual * ind;
		std::vector<stat::TestResult> results;
		
		calculator->reset();
		stat::TestResult tempResult;
		
		int missingInds=0;
		// when missing skip that ind
		for(unsigned int i=0; i < set->numInds(); i++){
				ind = (*set)[i];
							 
				ContinVariable::setInd(ind);
				GenotypeTerm::setInd(ind);       
				if(!useInd(ind, set)){
						missingInds++;
						continue;
				}
				
				tempResult.score = evaluateInd(ind);
				tempResult.status = ind->getStatus();
				results.push_back(tempResult);
		}
		float percentMissing=float(missingInds)/set->numInds() * 100.0;
		stringstream ss;
		ss << percentMissing;
		addOutputValues.clear();
		addOutputValues.push_back(ss.str() + "%");	
		calculator->evaluateAdditionalOutput(results);
}


///
/// Calculate/store scores for variables to use for cases
/// where there are no parents for the network.  Calculates overall
/// score for baseline network where no nodes are connected.
/// @param set Dataset
///
///
/// Calculate/store scores for variables to use for cases
/// where there are no parents for the network.  Calculates overall
/// score for baseline network where no nodes are connected.
/// @param set Dataset
///
void BayesSolutionCreator::setParentScores(Dataset* newSet){

	// check that all variables are genos (or any value that fits in a geno variable)
	if(newSet->numCovariates() > 0){
		throw AthenaExcept("Data set includes " + Stringmanip::numberToString(newSet->numCovariates()) + 
			"  continuous variables and GEBN requires that all values be in the data file.");
	}
	double total=0.0;	
	unsigned int numGenos = newSet->numGenos();
	parentScore.clear();
	
	TerminalSymbol * term;
	double nParams;
	
	for(unsigned int gIndex=0; gIndex < numGenos; gIndex++){
		string genoName = termHolder.getGenoName(gIndex);
		term = termHolder.getTerm(genoName);
 		calculator->reset();
		parentScore[term] = k2calcNoParent(gIndex, nParams);
		parentParams[term] = nParams;
		calculator->addIndScore(parentScore[term], nParams);
		// store with penalty - will be stored as the negative value (actual score)
		parentScore[term] = -calculator->getScore();
		total += parentScore[term];
	}
	
	term = termHolder.phenotype();
	// set value for Phenotype without parents
	calculator->reset();
	parentScore[term] = k2calcPhenoNoParent(nParams);
	parentParams[term] = nParams;
	calculator->addIndScore(parentScore[term], nParams);
	parentScore[term] = -calculator->getScore();
	total += parentScore[term];
	
	cout << "dataset total score=" << total << endl;
	newSet->setConstant(total);
	// base score will be stored as negative of the K2 (so it will be a positive number)
}


///
/// Calculate phenotype without parent score
/// @param nP (out) number of parameters in calculation
///
double BayesSolutionCreator::k2calcPhenoNoParent(double& nP){
	vector<unsigned int> totals(3,0);
	unsigned int nInds = currentSet->numInds();
	double result = 0.0;

	for(unsigned int i=0; i < nInds; i++){
		totals[int(currentSet->getInd(i)->getStatus())]++;
	}

	// alpha = 1
	// imaginary is simply number of levels of totals for K2
	double imaginary = double(totals.size());
	double alpha = 1.0;
	nP=-1;
	for(unsigned int i=0; i<totals.size(); i++){
		result += lgamma(totals[i]+alpha); // no need to get gamma of alpha (it is 0)
		// reduce imaginary when a cell is empty
		if(totals[i]==0)
			imaginary--;
		else
			nP++;
	}
	result += lgamma(imaginary) - lgamma(imaginary + nInds);
	return result;
	
}


///
/// Returns k2 score for node without a parent
/// @param gIndex genotype Index for Node
/// @param nP (out) number of parameters in calculation
/// @returns k2 score for node
///
double BayesSolutionCreator::k2calcNoParent(unsigned int gIndex, double& nP){
	
	vector<unsigned int> totals(3,0);
	unsigned int nInds = currentSet->numInds();
	double result = 0.0;
	
	for(unsigned int i=0; i < nInds; i++){
		totals[int(currentSet->getInd(i)->getGenotype(gIndex))]++;
	}
	
	nP =-1;
	// alpha = 1
	// imaginary is simply number of levels of totals for K2
	double imaginary = double(totals.size());
	double alpha = 1.0;
	for(unsigned int i=0; i<totals.size(); i++){
		result += lgamma(totals[i]+alpha); // no need to get gamma of alpha (it is 0)
		// reduce imaginary when a cell is empty
		if(totals[i]==0)
			imaginary--;
		else
			nP++;
	}
	result += lgamma(imaginary) - lgamma(imaginary + nInds);
	return result;
}

///
/// Create parent data combination 
/// @param parentValues
/// @param parents
/// @returns number of different levels(factors) in the parent combined values
///
int BayesSolutionCreator::configParentData(vector<int>& parentValues, vector<IndividualTerm*> &parents){
	// assume three levels (to hold SNP data)
	int nLevels = 3;
	vector<int> cumulativeLevels(parents.size(), 1);
	for(unsigned int i=1; i<cumulativeLevels.size(); i++){
		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels;
	}
	int nl = cumulativeLevels.back() * nLevels;
	deque<float> args;
	unsigned int nParents = parents.size();
	Individual* ind;
	
	for(unsigned int i=0; i < currentSet->numInds(); i++){
		ind = (*currentSet)[i];
		IndividualTerm::setInd(ind);

		int value = 0;	
		for(unsigned int j=0; j < nParents; j++){
			value += parents[j]->evaluate(args) * cumulativeLevels[j];
		}
		parentValues.push_back(value);
	}
	return nl;
}

///
/// Returns k2 score for node without a parent
/// @param childIndex genotype Index for Node
/// @param parents vector containing indices for parents
/// @nP (out) number of parameters used in calculation
/// @returns k2 score for node
///
double BayesSolutionCreator::k2calcWithParent(IndividualTerm* node, vector<IndividualTerm*> &parents,
	double& nP){

	// construct table with node values and the values for the parent combinations
	vector<int> parentValues;
	int parentLevels = configParentData(parentValues, parents);
	int nodeLevels = 3; // default value assuming SNP data 
	 
	int imaginary = 0;
	float alpha = 1.0;
	unsigned int numInds = currentSet->numInds();
	deque<float> args;
	
	vector<int> inner(parentLevels,0);
	vector<vector<int> > totals(nodeLevels, inner);
	vector<int> parentTotals(parentLevels, 0);
	vector<int> nodeTotals(nodeLevels, 0);
	Individual* ind;
	
	for(unsigned int i=0; i < numInds; i++){
		ind = (*currentSet)[i];
		IndividualTerm::setInd(ind);
		totals[node->evaluate(args)][parentValues[i]]++;
		parentTotals[parentValues[i]]++;
		nodeTotals[node->evaluate(args)]++;
	}
	
	// calculate conditional posterior probability
	double score = 0.0;
	for(int i=0; i<nodeLevels; i++){
		for(int j=0; j<parentLevels; j++){
			score += lgamma(totals[i][j] + alpha); // - lgamma(alpha) /* always zero for k2 */
		}
	}
	
	int finalParentLevels = parentLevels;
	for(int j=0; j < parentLevels; j++){
		if(parentTotals[j] == 0){
			finalParentLevels--;
		}
	}
	int finalNodeLevels=0;
	for(int j=0; j < nodeLevels; j++){
		if(nodeTotals[j] > 0)
			finalNodeLevels++;
	}

	imaginary = finalNodeLevels * finalParentLevels;
	nP = (finalNodeLevels-1) * finalParentLevels;
	
	for(int j=0; j < parentLevels; j++){
		if(parentTotals[j] > 0){
			score += lgamma(double(imaginary)/finalParentLevels) -
				lgamma(parentTotals[j] + double(imaginary)/finalParentLevels);
		}
	}
	return score;
}

void BayesSolutionCreator::setMapper(AthenaGrammarSI* m){
	mapper = m;
	mapper->setVariableSymbol("<v>");
}
