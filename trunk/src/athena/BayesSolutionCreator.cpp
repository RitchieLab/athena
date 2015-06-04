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
#include <KnuthComboGenerator.h>
#include <algorithm>
#include "MDR.h"

bool operator<(const OptScore& lhs, const OptScore& rhs){
      return lhs.key < rhs.key;
}

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
	currentSet = NULL;
	graphicExtension = ".dot";
	startOpt = "<v>";
// 	optSymbols.insert("<v>");
	
	leftOptBound = '(';
	rightOptBound = ')';
	numNodes=0;
	maximumSetSize = 500000;
}


///
/// Creates bayesian network for evaluation
/// @param pheno Phenotype from libGE that can be turned into a bayesian network
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
			scoreSet.clear();
			savedScores.clear();
// 			worstScore.key = "";
// 			worstScore.sc = 1.01;
		}
		
		establishSolution(symbols);
		repairNetwork(set);
}

///
/// Creates neural network for output as equation
/// sets up a postfix stack for evaluation of the neural network
/// @param pheno Phenotype from libGE that can be turned into a neural network
/// @param set Dataset
///
void BayesSolutionCreator::establishSolutionEquation(vector<string>& symbols){
	establishSolution(symbols);
	repairNetwork(currentSet);
}

///
/// Creates Bayesian network for evaluation
/// @param symbols Phenotype from libGE that can be turned into a Bayesian network
///
void BayesSolutionCreator::establishSolution(vector<string>& symbols){
	network.clearNodes();
	
	// convert symbols to DAGraph
// for(size_t i=0; i<symbols.size(); i++){
// cout << symbols[i] << " ";
// }
// cout << endl;
// throw AthenaExcept("Temporary placeholder");
// exit(1);
// stack based
// push individualterminals (contin or genos) on stack
// when reach a ')'
// pop off terminals until reach '('
// connect all those terminals as children of the next term on the stack
// push that parent term back on stack as it can be used as a child
// if the terminal is sc (means sibling connector) connect the next terminal
// with the previous one (both pushed back on stack)
	stack<TerminalSymbol*> terminals;
	TerminalSymbol* currTerm, *nextTerm, *stackTerm;
	// first terminal is a node that can be pushed on
	for(size_t symbIndex=0; symbIndex < symbols.size(); symbIndex++){
  	if(symbols[symbIndex]==" "){
// 			cout << "skip blank\n";
      continue;
  	}
		TerminalSymbol* currTerm = termHolder.getTerm(symbols[symbIndex]);
// cout << "currTerm is " << currTerm->getName() << endl;
		// when it is a connector connect the next terminal with the top
		// of the stack
		if(currTerm == termHolder.connector()){
			stackTerm = terminals.top(); // it will stay on the stack for future operations
			nextTerm = termHolder.getTerm(symbols[++symbIndex]);
			// add network connection from stackTerm to nextTerm
// cout << "add connection " << stackTerm->getName() << " to " << nextTerm->getName() << endl;
			network.addConnection(stackTerm, nextTerm);
			// add next Term to the stack
			terminals.push(nextTerm);
		}
		else if(currTerm == termHolder.rightParen()){
			// pop off terminals until reach left paren
			// then add connections for each 
			vector<TerminalSymbol*> toNodes;
			stackTerm = terminals.top();
			terminals.pop();
			while(stackTerm != termHolder.leftParen()){
				toNodes.push_back(stackTerm);
				stackTerm = terminals.top();
				terminals.pop();
			}
			// get from Node
			stackTerm = terminals.top();
			// add all connections
// cout << "add connection " << stackTerm->getName() << " to " << currTerm->getName() << endl;
// 			network.addConnection(stackTerm, currTerm);
			for(vector<TerminalSymbol*>::iterator iter=toNodes.begin(); iter != toNodes.end();
				iter++){
// cout << "add connection " << stackTerm->getName() << " to " << (*iter)->getName() << endl;
				network.addConnection(stackTerm, *iter);
			}
		}
		else if(currTerm->getTermType() == TerminalSymbol::Genotype ||
			currTerm->getTermType() == TerminalSymbol::Covariate ||
			currTerm == termHolder.leftParen()){
// cout << "push " << currTerm->getName() << endl;
				terminals.push(currTerm);
		} 
		// there is a final type of terminal 'nc' that denotes no connection.  it does not
		// need to be pushed onto the stack
	}
	
	// iterate through stack and store covars and genotypes for use in checking
	// for missing data
	genos.clear();
	covars.clear();

	for(GraphNodeIter iter=network.begin(); iter != network.end(); ++iter){
		if((*iter)->term->getTermType() == TerminalSymbol::Genotype){
			genos[(*iter)->term]=1;
		}
		else if((*iter)->term->getTermType() == TerminalSymbol::Covariate){
			covars[(*iter)->term]=1;
		}
	} 

// network.outputRstring(cout);
// cout << endl;

}



///
/// Creates Bayesian network for evaluation
/// @param symbols Phenotype from libGE that can be turned into a Bayesian network
///
// void BayesSolutionCreator::establishSolution(vector<string>& symbols){

// no longer needed with the repair network approach	
// 	unsigned int varCount=0;
// 	for(unsigned int i=0; i<symbols.size(); i++){
// 		if(symbols[i][0] == 'G' || symbols[i][0] == 'C'){
// 			varCount++;
// 		}
// 	}
// 	if(varCount > currentSet->numGenos() + currentSet->numCovariates()){
// 		throw AthenaExcept("not enough variables");
// 	}

// symbols.clear();
// symbols.push_back("(");
// symbols.push_back("G1");
// symbols.push_back("G5");
// symbols.push_back("G11");
// symbols.push_back(")");
// symbols.push_back("^");
// symbols.push_back("(");
// symbols.push_back("pheno");
// symbols.push_back(")");
// symbols.push_back("^");
// symbols.push_back("(");
// symbols.push_back("G11");
// symbols.push_back(")");

// 
// cout << "START:";
// for(unsigned int i=0; i<symbols.size(); i++){
// cout << symbols[i];
// }
// cout << endl;

	// make list of children and parents of pheno (can't be both) -- so switch any children
	// to a second variable
	// also can't have duplicate parents or children -- switch those to new variable
	// need to propagate changes back to the original genome -- tricky part
  // change vector of symbols	 
	 
// 	set<TerminalSymbol*> included, parentIncluded;
// 
// 	 stack<TerminalInfo> children;
// 	 
// 	 network.clearNodes();
	 
// testing the tarjan's implementation
// string g1="G1", g2="G2", g3="G3", g4="G4";
// TerminalSymbol* G1Term = termHolder.getTerm(g1);
// TerminalSymbol* G2Term = termHolder.getTerm(g2);
// TerminalSymbol* G3Term = termHolder.getTerm(g3);
// TerminalSymbol* G4Term = termHolder.getTerm(g4);
// network.addConnection(G2Term, G1Term);
// network.addConnection(G1Term, G3Term);
// network.addConnection(G3Term, G2Term);
// network.addConnection(G1Term, G4Term);
// network.SCC();

// cout << "fixed" << endl;
// network.removeConnection(G3Term, G2Term);
// network.SCC();

// exit(1);	 
	 
/*	 
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
	
	
			set<TerminalSymbol*> parentSet;
			for(vector<TerminalInfo>::iterator parIter=parents.begin(); parIter != parents.end();
				++parIter){
				parentSet.insert(parIter->term);
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
// cout << "selected newSymb=" << newSymb << endl;
									newTerm = termHolder.getTerm(newSymb);
								}while(phenoChildren.find(newTerm) != phenoChildren.end() ||
									parentSet.find(newTerm) != parentSet.end());
								phenoChildren.insert(parIter->term);
								parIter->term=newTerm;
								parIter->newValue = newCodon;
								parentSet.insert(parIter->term);
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
		else if((*iter)->term->getTermType() == TerminalSymbol::Covariate){
			covars[(*iter)->term]=1;
		}
	}
*/
// }

///
/// placeholder for optimization method for Bayes Networks
/// @param symbols vector of strings the can be converted into network
/// @param set Dataset for optimizing solution
/// @returns indicator of optimization
///
int BayesSolutionCreator::optimizeSolution(std::vector<std::string>& symbols, Dataset* set){
	// clear old symbols
	optValSymbols.clear();
	
// 	unsigned int varCount=0;
// 	for(unsigned int i=0; i<symbols.size(); i++){
// 		if(symbols[i][0] == 'G' || symbols[i][0] == 'C'){
// 			varCount++;
// 		}
// 	}
// 	if(varCount > currentSet->numGenos() + currentSet->numCovariates()){
// 		throw AthenaExcept("not enough variables");
// 	}
	

	vector<IndividualTerm*> variables;
	std::set<IndividualTerm*> phenoOrigParents;
	size_t combinationSize;
	
	for(size_t i=0; i<symbols.size(); i++){
// cout << symbols[i] << " ";
		TerminalSymbol* currTerm = termHolder.getTerm(symbols[i]);
		if(currTerm->getTermType() == TerminalSymbol::Genotype ||
			currTerm->getTermType() == TerminalSymbol::Covariate){
			variables.push_back( dynamic_cast<IndividualTerm*>(currTerm));
		}
		else if(currTerm->getTermType() == TerminalSymbol::Phenotype){
			combinationSize = variables.size();
			for(vector<IndividualTerm*>::iterator iter=variables.begin(); iter!=variables.end();
				++iter){
				phenoOrigParents.insert(*iter);
			}
		}
	}
// cout << endl;
	// no need to try combinations for balanced accuracy as no other combinations possible
	if(combinationSize == variables.size()){
		optimizedScore = 1e9;
		return 0;
	}
	
	// set up combination generator
	stat::KnuthComboGenerator generator;
	// set size of models to generate
  generator.ComboEnds(combinationSize, combinationSize);
  // set number of loci to use in generating the combinations
  generator.SetLoci(variables.size());
  generator.SetComboInterval(10000);
  generator.initialize_state();
  bool combosDone = false;
  vector<IndividualTerm*> modelTerms(combinationSize, NULL), bestTerms;
  float bestBalAcc=0.0;
	
  while(!combosDone){
  	combosDone = generator.GenerateCombinations();
  	for(vector<vector<unsigned int> >::iterator iter=generator.ComboList.begin(); iter != generator.ComboList.end();
  		++iter){
  		for(size_t i=0; i<iter->size(); i++){
  			modelTerms[i]=variables[(*iter)[i]-1];
  		} 		
//     float balAcc = MDR::calcBalAccuracy(currentSet, currentSet, modelTerms);
    float balAcc = getBalAccuracy(currentSet, modelTerms);
    if(balAcc > bestBalAcc){
    	bestBalAcc = balAcc;
    	bestTerms = modelTerms;
    }
// cout << "balAcc=" << balAcc << endl;
  	}
  }

	modelTerms.clear();
	for(vector<IndividualTerm*>::iterator modelIter=bestTerms.begin(); modelIter != bestTerms.end();
		++modelIter){
		std::set<IndividualTerm*>::iterator origIter = phenoOrigParents.find(*modelIter);
		if(origIter != phenoOrigParents.end()){
			phenoOrigParents.erase(origIter);
		}
		else{
			modelTerms.push_back(*modelIter);
		}
		
	}

// for(size_t i=0; i<bestTerms.size(); i++){
// cout << bestTerms[i]->getLabel() << " ";
// }
// cout << " bestBalAcc=" << bestBalAcc << endl;
	// if find one of best terms and it is not in original pheno Parent terms
	// swap it with one of the pheno parents in the variables array
	for(vector<IndividualTerm*>::iterator modelIter=modelTerms.begin(); modelIter != modelTerms.end();
		++modelIter){
// 		std::set<IndividualTerm*>::iterator origIter = phenoOrigParents.find(*modelIter);
		// model term not in original phenotype parents
// 		if(origIter == phenoOrigParents.end()){
			std::set<IndividualTerm*>::iterator swapIter = phenoOrigParents.begin();
			
// cout << "new var=" << (*modelIter)->getLabel() << " swap var=" << (*swapIter)->getLabel() << endl;
			// in variables find both new and old variables
			vector<IndividualTerm*>::iterator v1 =  find(variables.begin(), variables.end(), *modelIter);
			IndividualTerm* tempTerm = *v1;
			
			vector<vector<IndividualTerm*>::iterator> v1iters;
			vector<vector<IndividualTerm*>::iterator> v2iters;
			
			while(v1 != variables.end()){
// cout << (*v1)->getLabel() << endl;
// cout << "v1=" << (*v1)->getLabel() << " v2=" << (*v2)->getLabel() << endl;
// 				*v1 = *v2;
// 				*v2 = tempTerm;
// cout << "v1=" << (*v1)->getLabel() << " v2=" << (*v2)->getLabel() << endl;				
// cout << "look for another " << (*modelIter)->getLabel() << endl;
				v1iters.push_back(v1);
				v1 = find(v1+1, variables.end(), *modelIter);
			}
			
			vector<IndividualTerm*>::iterator v2 = find(variables.begin(), variables.end(), *swapIter);
			while(v2 != variables.end()){
// cout << "v1=" << (*v1)->getLabel() << " v2=" << (*v2)->getLabel() << endl;
// 				*v1 = *v2;
// 				*v2 = tempTerm;
// cout << "v1=" << (*v1)->getLabel() << " v2=" << (*v2)->getLabel() << endl;				
// cout << "look for another " << (*modelIter)->getLabel() << endl;
				v2iters.push_back(v2);
				v2 = find(v2+1, variables.end(), *swapIter);				
			}
			IndividualTerm* newv2 = *(v1iters[0]);
			IndividualTerm* newv1 = *(v2iters[0]);
			
		  for(vector<vector<IndividualTerm*>::iterator>::iterator iter=v1iters.begin();
		  	iter != v1iters.end(); ++iter){
		  	*(*iter) = newv1;
		  }
		  for(vector<vector<IndividualTerm*>::iterator>::iterator iter=v2iters.begin();
		  	iter != v2iters.end(); ++iter){
		  	*(*iter) = newv2;
		  }		  
			
			phenoOrigParents.erase(swapIter);
// 		}
	}

	// swap the best parents with variable position
	symbVector newSymbols;
	optSymbol tempSymbol;
	newSymbols.push_back(tempSymbol);
	tempSymbol.noNT=false;
	for(vector<IndividualTerm*>::iterator iter=variables.begin(); iter != variables.end();
		++iter){		
		newSymbols[0].symbol = (*iter)->getLabel();
		optValSymbols.push_back(newSymbols);
// cout << "newSymbols=>" << newSymbols[0].symbol << endl;
	}
// cout << "optValSymbols.size=" << optValSymbols.size() << endl;
	// force replacement in original network
	// switch this later to reflect need or lack thereof to change
	optimizedScore = 0.0;
	
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
/// Repairs network by removing loops
///
void BayesSolutionCreator::repairNetwork(Dataset* dSet){
	while(network.SCC()){
		// conduct repair
		std::vector<std::vector<GraphNode*> > loops = network.getLoops();
		
		TerminalSymbol* from, *to; 
		TerminalSymbol* fromBreak, *toBreak;
		double mutualInfo, worstMI;
		
		std::vector<std::vector<GraphNode*> >::iterator loopIter;
// cout << "number of loops=" << loops.size() << endl;
		for(loopIter = loops.begin(); loopIter != loops.end(); ++loopIter){
			worstMI = 1e10;
			set<TerminalSymbol*> loopTerminals;
// cout << "LOOP terminals: ";
			for(size_t i=0; i<(*loopIter).size(); ++i){
				// store all terminals that are part of loop
// cout << (*loopIter)[i]->term->getName() << " ";
				loopTerminals.insert((*loopIter)[i]->term);
			}
// cout << "\n";
				
			// for each terminal in the set calculate MI for each connection
			// to another terminal in the set -- future optimization is don't repeat calculations
			for(size_t i=0; i<(*loopIter).size(); ++i){
				for(set<GraphNode*>::iterator childIter=(*loopIter)[i]->children.begin(); childIter != (*loopIter)[i]->children.end();
					++childIter){
						
						// calculate only if child is in the loop
						if(loopTerminals.find((*childIter)->term)!=loopTerminals.end()){
							from = ((*loopIter)[i])->term;
							to = (*childIter)->term;
							mutualInfo = calcMI(from, to, dSet);
							if(mutualInfo < worstMI){
								toBreak = to;
								fromBreak = from;
								worstMI = mutualInfo;
							}
						}
					}
				}
				
				// the two nodes will be the from node and the to node before it
				// index 1 is the from node and index 0 is the to node for example
// 				from = ((*loopIter)[i])->term;
// 				to = ((*loopIter)[i-1])->term;
// 				mutualInfo = calcMI(from, to, dSet);
// cout << from->getLabel() << "->" << to->getLabel() << " ";
// // cout << " worstMI=" << worstMI << " mutualInfo=" << mutualInfo << endl;
// 				if(mutualInfo < worstMI){
// 					toBreak = to;
// 					fromBreak = from;
// 					worstMI = mutualInfo;
// 				}
// 			}
// cout << "removing " << toBreak->getLabel() << " from " << fromBreak->getLabel() << endl;
			network.removeConnection(fromBreak, toBreak);
// network.outputRstring(cout);
// cout << "\n";
		}
	}

// cout << "EXIT repairNetwork" << endl;
}


///
/// Calculate mutual information for a pair of variables
/// @param v1 
/// @param v2
///
double BayesSolutionCreator::calcMI(TerminalSymbol* v1, TerminalSymbol* v2, Dataset* dSet){
	double mutualInfo=0.0;
	
	vector<IndividualTerm*> terms(2,NULL);
	terms[0] = static_cast<IndividualTerm*>(v1);
	terms[1] = static_cast<IndividualTerm*>(v2);
// 	int parentLevels = configParentData(combinedValues, terms);
	
	// need total of each individual 
	// need total of combination
	// loop through each individuals and add the individual totals
	// add the combination total
	// assume three levels (to hold SNP data)
// 	int nLevels = 3; // this has to be changed also
	
	// set number of levels for each parent
	vector<int> nLevels(terms.size(), 0);
	int nl = 1;
	for(size_t i=0; i<terms.size(); i++){
		nLevels[i] = terms[i]->getNumLevels(currentSet);
		nl *= nLevels[i];
	}

	vector<int> combinedTotals(nl,0);
	size_t v1Levels = terms[0]->getNumLevels(currentSet);
	size_t v2Levels = terms[1]->getNumLevels(currentSet);
	vector<int> v1Totals(v1Levels, 0);
	vector<int> v2Totals(v2Levels, 0);

	
	vector<int> cumulativeLevels(terms.size(), 1);
	for(unsigned int i=1; i<cumulativeLevels.size(); i++){
		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels[i-1];
	}

	deque<float> args;
// 	unsigned int nTerms = terms.size();
	Individual* ind;
	unsigned int N = currentSet->numInds();
	

	
	for(unsigned int i=0; i < currentSet->numInds(); i++){
		ind = (*currentSet)[i];
		IndividualTerm::setInd(ind);

		int value = 0;
// 		for(unsigned int j=0; j < terms; j++){
		value += terms[0]->evaluate(args) * cumulativeLevels[0];
		v1Totals[terms[0]->evaluate(args)]++;
			
// 		}
		value += terms[1]->evaluate(args) * cumulativeLevels[1];
		v2Totals[terms[1]->evaluate(args)]++;
		
// 		parentValues.push_back(value);
		combinedTotals[value]++;
	}

// cout << "N=" << N << endl;
// cout << "\t0\t1\t2\n";
// cout << v1->getName();
// for(size_t i=0; i<v1Totals.size(); i++){
// 	cout << "\t" << v1Totals[i];
// }	
// cout << "\n";
// cout << v2->getName();
// for(size_t i=0; i<v2Totals.size(); i++){
// 	cout << "\t" << v2Totals[i];
// }
// cout << "\n";
// cout << "totals=" << N << endl;
// for(size_t i=0; i<v1Totals.size(); i++){
// 	for(size_t j=0; j<v2Totals.size(); j++){
// 		cout << "i=" << i << " j=" << j;
// 		int index = i * cumulativeLevels[0] + j * cumulativeLevels[1];
// 		cout << " tot=" << combinedTotals[index] << endl;
// 	}
// }
	
	// calculate mutual information
	// for each combination 
	int index;
	for(unsigned int i=0; i<v1Levels; i++){
		for(unsigned int j=0; j<v2Levels; j++){
			index = i * cumulativeLevels[0] + j * cumulativeLevels[1];
// cout << "v1 " << i << "=" << v1Totals[i] << endl;
// cout << "v2 " << j << "=" << v1Totals[j] << endl;
			if(combinedTotals[index] > 0){
				mutualInfo += combinedTotals[index] / double(N) * log(N*combinedTotals[index]/double(v1Totals[i]*v2Totals[j]));
			}
// cout << "i=" << i << " j=" << j << " tot=" << combinedTotals[index] << " mi=" << mutualInfo << endl;
// cout << "first part=" << combinedTotals[index] / double(N) << endl;
// cout <<  "second part without log=" << N*combinedTotals[index]/double(v1Totals[i]*v2Totals[j]) << endl;
// cout << "second part=" << log(N*combinedTotals[index]/double(v1Totals[i]*v2Totals[j])) << endl;
		}
	}
	return mutualInfo;
}



///
/// Evaluates the bayesian network and returns the K2 score
/// of the bayesian network.
/// @param set Dataset to use
///
float BayesSolutionCreator::evaluate(Dataset* set){

	 // when only one variable in network set to worst score
	 if(genos.size() + covars.size()==1){
// 	 cout << "only one-- worst possible" << endl;
		 return getWorst();
	 }
	 
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
// cout << "score for " << (*iter)->term->getName() << " with  parents: ";
// for(unsigned int i=0; i<parents.size(); i++){
// cout << parents[i]->getName() << " ";
// }
// cout << " is " << score << endl;
	 			calculator->addIndScore(score, nParams);
	 			cumulativeBaseScore += parentScore[(*iter)->term];
	 		}
	 }

// change score to be the - of the actual score
// include the total score 
// cout << "cumulativeBaseScore=" << cumulativeBaseScore << " actual score=" << calculator->getScore() << endl;
// cout <<  "total=" << set->getConstant() << " difference is " << cumulativeBaseScore + calculator->getScore() << endl;
// cout << "------------------------------------" << endl;
// so all I have to do is add the difference to the original and then take the - of it
// so that better scores are more positive
// at the end of the run convert the scores back to negative with simple -
// need to make the scores positive during fitness so that roulette wheel will work

/// ADD output to indicate any cases where all the models are worse than the original
/// IN the sum file

/// one check is to throw out variables without any variants
/// check for standard deviation of the column and throw out if zero
/// print warning and throw out

		// score will be the sum of total score plus the difference of the two
		// It will be scaled so that the better scores are greater for the fitness
		// The final score report will alter it so that the scores are the correct negative
		// values
		double diff = cumulativeBaseScore + calculator->getScore();	
		score = -(set->getConstant() + diff);
		if(score < 0.0){
			score = 0.0;
		}
// cout << "diff=" << diff << " final network score=" << score << endl;
// exit(1);
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
 		string genoPrefix="", continPrefix="";
 				
 		for(GraphNodeIter netIter = network.begin(); netIter != network.end(); ++netIter){
 			os << "[" << getLabel(netIter, holder,mapUsed, continMapUsed);
 			if(!(*netIter)->parents.empty()){
 				os << "|";
 				GraphNodeIter parentIter=(*netIter)->parents.begin();
 				os << getLabel(parentIter, holder, mapUsed, continMapUsed);
 				++parentIter;
				for(;parentIter != (*netIter)->parents.end(); ++parentIter){
					os << ":" << getLabel(parentIter, holder, mapUsed, continMapUsed);
				}
 			}
 			os << "]";
 		}	
}

// void BayesSolutionCreator::equationOutput(ostream& os, data_manage::Dataholder* holder,
// 			bool mapUsed, bool ottDummy, bool continMapUsed){
// 			
// 		// traverse the network and output it
// 		// find the phenotype
// 		GraphNodeIter phenoIter = network.begin();
// 		for(; phenoIter != network.end(); ++phenoIter){
// 			if((*phenoIter)->term == termHolder.phenotype()){
// 				break;
// 			}
// 		}
// 		
// 		GraphNodeIter phenoParentsEnd = (*phenoIter)->parents.end();
// 		bool first = true;
// 		string starter = "(", finish="";
// 		// first output all the phenotype parents
// 		for(GraphNodeIter iter=(*phenoIter)->parents.begin(); iter != phenoParentsEnd; 
// 			++iter){
// 			if(first){
// 				os << starter;
// 				first=false;
// 			}
// 			os << getLabel(iter, holder) << " ";
// 			finish = ")->";
// 		}
// 		os << finish << "pheno";
// 		if((*phenoIter)->children.size() > 0){
// 			os << "->(";
// 			for(GraphNodeIter childIter=(*phenoIter)->children.begin(); 
// 				childIter != (*phenoIter)->children.end(); ++childIter){			
// 				first = true;
// 				finish ="";
// 				if((*childIter)->parents.size() > 1){
// 					finish = ")->";
// 					for(GraphNodeIter childParentIter=(*childIter)->parents.begin(); 
// 						childParentIter != (*childIter)->parents.end(); ++childParentIter){	
// 						if((*childParentIter)->term == (*phenoIter)->term)
// 							continue;
// 						if(first){
// 							os << starter;
// 							first=false;
// 						}
// 						os << getLabel(childParentIter, holder) << " ";
// 					}
// 					os << finish;
// 				}
// 				os << getLabel(childIter, holder) << " ";
// 			}
// 			os << ")";
// 		}
// 		
// }

///
/// @param node GraphNodeIter
/// @returns variable index
///
string BayesSolutionCreator::getLabel(GraphNodeIter& node, Dataholder* holder,
	bool mapUsed, bool continMapUsed){ 
	string name = (*node)->term->getName();
	stringstream ss(name.substr(1,name.length()-1));	
	int num;
	ss >> num;
// cout << "holder=" << holder << " genoname="<< holder->getGenoName(num-1) << endl;
	if(holder != NULL)
		if((*node)->term->getTermType()==TerminalSymbol::Genotype)
			if(!mapUsed)
				return "G" +  holder->getGenoName(num-1);
			else
				return holder->getGenoName(num-1);
		else
			if(!continMapUsed)
				return "C" + holder->getCovarName(num-1);
			else
				return holder->getCovarName(num-1);
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
/// check whether the current solution is better than the default network
/// @param set Dataset
/// @param results TestResult vector to contain values
///
void BayesSolutionCreator::evaluateForOutput(Dataset* set){
	evaluateForOutput(set, set);
}

///
/// check whether the current solution is better than the default network
/// @param set Dataset
/// @param set Reference set that determines the conditional probabilities
/// @param results TestResult vector to contain values
///
void BayesSolutionCreator::evaluateForOutput(Dataset* set, Dataset* refSet){
	float score = evaluate(set);
	
	addOutputValues.clear();
	
	if(score + set->getConstant() > 0){	
		addOutputValues.push_back("+");
	} 
	else{
		addOutputValues.push_back("-");
	}

// equationOutput(cout,NULL,false,false,false);
// cout << endl;
/*
	// calculate balanced accuracy
	// find pheno
	vector<IndividualTerm*> parents;
	for(GraphNodeIter iter=network.begin(); iter != network.end(); ++iter){
		if((*iter)->term->getTermType() == TerminalSymbol::Phenotype){
			for(GraphNodeIter parIter=(*iter)->parents.begin(); parIter != (*iter)->parents.end();
				++parIter){
				parents.push_back(dynamic_cast<IndividualTerm*>((*parIter)->term));
			}
			break;
		}
	}

	float balAcc = MDR::calcBalAccuracy(currentSet, refSet, parents);
// cout << "in evaluateForOutput balAcc=" << balAcc << endl;
	addOutputValues.push_back(Stringmanip::numberToString(balAcc));
*/
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
// 	if(newSet->numCovariates() > 0){
// 		throw AthenaExcept("Data set includes " + Stringmanip::numberToString(newSet->numCovariates()) + 
// 			"  continuous variables and GEBN requires that all values be in the data file.");
// 	}
	
	double total=0.0;	
	unsigned int numGenos = newSet->numGenos(), numContins = newSet->numCovariates();
	parentScore.clear();
	
	TerminalSymbol * term;
	double nParams;
	
	// calculate genos scores when no parent
	for(unsigned int gIndex=0; gIndex < numGenos; gIndex++){
		string genoName = termHolder.getGenoName(gIndex);
		term = termHolder.getTerm(genoName);
 		calculator->reset();
		parentScore[term] = k2calcNoParent(gIndex, nParams);
		parentParams[term] = nParams;
		calculator->addIndScore(parentScore[term], nParams);
		// store with penalty - will be stored as the negative value (actual score)
		parentScore[term] = -calculator->getScore();
// cout << term->getName() << " score=" << parentScore[term] << endl;
		total += parentScore[term];
	}


	// calculate continuous variable scores when no parent
	for(unsigned int cIndex=0; cIndex < numContins; cIndex++){
		string continName = termHolder.getContinName(cIndex);
		term = termHolder.getTerm(continName);
		calculator->reset();
		parentScore[term]=k2calcNoParentContin(cIndex, nParams);
		parentParams[term]=nParams;
		calculator->addIndScore(parentScore[term], nParams);
		parentScore[term] = -calculator->getScore();
// cout << term->getName() << " score=" << parentScore[term] << endl;
		total += parentScore[term];	
	}
// cout << "total=" << total << endl;

// 	term = termHolder.phenotype();
	// set value for Phenotype without parents
// 	calculator->reset();
// 	parentScore[term] = k2calcPhenoNoParent(nParams);
// 	parentParams[term] = nParams;
// 	calculator->addIndScore(parentScore[term], nParams);
// 	parentScore[term] = -calculator->getScore();
// 	total += parentScore[term];
	newSet->setConstant(total);
}


///
/// Calculate phenotype without parent score
/// @param nP (out) number of parameters in calculation
///
double BayesSolutionCreator::k2calcPhenoNoParent(double& nP){
	vector<unsigned int> totals(currentSet->getNumStatusLevels(),0);
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
/// Returns k2 score for node without a parent
/// @param cIndex continuous variable Index for Node
/// @param nP (out) number of parameters in calculation
/// @returns k2 score for node
///
double BayesSolutionCreator::k2calcNoParentContin(unsigned int cIndex, double& nP){
	
	// create empty vector to hold totals
	vector<unsigned int> totals(currentSet->getNumLevels(cIndex),0);
	unsigned int nInds = currentSet->numInds();
	double result = 0.0;
	
	for(unsigned int i=0; i < nInds; i++){
		totals[int(currentSet->getInd(i)->getCovariate(cIndex))]++;
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
// 	int nLevels = 3; // this has to be changed also
	
	// set number of levels for each parent
	vector<int> nLevels(parents.size(), 0);
	int nl = 1;
	for(size_t i=0; i<parents.size(); i++){
		nLevels[i] = parents[i]->getNumLevels(currentSet);
		nl *= nLevels[i];
	}
	
	vector<int> cumulativeLevels(parents.size(), 1);
	for(unsigned int i=1; i<cumulativeLevels.size(); i++){
// 		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels;
		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels[i-1];
	}
// 	int nl = cumulativeLevels.back() * nLevels;
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
// 	int nodeLevels = 3; // default value assuming SNP data -- need to adjust calculation here
	int nodeLevels = node->getNumLevels(currentSet);
	 
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

///
/// Generates and return balanced accuracy output
///@return additional output 
///
vector<std::string> BayesSolutionCreator::getAdditionalFinalOutput(){
// 	vector<std::string> returnValues;
// 	if(notImproved){
// 		returnValues.push_back("*");
// 	}
// 	else{
// 		returnValues.push_back(" ");
// 	}
// 	
// 	// placeholder for balanced accuracy results
// 	returnValues.push_back(" ");		
	return addOutputValues;
}

// float BayesSolutionCreator::getBalAccuracy(Dataset* calcSet, vector<IndividualTerm*>& modelTerms){
// 
// 	string key;
// 	for(vector<IndividualTerm*>::iterator iter=modelTerms.begin(); iter != modelTerms.end();
// 		++iter){
// 			key += (*iter)->getLabel() + ":";
// 	}
// 	
// 	OptScore newScore(1.01, key);
// 	float balAcc;
// 	std::set<OptScore>::iterator scoreIter = scoreSet.find(newScore);
// 	if(scoreIter == scoreSet.end()){
// 	  newScore.sc = MDR::calcBalAccuracy(calcSet, calcSet, modelTerms);
// 	  if(scoreSet.size() < maximumSetSize){
// 	  	scoreSet.insert(newScore);
// 		}
// 	}
// 	else{
// 		newScore = (*scoreIter);
// 	}
// 	return newScore.sc;
// }

///
/// alternative implementation
///
float BayesSolutionCreator::getBalAccuracy(Dataset* calcSet, vector<IndividualTerm*> modelTerms){


// for(size_t i=0; i<modelTerms.size(); i++){
// cout << modelTerms[i]->getLabel() << " ";
// }

sort(modelTerms.begin(), modelTerms.end());	
	
	ScoreHolder::ScoreNode* scNode = savedScores.getScore(modelTerms);
// cout << "score=" << scNode->sc; //<< endl;
	if(scNode->sc == ScoreHolder::notFound()){
// cout << " NOT FOUND ";
		scNode->sc = MDR::calcBalAccuracy(calcSet, calcSet, modelTerms);
	}
// else{
// 	cout << " FOUND ";
// }
// 
// cout << scNode->sc << endl;
	return scNode->sc;
}



