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
#include "AthenaGrammarSI.h"

#include <iostream>
#include <sstream>
#include <map>
#include "AthenaExcept.h"

using namespace std;


///
/// Sets the rule to be used in recoding values during intialization
/// @param rule_string
///
void AthenaGrammarSI::setRestrictRule(std::string ruleString){
	Symbol searchSymbol(ruleString, NTSymbol);
	restrictRulePtr = findRule(searchSymbol);
}



///
/// Adds vector representing model by indexes in original genotype dataset array
/// @param indexes vector contain indexes with positions in original geno array of the dataset
/// 
void AthenaGrammarSI::addModel(std::vector<int> indexes){
	GrammarModel mod;
	mod.datasetIndexes = indexes;
	gramModels.push_back(mod);
}



///
/// After grammar is read, the models can be adjusted to have the codon values set 
/// for them rather than the original dataset indexes.  These codons can then be used
/// in initialization to replace the variables that were randomly selected.
/// @param dummyEncoded when true the codons can be either of two different choices
///
void AthenaGrammarSI::setModelCodons(bool dummyEncoded){
	
	// construct map with value being the codon value for the rule
	map<string, int> genoLocMap;
	for(unsigned int i=0; i<restrictRulePtr->size(); i++){
		genoLocMap[*((*restrictRulePtr)[i][0])] = i;
	}
	
	map<string, int>::iterator mapiter;
	
	for(vector<GrammarModel>::iterator iter=gramModels.begin(); iter != gramModels.end(); iter++){   
		for(unsigned int i=0; i < iter->datasetIndexes.size(); i++){
			// construct string to look for symbol (codon value) that matches
			if(dummyEncoded)
				iter->datasetIndexes[i] = iter->datasetIndexes[i] * 2;
			
			iter->datasetIndexes[i]++; // offset to start with V1 then V2 etc
			
			// 'coin flip' to see if take first or second entry when dummy encoded
			if(dummyEncoded && rand()/(RAND_MAX+1.0) > 0.5){
				iter->datasetIndexes[i]++;
			}
			
			stringstream ss;
			ss << "G" << iter->datasetIndexes[i];

			// now find codon for it by searching for match in productions of the restricted rule
			if((mapiter=genoLocMap.find(ss.str())) == genoLocMap.end()){
				throw AthenaExcept(ss.str() + " is not a valid genotype in grammar for initialization");
			}	
			
			iter->codonValues.push_back(mapiter->second);
		}
		
	}

}



///////////////////////////////////////////////////////////////////////////////
// Grow the derivation tree according to the grow or full method, up to the
// maximumDepth specified.
bool AthenaGrammarSI::growTree(DerivationTree &tree,const bool &growMethod,const unsigned int &maximumDepth){
	#if (DEBUG_LEVEL >= 2)
	cerr << "'bool GEGrammarSI::growTree(DerivationTree&, const bool&, const int&)' called\n";
	#endif

	// Stop conditions
	if(tree.getCurrentLevel()>maximumDepth){
		return false;
	}
	if(tree.getData()->getType()==TSymbol){
		return true;
	}

	Rule *rulePtr;
	
	if(!(rulePtr=findRule(*(tree.getData())))){// No definition for the current non-terminal found
		if((*(tree.getData())).substr(0,strlen("<GECodonValue"))=="<GECodonValue"){
			genotype.push_back(static_cast<CodonType>(genotype.getMaxCodonValue()*(static_cast<float>(rand())/RAND_MAX)));
			return true;
		}
		else if (*(tree.getData())=="<GEXOMarker>"){
			return true;
		}
		else{
			genotype.setValid(false);
			phenotype.setValid(false);
			return true;
		}
	}
	else{
		Rule::const_iterator prodIt=rulePtr->begin();
		int ii=0;
		vector<int> possibleRules;// The possible rules to choose from
		possibleRules.clear();
		bool recursiveRules=false;// Flags the presence of recursive rules on full method
		while(prodIt!=rulePtr->end()){
			// Choose from all rules growing the individual
			// up to maximumDepth
			if((tree.getCurrentLevel()+prodIt->getMinimumDepth())<=maximumDepth){
				if((!growMethod)&&(!recursiveRules)&&(prodIt->getRecursive())){
					// Choose only recursive rules from now on
					recursiveRules=true;
					possibleRules.clear();
				}
				if((growMethod)||(!recursiveRules)||((!growMethod)&&(recursiveRules)&&(prodIt->getRecursive()))){
					possibleRules.push_back(ii);
				}
			}
			prodIt++;
			ii++;
		}
		// possibleRules now contains all valid rules
		if(possibleRules.empty()){
			return false;
		}
		else if(rulePtr->size()>1){
		
		  // If this is a variable from the grammar (represented by <v>)
		  // and there is a model to be used from the biofilter 
		  // and all codons haven't been taken from the model
		  if( rulePtr == restrictRulePtr && currGramModel != gramModels.end() &&
		    currGramModel->codonValues.size() > currModelCodonIndex){
		    genotype.push_back(currGramModel->codonValues[currModelCodonIndex]);
		    currModelCodonIndex++;
			}
		  else{
			// Only choose production and insert it on genotype if there
			// is more than 1 production associated with current rule
			  genotype.push_back(possibleRules[static_cast<CodonType>(possibleRules.size()*(rand()/(RAND_MAX+1.0)))]);
			}
			
			// Save choice
			prodIt=rulePtr->begin()+genotype.back();
			// Perform "unmod" on choice
			genotype.back()+=static_cast<CodonType>((genotype.getMaxCodonValue()/rulePtr->size()*(rand()/(RAND_MAX+1.0))))*rulePtr->size();
		}
		else{
			// Otherwise set prodIt to point to the only production
			prodIt=rulePtr->begin();
		}
		// Insert symbols of chosen production on argument derivation tree,
		// and call grow tree for each symbol
		Production::const_iterator symbIt=prodIt->begin();
		bool result=true;
		unsigned int newMaxDepth=tree.getDepth();
		while((symbIt!=prodIt->end())&&result){
			tree.push_back(DerivationTree(*symbIt,tree.getCurrentLevel()+1,tree.getDepth()+1));
			result=growTree(tree.back(),growMethod,maximumDepth);
			// Update maximum depth of tree
			if(newMaxDepth<tree.back().getDepth()){
				newMaxDepth=tree.back().getDepth();
			}
			symbIt++;
		}
		genotype.setValid(result);
		phenotype.setValid(result);
		tree.setDepth(newMaxDepth);

		return result;
	}
}



///////////////////////////////////////////////////////////////////////////////
// Initialise the Genotype and Phenotype structures, according to the sensible
// initialisation technique for GE. If index is not set (or if it is set
// to UINT_MAX), initialise the structures as part of a series of calls to
// init(); if it is set, initialise the structures as being the index-th member
// (out of popSize) of the population.
// The next call to this routine will initialise the structures as being the
// next individual of a population of popSize individuals. There is an
// exception to this rule:
// * If a specific index is set, then the structures are initialised as being
// the index-th individual of a population of popSize, and the next call
// to this routine will initialise the index-th+1 individual (unless a
// specific index is set again).
bool AthenaGrammarSI::init(const unsigned int index){
	#if (DEBUG_LEVEL >= 2)
	cerr << "'bool GEGrammarSI::init(const unsigned int)' called\n";
	#endif
	if(index!=UINT_MAX){
		setIndex(index);
	}
	unsigned int maxDepth = getMaxDepth();

	// Check depth validity
	if(maxDepth<1){
		cerr << "Cannot initialize individual with maxDepth set to zero.\n";
		return false;
	}
	// Check for valid mapper
	if(!getValidGrammar()){
		cerr << "Invalid Mapper, cannot initialize individual.\n";
		return false;
	}
	// check if start symbol minimumDepth smaller or equal to newMaxDepth
	const Rule *startRule=getStartRule();
	if(startRule->getMinimumDepth()>=getMaxDepth()){// maxDepth is smaller
		cerr << "Current maxDepth (" << getMaxDepth() <<  ") is too small to initialize individual.\n";
		return false;
	}

	// Grow or Full?
	bool grow=(getIndex()<static_cast<unsigned int>(getPopSize()*getGrow()));
	// Clear genotype
	genotype.clear();
	// Clear derivation tree, and add start symbol
	derivationTree.clear();
	derivationTree.setData(getStartSymbol());
	derivationTree.setCurrentLevel(1);
	derivationTree.setDepth(1);

	currModelCodonIndex = 0;

	// Grow individual until reaching newDepth
	bool returnValue;
	returnValue=growTree(derivationTree,grow,getMaxDepth());
	if(returnValue){
		genotype.setValid(true);
		if (!genotype2Phenotype()){
			cerr << "WARNING: invalid phenotype structure produced with Sensible Initialisation\n";
		}
	}
	
	// Create tails if required
	unsigned int tailsize=0;
	if(getTailRatio()>0.0){
		tailsize=static_cast<unsigned int>(genotype.size()*getTailRatio());
	}
	else{
		tailsize=getTailSize();
	}
	if(tailsize){
		// Create tail of size tailsize
		for(unsigned int ii=0;ii<tailsize;ii++){
			genotype.push_back(static_cast<CodonType>(genotype.getMaxCodonValue()*(static_cast<float>(rand())/RAND_MAX)));
		}
	}
	setIndex(getIndex()+1);
	
	// increment currGramModel iterator if not already at end
	if(currGramModel != gramModels.end())
		++currGramModel;	
			
	return returnValue;
}



///
/// Establishes variable_codon_map for use in transferring genomes from 
/// one grammar to another.  The restrictRulePtr needs to be set first.
///
void AthenaGrammarSI::setVariableCodonMap(){
	for(unsigned int prod=0; prod < restrictRulePtr->size(); prod++){
		for(unsigned int symb=0; symb < (*restrictRulePtr)[prod].size(); symb++){
			variableCodonMap[(*(*restrictRulePtr)[prod][symb])] = prod;
		}
	}
}


///
/// Returns codon value associated with variable.  Codon is randomly chosen 
/// @param variable string matching terminal to find
/// @param maxCodonValue maximum possible codon value (get from genotype.getMaxCodonValue)
/// @return codon value
///
int AthenaGrammarSI::getCodonVarValue(std::string variable, int maxCodonValue){
	return variableCodonMap[variable] + 
		(static_cast<CodonType>((maxCodonValue/restrictRulePtr->size()*(rand()/(RAND_MAX+1.0))))*restrictRulePtr->size());
}


///
/// Takes a genome translates it using current rules and then alters the variables for new mapper rules
///
void AthenaGrammarSI::convertGenomeVariables(AthenaGrammarSI& newMapper, 
	const GA1DArrayGenome<int> &genome){

	setGenotype(genome);
	genotype2PhenotypeConvert(newMapper);
}


///////////////////////////////////////////////////////////////////////////////
// Updates the contents of the phenotype structure, based on the current
// genotype and the current grammar, and according to the standard GE
// mapping process. Returns true upon a successful mapping, and false
// otherwise, and also updates the valid field of the phenotype.
// With argument set to true, also updates derivationTree.
bool AthenaGrammarSI::genotype2PhenotypeConvert(AthenaGrammarSI& newMapper, const bool buildDerivationTree){

	#if (DEBUG_LEVEL >= 2)
	cerr << "'bool GEGrammar::genotype2Phenotype(const bool)' called\n";
	#endif
	bool returnValue=true;
	unsigned int newEffectiveSize=0;
	// Start by setting effectiveSize to 0
	genotype.setEffectiveSize(newEffectiveSize);

	phenotype.clear();
	if(buildDerivationTree){
		productions.clear();
	}
	// Quick safety checks
	//if((!getValidGrammar())||(!genotype.getValid())||(!getGenotype()->size())){
	if(!getValidGrammar()){
		phenotype.clear();
		phenotype.setValid(false);
		return false;
	}
	// Wraps counter and nonTerminals stack
	unsigned int wraps=0;
	stack<const Symbol*> nonTerminals;

	// Iterators
	iterator ruleIt;
	Rule::iterator prodIt;
	Genotype::iterator genoIt=genotype.begin();

	// Start with the start symbol
	nonTerminals.push(getStartSymbol());
	if(buildDerivationTree){
		// Use start symbol as the derivationTree node
		derivationTree.setData(getStartSymbol());
	}

	bool gotToUseWrap=false;
	// Get rid of all non-terminal symbols
	while((!nonTerminals.empty())&&(wraps<=getMaxWraps())){
		// Do a mapping step
		switch(genotype2PhenotypeStepConvert(newMapper,nonTerminals,genoIt,buildDerivationTree)){
			case -1:returnValue=false;
				break;
			case 0:	;
				break;
			case 1:	genoIt++;
				newEffectiveSize++;
				if(gotToUseWrap){
					wraps++;
					gotToUseWrap=false;
				}
		
				// Check if wrap is needed
				if(genoIt==genotype.end()){
					genoIt=genotype.begin();
					gotToUseWrap=true;
				}
				break;
			default:cerr << "Internal error in genotype2Phenotype().\n";
				cerr << "Execution aborted.\n";
				exit(0);
		}
	}

	// Was the mapping successful?
	if((wraps>getMaxWraps())||(!nonTerminals.empty())){
		returnValue=false;
		// Add remaining symbols in nonTerminals queue to phenotype
		while(!nonTerminals.empty()){
			phenotype.push_back(nonTerminals.top());
			nonTerminals.pop();
		}
	}

	phenotype.setValid(returnValue);
	genotype.setEffectiveSize(newEffectiveSize);
	genotype.setWraps(wraps);
	// Now build derivation tree, based on productions vector
	if(buildDerivationTree){
		derivationTree.clear();
		derivationTree.setData(getStartSymbol());
		derivationTree.setCurrentLevel(1);
		derivationTree.setDepth(1);
		vector<Production*>::iterator prodIt=productions.begin();
		buildDTree(derivationTree,prodIt);
	}

	return returnValue;
}



///
/// Performs one step of the mapping process, that is, maps the next
/// non-terminal symbol on the nonTerminals stack passed as argument, using the
/// codon at the position pointed by genoIt.
/// @param nonTerminals
/// @param genoIt
/// @param buildDerivationTree
/// @return number of codons consumed, -1 if not successful
///
int AthenaGrammarSI::genotype2PhenotypeStepConvert(AthenaGrammarSI& newMapper, 
	stack<const Symbol*> &nonTerminals, Genotype::iterator genoIt, bool buildDerivationTree){
	#if (DEBUG_LEVEL >= 2)
	cerr << "'bool GEGrammar::genotype2PhenotypeStep(stack<const Symbol*>&,Genotype::iterator&,bool)' called\n";
	#endif
	Rule::iterator prodIt;
	int returnValue;

	// Find the rule for the current non-terminal
	Rule *rulePtr=findRule(*(nonTerminals.top()));

	if(!rulePtr){// Undefined symbol - could be an extension symbol
		if(((*(nonTerminals.top())).substr(0,strlen("<GECodonValue"))==
			"<GECodonValue")&&(genoIt!=genotype.end())){
			// Insert codon value
			// Extract range for value from non-terminal specification
			int low=0,high=-1,pointer=strlen("<GECodonValue");
			// currentChar is the first character after "<GECodonValue"
			char currentChar=((*(nonTerminals.top())).substr(pointer,1))[0];
			// Look for range definitions
			while(currentChar!='>'){
				if(currentChar=='-'){
					// Low range specification
					currentChar=((*(nonTerminals.top())).substr(++pointer,1))[0];
					while((currentChar>='0')&&(currentChar<='9')){
						low=(low*10)+(currentChar-'0');
						currentChar=((*(nonTerminals.top())).substr(++pointer,1))[0];
					}
				}
				else if(currentChar=='+'){
					// High range specification
					currentChar=((*(nonTerminals.top())).substr(++pointer,1))[0];
					while((currentChar>='0')&&(currentChar<='9')){
						if(high==-1){
							high=0;
						}
						high=(high*10)+(currentChar-'0');
						currentChar=((*(nonTerminals.top())).substr(++pointer,1))[0];
					}
				}
				else{// Ignore errors
					currentChar=((*(nonTerminals.top())).substr(++pointer,1))[0];
				}
			}
			// High range was not specified, so set it to maximum
			if(high==-1){
				high=genotype.getMaxCodonValue();
			}
			// Remove non-terminal
			nonTerminals.pop();
			// Print value onto "codon"
			ostringstream codon;
			if(high==low){
				// Catch division by zero
				codon << low;
			}
			else{
				codon << (*genoIt%(high-low+1))+low;;
			}
			// Insert symbol with value onto phenotype
			phenotype.push_back(new Symbol(codon.str()));
			returnValue=1;
		}
		else{
			// Unknown symbol or special symbol that requires non-empty genotype
			// Include symbol on phenotype
			phenotype.push_back(nonTerminals.top());
			// Remove non-terminal
			nonTerminals.pop();
			// Invalidate mapping
			returnValue=-1;
		}
	} // end of !rulePtr

	// Allow recursive rules, but only if they consume a codon
	else if((rulePtr->getMinimumDepth()>=INT_MAX>>1)// Stuck on recursive rule
		&&(rulePtr->size()<=1)){// No codon will be consumed
		// Include symbol on phenotype
		phenotype.push_back(nonTerminals.top());
		// Remove non-terminal
		nonTerminals.pop();
		// Invalidate mapping
		returnValue=-1;
	}
	else{ // Usual progression of the mapping
		// Remove non-terminal
		nonTerminals.pop();
		// Choose production
		if((genoIt==genotype.end())&&(rulePtr->size()>1)){
			// Empty genotype, but symbol requires choice
			// Include symbol on phenotype
			phenotype.push_back(*(rulePtr->lhs.begin()));
			// Invalidate mapping
			returnValue=-1;
		}
		else{
			if(genoIt==genotype.end()){//Empty genotype
				prodIt=rulePtr->begin();
			}
			else{
				prodIt=rulePtr->begin()+*genoIt%(rulePtr->size());
			}
			// Place production on productions vector
			if(buildDerivationTree){
				productions.push_back(&*prodIt);
			}
			// Put all terminal symbols at start of production onto phenotype
			int s_start=0;
			int s_stop=(*prodIt).size();
			while((s_start<s_stop)&&((*prodIt)[s_start]->getType()==TSymbol)){
				phenotype.push_back((*prodIt)[s_start]);
				s_start++;
			}
			// Push all remaining symbols from production onto nonTerminals queue, backwards
			for(;s_stop>s_start;s_stop--){
				nonTerminals.push((*prodIt)[s_stop-1]);
			}
			// 0 or 1 choice for current rule, didn't consume genotype codon
			if(rulePtr->size()<=1){
				returnValue=0;
			}
			else{
				returnValue=1;
			}
		}
	}
	// Finally, pop all terminal symbols on top of stack and insert onto phenotype
	while((!nonTerminals.empty()) && (nonTerminals.top()->getType()==TSymbol)){
		phenotype.push_back(nonTerminals.top());
		nonTerminals.pop();
	}


	// when this rulePtr = restrictRulePtr need to adjust the value of the current genotype to reflect
	// the second mapper location for it
	if(rulePtr == restrictRulePtr){
		Phenotype::iterator phenoIter = phenotype.end();
		phenoIter--;
		phenoIter--;
	  *genoIt = newMapper.getCodonVarValue(*(*phenoIter), genotype.getMaxCodonValue());
	}
	
	return returnValue;
}


///
/// Determines size of genome block from starting codon passed
///
int AthenaGrammarSI::determineBlockLength(int startCodon){

	int blockSize = 0;
	
	stack<const Symbol*> nonTerminals;
	
	Symbol * startSymbol = new Symbol(codonVector[startCodon], NTSymbol);
	nonTerminals.push(startSymbol);
	
	unsigned int wraps=0;
	bool buildDerivationTree = false;
	
	// set genotype iterator to point 
	Genotype::iterator genoIt=genotype.begin();
	for(int i=0; i<startCodon; i++){
		genoIt++;
	}
	
	bool gotToUseWrap=false;
	// Get rid of all non-terminal symbols
	while((!nonTerminals.empty())&&(wraps<=getMaxWraps())){
		// Do a mapping step
		switch(genotype2PhenotypeStep(nonTerminals,genoIt,buildDerivationTree)){
			case -1:blockSize=-1;
				break;
			case 0:	; // no codon consumed here so don't update the multimap or vector
				break;
			case 1:	
			  // consumed a codon so need record block size increase
				blockSize++;
			  genoIt++;
				if(gotToUseWrap){
					wraps++;
					gotToUseWrap=false;
				}
		
				// Check if wrap is needed
				if(genoIt==genotype.end()){
					genoIt=genotype.begin();
					gotToUseWrap=true;
				}
				break;
			default:cerr << "Internal error in genotype2Phenotype().\n";
				cerr << "Execution aborted.\n";
				exit(0);
		}
	}  
	
	delete startSymbol;  
	return blockSize;
}



///
/// Constructs vector and multimap for use in crossovers
/// @param genome
///
void AthenaGrammarSI::establishCodons(const GA1DArrayGenome<int> &genome){
	// clear holding structures
	codonMap.clear();
	codonVector.clear();
	setGenotype(genome);
	fillCodons();
}



///
/// fills structures that identify codons -- No wrapping allowed in this case
///
int AthenaGrammarSI::fillCodons(){

	phenotype.clear();

	bool returnValue=true;
	unsigned int newEffectiveSize=0;

	bool buildDerivationTree = false;

	// Wraps counter and nonTerminals stack
	stack<const Symbol*> nonTerminals;

	// Iterators
	iterator ruleIt;
	Rule::iterator prodIt;
	Genotype::iterator genoIt=genotype.begin();

	// Start with the start symbol
	nonTerminals.push(getStartSymbol());
	unsigned int wraps=0;

	bool gotToUseWrap=false;
	// Get rid of all non-terminal symbols
	while((!nonTerminals.empty())&&(wraps<=getMaxWraps())){
	  const Symbol* currentNonTerminal = nonTerminals.top();
	  
		// Do a mapping step
		switch(genotype2PhenotypeStep(nonTerminals,genoIt,buildDerivationTree)){
			case -1:returnValue=false;
				break;
			case 0:	; // no codon consumed here so don't update the multimap or vector
				break;
			case 1:	
			  // consumed a codon so need to store the nonterminal
			  codonVector.push_back(*currentNonTerminal);
			  codonMap.insert(pair<string, int>(*currentNonTerminal, newEffectiveSize));
			  genoIt++;
				newEffectiveSize++;
				if(gotToUseWrap){
					wraps++;
					gotToUseWrap=false;
				}
							
				// Check if wrap is needed
				if(genoIt==genotype.end()){
					genoIt=genotype.begin();
					gotToUseWrap=true;
				}
				break;
			default:cerr << "Internal error in genotype2Phenotype().\n";
				cerr << "Execution aborted.\n";
				exit(0);
		}
	}

	// Was the mapping successful?
	if((wraps>getMaxWraps())||(!nonTerminals.empty())){
		returnValue=false;
		// Add remaining symbols in nonTerminals queue to phenotype
		while(!nonTerminals.empty()){
			phenotype.push_back(nonTerminals.top());
			nonTerminals.pop();
		}
	}
	phenotype.setValid(returnValue);
	genotype.setEffectiveSize(newEffectiveSize);
	genotype.setWraps(wraps);
	return returnValue;
}



///    
/// Returns codon matching by randomly selecting one of the 
/// codons associated with the string indicating the rule choice (nonterminal)
/// @param rule string
/// @return index of codon that matches
///
int AthenaGrammarSI::getMatchingCodon(string rule){

	int count = codonMap.count(rule);
	int returnIndex = -1;

	if(count > 0){
		pair<multimap<string, int>::iterator, multimap<string, int>::iterator> range =
				codonMap.equal_range(rule);
		// randomly select one of them and return that index for start 
		// of matching codon
		multimap<string, int>::iterator chosenIt = range.first;
		int chosenIterator = (rand() % count);
		int i=0;
		while(i++ < chosenIterator){
			chosenIt++;
		}
		returnIndex = chosenIt->second;
	}
		
	return returnIndex;
}



///
/// Sets genotype for use in optimization.  Marks start and end of 
/// every set of codons that specify a constant (weight) in the network.
/// @param genome GA1DArrayGenome to convert for optimization
/// @return vector with struct marking start and end of codons for
/// later conversion
///
vector<AthenaGrammarSI::codonBlocks> AthenaGrammarSI::setGenotypeOpt(const GA1DArrayGenome<int> &genome){

	// similar to original but tracks the location of blocks corresponding to the
	// optimized values (in GENN, those are the weights -- optimized by backpropagation)
	genotype.clear();
	
	if(genome.size()){
		// Copy elements of GAGenome onto genotype
		int ii=0;
		while(ii<genome.size()){
			genotype.push_back(genome.gene(ii++));
		}
		genotype.setValid(true);
	}

	vector<AthenaGrammarSI::codonBlocks> blocks;
	// call genotype2PhenotypeOpt -- sets block vector
	genotype2PhenotypeOpt(blocks);

	return blocks;  
}



///
/// Sets genotype for use in optimization.  Marks start and end of 
/// every set of codons that specify a constant (weight) in the network.
/// @param genome GA1DArrayGenome to convert for optimization
/// @param compressed_symbols Symbols that will be part of the 
/// compression.  So the end of the compressed codons will be the last one
/// before a codon that is not part of this set.
/// @param startSymbol symbol marking start of compression 
/// @return vector with struct marking start and end of codons for
/// later conversion
///
vector<AthenaGrammarSI::codonBlocks> AthenaGrammarSI::setGenotypeOpt(const GA1DArrayGenome<int> &genome,
	std::set<string> compressedSet, string startSymbol, bool singleOpt){

	setOptStartSymbol(startSymbol);
	setOptSymbolSet(compressedSet);
	optIsSingle = singleOpt;
	return setGenotypeOpt(genome);
}



///////////////////////////////////////////////////////////////////////////////
// Updates the contents of the phenotype structure, based on the current
// genotype and the current grammar, and according to the standard GE
// mapping process. Returns true upon a successful mapping, and false
// otherwise, and also updates the valid field of the phenotype.
// With argument set to true, also updates derivationTree.
bool AthenaGrammarSI::genotype2PhenotypeOpt(vector<AthenaGrammarSI::codonBlocks>& blocks){
	#if (DEBUG_LEVEL >= 2)
	cerr << "'bool GEGrammar::genotype2Phenotype(const bool)' called\n";
	#endif
	bool returnValue=true;
	unsigned int newEffectiveSize=0;
	// Start by setting effectiveSize to 0
	genotype.setEffectiveSize(newEffectiveSize);

	phenotype.clear();
	// Quick safety checks
	if(!getValidGrammar()){
		phenotype.clear();
		phenotype.setValid(false);
		return false;
	}
	// Wraps counter and nonTerminals stack
	unsigned int wraps=0;
	stack<const Symbol*> nonTerminals;

	// Iterators
	iterator ruleIt;
	Rule::iterator prodIt;
	Genotype::iterator genoIt=genotype.begin();

	// Start with the start symbol
	nonTerminals.push(getStartSymbol());
	const Symbol* startOptSymbol = startOptRulePtr->lhs.front();

	int startOptSite=0;
	bool inOpt=false;
	codonBlocks newBlock;

	bool gotToUseWrap=false;
	// Get rid of all non-terminal symbols
	while((!nonTerminals.empty())&&(wraps<=getMaxWraps())){
// cout << "non terminals top=" << *(nonTerminals.top()) << endl;
		// check top of nonTerminals to see if it is the start symbol for optimization
		if(!inOpt){
			if(startOptSymbol == nonTerminals.top()){
				startOptSite = newEffectiveSize;
				inOpt=true;
				if(optIsSingle){
				  inOpt=false;
				  newBlock.start = startOptSite;
				  newBlock.end = newEffectiveSize+1;
				  blocks.push_back(newBlock);
				}
			}
		}
		else{ // currently running through optimized section of solution
			// check to see if top of nonTerminals is still in optimized set
			string ntString = *(nonTerminals.top());
// cout << "ntString=" << ntString << endl;
			if(optSymbols.find(ntString) == optSymbols.end()){
				// found end of the current optimized 
				inOpt=false;
				newBlock.start = startOptSite;
				newBlock.end = newEffectiveSize;
				blocks.push_back(newBlock);
			}
		}

		switch(genotype2PhenotypeStepOpt(nonTerminals,genoIt)){
			case -1:returnValue=false;
				break;
			case 0:	;
				break;
			case 1:	genoIt++;
				newEffectiveSize++;
				if(gotToUseWrap){
					wraps++;
					gotToUseWrap=false;
				}
				// Check if wrap is needed
				if(genoIt==genotype.end()){
					//newEffectiveSize+=genotype.size();
					genoIt=genotype.begin();
					gotToUseWrap=true;
				}
				break;
			default:cerr << "Internal error in genotype2Phenotype().\n";
				cerr << "Execution aborted.\n";
				exit(0);
		}
	}

	// Was the mapping successful?
	if((wraps>getMaxWraps())||(!nonTerminals.empty())){
		returnValue=false;
		// Add remaining symbols in nonTerminals queue to phenotype
		while(!nonTerminals.empty()){
			phenotype.push_back(nonTerminals.top());
			nonTerminals.pop();
		}
	}
	phenotype.setValid(returnValue);
	genotype.setEffectiveSize(newEffectiveSize);
	genotype.setWraps(wraps);
	return returnValue;
}




///////////////////////////////////////////////////////////////////////////////
// Performs one step of the mapping process, that is, maps the next
// non-terminal symbol on the nonTerminals stack passed as argument, using the
// codon at the position pointed by genoIt.
// Returns number of codons consumed, -1 if not successful
int AthenaGrammarSI::genotype2PhenotypeStepOpt(stack<const Symbol*> &nonTerminals,
		Genotype::iterator genoIt){
		
	#if (DEBUG_LEVEL >= 2)
	cerr << "'bool GEGrammar::genotype2PhenotypeStep(stack<const Symbol*>&,Genotype::iterator&,bool)' called\n";
	#endif
	Rule::iterator prodIt;
	int returnValue;

	// Find the rule for the current non-terminal
	Rule *rulePtr=findRule(*(nonTerminals.top()));

	if(!rulePtr){// Undefined symbol - could be an extension symbol
		if(((*(nonTerminals.top())).substr(0,strlen("<GECodonValue"))==
			"<GECodonValue")&&(genoIt!=genotype.end())){
			// Insert codon value
			// Extract range for value from non-terminal specification
			int low=0,high=-1,pointer=strlen("<GECodonValue");
			// currentChar is the first character after "<GECodonValue"
			char currentChar=((*(nonTerminals.top())).substr(pointer,1))[0];
			// Look for range definitions
			while(currentChar!='>'){
				if(currentChar=='-'){
					// Low range specification
					currentChar=((*(nonTerminals.top())).substr(++pointer,1))[0];
					while((currentChar>='0')&&(currentChar<='9')){
						low=(low*10)+(currentChar-'0');
						currentChar=((*(nonTerminals.top())).substr(++pointer,1))[0];
					}
				}
				else if(currentChar=='+'){
					// High range specification
					currentChar=((*(nonTerminals.top())).substr(++pointer,1))[0];
					while((currentChar>='0')&&(currentChar<='9')){
						if(high==-1){
							high=0;
						}
						high=(high*10)+(currentChar-'0');
						currentChar=((*(nonTerminals.top())).substr(++pointer,1))[0];
					}
				}
				else{// Ignore errors
					currentChar=((*(nonTerminals.top())).substr(++pointer,1))[0];
				}
			}
			// High range was not specified, so set it to maximum
			if(high==-1){
				high=genotype.getMaxCodonValue();
			}
			// Remove non-terminal
			nonTerminals.pop();
			// Print value onto "codon"
			ostringstream codon;
			if(high==low){
				// Catch division by zero
				codon << low;
			}
			else{
				codon << (*genoIt%(high-low+1))+low;;
			}
			// Insert symbol with value onto phenotype
			phenotype.push_back(new Symbol(codon.str()));
			returnValue=1;
		}
		else{
			// Unknown symbol or special symbol that requires non-empty genotype
			// Include symbol on phenotype
			phenotype.push_back(nonTerminals.top());
			// Remove non-terminal
			nonTerminals.pop();
			// Invalidate mapping
			returnValue=-1;
		}
	}
	// Allow recursive rules, but only if they consume a codon
	else if((rulePtr->getMinimumDepth()>=INT_MAX>>1)// Stuck on recursive rule
		&&(rulePtr->size()<=1)){// No codon will be consumed
		// Include symbol on phenotype
		phenotype.push_back(nonTerminals.top());
		// Remove non-terminal
		nonTerminals.pop();
		// Invalidate mapping
		returnValue=-1;
	} /* end of check for when not an actual rule */
	else{
		// Remove non-terminal
		nonTerminals.pop();
		// Choose production
		if((genoIt==genotype.end())&&(rulePtr->size()>1)){
			// Empty genotype, but symbol requires choice
			// Include symbol on phenotype
			phenotype.push_back(*(rulePtr->lhs.begin()));
			// Invalidate mapping
			returnValue=-1;
		}
		else{
			if(genoIt==genotype.end()){//Empty genotype
				prodIt=rulePtr->begin();
			}
			else{ /* here is main sequence that builds the solution */
				prodIt=rulePtr->begin()+*genoIt%(rulePtr->size());
			}
			// Place production on productions vector
			// Put all terminal symbols at start of production onto phenotype
			int s_start=0;
			int s_stop=(*prodIt).size();
			while((s_start<s_stop)&&((*prodIt)[s_start]->getType()==TSymbol)){
				phenotype.push_back((*prodIt)[s_start]);
				s_start++;
			}
			// Push all remaining symbols from production onto nonTerminals queue, backwards
			for(;s_stop>s_start;s_stop--){
				nonTerminals.push((*prodIt)[s_stop-1]);
			}
			// 0 or 1 choice for current rule, didn't consume genotype codon
			if(rulePtr->size()<=1){
				returnValue=0;
			}
			else{
				returnValue=1;
			}
		}
	}
	// Finally, pop all terminal symbols on top of stack and insert onto phenotype
	while((!nonTerminals.empty()) && (nonTerminals.top()->getType()==TSymbol)){
		phenotype.push_back(nonTerminals.top());
		nonTerminals.pop();
	}

	return returnValue;
	
}



///
/// Returns vector of ints corresponding to block 
/// @param opt_value Value that needs to be converted into grammar representation
/// @return vector of int with values of codons needed
///
vector<int> AthenaGrammarSI::translateOptValue(symbVector& optimizedSymb){
	 
	if(optimizedSymb.size() == 1){
	  // find rhs that matches the value needed then return a codon for it
	  int codon;
	  string lhstr;
	  vector<int> codons(1,0);
// cout << "get codon for " << optimizedSymb[0].symbol << endl;
	  getRuleFromProduction(optimizedSymb[0].symbol, lhstr, codons[0]);
	  return codons;
	} 
	 
	stack<GramElement> terms;
	stack<GramElement> workingStack;
	GramElement tempGram;
	// set up stack of terminal characters
	for(vector<optSymbol>::reverse_iterator symbIter = optimizedSymb.rbegin(); 
		symbIter != optimizedSymb.rend(); ++symbIter){
// cout << symbIter->symbol << " " << symbIter->noNT << endl;		 
			tempGram.symbol = symbIter->symbol;
			tempGram.noNT = symbIter->noNT;
			terms.push(tempGram);
	}

	// have left marker and right marker set already (left_opt_bound, right_opt_bound)
	GramElement currTerm;
	GramElement tempTerm;
 
	// continue loop until working stack is empty (will be empty on the first pass)
	do{
	
	while(terms.size()){
		currTerm = terms.top();
		terms.pop();
		// when no nonterminal for this element
		if(currTerm.noNT){
			// now should check to see if the last one was a right_opt_bound character
			// -- ')' -- that marks the end of a functional set
			if(currTerm.symbol[0] == rightOptBound){
				// marks the end of a set of elements -- will pop off stack until reach the 
				// left hand element that marks the end -- '(' -- usually
				GramElement workingGram = workingStack.top();
				string workingItem = workingGram.symbol;
				string newProduction;
				tempTerm.codons.clear();
				tempTerm.codons.insert(tempTerm.codons.begin(), workingGram.codons.begin(), workingGram.codons.end());
				workingStack.pop();
				
				do{
					newProduction = workingItem + newProduction;
					workingGram = workingStack.top();
					workingStack.pop();
					workingItem = workingGram.symbol;
					// carry over codons 
					tempTerm.codons.insert(tempTerm.codons.begin(), workingGram.codons.begin(), workingGram.codons.end());
				}while(workingItem[0] != leftOptBound || workingItem.size() > 1);			
				// now check for the rule made out of these parameters
				string leftSide;			
				// check to see if it needs another element from working stack
				if(isArg.find(newProduction) != isArg.end()){
					workingGram = workingStack.top();
					workingStack.pop();
					workingItem = workingGram.symbol;
					tempTerm.codons.insert(tempTerm.codons.begin(), workingGram.codons.begin(), workingGram.codons.end());
					newProduction = workingItem + newProduction;
				}				
				int codonValue;
				bool ruleFound = getRuleFromProduction(newProduction, leftSide, codonValue);
				if(codonValue >= 0){
					// insert at beginning
					tempTerm.codons.insert(tempTerm.codons.begin(), codonValue);
				}
								
				if(leftSide.size() > 0)
					tempTerm.symbol = leftOptBound + leftSide + rightOptBound;
				else
					tempTerm.symbol = leftOptBound + newProduction + rightOptBound;
					
				// push onto terms so it will be evaluated next time through the loop
				tempTerm.noNT = false;
				terms.push(tempTerm);
			}
			else{
				// as this is marked as an element that doesn't get translated up to a nonterminal
				// just push on to the working stack
				workingStack.push(currTerm);
			}
		}
		else{
			// may have a non-terminal
			string leftSide, current=currTerm.symbol;	
			// if the element is an argument to another element -- (<num>) for Concat -- take that element 
			// from the working stack and check the new combined element to see if it can be 
			// defined by a rule
			if(isArg.find(current) != isArg.end()){
				currTerm.codons.insert(currTerm.codons.begin(), workingStack.top().codons.begin(), 
					workingStack.top().codons.end());
				current = workingStack.top().symbol + current;
				workingStack.pop();
			}
			
			int codonValue;
			bool rule_found = getRuleFromProduction(current, leftSide, codonValue);
			// repeat this loop as long as the path can be followed up the rules
			while(rule_found){
				if(codonValue >=0)
					currTerm.codons.insert(currTerm.codons.begin(), codonValue);
				currTerm.symbol = leftSide;
				current = leftSide;
				rule_found = getRuleFromProduction(current, leftSide, codonValue);
			}
			
			// place last value on to workingStack for later conversion
			workingStack.push(currTerm);
		}
	} // terms stack is now empty

	// now add to the term stack if any elements in the workingStack (take top two if more than one)
	tempTerm = workingStack.top();
	workingStack.pop();
	if(!workingStack.empty()){
		tempTerm.symbol = workingStack.top().symbol + tempTerm.symbol;
		tempTerm.codons.insert(tempTerm.codons.begin(), workingStack.top().codons.begin(),
			workingStack.top().codons.end());
		workingStack.pop();
	}
	tempTerm.noNT = false;
	terms.push(tempTerm);
	// continue running until the only terminal is the starting one for the optimized portion
	// of the grammar
	}while(terms.top().symbol.compare(optStartSymbol) != 0);

	return terms.top().codons;
}


/// 
/// Returns rule (left hand side) that gives the production_string passed
/// @param production_string
/// @param left_hand
/// @param codon
/// @return codon value
///
bool AthenaGrammarSI::getRuleFromProduction(string& productionString, string& leftHand, int&
	codon){

	map<string, reverseRule>::iterator revIter = reverseRules.find(productionString);
	if(revIter != reverseRules.end()){

		leftHand = revIter->second.rule;
		
		// adjust codon by unmodding it
		if(revIter->second.codon >= 0)
			 codon = revIter->second.codon + static_cast<CodonType>((genotype.getMaxCodonValue()/revIter->second.numCodons*(rand()/(RAND_MAX+1.0))))*revIter->second.numCodons;
		else
			codon = revIter->second.codon;
		return true;
	}
	else{
		leftHand = "";
		codon = -2;
		return false;
	}

}


///
/// Construct a map.  The key for the map is the right side of a grammar.  The left side
/// of the grammar are the values along with the codon value that results in the right side
/// of the rule (the key of the map).  It assumes that all right side portions of the grammar
/// are unique (when considering the entire line in the grammar).  When there is only one
/// choice for a production, the codon is set to -1.
///
void AthenaGrammarSI::constructReverseGrammar(){

	reverseRules.clear();
	reverseRule reverse;
	
	// iterate through the grammar
	for(vector<Rule>::iterator ruleIter = begin(); ruleIter != end(); ++ruleIter){
		string ruleString = *(ruleIter->lhs.front());
		for(unsigned int prod=0; prod < ruleIter->size(); prod++){
			// construct the right-side string that will be used to map back which codon to use
			string rightSideString;
			for(vector<Symbol*>::iterator symbIter=ruleIter->at(prod).begin(); 
				symbIter!=ruleIter->at(prod).end(); ++symbIter){
					rightSideString += **symbIter;
			}
			reverse.rule = ruleString;
			if(prod == 0 && ruleIter->size() == 1)
				reverse.codon=-1;
			else
				reverse.codon = prod;
					 
			reverse.numCodons = int(ruleIter->size());
			reverseRules[rightSideString] = reverse;
		}
	}
}



int AthenaGrammarSI::buildDerivationTree(){
		bool returnValue=true;
		unsigned int newEffectiveSize=0;
	// Start by setting effectiveSize to 0
	genotype.setEffectiveSize(newEffectiveSize);

	phenotype.clear();

		// have to fill the production vector
		productions.clear();
		
		derivationTree.setData(getStartSymbol());
	// Wraps counter and nonTerminals stack
	unsigned int wraps=0;
	stack<const Symbol*> nonTerminals;

	// Iterators
	iterator ruleIt;
	Rule::iterator prodIt;
	Genotype::iterator genoIt=genotype.begin();
	
		// Start with the start symbol
	nonTerminals.push(getStartSymbol());
	bool gotToUseWrap=false;
	
	derivationTree.setData(getStartSymbol());
	
	// Get rid of all non-terminal symbols
	while((!nonTerminals.empty())&&(wraps<=getMaxWraps())){
		// Do a mapping step
		// call with buildDerivationTree=true
		switch(genotype2PhenotypeStep(nonTerminals,genoIt,true)){
			case -1:returnValue=false;
				break;
			case 0:	;
				break;
			case 1:	genoIt++;
				newEffectiveSize++;
				if(gotToUseWrap){
					wraps++;
					gotToUseWrap=false;
				}
				// Check if wrap is needed
				if(genoIt==genotype.end()){
					genoIt=genotype.begin();
					gotToUseWrap=true;
				}
				break;
			default:cerr << "Internal error in genotype2Phenotype().\n";
				cerr << "Execution aborted.\n";
				exit(0);
		}
	}	
	
		if((wraps>getMaxWraps())||(!nonTerminals.empty())){
		returnValue=false;
		// Add remaining symbols in noerminals queue to phenotype
		while(!nonTerminals.empty()){
			phenotype.push_back(nonTerminals.top());
			nonTerminals.pop();
		}
	}
	
		phenotype.setValid(returnValue);
	  genotype.setEffectiveSize(newEffectiveSize);
	  genotype.setWraps(wraps); 
		
		derivationTree.clear();
		derivationTree.setData(getStartSymbol());
		derivationTree.setCurrentLevel(1);
		derivationTree.setDepth(1);
		vector<Production*>::iterator prodIterator=productions.begin();
		buildDTree(derivationTree,prodIterator);
		
		return getMax(derivationTree);
}



///
/// Recursively search for the maximum level in the tree
///
int AthenaGrammarSI::getMax(DerivationTree& tree){
		DerivationTree::iterator treeIt = tree.begin();
		int maxDepth=0, depth=0;
		while(treeIt != tree.end()){
				if(treeIt->empty()){
						depth = treeIt->getCurrentLevel();
				}
				else{
						depth = getMax(*treeIt);
				}
				if(depth > maxDepth){
						maxDepth=depth;
				}
				treeIt++;
		}
		return maxDepth;
}

///
/// Return a new variable by selecting a random variable using the variable rule pointer
/// @param Stores the new codon value
///
std::string AthenaGrammarSI::getNewVariable(int& newCodon){
  // return the string associated with this codon
  newCodon = GARandomInt(0,RAND_MAX);
  Rule::iterator prodIt=varRulePtr->begin()+newCodon%(varRulePtr->size());
  return (*(*prodIt)[0]);
}



///////////////////////////////////////////////////////////////////////////////
// Updates the contents of the phenotype structure, based on the current
// genotype and the current grammar, and according to the standard GE
// mapping process. Returns true upon a successful mapping, and false
// otherwise, and also updates the valid field of the phenotype.
// With argument set to true, also updates derivationTree.
void AthenaGrammarSI::changeVariables(GA1DArrayGenome<int>& genome, map<int, SolutionCreator::TerminalInfo> variableMap){

	bool returnValue=true;
	unsigned int newEffectiveSize=0;

	phenotype.clear();

	// Wraps counter and nonterminals stack
	unsigned int wraps=0;
	stack<const Symbol*> nonterminals;

	// Iterators
	iterator ruleIt;
	Rule::iterator prodIt;
	int genoIt=0;
	
	// Start with the start symbol
	nonterminals.push(getStartSymbol());

	bool gotToUseWrap=false;
	// Get rid of all non-terminal symbols
	while((!nonterminals.empty())&&(wraps<=getMaxWraps())){
		// Do a mapping step
		switch(genotype2PhenotypeStepReplaceVar(nonterminals,genoIt,genome, variableMap)){
			case -1:returnValue=false;
				break;
			case 0:	;
				break;
			case 1:	genoIt++;
				newEffectiveSize++;
				if(gotToUseWrap){
					wraps++;
					gotToUseWrap=false;
				}
				// Check if wrap is needed
				if(genoIt >= genome.length()){
					//newEffectiveSize+=genotype.size();
					genoIt=0;
					gotToUseWrap=true;
				}
				break;
			default:cerr << "Internal error in genotype2Phenotype().\n";
				cerr << "Execution aborted.\n";
				exit(0);
		}
	}


}



///////////////////////////////////////////////////////////////////////////////
// Performs one step of the mapping process, that is, maps the next
// non-terminal symbol on the nonterminals stack passed as argument, using the
// codon at the position pointed by genoIt.
// Returns number of codons consumed, -1 if not successful
int AthenaGrammarSI::genotype2PhenotypeStepReplaceVar(stack<const Symbol*> &nonterminals,
		int& genoIt, GA1DArrayGenome<int>& genome, map<int, SolutionCreator::TerminalInfo>& variableMap){

	Rule::iterator prodIt;
	int returnValue;

	// Find the rule for the current non-terminal
	Rule *rulePtr=findRule(*(nonterminals.top()));

	if(!rulePtr){// Undefined symbol - could be an extension symbol
		if(((*(nonterminals.top())).substr(0,strlen("<GECodonValue"))==
			"<GECodonValue")&&(genoIt<genome.length())){
			// Insert codon value
			// Extract range for value from non-terminal specification
			int low=0,high=-1,pointer=strlen("<GECodonValue");
			// currentChar is the first character after "<GECodonValue"
			char currentChar=((*(nonterminals.top())).substr(pointer,1))[0];
			// Look for range definitions
			while(currentChar!='>'){
				if(currentChar=='-'){
					// Low range specification
					currentChar=((*(nonterminals.top())).substr(++pointer,1))[0];
					while((currentChar>='0')&&(currentChar<='9')){
						low=(low*10)+(currentChar-'0');
						currentChar=((*(nonterminals.top())).substr(++pointer,1))[0];
					}
				}
				else if(currentChar=='+'){
					// High range specification
					currentChar=((*(nonterminals.top())).substr(++pointer,1))[0];
					while((currentChar>='0')&&(currentChar<='9')){
						if(high==-1){
							high=0;
						}
						high=(high*10)+(currentChar-'0');
						currentChar=((*(nonterminals.top())).substr(++pointer,1))[0];
					}
				}
				else{// Ignore errors
					currentChar=((*(nonterminals.top())).substr(++pointer,1))[0];
				}
			}
			// High range was not specified, so set it to maximum
			if(high==-1){
				high=genotype.getMaxCodonValue();
			}
			// Remove non-terminal
			nonterminals.pop();
			// Print value onto "codon"
			ostringstream codon;
			if(high==low){
				// Catch division by zero
				codon << low;
			}
			else{
				codon << (genome.gene(genoIt)%(high-low+1))+low;;
			}
			// Insert symbol with value onto phenotype
			phenotype.push_back(new Symbol(codon.str()));
			int phenoIndex = int(phenotype.size())-1;
			if(variableMap.find(phenoIndex)!=variableMap.end()){
			  // change genotype
				genome.gene(genoIt, variableMap[phenoIndex].newValue);
			}
			returnValue=1;
		}
		else{
			// Unknown symbol or special symbol that requires non-empty genotype
			// Include symbol on phenotype
			phenotype.push_back(nonterminals.top());
			// Remove non-terminal
			nonterminals.pop();
			// Invalidate mapping
			returnValue=-1;
		}
	}
	//else if(rulePtr->getMinimumDepth()>=INT_MAX>>1){// Stuck on recursive rule
	// Allow recursive rules, but only if they consume a codon
	else if((rulePtr->getMinimumDepth()>=INT_MAX>>1)// Stuck on recursive rule
		&&(rulePtr->size()<=1)){// No codon will be consumed
		// Include symbol on phenotype
		phenotype.push_back(nonterminals.top());
		// Remove non-terminal
		nonterminals.pop();
		// Invalidate mapping
		returnValue=-1;
	}
	else{
		// Remove non-terminal
		nonterminals.pop();
		// Choose production
		if(genoIt == genome.length()&&(rulePtr->size()>1)){
			// Empty genotype, but symbol requires choice
			// Include symbol on phenotype
			phenotype.push_back(*(rulePtr->lhs.begin()));
			// Invalidate mapping
			returnValue=-1;
		}
		else{
			if(genoIt>genome.length()){//Empty genotype
				prodIt=rulePtr->begin();
			}
			else{
				prodIt=rulePtr->begin()+(genome.gene(genoIt))%(rulePtr->size());
			}
			// Place production on productions vector
			// Put all terminal symbols at start of production onto phenotype
			int s_start=0;
			int s_stop=(*prodIt).size();
			while((s_start<s_stop)&&((*prodIt)[s_start]->getType()==TSymbol)){
				phenotype.push_back((*prodIt)[s_start]);
				// check to see if need to convert this genotype
				int phenoIndex = int(phenotype.size())-1;
				if(variableMap.find(phenoIndex)!=variableMap.end()){
				  // change genotype
				  genome.gene(genoIt, variableMap[phenoIndex].newValue);
				}
				
				s_start++;
			}
			// Push all remaining symbols from production onto nonterminals queue, backwards
			for(;s_stop>s_start;s_stop--){
				nonterminals.push((*prodIt)[s_stop-1]);
			}
			// 0 or 1 choice for current rule, didn't consume genotype codon
			if(rulePtr->size()<=1){
				returnValue=0;
			}
			else{
				returnValue=1;
			}
		}
	}
	// Finally, pop all terminal symbols on top of stack and insert onto phenotype
	while((!nonterminals.empty()) && (nonterminals.top()->getType()==TSymbol)){
		phenotype.push_back(nonterminals.top());
		nonterminals.pop();
	}
	return returnValue;
}

