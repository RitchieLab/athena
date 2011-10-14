#include "HemannGrammarSI.h"

#include <iostream>
#include <sstream>
#include <map>
#include "HemannExcept.h"

using namespace std;


///
/// Sets the rule to be used in recoding values during intialization
/// @param rule_string
///
void HemannGrammarSI::setRestrictRule(std::string rule_string){
  Symbol searchSymbol(rule_string, NTSymbol);
  restrictRulePtr = findRule(searchSymbol);
// cout << "number of productions=" << restrictRulePtr->size() << endl;
// for(unsigned int i=0; i<restrictRulePtr->size(); i++){
//   cout << i << "=" << *((*restrictRulePtr)[i][0]) << endl;
// }
// exit(1);
}

///
/// Adds vector representing model by indexes in original genotype dataset array
/// @param indexes vector contain indexes with positions in original geno array of the dataset
/// 
void HemannGrammarSI::addModel(std::vector<int> indexes){
  GrammarModel mod;
  mod.dataset_indexes = indexes;
  gramModels.push_back(mod);
}


///
/// After grammar is read, the models can be adjusted to have the codon values set 
/// for them rather than the original dataset indexes.  These codons can then be used
/// in initialization to replace the variables that were randomly selected.
/// @param dummyEncoded when true the codons can be either of two different choices
///
void HemannGrammarSI::setModelCodons(bool dummyEncoded){
  
  // construct map with value being the codon value for the rule
  map<string, int> geno_loc_map;
  for(unsigned int i=0; i<restrictRulePtr->size(); i++){
    geno_loc_map[*((*restrictRulePtr)[i][0])] = i;
  }  

// cout << "map size = " << geno_loc_map.size() << endl;
// cout << "gramModels size=" << gramModels.size() << endl;

  map<string, int>::iterator mapiter;
  
  for(vector<GrammarModel>::iterator iter=gramModels.begin(); iter != gramModels.end(); iter++){
// cout << iter->dataset_indexes.size() << endl;   
    for(unsigned int i=0; i < iter->dataset_indexes.size(); i++){
 
//  cout << "start index = " << iter->dataset_indexes[i] << endl;
 
      // construct string to look for symbol (codon value) that matches
      if(dummyEncoded)
        iter->dataset_indexes[i] = iter->dataset_indexes[i] * 2;
      
      iter->dataset_indexes[i]++; // offset to start with V1 then V2 etc
      
      // 'coin flip' to see if take first or second entry when dummy encoded
      if(dummyEncoded && rand()/(RAND_MAX+1.0) > 0.5){
// cout << "USE SECOND" << endl;
        iter->dataset_indexes[i]++;
      }
      
      stringstream ss;
      ss << "G" << iter->dataset_indexes[i];

// cout << "looking for " << ss.str() << endl;
      
      // now find codon for it by searching for match in productions of the restricted rule
      if((mapiter=geno_loc_map.find(ss.str())) == geno_loc_map.end()){
        throw HemannExcept(ss.str() + " is not a valid genotype in grammar for initialization");
      }
// cout << "iter codon value is " << mapiter->second << endl;
      iter->codon_values.push_back(mapiter->second);
    }
    
  }

}


///////////////////////////////////////////////////////////////////////////////
// Grow the derivation tree according to the grow or full method, up to the
// maximumDepth specified.
bool HemannGrammarSI::growTree(DerivationTree &tree,const bool &growMethod,const unsigned int &maximumDepth){
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
	
// cout << "$" << *(tree.getData()) << "$" << gramModels.size() << endl;
	
	if(!(rulePtr=findRule(*(tree.getData())))){// No definition for the current non-terminal found
		if((*(tree.getData())).substr(0,strlen("<GECodonValue"))=="<GECodonValue"){
		//if (*(tree.getData())=="<GECodonValue>"){
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
		    currGramModel->codon_values.size() > currModelCodonIndex){
// cout << "adding codon to model " << currGramModel->codon_values[currModelCodonIndex] << endl;
		    genotype.push_back(currGramModel->codon_values[currModelCodonIndex]);
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

// cout << "reached end of growTree with result=" << result << endl;
		
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
bool HemannGrammarSI::init(const unsigned int index){
	#if (DEBUG_LEVEL >= 2)
	cerr << "'bool GEGrammarSI::init(const unsigned int)' called\n";
	#endif
	if(index!=UINT_MAX){
		setIndex(index);
	}
unsigned int maxDepth = getMaxDepth();

	// Check depth validity
	if(maxDepth<1){
		cerr << "Cannot initialise individual with maxDepth set to zero.\n";
		return false;
	}
	// Check for valid mapper
	if(!getValidGrammar()){
		cerr << "Invalid Mapper, cannot initialise individual.\n";
		return false;
	}
	// check if start symbol minimumDepth smaller or equal to newMaxDepth
	const Rule *startRule=getStartRule();
	if(startRule->getMinimumDepth()>=getMaxDepth()){// maxDepth is smaller
		cerr << "Current maxDepth (" << getMaxDepth() <<  ") is too small to initialise individual.\n";
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

//cout << "maxDepth=" << maxDepth << endl;
	// Grow individual until reaching newDepth
	bool returnValue;
	returnValue=growTree(derivationTree,grow,getMaxDepth());
	if(returnValue){
		genotype.setValid(true);
		if (!genotype2Phenotype()){
			cerr << "WARNING: invalid phenotype structure produced with Sensible Initialisation\n";
		}
	}
	
// output the new one here
// Phenotype const *phenotype=getPhenotype();
// cout << "pheno_size=" << phenotype->size() << " ";
//   for(unsigned int i=0; i<phenotype->size(); i++){
// ////       cout << "|" << *((*phenotype)[i]) << "|";
//       cout << *((*phenotype)[i]) << "|";
//   }
//   cout << endl << endl;	
/// debugging
	
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
void HemannGrammarSI::setVariableCodonMap(){
  for(unsigned int prod=0; prod < restrictRulePtr->size(); prod++){
// cout << "prod=" << prod << endl;
    for(unsigned int symb=0; symb < (*restrictRulePtr)[prod].size(); symb++){
      variable_codon_map[(*(*restrictRulePtr)[prod][symb])] = prod;
    }
  }
  
// check codon values
// std::map<std::string, int>::iterator varIter;
// for(varIter = variable_codon_map.begin(); varIter != variable_codon_map.end(); varIter++){
//   cout << varIter->first << " <=> " << varIter->second << endl;
// }
// cout << "++++++++" << endl;
  
}


///
/// Returns codon value associated with variable.  Codon is randomly chosen 
/// @param variable string matching terminal to find
/// @param maxCodonValue maximum possible codon value (get from genotype.getMaxCodonValue)
/// @return codon value
///
int HemannGrammarSI::getCodonVarValue(std::string variable, int maxCodonValue){
// static_cast<CodonType>((genotype.getMaxCodonValue()/restrictRulePtr->size()*(rand()/(RAND_MAX+1.0))))*rulePtr->size();
// static_cast<CodonType>(possibleRules.size()*(rand()/(RAND_MAX+1.0)))]
  // return is unmodded 
  
//   int value = variable_codon_map[variable] + 
//     (static_cast<CodonType>((maxCodonValue/restrictRulePtr->size()*(rand()/(RAND_MAX+1.0))))*restrictRulePtr->size());
//   cout << " converted to " << value << endl;
//   cout << "altering codon for " << variable << endl;
  
  return variable_codon_map[variable] + 
    (static_cast<CodonType>((maxCodonValue/restrictRulePtr->size()*(rand()/(RAND_MAX+1.0))))*restrictRulePtr->size());
    
}


///
/// Takes a genome translates it using current rules and then alters the variables for new mapper rules
///
void HemannGrammarSI::convertGenomeVariables(HemannGrammarSI& newMapper, 
  const GA1DArrayGenome<int> &genome){

  setGenotype(genome);


// cout << "Genotype in HemannGrammarSI before any conversion done" << endl;
// Genotype::const_iterator genIt=getGenotype()->begin();
// int i=0;
// while(genIt!=getGenotype()->end()){
// cout << int(genome.gene(i)) << "=" << *genIt << " ";
// //     genome.gene(i, *genIt);
// // cout << int(genome.gene(i)) << " ";
//  	  genIt++;
//     i++;
// }
// cout << endl;
  
  genotype2PhenotypeConvert(newMapper);

// cout << "Genotype in HemannGrammarSI after any conversion done" << endl;
// Genotype::const_iterator genIt2=getGenotype()->begin();
// // int i=0;
// i=0;
// while(genIt2!=getGenotype()->end()){
// cout << int(genome.gene(i)) << "=" << *genIt2 << " ";
// //     genome.gene(i, *genIt);
// // cout << int(genome.gene(i)) << " ";
//  	  genIt2++;
//     i++;
// }
// cout << endl << endl;

}


///////////////////////////////////////////////////////////////////////////////
// Updates the contents of the phenotype structure, based on the current
// genotype and the current grammar, and according to the standard GE
// mapping process. Returns true upon a successful mapping, and false
// otherwise, and also updates the valid field of the phenotype.
// With argument set to true, also updates derivationTree.
bool HemannGrammarSI::genotype2PhenotypeConvert(HemannGrammarSI& newMapper, const bool buildDerivationTree){

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

//cout << "in genotype2Phenotype " << endl;

	// Wraps counter and nonterminals stack
	unsigned int wraps=0;
	stack<const Symbol*> nonterminals;

	// Iterators
	iterator ruleIt;
	Rule::iterator prodIt;
	Genotype::iterator genoIt=genotype.begin();

	// Start with the start symbol
	nonterminals.push(getStartSymbol());
	if(buildDerivationTree){
		// Use start symbol as the derivationTree node
		derivationTree.setData(getStartSymbol());
	}

	bool gotToUseWrap=false;
	// Get rid of all non-terminal symbols
	while((!nonterminals.empty())&&(wraps<=getMaxWraps())){
		// Do a mapping step
		switch(genotype2PhenotypeStepConvert(newMapper,nonterminals,genoIt,buildDerivationTree)){
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
				
// check nonterminals size -- may indicate when finished ?
// cout << "effsize=" << newEffectiveSize << " nonterm.size = " << nonterminals.size() << endl;
				
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
	//newEffectiveSize+=(genoIt-genotype.begin());
	// Was the mapping successful?
	if((wraps>getMaxWraps())||(!nonterminals.empty())){
		returnValue=false;
		// Add remaining symbols in nonterminals queue to phenotype
		while(!nonterminals.empty()){
			phenotype.push_back(nonterminals.top());
			nonterminals.pop();
		}
	}
// cout << "VALID = " << returnValue << endl;
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
// cout << "COMPLETED GENOTYPE TRANSLATION" << endl;	
	
	return returnValue;
}


///
/// Performs one step of the mapping process, that is, maps the next
/// non-terminal symbol on the nonterminals stack passed as argument, using the
/// codon at the position pointed by genoIt.
/// @param nonterminals
/// @param genoIt
/// @param buildDerivationTree
/// @return number of codons consumed, -1 if not successful
///
int HemannGrammarSI::genotype2PhenotypeStepConvert(HemannGrammarSI& newMapper, 
  stack<const Symbol*> &nonterminals, Genotype::iterator genoIt, bool buildDerivationTree){
	#if (DEBUG_LEVEL >= 2)
	cerr << "'bool GEGrammar::genotype2PhenotypeStep(stack<const Symbol*>&,Genotype::iterator&,bool)' called\n";
	#endif
	Rule::iterator prodIt;
	int returnValue;

	// Find the rule for the current non-terminal
	Rule *rulePtr=findRule(*(nonterminals.top()));

	//cerr << "mapping " << *(nonterminals.top()) << " with " << *genoIt << "\n";
/*
bool output = false;
if((*(nonterminals.top())).compare("<v>")==0){
cout << "mapping " << *(nonterminals.top()) << endl;
cout << *genoIt << endl;
output = true;
}
*/

	if(!rulePtr){// Undefined symbol - could be an extension symbol
		if(((*(nonterminals.top())).substr(0,strlen("<GECodonValue"))==
			"<GECodonValue")&&(genoIt!=genotype.end())){
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
				codon << (*genoIt%(high-low+1))+low;;
			}
			// Insert symbol with value onto phenotype
			phenotype.push_back(new Symbol(codon.str()));
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
	} // end of !rulePtr
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
	else{ // Usual progression of the mapping
		// Remove non-terminal
// cout << "consumed nonterminal " << *(nonterminals.top()) << endl;
		nonterminals.pop();
//if(output) cout << "nonterminals popped " << endl;
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
//if(output) cout << "rulePtr->size()=" << rulePtr->size() <<  " and *genoIt= " << *genoIt
//  << " so rule number is " << *genoIt%(rulePtr->size()) << endl;
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


	// when this rulePtr = restrictRulePtr need to adjust the value of the current genotype to reflect
	// the second mapper location for it
	if(rulePtr == restrictRulePtr){
	  // update genoIt
// 
// for(Phenotype::iterator phenoIter = phenotype.begin(); phenoIter != phenotype.end(); phenoIter++){
//   cout << *(*phenoIter) << " ";
// }
// cout << endl;

    Phenotype::iterator phenoIter = phenotype.end();
    phenoIter--;
    phenoIter--;
//     cout << "updated codon " << *genoIt << " for variable " << **phenoIter << endl;
	  *genoIt = newMapper.getCodonVarValue(*(*phenoIter), genotype.getMaxCodonValue());
	}
	
// cout << "nonterminals now have size=" << nonterminals.size() << endl;
	return returnValue;
}


///
/// Determines size of genome block from starting codon passed
///
int HemannGrammarSI::determineBlockLength(int startCodon){

  int blocksize = 0;
  
  stack<const Symbol*> nonterminals;
  
  Symbol * startSymbol = new Symbol(codon_vector[startCodon], NTSymbol);
  nonterminals.push(startSymbol);
  
  unsigned int wraps=0;
  bool buildDerivationTree = false;
  
  // set genotype iterator to point 
  Genotype::iterator genoIt=genotype.begin();
  for(int i=0; i<startCodon; i++){
    genoIt++;
  }
  
	bool gotToUseWrap=false;
	// Get rid of all non-terminal symbols
	while((!nonterminals.empty())&&(wraps<=getMaxWraps())){
		// Do a mapping step
		switch(genotype2PhenotypeStep(nonterminals,genoIt,buildDerivationTree)){
			case -1:blocksize=-1;
				break;
			case 0:	; // no codon consumed here so don't update the multimap or vector
				break;
			case 1:	
			  // consumed a codon so need record block size increase
        blocksize++;
			  genoIt++;
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
  
  delete startSymbol;  
  return blocksize;
}



///
/// Constructs vector and multimap for use in crossovers
/// @param genome
///
void HemannGrammarSI::establishCodons(const GA1DArrayGenome<int> &genome){
  
  // clear holding structures
  codon_map.clear();
  codon_vector.clear();
  
  setGenotype(genome);
  
  fillCodons();
  
}

///
/// fills structures that identify codons -- No wrapping allowed in this case
///
int HemannGrammarSI::fillCodons(){

	phenotype.clear();

	bool returnValue=true;
	unsigned int newEffectiveSize=0;

  bool buildDerivationTree = false;

	// Wraps counter and nonterminals stack
	stack<const Symbol*> nonterminals;

	// Iterators
	iterator ruleIt;
	Rule::iterator prodIt;
	Genotype::iterator genoIt=genotype.begin();

	// Start with the start symbol
	nonterminals.push(getStartSymbol());
  unsigned int wraps=0;

	bool gotToUseWrap=false;
	// Get rid of all non-terminal symbols
	while((!nonterminals.empty())&&(wraps<=getMaxWraps())){
	  const Symbol* currentNonTerminal = nonterminals.top();
	  
		// Do a mapping step
		switch(genotype2PhenotypeStep(nonterminals,genoIt,buildDerivationTree)){
			case -1:returnValue=false;
				break;
			case 0:	; // no codon consumed here so don't update the multimap or vector
				break;
			case 1:	
			  // consumed a codon so need to store the nonterminal
			  codon_vector.push_back(*currentNonTerminal);
			  codon_map.insert(pair<string, int>(*currentNonTerminal, newEffectiveSize));
			  genoIt++;
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
	if((wraps>getMaxWraps())||(!nonterminals.empty())){
		returnValue=false;
		// Add remaining symbols in nonterminals queue to phenotype
		while(!nonterminals.empty()){
			phenotype.push_back(nonterminals.top());
			nonterminals.pop();
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
int HemannGrammarSI::getMatchingCodon(string rule){

  int count = codon_map.count(rule);
  int returnIndex = -1;

  if(count > 0){
    pair<multimap<string, int>::iterator, multimap<string, int>::iterator> range =
        codon_map.equal_range(rule);
    // randomly select one of them and return that index for start 
    // of matching codon
    multimap<string, int>::iterator chosenIt = range.first;
    int chosen_iterator = (rand() % count);
    int i=0;
    while(i++ < chosen_iterator){
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
vector<HemannGrammarSI::codonBlocks> HemannGrammarSI::setGenotypeOpt(const GA1DArrayGenome<int> &genome){

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
  
  vector<HemannGrammarSI::codonBlocks> blocks;
  
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
vector<HemannGrammarSI::codonBlocks> HemannGrammarSI::setGenotypeOpt(const GA1DArrayGenome<int> &genome,
  std::set<string> compressed_set, string startSymbol){

  setOptStartSymbol(startSymbol);
  setOptSymbolSet(compressed_set);
  
  return setGenotypeOpt(genome);

}


///////////////////////////////////////////////////////////////////////////////
// Updates the contents of the phenotype structure, based on the current
// genotype and the current grammar, and according to the standard GE
// mapping process. Returns true upon a successful mapping, and false
// otherwise, and also updates the valid field of the phenotype.
// With argument set to true, also updates derivationTree.
bool HemannGrammarSI::genotype2PhenotypeOpt(vector<HemannGrammarSI::codonBlocks>& blocks){
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
	// Wraps counter and nonterminals stack
	unsigned int wraps=0;
	stack<const Symbol*> nonterminals;

	// Iterators
	iterator ruleIt;
	Rule::iterator prodIt;
	Genotype::iterator genoIt=genotype.begin();

	// Start with the start symbol
	nonterminals.push(getStartSymbol());

  const Symbol* startOptSymbol = startOptRulePtr->lhs.front();

  int startOptSite=0;
  bool inOpt=false;
  codonBlocks new_block;


	bool gotToUseWrap=false;
	// Get rid of all non-terminal symbols
	while((!nonterminals.empty())&&(wraps<=getMaxWraps())){
    // check top of nonterminals to see if it is the start symbol for optimization
    if(!inOpt){
      if(startOptSymbol == nonterminals.top()){
        startOptSite = newEffectiveSize;
        inOpt=true;
      }
    }
    else{ // currently running through optimized section of solution
      // check to see if top of nonterminals is still in optimized set
      string ntString = *(nonterminals.top());
      if(opt_symbols.find(ntString) == opt_symbols.end()){
        // found end of the current optimized 
        inOpt=false;
        new_block.start = startOptSite;
        new_block.end = newEffectiveSize;
        blocks.push_back(new_block);
      }
    }

		switch(genotype2PhenotypeStepOpt(nonterminals,genoIt)){
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
	//newEffectiveSize+=(genoIt-genotype.begin());
	// Was the mapping successful?
	if((wraps>getMaxWraps())||(!nonterminals.empty())){
		returnValue=false;
		// Add remaining symbols in nonterminals queue to phenotype
		while(!nonterminals.empty()){
			phenotype.push_back(nonterminals.top());
			nonterminals.pop();
		}
	}
	phenotype.setValid(returnValue);
	genotype.setEffectiveSize(newEffectiveSize);
	genotype.setWraps(wraps);
	return returnValue;
}


///////////////////////////////////////////////////////////////////////////////
// Performs one step of the mapping process, that is, maps the next
// non-terminal symbol on the nonterminals stack passed as argument, using the
// codon at the position pointed by genoIt.
// Returns number of codons consumed, -1 if not successful
int HemannGrammarSI::genotype2PhenotypeStepOpt(stack<const Symbol*> &nonterminals,
		Genotype::iterator genoIt){
		
	#if (DEBUG_LEVEL >= 2)
	cerr << "'bool GEGrammar::genotype2PhenotypeStep(stack<const Symbol*>&,Genotype::iterator&,bool)' called\n";
	#endif
	Rule::iterator prodIt;
	int returnValue;

	// Find the rule for the current non-terminal
	Rule *rulePtr=findRule(*(nonterminals.top()));

	if(!rulePtr){// Undefined symbol - could be an extension symbol
		if(((*(nonterminals.top())).substr(0,strlen("<GECodonValue"))==
			"<GECodonValue")&&(genoIt!=genotype.end())){
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
				codon << (*genoIt%(high-low+1))+low;;
			}
			// Insert symbol with value onto phenotype
			phenotype.push_back(new Symbol(codon.str()));
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
	// Allow recursive rules, but only if they consume a codon
	else if((rulePtr->getMinimumDepth()>=INT_MAX>>1)// Stuck on recursive rule
		&&(rulePtr->size()<=1)){// No codon will be consumed
		// Include symbol on phenotype
		phenotype.push_back(nonterminals.top());
		// Remove non-terminal
		nonterminals.pop();
		// Invalidate mapping
		returnValue=-1;
	} /* end of check for when not an actual rule */
	else{
		// Remove non-terminal
		nonterminals.pop();
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


///
/// Returns vector of ints corresponding to block 
/// @param opt_value Value that needs to be converted into grammar representation
/// @return vector of int with values of codons needed
///
vector<int> HemannGrammarSI::translateOptValue(symbVector& optimized_symb){
   
  stack<gramelement> terms;
  stack<gramelement> workingStack;
  gramelement tempGram;
  
  // set up stack of terminal characters
  for(vector<optSymbol>::reverse_iterator symbIter = optimized_symb.rbegin(); 
    symbIter != optimized_symb.rend(); ++symbIter){
     
      tempGram.symbol = symbIter->symbol;
      tempGram.noNT = symbIter->noNT;
      terms.push(tempGram);
  }
  
  // have left marker and right marker set already (left_opt_bound, right_opt_bound)
  gramelement currTerm;
  gramelement tempTerm;
 
  // continue loop until working stack is empty (will be empty on the first pass)
  do{
  
  while(terms.size()){
    currTerm = terms.top();
    terms.pop();
    // when no nonterminal for this element
    if(currTerm.noNT){
      // now should check to see if the last one was a right_opt_bound character
      // -- ')' -- that marks the end of a functional set
      if(currTerm.symbol[0] == right_opt_bound){
        // marks the end of a set of elements -- will pop off stack until reach the 
        // left hand element that marks the end -- '(' -- usually
        gramelement workingGram = workingStack.top();
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
        }while(workingItem[0] != left_opt_bound || workingItem.size() > 1);
        
        // now check for the rule made out of these parameters
        string leftSide;
        
        // check to see if it needs another element from working stack
        if(is_arg.find(newProduction) != is_arg.end()){
          workingGram = workingStack.top();
          workingStack.pop();
          workingItem = workingGram.symbol;
          tempTerm.codons.insert(tempTerm.codons.begin(), workingGram.codons.begin(), workingGram.codons.end());
          newProduction = workingItem + newProduction;
        }
        
        int codon_value;
        bool rule_found = getRuleFromProduction(newProduction, leftSide, codon_value);
        if(codon_value >= 0){
          // codons.push_back(codon_value);
          // insert at beginning
          tempTerm.codons.insert(tempTerm.codons.begin(), codon_value);
        }
        
        if(leftSide.size() > 0)
          tempTerm.symbol = left_opt_bound + leftSide + right_opt_bound;
        else
          tempTerm.symbol = left_opt_bound + newProduction + right_opt_bound;
          
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
      if(is_arg.find(current) != is_arg.end()){
        currTerm.codons.insert(currTerm.codons.begin(), workingStack.top().codons.begin(), workingStack.top().codons.end());
        current = workingStack.top().symbol + current;
        workingStack.pop();
// cout << "currTerm codons are now ";
// for(unsigned int i=0; i<currTerm.codons.size(); i++){
// cout << currTerm.codons[i] << " ";
// }
// cout << endl;
      }
      
      int codon_value;
      bool rule_found = getRuleFromProduction(current, leftSide, codon_value);
      // repeat this loop as long as the path can be followed up the rules
      while(rule_found){
//         codons.push_back(codon_value);
        if(codon_value >=0)
          currTerm.codons.insert(currTerm.codons.begin(), codon_value);
        currTerm.symbol = leftSide;
        current = leftSide;
        rule_found = getRuleFromProduction(current, leftSide, codon_value);
      }
      
// cout << "currTerm codons are now ";
// for(unsigned int i=0; i<currTerm.codons.size(); i++){
// cout << currTerm.codons[i] << " ";
// }
// cout << endl;
      
      // place last value on to workingStack for later conversion
      workingStack.push(currTerm);
// cout << "pushed on workingStack " << currTerm.symbol << " size of stack is " << workingStack.size() << " => ";
// stack<gramelement> tempstack = workingStack;
// while(tempstack.size() > 0){
// cout << tempstack.top().symbol << " [";
// for(unsigned int i=0; i<tempstack.top().codons.size(); i++){
//   cout << tempstack.top().codons[i] << " ";
// }
// cout << "] ";
// tempstack.pop();
// }
// cout << endl;
    }   
    
  } // terms stack is now empty
  
  
// cout << "cleared the terms stack" << endl;
  
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
// cout << "terms top is " << tempTerm.symbol << endl;
// cout << "start symbol is " << opt_start_symbol << endl;
  // continue running until the only terminal is the starting one for the optimized portion
  // of the grammar
  }while(terms.top().symbol.compare(opt_start_symbol) != 0);


// cout << "top term for return is " << terms.top().symbol << endl;
// for(vector<int>::iterator codoniter=terms.top().codons.begin();  
//   codoniter!=terms.top().codons.end(); ++codoniter){
//     
//   cout << *codoniter << " ";
// }
// cout << endl;

  return terms.top().codons;
  
}


/// 
/// Returns rule (left hand side) that gives the production_string passed
/// @param production_string
/// @param left_hand
/// @param codon
/// @return codon value
///
bool HemannGrammarSI::getRuleFromProduction(string& production_string, string& left_hand, int&
  codon){

  map<string, reverseRule>::iterator revIter = reverseRules.find(production_string);
  if(revIter != reverseRules.end()){

    left_hand = revIter->second.rule;
    
    // adjust codon by unmodding it
    if(revIter->second.codon >= 0)
       codon = revIter->second.codon + static_cast<CodonType>((genotype.getMaxCodonValue()/revIter->second.numcodons*(rand()/(RAND_MAX+1.0))))*revIter->second.numcodons;
    else
      codon = revIter->second.codon;
    return true;
  }
  else{
    left_hand = "";
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
void HemannGrammarSI::constructReverseGrammar(){

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
           
      
      reverse.numcodons = int(ruleIter->size());
        
      reverseRules[rightSideString] = reverse;
    }
    
  }
}



