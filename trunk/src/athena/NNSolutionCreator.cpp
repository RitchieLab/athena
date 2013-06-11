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
#include "NNSolutionCreator.h"
#include "TerminalSymbol.h"
#include "Terminals.h"
#include "ExpressionTree.h"
#include "BackPropAnnie.h"

///
/// Constructor
///
NNSolutionCreator::NNSolutionCreator(){
		initialize();
}

///
/// Alternative constructor
///
NNSolutionCreator::NNSolutionCreator(vector<string>& symbols){
	initialize();
}


///
/// Initializes variables
///
void NNSolutionCreator::initialize(){
	terminalsSet = false;
	calculator = NULL;
	graphicExtension = ".dot";
	startOpt = "<cop>";
	
	string symbols[] = {"<cop>","<op>","<Concat>","(",")", "<num>","Concat","1","2","3","4",
		"5","6","7","8","9","0","<dig>","+","-","*","/"};
	
	// set up list of symbols included in setting weights in neural network
	optSymbols.insert(symbols, symbols+22);
	
	// set up symbols that take an argument  
	optArgSymbols.insert("(<num>)");
	
	leftOptBound = '(';
	rightOptBound = ')';
}


///
/// Creates neural network for evaluation
/// sets up a postfix stack for evaluation of the neural network
/// @param pheno Phenotype from libGE that can be turned into a neural network
/// @param set Dataset
///
void NNSolutionCreator::establishSolution(vector<string>& symbols, Dataset* set){
		if(!terminalsSet)
				setVariables(set);
		establishSolution(symbols);
}


///
/// Creates neural network for evaluation
/// sets up a postfix stack for evaluation of the neural network
/// @param symbols Phenotype from libGE that can be turned into a neural network
///
void NNSolutionCreator::establishSolution(vector<string>& symbols){
	 
	 nnDepth = 0;
	 nnTerminalSize= symbols.size();
	 
		// as part of stack look for variables and keep a list of those
		// to use in evaluation for skipping missing dataPq
		
		// Anything except parentheses with a priority of zero is pushed to the postfix stack
		// Left paren always pushed on stack
		// When get to right paren, pop off stack until left paren encountered
		// if the symbol being scanned has a higher precedence than the symbol at the
		//     top of the stack, the symbol scanned is pushed on the stack
		// If the precedence of the symbol being scanned is lower than or equal to the
		//     precedence of the symbol at the top of the stack, the stack is popped to the output
		// When the terminating symbol is reached on the input scan, the stack is popped to the
		//     output until the terminating symbol is also reached on the stack. Then the algorithm terminates.
	 
		vector<TerminalSymbol *> tempPostFix;
	 
		postFixStack.clear();
		
		vector<TerminalSymbol *> symbol_stack;
		
		vector<TerminalSymbol *> stack;
		
		deque<float> args, blank;
		TerminalSymbol* constantTerminalPtr;

		unsigned int inputSize = symbols.size();
		for(unsigned int inputCount=0; inputCount < inputSize; inputCount++){
				// priority is zero when element is a number, variable, comma or parentheses
				// may need a dynamic cast
				TerminalSymbol* currTerm = termHolder.getTerm(symbols[inputCount]);
				
				// add check for concatenation operator here and replace with constant if necessary
				if(currTerm == termHolder.concaten()){
					args.clear();        
					inputCount+=2; //skip next one which will be a '('
					while(termHolder.getTerm(symbols[inputCount]) != termHolder.rightParen()){
						args.push_back(termHolder.getTerm(symbols[inputCount++])->evaluate(blank));
					}
					// have to pop off last one which is variable value
					args.pop_back();
					constantTerminalPtr = new Constant(termHolder.concaten()->evaluate(args));
					tempPostFix.push_back(constantTerminalPtr);
					constants.push_back(constantTerminalPtr);
					continue; // go to next string
				}
				
				if(!currTerm->getPriority()){
				 if(currTerm == termHolder.leftParen()){
						stack.push_back(currTerm);
				 }
				 else if(currTerm == termHolder.rightParen()){
					while(stack.back() != termHolder.leftParen()){
						tempPostFix.push_back(stack.back());            
						stack.pop_back();
					}
					// pop off left paren
					stack.pop_back();
				}
				else if(currTerm == termHolder.comma()){
					while(stack.back() != termHolder.leftParen()){
						tempPostFix.push_back(stack.back());
						stack.pop_back();
					}
					// in this case do not pop off left paren as
					// it is being used to mark off arguments
					// to function calls
				}
				else{
					//  numbers and variables are pushed to output
					tempPostFix.push_back(currTerm);
				}
			}
			else{ // element is an operator
				while(stack.size() > 0 && 
							(currTerm->getPriority() <= stack.back()->getPriority())){
						tempPostFix.push_back(stack.back());
						stack.pop_back();
					}
				stack.push_back(currTerm);
			}
	}
 
	// after last terminal push anything on the stack to the postfix stack
	for(vector<TerminalSymbol *>::iterator stackIter=stack.begin(); stackIter != stack.end(); stackIter++)
		tempPostFix.push_back(*stackIter);
		
	// now convert any constants and store in postFixStack
	compressOperator(tempPostFix, postFixStack);
		
	// iterate through stack and store covars and genotypes for use in checking
	// for missing data
	genos.clear();
	covars.clear();
 
	for(unsigned int i=0; i<postFixStack.size(); i++){

			if(postFixStack[i]->getTermType() == TerminalSymbol::Genotype){
				 genos[postFixStack[i]]++;
			}
			else if (postFixStack[i]->getTermType() == TerminalSymbol::Covariate){
				 covars[postFixStack[i]]++;
			}
	}
}



///
/// optimize solution by running back propagation
/// @param symbols vector of strings the can be converted into network
/// @param set Dataset for optimizing solution
/// @returns number of epochs trained
///
int NNSolutionCreator::optimizeSolution(std::vector<std::string>& symbols, Dataset* set){

	BackPropAnnie bp;
	// clear old symbols
	optValSymbols.clear();

	// this creates the postFixStack
	establishSolution(symbols, set);
 
	int numEpochsTrained = bp.runBackProp(postFixStack, set);
	
	optValues = bp.getWeights();
	optimizedScore = bp.getOptimizedScore();
	
	vector<float>::iterator iter;
	
	symbVector newSymbols;
	for(iter = optValues.begin(); iter != optValues.end(); ++iter){
		termHolder.terminalsFromConstant(*iter, newSymbols);
		optValSymbols.push_back(newSymbols);
	}
	return numEpochsTrained;
}


///
/// Compresses the operator calculations for generating 
///
void NNSolutionCreator::compressOperator(vector<TerminalSymbol*> & postFixStack,
	vector<TerminalSymbol*>& new_stack){
	// for any operator compresses them so that stack will not have redundant information
	// any operator can be compressed into a constant value
	
	// 1.  for any non-constant / non-operator push on to new stack
	// 2.  when find constant evaluate until find non-constant and then take that value from
	// stack and create a new constant for the new stack
	
	TerminalSymbol* newConstant;
	
	vector<float> stack;
	deque<float> args;
	int numArgs;
	
	vector<float>::iterator iter;
	
	for(unsigned int i=0; i < postFixStack.size(); i++){
		
		TerminalSymbol::TerminalType termType = postFixStack[i]->getTermType();
		
		if(termType == TerminalSymbol::Operator){
			// when it is an operator evaluate and push back on to stack
			numArgs = postFixStack[i]->getNumArgs();
			for(int k=0; k<numArgs; k++){
				args.push_front(stack.back());
				stack.pop_back();
			}
		 stack.push_back(postFixStack[i]->evaluate(args));
		}
		else if(termType == TerminalSymbol::Constant){
			args.clear();
			stack.push_back(postFixStack[i]->evaluate(args));
		}
		else{
			// when not an operator or constant take the argument from the stack and
			// convert it to a Constant and then put on the new_stack followed by the
			// original terminal
			// start at bottom of stack so that you get any constants that need to be
			// carried over for later evaluation
			for(iter = stack.begin(); iter != stack.end(); ++iter){
				newConstant = new Constant(*iter);
				constants.push_back(newConstant);
				new_stack.push_back(newConstant);
			}
			stack.clear();
			
			new_stack.push_back(postFixStack[i]);
			
		}
		
	} 
}



///
/// Returns vector with indexes of SNPs in the latest solution
/// @return vector of indexes
///
vector<int> NNSolutionCreator::getGeneIndexes(){
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
inline vector<int> NNSolutionCreator::getCovarIndexes(){
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
void NNSolutionCreator::setVariables(Dataset* set){
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
float NNSolutionCreator::evaluateInd(Individual* ind){

		termHolder.setInd(ind);
		
		vector<float> stack;
		deque<float> args;
		
		int numArgs=0, k;
		unsigned int postFixSize = postFixStack.size();
		
		for(unsigned int j=0; j<postFixSize; j++){
				args.clear();
				numArgs = postFixStack[j]->getNumArgs();
				// evaluate and push any elements that don't take arguments
				if(!(numArgs=postFixStack[j]->getNumArgs())){
						stack.push_back(postFixStack[j]->evaluate(args));
				}
				 // need to determine number of arguments
				 // pop off arguments and then push back new result
				 // when number of args for the element is < 0 
				 // top number on stack should be the number of 
				 // arguments for the function
				else{
						if(numArgs<0){
								numArgs = (unsigned int)stack.back();
								stack.pop_back();
						}
						// construct deque of arguments
						for(k=0; k<numArgs; k++){
								args.push_front(stack.back());
								stack.pop_back();
						}
						// evaluate current operator and push onto stack
						stack.push_back(postFixStack[j]->evaluate(args));
				}
		}
 
		// only value on stack should be final value
		return stack.back();
}


///
/// Checks to see if the individual has any missing data
/// @param ind
/// @param set
/// @return true if the individual should be used
///
bool NNSolutionCreator::useInd(Individual* ind, Dataset* set){
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
/// Evaluates the neural network and returns the balanced accuracy
/// of the neural network.
/// @param set Dataset to use
///
float NNSolutionCreator::evaluate(Dataset* set){
	 
		// to check for missing data use the list of variables gotten from the
		// creation of the stack and then check each ind to  make sure
		// there is a value in each
		Individual * ind;
		nIndsEvaluated = 0;

		calculator->reset();

		// when missing skip that ind
		for(unsigned int i=0; i < set->numInds(); i++){
				ind = (*set)[i];
		
				ContinVariable::setInd(ind);
				GenotypeTerm::setInd(ind);
			
				if(!useInd(ind, set))
						continue;

				nIndsEvaluated++;
				
				calculator->addIndScore(evaluateInd(ind), ind->getStatus());
		}
		return calculator->getScore();  
}



///
/// For Networks, the detailed logging determines the depth in nodes of the
/// model.
///
void NNSolutionCreator::detailedLogging(){
	ExpressionTree extree;
	extree.convertPostFix(postFixStack);
	nnDepth = extree.getMaxDepth();
}

///
/// Return the network depth
///
unsigned int NNSolutionCreator::getDetailedLog(){
		return nnDepth;
}

///
/// writes a dot compatible text file representing the network
/// @param os ostream to write to
///
void NNSolutionCreator::graphicalOutput(ostream& os, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed){
	
	// network will be part of this tree
	ExpressionTree extree;
	
	extree.convertPostFix(postFixStack);
	extree.outputDot(os,holder, mapUsed, ottDummy, continMapUsed);
	extree.clearConstants();
}


/// 
/// Performs evaluation on the dataset.  Each individual has their result passed
/// to the output stream provided.
/// @param set Dataset 
/// @param os ostream 
/// @param results TestResult vector for holding result values
///
float NNSolutionCreator::evaluateWithOutput(Dataset* set, ostream& os){

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
void NNSolutionCreator::evaluateForOutput(Dataset* set){
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
				
				tempResult.score = evaluateInd(ind);;
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


