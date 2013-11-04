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
#include "SymRegressSolutionCreator.h"
#include "TerminalSymbol.h"
#include "Terminals.h"
#include "ExpressionTree.h"

///
/// sets up a postfix stack for evaluation of the symbolic regression equation
/// @param pheno Phenotype from libGE that can be turned into a symbolic regression equation
/// @param set Dataset
///
void SymRegressSolutionCreator::establishSolution(vector<string>& symbols, Dataset* set){
		if(!terminalsSet)
				setVariables(set);
		establishSolution(symbols);
}

///
/// Creates symbolic regression equation for evaluation
/// sets up a postfix stack for evaluation of the symbolic regression equation
/// @param symbols Phenotype from libGE that can be turned into a symbolic regression equation
///
void SymRegressSolutionCreator::establishSolution(vector<string>& symbols){
	 
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
	 
		postFixStack.clear();
		
		vector<TerminalSymbol *> symbolStack;
		
		vector<TerminalSymbol *> stack;
		
		deque<float> args, blank;
		TerminalSymbol* constantTerminalPtr;

		unsigned int inputSize = symbols.size();
		for(unsigned int inputCount=0; inputCount < inputSize; inputCount++){
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
					postFixStack.push_back(constantTerminalPtr);
					constants.push_back(constantTerminalPtr);
					continue; // go to next string
				}
				
				if(!currTerm->getPriority()){
				 if(currTerm == termHolder.leftParen()){
						stack.push_back(currTerm);
				 }
				 else if(currTerm == termHolder.rightParen()){
					while(stack.back() != termHolder.leftParen()){
						postFixStack.push_back(stack.back());            
						stack.pop_back();
					}
					// pop off left paren
					stack.pop_back();
				}
				else if(currTerm == termHolder.comma()){
					while(stack.back() != termHolder.leftParen()){
						postFixStack.push_back(stack.back());
						stack.pop_back();
					}
					// in this case do not pop off left paren as
					// it is being used to mark off arguments
					// to function calls
				}
				else{
					//  numbers and variables are pushed to output
					postFixStack.push_back(currTerm);
				}
			}
			else{ // element is an operator
				while(stack.size() > 0 && 
							(currTerm->getPriority() <= stack.back()->getPriority())){
						postFixStack.push_back(stack.back());
						stack.pop_back();
					}
				stack.push_back(currTerm);
			}
	}
 
 
	// after last terminal push anything on the stack to the postfix stack
	for(vector<TerminalSymbol *>::iterator stackIter=stack.begin(); stackIter != stack.end(); stackIter++)
		postFixStack.push_back(*stackIter);

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
/// evaluates single individual and returns value for that individual
///
float SymRegressSolutionCreator::evaluateInd(Individual* ind){
	return TerminalSymbol::ActivateSigmoid(NNSolutionCreator::evaluateInd(ind));
}



///
/// Creates neural network for output as equation
/// sets up a postfix stack for evaluation of the neural network
/// @param pheno Phenotype from libGE that can be turned into a neural network
/// @param set Dataset
///
void SymRegressSolutionCreator::establishSolutionEquation(vector<string>& symbols){
	equationStack.clear();
	for(vector<string>::iterator iter=symbols.begin(); iter != symbols.end(); ++iter){
		if(iter->compare("Concat")==0){
			// first will be "("
			iter+=2;
			vector<string> constTerms;
			while(iter->compare(")")!=0){
				constTerms.push_back(*iter);
				++iter;
			}
			string constantTerm;
			for(unsigned int i=0; i<constTerms.size()-1;i++){	
				constantTerm += constTerms[i];
			}
			equationStack.push_back(constantTerm);
		}
		else{
			equationStack.push_back(*iter);
		}
	}
}


///
/// writes output as an equation	
///
void SymRegressSolutionCreator::equationOutput(ostream& os, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed){		
	for(vector<string>::iterator iter=equationStack.begin(); iter != equationStack.end(); ++iter){
		os << *iter;
	}
}

