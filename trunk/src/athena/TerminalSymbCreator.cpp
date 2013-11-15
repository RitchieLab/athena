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

#include "TerminalSymbCreator.h"
#include "Terminals.h"
#include <Stringmanip.h>
#include <sstream>

#include <iostream>
#include <iomanip>

#include "GEObjective.h"

using namespace data_manage;

///
/// Terminal creator class
///
TerminalSymbCreator::TerminalSymbCreator(){
	grammarOptimizer=OptGrammarFactory::getFactory().create("COMPLEX");
}


///
/// Terminal destructor
///
TerminalSymbCreator::~TerminalSymbCreator(){
	if(grammarOptimizer != NULL){
		delete grammarOptimizer;
	}
}

///
/// Sets optimization grammar type
///
void TerminalSymbCreator::setGrammerOptimization(string optName){
	if(grammarOptimizer != NULL){
		delete grammarOptimizer;
	}
	grammarOptimizer=OptGrammarFactory::getFactory().create(optName);
}


///
/// Returns a terminal associated with the string passed.
/// @param symbol string with name of terminal
/// @return TerminalSymbol pointer
/// @throws AthenaExcept when no match
///
TerminalSymbol* TerminalSymbCreator::getTerm(string& symbol){
		if(terminalMap.find(symbol) == terminalMap.end())
				throw AthenaExcept("No terminal match for " + symbol);
		
		return terminalMap[symbol]; 
}


///
/// Adds variable terminals to map for use in constructing 
/// neural networks.
/// @param num_variables
/// 
void TerminalSymbCreator::addGenotypeVariables(int num_variables){
		int g;
		string name;
		for(g=1; g<=num_variables; g++){
				name = "G" + Stringmanip::numberToString(g);
				terminalMap[name] = new GenotypeTerm(name, g);
				
		}
		
}


///
/// Adds continuous variable terminals to map for use in constructing
/// neural networks
/// @param numVariables
///
void TerminalSymbCreator::addContinVariables(int numVariables){ 
		int c;
		string name;
 
		for(c=1; c<=numVariables; c++){
				name = "C" + Stringmanip::numberToString(c);
				terminalMap[name] = new ContinVariable(name, c);
		}
		
}


///
/// Sets the individual in the Genotype and Variable classes
/// @param ind Individual 
///
void TerminalSymbCreator::setInd(Individual* ind){
		GenotypeTerm::setInd(ind);
		ContinVariable::setInd(ind);
}


///
/// Function for returning appropriate pointer based
/// on symbol passed in
/// @param symbol for selecting appropriate terminal
/// @param num_args number of arguments the terminal accepts
/// @return new Terminal 
///
void TerminalSymbCreator::createTerminals(int numGenotypes, 
				int numCovariates){
		
	// create terminals
	terminalMap["+"] = new Addition("+", 2);
	terminalMap["-"] = new Subtraction("-", 2);
	terminalMap["*"] = new Multiplication("*", 2);
	terminalMap["/"] = new Division("/", 2);
	terminalMap["^"] = new Power("^", 2);
	terminalMap["PA"] = new pAdd("PA", -1);
	terminalMap["PS"] = new pSub("PS", -1);
	terminalMap["PM"] = new pMult("PM", -1);
	terminalMap["PD"] = new pDiv("PD", -1);
	terminalMap["W"] = new Weight("W", 2);
	terminalMap["."] = new Dot (".", 0);
	terminalMap["Concat"] = new ConCat("Concat", -1);
	terminalMap["PAND"] = new pAnd("PAND", -1);
	terminalMap["PNAND"] = new pNand("PNAND", -1);
	terminalMap["POR"] = new pOr("POR", -1);
	terminalMap["PNOR"] = new pNor("PNOR",-1);
	terminalMap["PXOR"] = new pXor("PXOR", -1);
	terminalMap["AND"] = new And("AND",2);
	terminalMap["NAND"] = new Nand("NAND",2);
	terminalMap["OR"] = new Or("OR",2);
	terminalMap["NOR"] = new Nor("NOR",2);
	terminalMap["XOR"] = new Xor("XOR",2);
	terminalMap["("] = new TerminalSymbol("(", 0);
	terminalMap[")"] = new TerminalSymbol(")", 0);
	terminalMap[","] = new TerminalSymbol(",", 0);
	terminalMap["sin"] = new Sine("Sine", 1);
	terminalMap["log"] = new LogF("Log", 1);
	terminalMap["cosin"] = new Cosine("Cosine", 1);
	terminalMap["tan"] = new Tangent("Tangent", 1);
	terminalMap["Bias"] = new BiasTerm("Bias", 0);
	
	// set the special terminals for quick access
	rParen = terminalMap[")"];
	lParen = terminalMap["("];
	commaPtr = terminalMap[","];
	concat = terminalMap["Concat"];
	
	// create the constants
	terminalMap["-1"] = createConstant("-1");
	terminalMap["0"] = createConstant("0");
	terminalMap["1"] = createConstant("1");
	terminalMap["2"] = createConstant("2");
	terminalMap["3"] = createConstant("3");
	terminalMap["4"] = createConstant("4");
	terminalMap["5"] = createConstant("5");
	terminalMap["6"] = createConstant("6");
	terminalMap["7"] = createConstant("7");
	terminalMap["8"] = createConstant("8");
	terminalMap["9"] = createConstant("9");
	terminalMap["10"] = createConstant("10");
	terminalMap["11"] = createConstant("11");
	terminalMap["12"] = createConstant("12");
	terminalMap["13"] = createConstant("13");
	terminalMap["14"] = createConstant("14");
	terminalMap["15"] = createConstant("15");
	terminalMap["16"] = createConstant("16");
	terminalMap["17"] = createConstant("17");
	terminalMap["18"] = createConstant("18");
	terminalMap["19"] = createConstant("19");
	terminalMap["20"] = createConstant("20");
	
	addGenotypeVariables(numGenotypes);
	addContinVariables(numCovariates);
	
	// create optSymbol structs for use in creating terminals from constant value
// 	leftParenSymb.symbol = "(";
// 	leftParenSymb.noNT = true;
// 	rightParenSymb.symbol = ")";
// 	rightParenSymb.noNT = true;
// 	concatSymb.symbol = "Concat";
// 	concatSymb.noNT = false;
// 	periodSymb.symbol = ".";
// 	periodSymb.noNT = true;
}

/// 
/// Creates a constant TerminalSymbol
/// @param symbol  
///
TerminalSymbol * TerminalSymbCreator::createConstant(const string& symbol){
	return new Constant(symbol, 0);
}


void TerminalSymbCreator::addConstant(const string& symbol){
	terminalMap[symbol]=createConstant(symbol);
	
// 	constantSet.insert(terminalMap[symbol]);
	grammarOptimizer->addNumberConstant(terminalMap[symbol]);
}


// TerminalSymbol* TerminalSymbCreator::getClosestConstant(float value){
// 	ConstantPointer temp;
// 	temp.value = value;
// 	std::set<ConstantPointer>::iterator upperIter = constantSet.upper_bound(temp);
// 	std::set<ConstantPointer>::iterator lowerIter = constantSet.lower_bound(temp);
// 	if((value-lowerIter->value) < (upperIter->value-value))
// 		return lowerIter->cons;
// 	else
// 		return upperIter->cons;
// }


/// 
/// Creates a string of symbols for grammar use based on a constant value
/// @param value float version of the number to create in the grammar
///
void TerminalSymbCreator::terminalsFromConstant(float value, 
	symbVector& optSymbols){
	
	grammarOptimizer->terminalsFromConstant(value, optSymbols);
	
// 	optSymbols.clear();
// 	
// 	stringstream ss;
// 
// 	bool negative=false;
// 	
// 	if(value < 0){
// 		value = fabs(value);
// 		negative = true;
// 	}
// 	
// 	optSymbol tempSymb;
// 	
// 	// when too small just make the weight be zero and return
// 	if(value < 0.0001){
// 		optSymbols.push_back(concatSymb);
// 		optSymbols.push_back(leftParenSymb);
// 		tempSymb.symbol = "0";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);    
// 		tempSymb.symbol = "1";
// 		tempSymb.noNT = true;
// 		optSymbols.push_back(tempSymb);
// 		optSymbols.push_back(rightParenSymb);
// 		return;
// 	}
// 	
// 	// need to subtract from 0 for a negative number
// 	if(negative){
// 		optSymbols.push_back(leftParenSymb);
// 		optSymbols.push_back(concatSymb);
// 		optSymbols.push_back(leftParenSymb);
// 		tempSymb.symbol = "0";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);
// 		tempSymb.symbol = "1";
// 		tempSymb.noNT = true;
// 		optSymbols.push_back(tempSymb);
// 		optSymbols.push_back(rightParenSymb);   
// 		tempSymb.symbol = "-";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);
// 
// 	}
// 	
// 	// check size of value to see if it is too big 
// 	// when too big need to multiply by appropriate value and
// 	// when too small need to divide by appropriate value
// 	// -- as of now range is from 99.99 * 99.99 to 0.01 * 0.01
// 	// or (9998.0001 to 0.0001) 
// 	if(value > 9998.0001)
// 		value = 9998.0001;
// 	
// 	if(value <= 99.99 and value >= 0.01){
// 		optSymbols.push_back(concatSymb);
// 		optSymbols.push_back(leftParenSymb); 
// 		getStringFromNum(value, optSymbols);
// 		optSymbols.push_back(rightParenSymb);
// 	}
// 	else if(value > 99.99){
// 		optSymbols.push_back(leftParenSymb);
// 		float newvalue = value / 99.99;
// 		optSymbols.push_back(concatSymb);
// 		optSymbols.push_back(leftParenSymb);
// 		
// 		getStringFromNum(newvalue, optSymbols);
// 		optSymbols.push_back(rightParenSymb);
// 	
// 		tempSymb.symbol = "*";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);
// 		
// 		optSymbols.push_back(concatSymb);
// 		optSymbols.push_back(leftParenSymb);
// 
// 		tempSymb.symbol = "9";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);
// 		
// 		tempSymb.symbol = "9";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);
// 		
// 		optSymbols.push_back(periodSymb);
// 		tempSymb.symbol = "9";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);    
// 		
// 		tempSymb.symbol = "9";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);
// 		
// 		tempSymb.symbol = "5";
// 		tempSymb.noNT = true;
// 		optSymbols.push_back(tempSymb);
// 		
// 		optSymbols.push_back(rightParenSymb);
// 		optSymbols.push_back(rightParenSymb);
// 	}
// 	else if(value < 0.01){
// 		optSymbols.push_back(leftParenSymb);
// 
// 		float newvalue = value / 0.01;
// 		optSymbols.push_back(concatSymb);
// 	
// 		optSymbols.push_back(leftParenSymb);
// 		getStringFromNum(newvalue, optSymbols);
// 		optSymbols.push_back(rightParenSymb);
// 
// 		tempSymb.symbol = "*";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);    
// 		optSymbols.push_back(concatSymb);
// 		optSymbols.push_back(leftParenSymb);
// 		optSymbols.push_back(periodSymb);
// 
// 		tempSymb.symbol = "0";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb); 
// 
// 		tempSymb.symbol = "1";
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);     
// 
// 		tempSymb.symbol = "3";
// 		tempSymb.noNT = true;
// 		optSymbols.push_back(tempSymb);     
// 
// 		optSymbols.push_back(rightParenSymb);
// 		
// 		optSymbols.push_back(rightParenSymb);
// 	}
// 
// 	if(negative)
// 		optSymbols.push_back(rightParenSymb);

}


///
/// converts the value into a string -- assumes values are no greater than 99.99
/// @param value
/// @param symbols to add to 
///
// void TerminalSymbCreator::getStringFromNum(float value, symbVector& optSymbols){
// 	
// 	string numString;
// 	stringstream ss;
// 	ss << value;
// 	numString = ss.str(); 
// 	// no decimal point
// 	if(int(value) <= 9 && numString.find(".")==string::npos){ 
// 		ss.str("");
// 		ss << int(value);
// 		numString = ss.str();
// 	}
// 	else{ // make sure there are 2 numbers after decimal point
// 		ss.str("");
// 	 ss << setiosflags(ios::fixed) << setprecision(2) << value;
// 		numString = ss.str();
// 	}
// 	
// 	optSymbol tempSymb;
// 	// add symbols to the symbVector (each digit and the decimal)
// 	for(string::iterator striter = numString.begin(); striter != numString.end();
// 		striter++){
// 		tempSymb.symbol = *striter;
// 		tempSymb.noNT = false;
// 		optSymbols.push_back(tempSymb);
// 	}
// 	
// 	ss.str("");
// 	ss << numString.length();
// 	tempSymb.symbol = ss.str();
// 	tempSymb.noNT = true;
// 	optSymbols.push_back(tempSymb);   
// }


