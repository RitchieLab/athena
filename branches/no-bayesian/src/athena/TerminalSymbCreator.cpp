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
/// Creates name for genotype variable
/// @param varIndex variable index
/// @returns Name for use in mapping
///
std::string TerminalSymbCreator::getGenoName(int varIndex){
	return "G" + Stringmanip::numberToString(varIndex+1);
}

///
/// Adds variable terminals to map for use in constructing
/// neural networks.
/// @param num_variables
///
void TerminalSymbCreator::addGenotypeVariables(int numVariables){
		int g;
		string name;
		for(g=0; g<numVariables; g++){
				name = getGenoName(g);
				terminalMap[name] = new GenotypeTerm(name, g);

		}

}


///
/// Creates name for continuous variable
/// @param varIndex variable index
/// @returns Name for use in mapping
///
std::string TerminalSymbCreator::getContinName(int varIndex){
	return "C" + Stringmanip::numberToString(varIndex+1);
}


///
/// Adds continuous variable terminals to map for use in constructing
/// neural networks
/// @param numVariables
///
void TerminalSymbCreator::addContinVariables(int numVariables){
		int c;
		string name;

		for(c=0; c<numVariables; c++){
				name = getContinName(c);
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
	terminalMap["->"] = new TerminalSymbol("^", 0);
	terminalMap["nc"] = new TerminalSymbol("nc",0);
	terminalMap["pheno"] = new PhenotypeTerm("Pheno");

	// set the special terminals for quick access
	rParen = terminalMap[")"];
	lParen = terminalMap["("];
	commaPtr = terminalMap[","];
	concat = terminalMap["Concat"];
	connect = terminalMap["^"];
	phenoPtr = terminalMap["pheno"];

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

	grammarOptimizer->addNumberConstant(terminalMap[symbol]);
}


///
/// Creates a string of symbols for grammar use based on a constant value
/// @param value float version of the number to create in the grammar
///
void TerminalSymbCreator::terminalsFromConstant(float value,
	symbVector& optSymbols){

	grammarOptimizer->terminalsFromConstant(value, optSymbols);
}



