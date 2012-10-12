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
}

///
/// Returns a terminal associated with the string passed.
/// @param symbol string with name of terminal
/// @return TerminalSymbol pointer
/// @throws AthenaExcept when no match
///
TerminalSymbol* TerminalSymbCreator::get_term(string& symbol){
    if(TerminalMap.find(symbol) == TerminalMap.end())
        throw AthenaExcept("No terminal match for " + symbol);
    
    return TerminalMap[symbol]; 
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
        name = "G" + Stringmanip::itos(g);
        TerminalMap[name] = new GenotypeTerm(name, g);
        
    }
    
}


///
/// Adds continuous variable terminals to map for use in constructing
/// neural networks
/// @param num_variables
///
void TerminalSymbCreator::addContinVariables(int num_variables){
    
    int c;
    string name;
 
    for(c=1; c<=num_variables; c++){
        name = "C" + Stringmanip::itos(c);
        TerminalMap[name] = new ContinVariable(name, c);
    }
    
}


///
/// Sets the individual in the Genotype and Variable classes
/// @param ind Individual 
///
void TerminalSymbCreator::set_ind(Individual* ind){
    GenotypeTerm::set_ind(ind);
    ContinVariable::set_ind(ind);
}


///
/// Function for returning appropriate pointer based
/// on symbol passed in
/// @param symbol for selecting appropriate terminal
/// @param num_args number of arguments the terminal accepts
/// @return new Terminal 
///
void TerminalSymbCreator::create_terminals(int num_genotypes, 
        int num_covariates){
    
  // create terminals
  TerminalMap["+"] = new Addition("+", 2);
  TerminalMap["-"] = new Subtraction("-", 2);
  TerminalMap["*"] = new Multiplication("*", 2);
  TerminalMap["/"] = new Division("/", 2);
  TerminalMap["^"] = new Power("^", 2);
  TerminalMap["PA"] = new pAdd("PA", -1);
  TerminalMap["PS"] = new pSub("PS", -1);
  TerminalMap["PM"] = new pMult("PM", -1);
  TerminalMap["PD"] = new pDiv("PD", -1);
  TerminalMap["W"] = new Weight("W", 2);
  TerminalMap["."] = new Dot (".", 0);
  TerminalMap["Concat"] = new ConCat("Concat", -1);
  TerminalMap["PAND"] = new pAnd("PAND", -1);
  TerminalMap["PNAND"] = new pNand("PNAND", -1);
  TerminalMap["POR"] = new pOr("POR", -1);
  TerminalMap["PNOR"] = new pNor("PNOR",-1);
  TerminalMap["PXOR"] = new pXor("PXOR", -1);
  TerminalMap["AND"] = new And("AND",2);
  TerminalMap["NAND"] = new Nand("NAND",2);
  TerminalMap["OR"] = new Or("OR",2);
  TerminalMap["NOR"] = new Nor("NOR",2);
  TerminalMap["XOR"] = new Xor("XOR",2);
  TerminalMap["("] = new TerminalSymbol("(", 0);
  TerminalMap[")"] = new TerminalSymbol(")", 0);
  TerminalMap[","] = new TerminalSymbol(",", 0);
  TerminalMap["sin"] = new Sine("Sine", 1);
  TerminalMap["log"] = new LogF("Log", 1);
  TerminalMap["cosin"] = new Cosine("Cosine", 1);
  TerminalMap["tan"] = new Tangent("Tangent", 1);
  
  // set the special terminals for quick access
  rparen = TerminalMap[")"];
  lparen = TerminalMap["("];
  commaptr = TerminalMap[","];
  concat = TerminalMap["Concat"];
  
  // create the constants
  TerminalMap["-1"] = create_constant("-1");
  TerminalMap["0"] = create_constant("0");
  TerminalMap["1"] = create_constant("1");
  TerminalMap["2"] = create_constant("2");
  TerminalMap["3"] = create_constant("3");
  TerminalMap["4"] = create_constant("4");
  TerminalMap["5"] = create_constant("5");
  TerminalMap["6"] = create_constant("6");
  TerminalMap["7"] = create_constant("7");
  TerminalMap["8"] = create_constant("8");
  TerminalMap["9"] = create_constant("9");
  TerminalMap["10"] = create_constant("10");
  TerminalMap["11"] = create_constant("11");
  TerminalMap["12"] = create_constant("12");
  TerminalMap["13"] = create_constant("13");
  TerminalMap["14"] = create_constant("14");
  TerminalMap["15"] = create_constant("15");
  TerminalMap["16"] = create_constant("16");
  TerminalMap["17"] = create_constant("17");
  TerminalMap["18"] = create_constant("18");
  TerminalMap["19"] = create_constant("19");
  TerminalMap["20"] = create_constant("20");
  
  addGenotypeVariables(num_genotypes);
  addContinVariables(num_covariates);
  
  // create optSymbol structs for use in creating terminals from constant value
  left_paren_symb.symbol = "(";
  left_paren_symb.noNT = true;
  right_paren_symb.symbol = ")";
  right_paren_symb.noNT = true;
  concat_symb.symbol = "Concat";
  concat_symb.noNT = false;
  period_symb.symbol = ".";
  period_symb.noNT = true;
}

/// 
/// Creates a constant TerminalSymbol
/// @param symbol  
///
TerminalSymbol * TerminalSymbCreator::create_constant(const string& symbol){
  return new Constant(symbol, 0);
}

///
/// Creates a string of symbols for grammar use based on a constant value
/// @param value float version of the number to create in the grammar
///
void TerminalSymbCreator::terminalsFromConstant(float value, 
  symbVector& opt_symbols){
  
  opt_symbols.clear();
  
  stringstream ss;

  bool negative=false;
  
  if(value < 0){
    value = fabs(value);
    negative = true;
  }
  
  optSymbol tempSymb;
  
  // when too small just make the weight be zero and return
  if(value < 0.0001){
    opt_symbols.push_back(concat_symb);
    opt_symbols.push_back(left_paren_symb);
    tempSymb.symbol = "0";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);    
    tempSymb.symbol = "1";
    tempSymb.noNT = true;
    opt_symbols.push_back(tempSymb);
    opt_symbols.push_back(right_paren_symb);
    return;
  }
  
  // need to subtract from 0 for a negative number
  if(negative){
    opt_symbols.push_back(left_paren_symb);
    opt_symbols.push_back(concat_symb);
    opt_symbols.push_back(left_paren_symb);
    tempSymb.symbol = "0";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);
    tempSymb.symbol = "1";
    tempSymb.noNT = true;
    opt_symbols.push_back(tempSymb);
    opt_symbols.push_back(right_paren_symb);   
    tempSymb.symbol = "-";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);

  }
  
  // check size of value to see if it is too big 
  // when too big need to multiply by appropriate value and
  // when too small need to divide by appropriate value
  // -- as of now range is from 99.99 * 99.99 to 0.01 * 0.01
  // or (9998.0001 to 0.0001) 
  if(value > 9998.0001)
    value = 9998.0001;
  
  if(value <= 99.99 and value >= 0.01){
    opt_symbols.push_back(concat_symb);
    opt_symbols.push_back(left_paren_symb); 
    getStringFromNum(value, opt_symbols);
    opt_symbols.push_back(right_paren_symb);
  }
  else if(value > 99.99){
    opt_symbols.push_back(left_paren_symb);
    float newvalue = value / 99.99;
    opt_symbols.push_back(concat_symb);
    opt_symbols.push_back(left_paren_symb);
    
    getStringFromNum(newvalue, opt_symbols);
    opt_symbols.push_back(right_paren_symb);
  
    tempSymb.symbol = "*";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);
    
    opt_symbols.push_back(concat_symb);
    opt_symbols.push_back(left_paren_symb);

    tempSymb.symbol = "9";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);
    
    tempSymb.symbol = "9";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);
    
    opt_symbols.push_back(period_symb);
    tempSymb.symbol = "9";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);    
    
    tempSymb.symbol = "9";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);
    
    tempSymb.symbol = "5";
    tempSymb.noNT = true;
    opt_symbols.push_back(tempSymb);
    
    opt_symbols.push_back(right_paren_symb);
    opt_symbols.push_back(right_paren_symb);
  }
  else if(value < 0.01){
    opt_symbols.push_back(left_paren_symb);

    float newvalue = value / 0.01;
    opt_symbols.push_back(concat_symb);
  
    opt_symbols.push_back(left_paren_symb);
    getStringFromNum(newvalue, opt_symbols);
    opt_symbols.push_back(right_paren_symb);

    tempSymb.symbol = "*";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);    
    opt_symbols.push_back(concat_symb);
    opt_symbols.push_back(left_paren_symb);
    opt_symbols.push_back(period_symb);

    tempSymb.symbol = "0";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb); 

    tempSymb.symbol = "1";
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);     

    tempSymb.symbol = "3";
    tempSymb.noNT = true;
    opt_symbols.push_back(tempSymb);     

    opt_symbols.push_back(right_paren_symb);
    
    opt_symbols.push_back(right_paren_symb);
  }

  if(negative)
    opt_symbols.push_back(right_paren_symb);

}


///
/// converts the value into a string -- assumes values are no greater than 99.99
/// @param value
/// @param symbols to add to 
///
void TerminalSymbCreator::getStringFromNum(float value, symbVector& opt_symbols){
  
  string num_string;
  stringstream ss;
  ss << value;
  num_string = ss.str(); 
  // no decimal point
  if(int(value) <= 9 && num_string.find(".")==string::npos){ 
    ss.str("");
    ss << int(value);
    num_string = ss.str();
  }
  else{ // make sure there are 2 numbers after decimal point
    ss.str("");
   ss << setiosflags(ios::fixed) << setprecision(2) << value;
    num_string = ss.str();
  }
  
  optSymbol tempSymb;
  // add symbols to the symbVector (each digit and the decimal)
  for(string::iterator striter = num_string.begin(); striter != num_string.end();
    striter++){
    tempSymb.symbol = *striter;
    tempSymb.noNT = false;
    opt_symbols.push_back(tempSymb);
  }
  
  ss.str("");
  ss << num_string.length();
  tempSymb.symbol = ss.str();
  tempSymb.noNT = true;
  opt_symbols.push_back(tempSymb);   
}


