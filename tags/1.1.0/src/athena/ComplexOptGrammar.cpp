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
#include "ComplexOptGrammar.h"
#include<sstream>
#include <iostream>
#include <iomanip>

string const ComplexOptGrammar::optGrammarName = ComplexOptGrammar::registerCalc("COMPLEX");


ComplexOptGrammar::ComplexOptGrammar(){
	
		// create optSymbol structs for use in creating terminals from constant value
	leftParenSymb.symbol = "(";
	leftParenSymb.noNT = true;
	rightParenSymb.symbol = ")";
	rightParenSymb.noNT = true;
	concatSymb.symbol = "Concat";
	concatSymb.noNT = false;
	periodSymb.symbol = ".";
	periodSymb.noNT = true;

}


/// 
/// Creates a string of symbols for grammar use based on a constant value
/// @param value float version of the number to create in the grammar
///
void ComplexOptGrammar::terminalsFromConstant(float value, 
	symbVector& optSymbols){
	
	optSymbols.clear();
	
	stringstream ss;

	bool negative=false;
	
	if(value < 0){
		value = fabs(value);
		negative = true;
	}
	
	optSymbol tempSymb;
	
	// when too small just make the weight be zero and return
	if(value < 0.0001){
		optSymbols.push_back(concatSymb);
		optSymbols.push_back(leftParenSymb);
		tempSymb.symbol = "0";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);    
		tempSymb.symbol = "1";
		tempSymb.noNT = true;
		optSymbols.push_back(tempSymb);
		optSymbols.push_back(rightParenSymb);
		return;
	}
	
	// need to subtract from 0 for a negative number
	if(negative){
		optSymbols.push_back(leftParenSymb);
		optSymbols.push_back(concatSymb);
		optSymbols.push_back(leftParenSymb);
		tempSymb.symbol = "0";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);
		tempSymb.symbol = "1";
		tempSymb.noNT = true;
		optSymbols.push_back(tempSymb);
		optSymbols.push_back(rightParenSymb);   
		tempSymb.symbol = "-";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);

	}
	
	// check size of value to see if it is too big 
	// when too big need to multiply by appropriate value and
	// when too small need to divide by appropriate value
	// -- as of now range is from 99.99 * 99.99 to 0.01 * 0.01
	// or (9998.0001 to 0.0001) 
	if(value > 9998.0001)
		value = 9998.0001;
	
	if(value <= 99.99 and value >= 0.01){
		optSymbols.push_back(concatSymb);
		optSymbols.push_back(leftParenSymb); 
		getStringFromNum(value, optSymbols);
		optSymbols.push_back(rightParenSymb);
	}
	else if(value > 99.99){
		optSymbols.push_back(leftParenSymb);
		float newvalue = value / 99.99;
		optSymbols.push_back(concatSymb);
		optSymbols.push_back(leftParenSymb);
		
		getStringFromNum(newvalue, optSymbols);
		optSymbols.push_back(rightParenSymb);
	
		tempSymb.symbol = "*";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);
		
		optSymbols.push_back(concatSymb);
		optSymbols.push_back(leftParenSymb);

		tempSymb.symbol = "9";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);
		
		tempSymb.symbol = "9";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);
		
		optSymbols.push_back(periodSymb);
		tempSymb.symbol = "9";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);    
		
		tempSymb.symbol = "9";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);
		
		tempSymb.symbol = "5";
		tempSymb.noNT = true;
		optSymbols.push_back(tempSymb);
		
		optSymbols.push_back(rightParenSymb);
		optSymbols.push_back(rightParenSymb);
	}
	else if(value < 0.01){
		optSymbols.push_back(leftParenSymb);

		float newvalue = value / 0.01;
		optSymbols.push_back(concatSymb);
	
		optSymbols.push_back(leftParenSymb);
		getStringFromNum(newvalue, optSymbols);
		optSymbols.push_back(rightParenSymb);

		tempSymb.symbol = "*";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);    
		optSymbols.push_back(concatSymb);
		optSymbols.push_back(leftParenSymb);
		optSymbols.push_back(periodSymb);

		tempSymb.symbol = "0";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb); 

		tempSymb.symbol = "1";
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);     

		tempSymb.symbol = "3";
		tempSymb.noNT = true;
		optSymbols.push_back(tempSymb);     

		optSymbols.push_back(rightParenSymb);
		
		optSymbols.push_back(rightParenSymb);
	}

	if(negative)
		optSymbols.push_back(rightParenSymb);

}


///
/// converts the value into a string -- assumes values are no greater than 99.99
/// @param value
/// @param symbols to add to 
///
void ComplexOptGrammar::getStringFromNum(float value, symbVector& optSymbols){
	
	string numString;
	stringstream ss;
	ss << value;
	numString = ss.str(); 
	// no decimal point
	if(int(value) <= 9 && numString.find(".")==string::npos){ 
		ss.str("");
		ss << int(value);
		numString = ss.str();
	}
	else{ // make sure there are 2 numbers after decimal point
		ss.str("");
	 ss << setiosflags(ios::fixed) << setprecision(2) << value;
		numString = ss.str();
	}
	
	optSymbol tempSymb;
	// add symbols to the symbVector (each digit and the decimal)
	for(string::iterator striter = numString.begin(); striter != numString.end();
		striter++){
		tempSymb.symbol = *striter;
		tempSymb.noNT = false;
		optSymbols.push_back(tempSymb);
	}
	
	ss.str("");
	ss << numString.length();
	tempSymb.symbol = ss.str();
	tempSymb.noNT = true;
	optSymbols.push_back(tempSymb);   
}
