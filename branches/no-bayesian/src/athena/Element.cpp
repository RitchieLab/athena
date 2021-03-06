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
//Element.cpp

#include "Element.h"

///
/// Constructor for types that contain a value
/// @param symbol string containing symbol used for the element
/// @param terminalStatus bool that indicates whether this is a terminal element
/// @param priorityLevel int specifying how to prioritize evaluation for this operator
///
Element::Element(string symbol, bool terminalStatus, int priorityLevel){
	name = symbol;
	isTerminal = terminalStatus;
	varArgs = false;
	priority = priorityLevel;
	label = "nonterminal";
	shape = "box";
	style = "bold";
	type = label;
}



///
/// Constructor for types that contain a value
/// @param value 
/// @param terminal_status bool that indicates whether this is a terminal element
/// @param priority_level int specifying how to priority evaluation for this operator
///
Element::Element(DATATYPE value, bool terminalStatus, int priorityLevel){
	stringstream ss;
	ss << value;
	name = ss.str();
	isterminal = terminalStatus;
	var_args = false;
	priority = priorityLevel;
	label = "nonterminal";
	shape = "box";
	style = "bold";
	type = label;
}


///
/// Tracks whether Element is a terminal
/// If not a terminal then element must expand in the
/// grammar
/// @return terminal status
///
bool Element::getIsTerminal() const{
	return isTerminal;
}


///
/// Returns number of elements needed
/// for this element
/// @return number of elements
///
int Element::getNumArgs() const{
	return numArgs;
}


///
/// Sets the number of arguments
/// @param nargs number of arguments needed by the argument
///
void Element::set_num_args(int nArgs){
	numArgs = nArgs;
}


///
/// Returns the symbol of the element
/// @return symbol 
///
string Element::getName() const{
	return name;
}


///
/// Returns priority level for the element
/// @returns priority
///
int Element::getPriority() const{
	return priority;
}


///
/// Results are scaled from -1 to 1
/// @param x DATATYPE with result
/// @return scaled result
///
DATATYPE Element::ActivateSigmoid(DATATYPE x)
{
				if(x < -709) return -1.0;
				if(x > 709)  return 1.0;
				return(1.0 / (1.0 + exp(-x)));
}


///
/// Adjusts the result if it is infinite 
/// or not a number
/// @param x DATATYPE representing result
/// @return adjusted result 
///
DATATYPE Element::AdjustResult(DATATYPE x)
{
				if(isinf(x) == 1) {
								return 1.0;
				}
				if(isinf(x) == -1) {
								return -1.0;
				}
				if(isnan(x)) {
								return 0.0;
				}

				return x;
}


///
/// Overloads output for element
/// @param os ostream to write to
/// @param el Element to output
/// @returns output stream
///
ostream & operator << (ostream & os, const Element & el){
				os << el.name;
				return os;
}

