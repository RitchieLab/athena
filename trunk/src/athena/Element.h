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
//Element.h

///
/// This is the base class for all elements of the grammar
/// Any element that is a nonterminal uses this base class
/// Also sytactical elements such as commas can use this class
/// Elements that perform operations or return values should inherit
/// from this class (see Terminals.h)
///


#ifndef __ELEMENT_H__
#define __ELEMENT_H__

#include<string>
#include "Stringmanip.h"
#include<vector>
#include<math.h>
#include<deque>
#include "defines.h"

#include<iostream> // remove later

using namespace std;

class Element{
	
	public:
 
		/// constructor for nonterminal
		Element(string symbol, bool terminalStatus = false, int priorityLevel=0);
		/// constructor for terminal
		Element(DATATYPE value, bool terminalStatus, int priorityLevel=0);
		virtual ~Element(){};
		
		/// Returns true for is a terminal
		bool getIsTerminal() const;

		/// Returns number of arguments needed by this element
		int getNumArgs() const;

		/// Sets the number of arguments
		void setNumArgs(int nArgs);

		/// Returns symbol showing which element this is
		string getName() const;

		/// Returns priority level for evaluating this element
		int getPriority() const;
		
		/// scale result using sigmoid function
		static DATATYPE ActivateSigmoid(DATATYPE x);    
		/// adjust results for infinite or nan
		static DATATYPE AdjustResult(DATATYPE x);
		
		/// nonterminals return 0 for evaluation
		virtual DATATYPE evaluate(deque<DATATYPE> & elements){return 0;};

		friend ostream & operator << (ostream & os, const Element & el);
		
		/// returns label for dot file production
		virtual string getLabel(){return label;}
		
		/// returns style for dot file production
		string getStyle(){return style;}
		
		/// returns shape for dot file production
		string getShape(){return shape;}
		
		/// returns type (used for dot)
		string getType(){return type;}
	
	protected:
		 string name, label, style, shape, type;
		 bool isTerminal, varArgs;
		 int numArgs, priority;

};


#endif   
