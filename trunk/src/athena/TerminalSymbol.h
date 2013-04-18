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
/* 
 * File:   TerminalSymbol.h
 * Author: dudeksm
 *
 * Created on November 14, 2008, 3:46 PM
 */

#ifndef _TERMINALSYMBOL_H
#define	_TERMINALSYMBOL_H

#include <math.h>
#include <deque>
#include <string>

using namespace std;

///
/// Base class for symbols used in neural networks
///

class TerminalSymbol{
		
public:
		
		enum TerminalType{
				Covariate,
				Genotype,
				Operator,
				Constant,
				NotVariable,
				Neuron,
				Weight,
				PreOperator,
				Bias
		};
		
		TerminalSymbol(){name = ""; priority = 0;}
		
		virtual ~TerminalSymbol(){}
		
		TerminalSymbol(std::string termName, int priorityLevel, TerminalType tType=NotVariable)
			{name = termName; termType = tType; priority=priorityLevel;}
		
		/// Returns number of arguments needed by this element
		int getNumArgs() const {return numArgs;}
		
		/// Sets the number of arguments
		void setNumArgs(int nargs){numArgs = nargs;}
		
		/// Returns priority level for evaluating this element
		int getPriority() const {return priority;}
		
		/// Sets the priority level for evaluating the element
		void setPriority(int p){priority = p;}

		/// scale result using sigmoid function
		static float ActivateSigmoid(float x);    
		/// adjust results for infinite or nan
		static float AdjustResult(float x);
		
		/// nonterminals return 0 for evaluation
		virtual float evaluate(std::deque<float> & elements){return 0;}
			 
		/// returns name of terminal symbol
		std::string getName(){return name;}
		
		/// returns correct indicator for covariate or genotype or not
		TerminalType getTermType(){return termType;}
		
		/// returns label for dot file production
		virtual string getLabel(){return label;}
		
		/// returns style for dot file production
		string getStyle(){return style;}
		
		/// returns shape for dot file production
		string getShape(){return shape;}
		
		/// returns type (used for dot)
		string getType(){return type;}
		

		
	protected:
		 std::string name, label, style, shape, type;
		 bool varArgs;
		 int numArgs, priority;    
		 TerminalType termType;
		
};



#endif	/* _TERMINALSYMBOL_H */


