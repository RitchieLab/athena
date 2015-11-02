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
//ExpressionTree.h

// converts postfix stack
// into an expression tree

#ifndef __EXPRESSIONTREE_H__
#define __EXPRESSIONTREE_H__

#include "Tree.hh"
#include "Terminals.h"
#include <iostream>
#include <Dataholder.h>

class ExpressionTree{

	public:
	
		ExpressionTree();
		~ExpressionTree();
	
		/// returns prefix stack after conversion from the post fix stack
		void convertPostFix(vector<TerminalSymbol*> & postFixStack);
		
		/// output tree in dot language
		void outputDot(ostream & out, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed);

		/// clears constants
		void clearConstants();
		
		struct ElementNode{
			TerminalSymbol* el;
			string id;
		};
		
		unsigned int getMaxDepth();
		
		/// returns iterator to beginning of tree
		inline tree<ElementNode>::iterator begin(){return expressTree.begin();}
		
		/// returns iterator to end of tree
		inline tree<ElementNode>::iterator end(){return expressTree.end();}
		
		/// returns indexed child iterator 
		inline tree<ElementNode>::iterator child(tree<ElementNode>::iterator& iter,
			int childIndex){return expressTree.child(iter, childIndex);}
		
		/// returns number of children for current iterator
		inline int numberOfChildren(tree<ElementNode>::iterator& iter){
			return expressTree.number_of_children(iter);}
		
		void outputEquation(ostream& out, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapused);
			
		std::string getEquation();
		
	private:
		unsigned int incrementDepth(tree<ElementNode>::iterator baseIter, 
				unsigned int currDepth);
	
		void setOperators(string& label, string& op, string& prefix);
	
		void outputEquationTree(ostream& out, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapused, tree<ElementNode>::iterator baseIter);
			
		void getEquationTree(ostream& out, tree<ElementNode>::iterator baseIter);	
	
		std::string alterLabel(data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, std::string label, bool continMapUsed, bool equationOut);
			
		void compressOperator(std::vector<TerminalSymbol*> & postFixStack,
			std::vector<TerminalSymbol*>& newStack);
	
		std::map<std::string, std::string> operatorMap, prefixMap;
	
		tree<ElementNode> expressTree;
		tree<ElementNode>::iterator extreeIter;
		
		vector<TerminalSymbol*> constants;
		
};

typedef tree<ExpressionTree::ElementNode>::iterator ExpressTreeIter;


#endif
