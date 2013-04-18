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
//ExpressionTree.cpp

#include "ExpressionTree.h"
#include "TerminalSymbol.h"
#include <map>
#include <sstream>


///
/// returns expression tree created from the postfix stack passed to it
///
void ExpressionTree::convertPostFix(vector<TerminalSymbol*> & postFixStack){

	vector<ElementNode> prefixStack;
	tree<ElementNode> argTree;
	vector<TerminalSymbol*> newStack;

	compressOperator(postFixStack, newStack);
	
	vector<tree<ElementNode> > stack;
	tree<ElementNode>::iterator top;
	deque<float> args;
	
	TerminalSymbol* currElement;
	int numArgs;
	
	for(unsigned int i=0; i<newStack.size(); i++){
		currElement = newStack[i];
		numArgs = currElement->getNumArgs();
		tree<ElementNode> currTree;
		
		top = currTree.begin();
		ElementNode currNode;
		currNode.el = currElement;
		currTree.insert(top, currNode);    
		// when need to add elements get them from stack
		if(numArgs != 0){
			// when a variable operator pop off top element from stack
			// and set the number of arguments to get using that value
			if(numArgs < 0){
				argTree = stack.back(); // get last element
				stack.pop_back(); // remove last element
				top = argTree.begin(); // variable argument amount will be at top of tree
				numArgs = int(top->el->evaluate(args));
			}
			
			vector<tree<ElementNode> > tempStack;
			
			// now pop off number of trees from stack and add to current tree
			// as children -- add in reverse order to match the network
			for(int currArg=0; currArg < numArgs; currArg++){
				tempStack.push_back(stack.back());
				stack.pop_back();
			}
			
			// now put in correct order onto original tree
			for(int currArg=0; currArg < numArgs; currArg++){
				argTree = tempStack.back();
				tempStack.pop_back();
				currTree.append_child(currTree.begin(), argTree.begin());
			}
			
		}
		
		stack.push_back(currTree);
	}
	
	// stack should now contain a single tree that is an expression
	// tree of the network 
	expressTree = stack.back();
}


///
/// Compresses the operator calculations for generating 
///
void ExpressionTree::compressOperator(vector<TerminalSymbol*> & postFixStack,
	vector<TerminalSymbol*>& newStack){
	
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
			// start at bottom of stack so that you get any constants that need to be
			// carried over for later evaluation
			for(iter = stack.begin(); iter != stack.end(); ++iter){
				newConstant = new Constant(*iter);
				constants.push_back(newConstant);
				newStack.push_back(newConstant);
			}
			stack.clear();
			
			newStack.push_back(postFixStack[i]);
		}
	}
	
}



///
/// Adjust label for genotypes based on map file contents
///
string ExpressionTree::alterLabel(data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, string label, bool continMapUsed){
			
	if(mapUsed && label[0] == 'G'){
		 stringstream ss(label.substr(1,label.length()-1));
		 int num;
		 ss >> num;
		 if(ottDummy)
				num = (num-1)/2;
			else
				num -= 1;
			label = holder->getGenoName(num);
	}
	else if(continMapUsed && label[0] == 'C'){
		stringstream ss(label.substr(1,label.length()-1));
		int num;
		ss >> num;
		num -= 1;
		label = holder->getCovarName(num);
	}        
	else{
		if(label[0]=='G'){
	  stringstream ss(label.substr(1,label.length()-1));
			int num;
			ss >> num;
			if(ottDummy)
				num = (num-1)/2;
			else
				num -= 1;
			label = "G" + holder->getGenoName(num);
		 }
	}
		 
	return label;
}


///
/// Clears constants created during the construction of the tree
///
void ExpressionTree::clearConstants(){
	for(unsigned int i=0; i < constants.size(); i++){
		delete constants[i];
	}
	constants.clear();
}


///
/// Free any left over memory
///
ExpressionTree::~ExpressionTree(){
	clearConstants();
}

///
/// Returns maximum depth of the tree.  For neural networks the depth of the nodes will
/// be one less than the maximum depth of the tree.
/// @return maximum depth of the tree
///
unsigned int ExpressionTree::getMaxDepth(){
		return incrementDepth(expressTree.begin(), 0);
}

///
/// Recursively traverse tree and record deepest depth of a neural network node
/// @param iter
/// @param currdepth
/// @return maximum depth found
///
unsigned int ExpressionTree::incrementDepth(tree<ElementNode>::iterator baseIter, unsigned int currDepth){
		tree<ElementNode>::iterator childIter;
		
		unsigned int maxDepth=0, depth;
		
		for(int child=0; child < int(expressTree.number_of_children(baseIter)); child++){  
				childIter = expressTree.child(baseIter, child);  
				depth=incrementDepth(childIter, currDepth);
				if(depth > maxDepth){
						maxDepth=depth;
				}
		}
		if(baseIter->el->getType()[0] == 'P')
			maxDepth+=1;
			
		return maxDepth;
}


///
/// Output expression tree in dot language for use by Graphviz to
/// create an image file
///
void ExpressionTree::outputDot(ostream & out, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed){

	tree<ElementNode>::iterator iter;
	// need counter to set the title for the node in the graph
	map<string, int> typeCount;
	string type;

	for(iter=expressTree.begin(); iter != expressTree.end(); iter++){
		type = iter->el->getType();
		if(typeCount.find(type) == typeCount.end())
			typeCount[type] = 0;
		typeCount[type]++;
 
		stringstream ss;
		ss << typeCount[type];
		string number;
		ss >> number;
		
		iter->id = iter->el->getType() + number;
	}
 
	out << "digraph G{\n";
	out << "\tgraph [ dpi = 300 ];\n";
	out << "\tsize=\"7.5,11.0\";\n";
	out << "\tdir=\"none\";\n";
	out << "\trankdir=\"LR\";\n";
	out << "\torientation=\"landscape\";\n";
	
	// each node needs to point to its parent
	tree<ElementNode>::iterator parent;
	iter = expressTree.begin();
 
	string label = iter->el->getLabel();

	label = alterLabel(holder, mapUsed, ottDummy, label, continMapUsed);

	out << "\t" <<  iter->id << " [shape=\"" << iter->el->getShape() << "\" style=\"" << 
		iter->el->getStyle() << "\" label=\"" << label << "\"];" << endl;
	iter++;
	
	for(; iter != expressTree.end(); iter++){
		parent = expressTree.parent(iter);
		out << "\t" << iter->id << "->" << parent->id << ";" << endl;
 
		string label = iter->el->getLabel();
		label = alterLabel(holder, mapUsed, ottDummy, label, continMapUsed);
 
		out << "\t" << iter->id << " [shape=\"" << iter->el->getShape() << "\" style=\"" << 
			iter->el->getStyle() << "\" label=\"" << label << "\"];" << endl;  
	}
	out << "}" << endl; 
}
