#include "DAGraph.h"
#include <sstream>
using namespace std;
#include <iostream>

DAGraph::DAGraph(){
}
	
///
/// Destructor deletes all the nodes
///
DAGraph::~DAGraph(){
	clearNodes();
}

///
/// Free memory from nodes
///
void DAGraph::clearNodes(){
	map<TerminalSymbol*, GraphNode*>::iterator iter;

	// delete node pointers
	for(iter = termMap.begin(); iter != termMap.end(); ++iter){
		delete iter->second;
	}
	// empty map
	termMap.clear();
	nodes.clear();
}


/// connects graphnode from with node containing to
void DAGraph::addConnection(TerminalSymbol* from, TerminalSymbol* to){

	GraphNode* fromNode, *toNode;

	map<TerminalSymbol*, GraphNode*>::iterator fromIter = termMap.find(from);
	if(fromIter == termMap.end()){
		fromNode = addNode(from);
	}
	else{
		fromNode = fromIter->second;
	}
	
	map<TerminalSymbol*, GraphNode*>::iterator toIter = termMap.find(to);
	if(toIter == termMap.end()){
		toNode = addNode(to);
	}
	else{
		toNode = toIter->second;
	}
	
	//don't allow nodes to point to each other
	if(fromNode->parents.find(toNode) == fromNode->parents.end()){
		fromNode->children.insert(toNode);
		toNode->parents.insert(fromNode);
	}

}
	
/// adds new GraphNode to graph
GraphNode* DAGraph::addNode(TerminalSymbol* term){
	GraphNode * newNode = new GraphNode;
	newNode->term = term;
	termMap[term]=newNode;
	nodes.insert(newNode);
	return newNode;
}


string DAGraph::alterLabel(data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, string label){
	
	string newLabel;
	
	if(label[0]=='G'){
		stringstream ss(label.substr(1,label.length()-1));
		int num, numOrig;
		ss >> numOrig;
		if(ottDummy)
			num = (numOrig-1)/2;
		else
			num = numOrig-1;
		newLabel = holder->getGenoName(num);
		if(!mapUsed)
			newLabel = "G" + newLabel;
	}
	else{
		newLabel = label;
	}
	
	return newLabel;
}


///
/// Output expression tree in dot language for use by Graphviz to
/// create an image file
///
void DAGraph::outputDot(ostream & out, data_manage::Dataholder* holder,
			bool mapUsed, bool ottDummy, bool continMapUsed){
			
	out << "digraph G{\n";
	out << "\tgraph [ dpi = 300 ];\n";
	out << "\tsize=\"8.5,11.0\";\n";
	out << "\tdir=\"none\";\n";
	out << "\trankdir=\"LR\";\n";
	out << "\torientation=\"portrait\";\n";
	
	// list all nodes and have parents point to children
	for(GraphNodeIter iter = begin(); iter != nodes.end(); ++iter){
		string label = alterLabel(holder, mapUsed, ottDummy, (*iter)->term->getLabel());
		out << "\t" <<  label << " [shape=\"" << (*iter)->term->getShape() << "\" style=\"" << 
			(*iter)->term->getStyle() << "\" label=\"" << label << "\"];" << endl;
		for(GraphNodeIter childIter = (*iter)->children.begin(); childIter != (*iter)->children.end();
			++childIter){
			string childLabel = alterLabel(holder, mapUsed, ottDummy, (*childIter)->term->getLabel());
			out << "\t" << label << "->" << childLabel << ";" << endl;
		}	
			
	}
	out << "}" << endl;
			
}



