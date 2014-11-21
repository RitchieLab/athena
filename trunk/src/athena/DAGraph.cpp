#include "DAGraph.h"
#include <sstream>
using namespace std;
#include <iostream>
#define NIL -1

#include <stdlib.h> 

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
// cout << "add connection " << from->getName() << " to " << to->getName() << endl;	
	// check that node isn't same
	if(fromNode == toNode){
// cout <<  "skip connection" << endl;
		return;
	}
	
	fromNode->children.insert(toNode);
	toNode->parents.insert(fromNode);
	
	//don't allow nodes to point to each other -- now have repair operator so can
	// allow this situation
// 	if(fromNode->parents.find(toNode) == fromNode->parents.end()){
// 		fromNode->children.insert(toNode);
// 		toNode->parents.insert(fromNode);
// 	}

}
	
///
/// Removes an edge from the graph
/// @param from node 
/// @param to node
///
void DAGraph::removeConnection(TerminalSymbol* from, TerminalSymbol* to){
	GraphNode* fromNode, *toNode;
	fromNode = termMap.find(from)->second;
	toNode = termMap.find(to)->second;

// if(fromNode->children.find(toNode) == fromNode->children.end()){
// cout << "can't find child " << toNode->term->getName() <<  endl;
// }
// if(toNode->parents.find(fromNode) == toNode->parents.end()){
// cout << "can't find parent " << fromNode->term->getName() <<  endl;
// exit(1);
// }
// 	
	fromNode->children.erase(toNode);
	toNode->parents.erase(fromNode);
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
	else if(label[0]=='C'){
		stringstream ss(label.substr(1,label.length()-1));
		int num;
		ss >> num;
		num--;
		newLabel = holder->getCovarName(num);
		if(!mapUsed)
			newLabel = "C" + newLabel;
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


void DAGraph::outputRstring(ostream & out){
	for(GraphNodeIter iter = begin(); iter != nodes.end(); ++iter){
		out <<  "[" << (*iter)->term->getName();
		if(!(*iter)->parents.empty()){
			out << "|";
			std::set<GraphNode*>::iterator parentIter=(*iter)->parents.begin();
			out << (*parentIter)->term->getName();
			++parentIter;
			for(;parentIter != (*iter)->parents.end(); ++parentIter){
				out << ":" << (*parentIter)->term->getName();
			}
		}
		out << "]";
	}	
}


///
/// Tarjan's strongly connected components -- adapted from http://www.geeksforgeeks.org/tarjan-algorithm-find-strongly-connected-components/
/// @returns true if any components strongly connected (cycles)
///
bool DAGraph::SCC(){
	loops.clear();
	
// temporary show all parents and all children of the nodes
// for(GraphNodeIter iter = begin(); iter != nodes.end(); ++iter){
// cout << (*iter)->term->getName() << " children:";
// for(GraphNodeIter cIter = (*iter)->children.begin(); cIter!=(*iter)->children.end(); ++cIter){
// cout << " " << (*cIter)->term->getName();
// }
// cout << " parents:";
// for(GraphNodeIter pIter = (*iter)->parents.begin(); pIter!=(*iter)->parents.end(); ++pIter){
// cout << " " << (*pIter)->term->getName();
// }
// cout << "\n";
// }
	
	
	size_t nNodes = nodes.size();
// cout << "in DAGraph SCC " << nodes.size() << endl;	
	vector<GraphNode*> nodeVec(nNodes,NULL);
	map<GraphNode*, int> nodeMap;
	size_t i=0;
	for(GraphNodeIter iter = begin(); iter != nodes.end(); ++iter){
		nodeVec[i] = *iter;
		nodeMap[*iter] = i;
// cout << "graphnode " << (*iter)->term->getName() << " has index=" << i << endl;
		i++;
	}	
	
// temporary show all parents and all children of the nodes as indexes
// for(GraphNodeIter iter = begin(); iter != nodes.end(); ++iter){
// cout << (*iter)->term->getName() << " children:";
// for(GraphNodeIter cIter = (*iter)->children.begin(); cIter!=(*iter)->children.end(); ++cIter){
// cout << " " << (*cIter)->term->getName() << "-" << nodeMap[*cIter];
// }
// cout << " parents:";
// for(GraphNodeIter pIter = (*iter)->parents.begin(); pIter!=(*iter)->parents.end(); ++pIter){
// cout << " " << (*pIter)->term->getName() << "-" << nodeMap[*pIter];
// }
// cout << "\n";
// }	
	
	
	vector<int> disc(nNodes,NIL);
	vector<int> low(nNodes,NIL);
	vector<bool> stackMember(nNodes, false);
	stack<int> *st = new stack<int>();
	int time = 0;
	for(i=0; i<nNodes; i++){
		if(disc[i] == NIL){
			SCCUtil(i, disc, low, st, stackMember, nodeVec, &time, nodeMap);
		}
	}
 	delete st;
 	
 	// report any loops
//  cout << "in DAGraph SCC loops.size=" << loops.size() << endl;
 	if(loops.empty()){
//  		cout << "NO LOOPS" << endl;
 		return false;
 	}
 	else{
//   cout << loops.back().size() << endl;
 		return true;
 	}
 	
}

///
/// A recursive function that finds and prints strongly connected
/// components using DFS traversal
/// @param u The vertex to be visited next (as listed in nodeVec)
/// @param disc Stores discovery times of visited vertices
/// @param low Earliest visited vertex (the vertex with minimum
///             discovery time) that can be reached from subtree
///             rooted with current vertex
/// @param st To store all the connected ancestors (could be part
///           of SCC)
/// @param stackMember bit/index vector for faster check whether
///                  a node is in stack
/// @param nodeVec vector holding pointers to graphnode
/// @param time
///
void DAGraph::SCCUtil(int u, vector<int>& disc, vector<int>& low, stack<int>* st,
		vector<bool>& stackMember, vector<GraphNode*>& nodeVec, int* time,
		map<GraphNode*, int>& nodeMap){
		
		// Initialize discovery time and low value
		disc[u] = low[u] = ++(*time);
		st->push(u);
    stackMember[u] = true;

    // Go through all vertices adjacent to this
    // Go through all children of the current node
//     list<int>::iterator i;
    std::set<GraphNode*>::iterator childIter;
// cout << "parent is " << nodeVec[u]->term->getName() << " u=" << u << "\n";
    for(childIter = nodeVec[u]->children.begin(); childIter != nodeVec[u]->children.end();
    	++childIter){
    	int v = nodeMap[*childIter];
// cout << "child is " << nodeciutVec[v]->term->getName() << " v=" << v << " u(parent)=" << u << "\n";
    	if(disc[v] == -1){
    		SCCUtil(v,disc,low,st,stackMember,nodeVec,time,nodeMap);
    		// Check if the subtree rooted with 'v' has a connection
    		// to one of the ancestors of 'u'
    		low[u] = min(low[u], low[v]);
    	}
    	// Update low value of 'u' only if 'v' is still in stack
    	// (i.e. it's a back edge, not cross edge
    	else if(stackMember[v]==true){
    		low[u] = min(low[u], disc[v]);
    	}
    	
    }
    
//     for (i = adj[u].begin(); i != adj[u].end(); ++i)
//     {
//         int v = *i;  // v is current adjacent of 'u'
//  
//         // If v is not visited yet, then recur for it
//         if (disc[v] == -1)
//         {
//             SCCUtil(v, disc, low, st, stackMember);
//  
//             // Check if the subtree rooted with 'v' has a
//             // connection to one of the ancestors of 'u'
//             low[u]  = min(low[u], low[v]);
//         }
//  
//         // Update low value of 'u' only of 'v' is still in stack
//         // (i.e. it's a back edge, not cross edge).
//         else if (stackMember[v] == true)
//             low[u]  = min(low[u], disc[v]);
//     }
    
    // head node found, pop the stack and print an SCC
    int w = 0;  // To store stack extracted vertices
    if (low[u] == disc[u])
    {
    		vector<GraphNode*> loopHolder;
        while (st->top() != u)
        {
            w = (int) st->top();
//               cout << w << "-" <<  nodeVec[w]->term->getLabel() << " ";
            loopHolder.push_back(nodeVec[w]);
            stackMember[w] = false;
            st->pop();
        }
        w = (int) st->top();
//           cout << w <<  "-" <<  nodeVec[w]->term->getLabel()<< "\n";
        if(loopHolder.size()>0){
        	loopHolder.push_back(nodeVec[w]);
        	loops.push_back(loopHolder);
        }
        stackMember[w] = false;
        st->pop();
    }
		
}


