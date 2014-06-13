/*
Copyright Marylyn Ritchie 2014

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
 * File:   DAGraph.h
 * Author: dudeksm
 *
 * Created on March 5, 2014, 11:10 AM
 */

#ifndef _DAGRAPH_H
#define	_DAGRAPH_H

#include <set>
#include <map>
#include "Terminals.h"
#include <Dataholder.h>

struct GraphNode{
	std::set<GraphNode*> children, parents;
	TerminalSymbol* term;
	int genomeIndex;
};

typedef std::set<GraphNode*>::iterator GraphNodeIter;

class DAGraph{

public:
	DAGraph();
	~DAGraph();
	
	/// connects graphnode from with node containing to
	void addConnection(TerminalSymbol* from, TerminalSymbol* to);
	
	/// adds new GraphNode to graph
	GraphNode* addNode(TerminalSymbol* term);
	
	/// output tree in dot language
	void outputDot(ostream & out, data_manage::Dataholder* holder,
		bool mapUsed, bool ottDummy, bool continMapUsed);
	
	void clearNodes();
	
	/// return number of nodes
	unsigned int numNodes(){return nodes.size();}
	
	GraphNodeIter begin(){return nodes.begin();}
	GraphNodeIter end(){return nodes.end();}
	
private:

	std::string alterLabel(data_manage::Dataholder* holder,
		bool mapUsed, bool ottDummy, std::string label);

	std::map<TerminalSymbol*, GraphNode*> termMap;
	std::set<GraphNode*> nodes;
	
};


#endif
