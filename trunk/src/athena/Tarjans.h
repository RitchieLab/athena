/*
Copyright Marylyn Ritchie 2015

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

#ifndef _TARJANS_H
#define	_TARJANS_H

// #include <ga/GA2DBinStrGenome.h>
#include <vector>
#include <list>
#include <stack>

///
/// Provides objective function for use with GALIb and GE library
///
class Tarjans{

public:
    Tarjans(int V);   // Constructor
    ~Tarjans();
    void addEdge(int v, int w);   // function to add an edge to graph
    void removeEdge(int v, int w);
    bool SCC();    // prints strongly connected components

    std::vector<std::vector<int> > getLoops(){return loops;}

private:
    int V;    // No. of vertices
    std::list<int> *adj;    // A dynamic array of adjacency lists

    // A Recursive DFS based function used by SCC()
    void SCCUtil(int u, int disc[], int low[],
                 std::stack<int> *st, bool stackMember[]);

    std::vector<std::vector<int> > loops;

};

#endif