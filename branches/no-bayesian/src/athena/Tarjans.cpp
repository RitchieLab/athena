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

#include "Tarjans.h"
#define NIL -1

using namespace std;

Tarjans::Tarjans(int V)
{
    this->V = V;
    adj = new list<int>[V];
}

Tarjans::~Tarjans(){
	delete [] adj;
}

void Tarjans::addEdge(int v, int w)
{
    adj[v].push_back(w);
}

///
/// Remove connection v->w
///
void Tarjans::removeEdge(int v, int w){
	adj[v].remove(w);
}

// A recursive function that finds and prints strongly connected
// components using DFS traversal
// u --> The vertex to be visited next
// disc[] --> Stores discovery times of visited vertices
// low[] -- >> earliest visited vertex (the vertex with minimum
//             discovery time) that can be reached from subtree
//             rooted with current vertex
// *st -- >> To store all the connected ancestors (could be part
//           of SCC)
// stackMember[] --> bit/index array for faster check whether
//                  a node is in stack
void Tarjans::SCCUtil(int u, int disc[], int low[], stack<int> *st,
                    bool stackMember[])
{
    // A static variable is used for simplicity, we can avoid use
    // of static variable by passing a pointer.
    static int time = 0;

    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
    st->push(u);
    stackMember[u] = true;

    // Go through all vertices adjacent to this
    list<int>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i)
    {
        int v = *i;  // v is current adjacent of 'u'

        // If v is not visited yet, then recur for it
        if (disc[v] == -1)
        {
            SCCUtil(v, disc, low, st, stackMember);

            // Check if the subtree rooted with 'v' has a
            // connection to one of the ancestors of 'u'
            // Case 1 (per above discussion on Disc and Low value)
            low[u]  = min(low[u], low[v]);
        }

        // Update low value of 'u' only of 'v' is still in stack
        // (i.e. it's a back edge, not cross edge).
        // Case 2 (per above discussion on Disc and Low value)
        else if (stackMember[v] == true)
            low[u]  = min(low[u], disc[v]);
    }

    // head node found, pop the stack and print an SCC
    int w = 0;  // To store stack extracted vertices
    if (low[u] == disc[u])
    {
    		vector<int> loopHolder;
        while (st->top() != u)
        {
            w = (int) st->top();
//             cout << w << " ";
            loopHolder.push_back(w);
            stackMember[w] = false;
            st->pop();
        }
        w = (int) st->top();
        // store any loops that occur
        if(!loopHolder.empty()){
        	loopHolder.push_back(w);
        	loops.push_back(loopHolder);
        }
//        cout << w << "\n";
        stackMember[w] = false;
        st->pop();
    }
}

// The function to do DFS traversal. It uses SCCUtil()
bool Tarjans::SCC()
{
		loops.clear();
    int *disc = new int[V];
    int *low = new int[V];
    bool *stackMember = new bool[V];
    stack<int> *st = new stack<int>();

    // Initialize disc and low, and stackMember arrays
    for (int i = 0; i < V; i++)
    {
        disc[i] = NIL;
        low[i] = NIL;
        stackMember[i] = false;
    }

    // Call the recursive helper function to find strongly
    // connected components in DFS tree with vertex 'i'
    for (int i = 0; i < V; i++)
        if (disc[i] == NIL)
            SCCUtil(i, disc, low, st, stackMember);

		// clean up
    delete st;
    delete [] disc;
    delete [] low;
    delete [] stackMember;

     	if(loops.empty()){
//  		cout << "NO LOOPS" << endl;
 		return false;
 	}
 	else{
//   cout << loops.back().size() << endl;
 		return true;
 	}

}