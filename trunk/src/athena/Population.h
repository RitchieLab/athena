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
 * File:   Population.h
 * Author: dudeksm
 *
 * Created on November 10, 2008, 1:26 PM
 */

#ifndef _POPULATION_H
#define	_POPULATION_H

#include "rbtree.h"
#include "Solution.h"

// has static variable that allows for sorting either in descending (0) or ascending order (2)
struct floatLT{
	 static int orderAdjust;
	 int operator()(const float l, const float r) const {
			if (l<r) return 1 - orderAdjust;
			if (l>r) return -1 + orderAdjust;
			return 0;
	 }
};

typedef Utility::RBTree<float, Solution*, floatLT> SolutionTree;
typedef Utility::RBTreeNode<float, Solution*, floatLT> SolutionTreeNode;


///
/// Abstract base class for holding solutions (neural networks, svm, etc.)
///
class Population{
		
public:
		
		Population(){floatLT::orderAdjust=0; currNode = NULL;convScores=false;}
		
		~Population();
		
		/// returns number of solutions
		unsigned int numSolutions(){return solutions.GetCount();}
		
		/// sets number of solutions
		void setPopSize(int popSize){solutions.SetMaxSize(popSize);}
		
		/// gets number of solutions
		int getPopSize(){return solutions.GetCount();}
		
		/// returns best solution
		Solution* best(){currNode = solutions.GetFirst(); return solutions.GetFirst()->GetData();}
		
		Solution* operator[](int index){
			SolutionTreeNode* currentNode;
			int count=0;
			for(currentNode = solutions.GetFirst(); count < index;
				currentNode = currentNode->GetNext()){       
				count++;
			}
			return currentNode->GetData();
		}
		
		/// returns NULL when no more to get and then resets the iterator to first position
		Solution* GetNext(){
			if(currNode == solutions.GetLast()){
				currNode = solutions.GetFirst();
				return NULL;
			}
			else{
				currNode = currNode->GetNext();
				Solution * ret = currNode->GetData();
				return ret;
			}
		}
		
		/// assignment operator
		Population& operator=(const Population& other);
		
		/// copy constructor
		Population(const Population& p);
		
		/// copy constructor
		void copy(const Population& other);
		
		/// clears the population
		void clear();
		
		/// Add a solution to the tree
		void insert(Solution* sol){solutions.Insert(sol->fitness(), sol);currNode=solutions.GetFirst();}
				
		/// outputs all solutions in tree to  check order
		void outputTree();
		
		/// sets order to be ascending
		void sortAscending(){floatLT::orderAdjust=2;}

	 /// if true scores stored should be converted for output(mean square to r-squared)
	 void setConvertScores(bool val){convScores=val;}
	 bool getConvertScores(){return convScores;} 
	 
	 /// converts scores of population
	 void convertScores(Dataset* train, Dataset* test);
		
	 /// convert scores of population when no testing set used
	 void convertScores(Dataset* train);
		
private:
		
		SolutionTree solutions;
		SolutionTreeNode* currNode;
		bool convScores;
};

#endif	/* _POPULATION_H */

