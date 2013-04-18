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
#include "Population.h"
#include "Terminals.h"

int floatLT::orderAdjust=0;

///
/// Assignment operator
///
Population& Population::operator=(const Population& other){
		if(this != &other){
			copy(other);   
		}
		return *this;
}


///
/// copy population
///
void Population::copy(const Population& other){
		SolutionTreeNode* currentNode;
		for(currentNode = other.solutions.GetFirst(); currentNode != other.solutions.GetLast();
			currentNode = currentNode->GetNext()){
				// need to create a clone operator  
				Solution* sol = currentNode->GetData()->clone();
				solutions.Insert(currentNode->GetKey(), sol);
		}
		currNode = solutions.GetFirst();
		convScores = other.convScores;
}



///
/// Outputs fitness of the solutions in tree in order
/// 
void Population::outputTree(){
		SolutionTreeNode* currentNode;
		for(currentNode = solutions.GetFirst(); currentNode != solutions.GetLast();
			currentNode = currentNode->GetNext()){
				// need to create a clone operator 
			 cout << "tree solution fitness=" << 
					 currentNode->GetData()->fitness() << endl;
		}
}


///
/// Copy constructor
///
Population::Population(const Population& p){
		copy(p);
}


///
/// Destructor frees any Solutions
///
Population::~Population(){
		clear();
}



///
/// Clears the tree and deletes the Solution pointers
///
void Population::clear(){
		
		if(solutions.GetCount() == 0)
				return;

		SolutionTreeNode* currentNode;
		for(currentNode = solutions.GetFirst(); currentNode != solutions.GetLast();
			currentNode = currentNode->GetNext()){
				delete currentNode->GetData();
		}
		delete currentNode->GetData();
		solutions.Clear();
}



///
/// converts scores of population
/// @param train Dataset
/// @param test Dataset
///
void Population::convertScores(Dataset* train, Dataset* test){

		SolutionTreeNode* node;
		for(node = solutions.GetFirst(); node != solutions.GetLast();
			node = node->GetNext()){
			node->GetData()->adjustScoreOut(train, test);
		}
}


///
/// converts scores of population
/// @param train Dataset
///
void Population::convertScores(Dataset* train){

		SolutionTreeNode* node;
		for(node = solutions.GetFirst(); node != solutions.GetLast();
			node = node->GetNext()){
			node->GetData()->adjustScoreOut(train);
		}

}

