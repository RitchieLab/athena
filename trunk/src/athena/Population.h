/* 
 * File:   Population.h
 * Author: dudeksm
 *
 * Created on November 10, 2008, 1:26 PM
 */
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

#ifndef _POPULATION_H
#define	_POPULATION_H

#include "rbtree.h"
#include "Solution.h"

// has static variable that allows for sorting either in descending (0) or ascending order (2)
struct floatLT{
   static int order_adjust;
   int operator()(const float l, const float r) const {
      if (l<r) return 1 - order_adjust;
      if (l>r) return -1 + order_adjust;
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
    
    Population(){floatLT::order_adjust=0; currnode = NULL;convertScores=false;}
    
    ~Population();
    
    /// returns number of solutions
    unsigned int num_solutions(){return solutions.GetCount();}
    
    /// sets number of solutions
    void set_pop_size(int pop_size){solutions.SetMaxSize(pop_size);}
    
    /// gets number of solutions
    int get_pop_size(){return solutions.GetCount();}
    
    /// returns best solution
    Solution* best(){currnode = solutions.GetFirst(); return solutions.GetFirst()->GetData();}
    
    Solution* operator[](int index){
      SolutionTreeNode* current_node;
      int count=0;
      for(current_node = solutions.GetFirst(); count < index;
        current_node = current_node->GetNext()){       
        count++;
      }
      return current_node->GetData();
    }
    
    /// returns NULL when no more to get and then resets the iterator to first position
    Solution* GetNext(){
      if(currnode == solutions.GetLast()){
        currnode = solutions.GetFirst();
        return NULL;
      }
      else{
        currnode = currnode->GetNext();
        Solution * ret = currnode->GetData();
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
    void insert(Solution* sol){solutions.Insert(sol->fitness(), sol);currnode=solutions.GetFirst();}
        
    /// outputs all solutions in tree to  check order
    void output_tree();
    
    /// sets order to be ascending
    void sort_ascending(){floatLT::order_adjust=2;}

   /// if true scores stored should be converted for output(mean square to r-squared)
   void setConvertScores(bool val){convertScores=val;}
   bool getConvertScores(){return convertScores;} 
   
   /// converts scores of population
   void convert_scores(Dataset* train, Dataset* test);
    
   /// convert scores of population when no testing set used
   void convert_scores(Dataset* train);
    
private:
    
    SolutionTree solutions;
    SolutionTreeNode* currnode;
    bool convertScores;
};

#endif	/* _POPULATION_H */

