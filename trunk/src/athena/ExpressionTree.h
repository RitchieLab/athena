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
  
  	~ExpressionTree();
  
    /// returns prefix stack after conversion from the post fix stack
    void convert_postfix(vector<TerminalSymbol*> & postfix_stack);
    
    /// output tree in dot language
    void output_dot(ostream & out, data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy, bool continmap_used);

    /// clears constants
    void clear_constants();
    
    struct Element_node{
      TerminalSymbol* el;
      string id;
    };
    
    unsigned int get_max_depth();
    
    /// returns iterator to beginning of tree
    inline tree<Element_node>::iterator begin(){return express_tree.begin();}
    
    /// returns iterator to end of tree
    inline tree<Element_node>::iterator end(){return express_tree.end();}
    
    /// returns indexed child iterator 
    inline tree<Element_node>::iterator child(tree<Element_node>::iterator& iter,
      int childIndex){return express_tree.child(iter, childIndex);}
    
    /// returns number of children for current iterator
    inline int number_of_children(tree<Element_node>::iterator& iter){
      return express_tree.number_of_children(iter);}
      
  private:
    unsigned int increment_depth(tree<Element_node>::iterator baseIter, 
        unsigned int currdepth);
  
    std::string alter_label(data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy, std::string label, bool continmap_used);
      
    void compress_operator(std::vector<TerminalSymbol*> & postfix_stack,
      std::vector<TerminalSymbol*>& new_stack);
  
    tree<Element_node> express_tree;
    tree<Element_node>::iterator extree_iter;
    
    vector<TerminalSymbol*> constants;
    
};

typedef tree<ExpressionTree::Element_node>::iterator ExpressTreeIter;


#endif
