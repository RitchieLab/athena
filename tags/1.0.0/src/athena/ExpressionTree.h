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
  
    /// returns prefix stack after conversion from the post fix stack
    void convert_postfix(vector<TerminalSymbol*> & postfix_stack);
    
    /// output tree in dot language
    void output_dot(ostream & out, data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy);

    /// clears constants
    void clear_constants();
    
    struct Element_node{
      TerminalSymbol* el;
      string id;
    };
    
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
  
    std::string alter_label(data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy, std::string label);
      
    void compress_operator(std::vector<TerminalSymbol*> & postfix_stack,
      std::vector<TerminalSymbol*>& new_stack);
  
    tree<Element_node> express_tree;
    tree<Element_node>::iterator extree_iter;
    
    vector<TerminalSymbol*> constants;
    
};

typedef tree<ExpressionTree::Element_node>::iterator ExpressTreeIter;


#endif
