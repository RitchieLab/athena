//ExpressionTree.cpp

#include "ExpressionTree.h"
#include "TerminalSymbol.h"
#include <map>
#include <sstream>


///
/// returns expression tree created from the postfix stack passed to it
///
void ExpressionTree::convert_postfix(vector<TerminalSymbol*> & postfix_stack){

  vector<Element_node> prefix_stack;
  tree<Element_node> arg_tree;
  vector<TerminalSymbol*> new_stack;

  compress_operator(postfix_stack, new_stack);
  
  vector<tree<Element_node> > stack;
  tree<Element_node>::iterator top;
  deque<float> args;
  
  TerminalSymbol* currElement;
  int num_args;
  
  for(unsigned int i=0; i<new_stack.size(); i++){
    currElement = new_stack[i];
    num_args = currElement->get_num_args();
    tree<Element_node> curr_tree;
    
    top = curr_tree.begin();
    Element_node curr_node;
    curr_node.el = currElement;
    curr_tree.insert(top, curr_node);    
    // when need to add elements get them from stack
    if(num_args != 0){
      // when a variable operator pop off top element from stack
      // and set the number of arguments to get using that value
      if(num_args < 0){
        arg_tree = stack.back(); // get last element
        stack.pop_back(); // remove last element
        top = arg_tree.begin(); // variable argument amount will be at top of tree
        num_args = int(top->el->evaluate(args));
      }
      
      vector<tree<Element_node> > temp_stack;
      
      // now pop off number of trees from stack and add to current tree
      // as children -- add in reverse order to match the network
      for(int curr_arg=0; curr_arg < num_args; curr_arg++){
        temp_stack.push_back(stack.back());
        stack.pop_back();
      }      
      
      // now put in correct order onto original tree
      for(int curr_arg=0; curr_arg < num_args; curr_arg++){
        arg_tree = temp_stack.back();
        temp_stack.pop_back();
        curr_tree.append_child(curr_tree.begin(), arg_tree.begin());
      }
      
    }
    
    stack.push_back(curr_tree);
  }
  
  // stack should now contain a single tree that is an expression
  // tree of the network 
  express_tree = stack.back();
}


///
/// Compresses the operator calculations for generating 
///
void ExpressionTree::compress_operator(vector<TerminalSymbol*> & postfix_stack,
  vector<TerminalSymbol*>& new_stack){
  
  // for any operator compresses them so that stack will not have redundant information
  // any operator can be compressed into a constant value
  
  // 1.  for any non-constant / non-operator push on to new stack
  // 2.  when find constant evaluate until find non-constant and then take that value from
  // stack and create a new constant for the new stack
  
  TerminalSymbol* newConstant;
  
  vector<float> stack;
  deque<float> args;
  int num_args;
  vector<float>::iterator iter;
  
  for(unsigned int i=0; i < postfix_stack.size(); i++){
    
    TerminalSymbol::TerminalType termType = postfix_stack[i]->get_term_type();
    
    if(termType == TerminalSymbol::Operator){
      // when it is an operator evaluate and push back on to stack
      num_args = postfix_stack[i]->get_num_args();
      for(int k=0; k<num_args; k++){
        args.push_front(stack.back());
        stack.pop_back();
      }
     stack.push_back(postfix_stack[i]->evaluate(args));
    }
    else if(termType == TerminalSymbol::Constant){
      args.clear();
      stack.push_back(postfix_stack[i]->evaluate(args));
    }
    else{
      // start at bottom of stack so that you get any constants that need to be
      // carried over for later evaluation
      for(iter = stack.begin(); iter != stack.end(); ++iter){
        newConstant = new Constant(*iter);
        constants.push_back(newConstant);
        new_stack.push_back(newConstant);
      }
      stack.clear();
      
      new_stack.push_back(postfix_stack[i]);
    }
  }
  
}


///
/// Adjust label for genotypes based on map file contents
///
string ExpressionTree::alter_label(data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy, string label){
      
  if(map_used && label[0] == 'G'){
     stringstream ss(label.substr(1,label.length()-1));
     int num;
     ss >> num;
     if(ott_dummy)
        num = (num-1)/2;
      else
        num -= 1;
      label = holder->get_geno_name(num);
  }        
   
  return label;
}

///
/// Clears constants created during the construction of the tree
///
void ExpressionTree::clear_constants(){
  for(unsigned int i=0; i < constants.size(); i++){
    delete constants[i];
  }
}

///
/// Output expression tree in dot language for use by Graphviz to
/// create an image file
///
void ExpressionTree::output_dot(ostream & out, data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy){

  tree<Element_node>::iterator iter;
  // need counter to set the title for the node in the graph
  map<string, int> type_count;
  string type;

  for(iter=express_tree.begin(); iter != express_tree.end(); iter++){
    type = iter->el->get_type();
    if(type_count.find(type) == type_count.end())
      type_count[type] = 0;
    type_count[type]++;
 
    stringstream ss;
    ss << type_count[type];
    string number;
    ss >> number;
    
    iter->id = iter->el->get_type() + number;
  }
 
  out << "digraph G{\n";
  out << "\tsize=\"7.5,11.0\";\n";
  out << "\tdir=\"none\";\n";
  out << "\trankdir=\"LR\";\n";
  out << "\torientation=\"landscape\";\n";
  
  // each node needs to point to its parent
  tree<Element_node>::iterator parent;
  iter = express_tree.begin();
 
  string label = iter->el->get_label();

  label = alter_label(holder, map_used, ott_dummy, label);

  out << "\t" <<  iter->id << " [shape=\"" << iter->el->get_shape() << "\" style=\"" << 
    iter->el->get_style() << "\" label=\"" << label << "\"];" << endl;
  iter++;
  
  for(; iter != express_tree.end(); iter++){
    parent = express_tree.parent(iter);
    out << "\t" << iter->id << "->" << parent->id << ";" << endl;
 
    string label = iter->el->get_label();
    label = alter_label(holder, map_used, ott_dummy, label);
 
    out << "\t" << iter->id << " [shape=\"" << iter->el->get_shape() << "\" style=\"" << 
      iter->el->get_style() << "\" label=\"" << label << "\"];" << endl;  
  }
  out << "}" << endl; 
}
