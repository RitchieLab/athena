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
#include "NNSolutionCreator.h"
#include "TerminalSymbol.h"
#include "Terminals.h"
#include "ExpressionTree.h"
#include "BackPropAnnie.h"
#include "BackPropSMD.h"

// #include <iostream>

///
/// Constructor
///
NNSolutionCreator::NNSolutionCreator(){
    initialize();
}

///
/// Alternative constructor
///
NNSolutionCreator::NNSolutionCreator(vector<string>& symbols){
  initialize();
}


///
/// Initializes variables
///
void NNSolutionCreator::initialize(){
  terminals_set = false;
  calculator = NULL;
  graphic_extension = ".dot";
  startopt = "<cop>";
  
  string symbols[] = {"<cop>","<op>","<Concat>","(",")", "<num>","Concat","1","2","3","4",
    "5","6","7","8","9","0","<dig>","+","-","*","/"};
  
  // set up list of symbols included in setting weights in neural network
  optsymbols.insert(symbols, symbols+22);
  
  // set up symbols that take an argument  
  optargsymbols.insert("(<num>)");
  
  left_opt_bound = '(';
  right_opt_bound = ')';
  
//   nnlog = NULL;
}


///
/// Creates neural network for evaluation
/// sets up a postfix stack for evaluation of the neural network
/// @param pheno Phenotype from libGE that can be turned into a neural network
/// @param set Dataset
///
void NNSolutionCreator::establish_solution(vector<string>& symbols, Dataset* set){
    if(!terminals_set)
        set_variables(set);
    establish_solution(symbols);
}


///
/// Creates neural network for evaluation
/// sets up a postfix stack for evaluation of the neural network
/// @param symbols Phenotype from libGE that can be turned into a neural network
///
void NNSolutionCreator::establish_solution(vector<string>& symbols){
   
   nn_depth = 0;
   nn_terminal_size= symbols.size();
   
    // as part of stack look for variables and keep a list of those
    // to use in evaluation for skipping missing dataPq
    
    // Anything except parentheses with a priority of zero is pushed to the postfix stack
    // Left paren always pushed on stack
    // When get to right paren, pop off stack until left paren encountered
    // if the symbol being scanned has a higher precedence than the symbol at the
    //     top of the stack, the symbol scanned is pushed on the stack
    // If the precedence of the symbol being scanned is lower than or equal to the
    //     precedence of the symbol at the top of the stack, the stack is popped to the output
    // When the terminating symbol is reached on the input scan, the stack is popped to the
    //     output until the terminating symbol is also reached on the stack. Then the algorithm terminates.
   
    vector<TerminalSymbol *> temp_postfix;
   
    postfix_stack.clear();
    
    vector<TerminalSymbol *> symbol_stack;
    
    vector<TerminalSymbol *> stack;
    
    deque<float> args, blank;
    TerminalSymbol* constant_terminal_ptr;

    unsigned int input_size = symbols.size();
    for(unsigned int input_count=0; input_count < input_size; input_count++){
        // priority is zero when element is a number, variable, comma or parentheses
        // may need a dynamic cast
        TerminalSymbol* curr_term = term_holder.get_term(symbols[input_count]);
        
        // add check for concatenation operator here and replace with constant if necessary
        if(curr_term == term_holder.concaten()){
          args.clear();        
          input_count+=2; //skip next one which will be a '('
          while(term_holder.get_term(symbols[input_count]) != term_holder.right_paren()){
            args.push_back(term_holder.get_term(symbols[input_count++])->evaluate(blank));
          }
          // have to pop off last one which is variable value
          args.pop_back();
          constant_terminal_ptr = new Constant(term_holder.concaten()->evaluate(args));
          temp_postfix.push_back(constant_terminal_ptr);
          constants.push_back(constant_terminal_ptr);
          continue; // go to next string
        }
        
        if(!curr_term->get_priority()){
         if(curr_term == term_holder.left_paren()){
            stack.push_back(curr_term);
         }
         else if(curr_term == term_holder.right_paren()){
          while(stack.back() != term_holder.left_paren()){
            temp_postfix.push_back(stack.back());            
            stack.pop_back();
          }
          // pop off left paren
          stack.pop_back();
        }
        else if(curr_term == term_holder.comma()){
          while(stack.back() != term_holder.left_paren()){
            temp_postfix.push_back(stack.back());
            stack.pop_back();
          }
          // in this case do not pop off left paren as
          // it is being used to mark off arguments
          // to function calls
        }
        else{
          //  numbers and variables are pushed to output
          temp_postfix.push_back(curr_term);
        }
      }
      else{ // element is an operator
        while(stack.size() > 0 && 
              (curr_term->get_priority() <= stack.back()->get_priority())){
            temp_postfix.push_back(stack.back());
            stack.pop_back();
          }
        stack.push_back(curr_term);
      }
  }
 
 
  // after last terminal push anything on the stack to the postfix stack
  for(vector<TerminalSymbol *>::iterator stack_iter=stack.begin(); stack_iter != stack.end(); stack_iter++)
    temp_postfix.push_back(*stack_iter);
    
    
  // now convert any constants and store in postfix_stack
  compress_operator(temp_postfix, postfix_stack);
    
    
  // iterate through stack and store covars and genotypes for use in checking
  // for missing data
  genos.clear();
  covars.clear();
 
  for(unsigned int i=0; i<postfix_stack.size(); i++){

      if(postfix_stack[i]->get_term_type() == TerminalSymbol::Genotype){
         genos[postfix_stack[i]]++;
      }
      else if (postfix_stack[i]->get_term_type() == TerminalSymbol::Covariate){
         covars[postfix_stack[i]]++;
      }
  }
}


///
/// optimize solution by running back propagation
/// @param symbols vector of strings the can be converted into network
/// @param set Dataset for optimizing solution
/// @returns number of epochs trained
///
int NNSolutionCreator::optimizeSolution(std::vector<std::string>& symbols, Dataset* set){

  BackPropAnnie bp;
  // clear old symbols
  opt_val_symbols.clear();

  // this creates the postfix_stack
  establish_solution(symbols, set);
 
  int numEpochsTrained = bp.runBackProp(postfix_stack, set);
  
  opt_values = bp.getWeights();
  optimized_score = bp.getOptimizedScore();
  
  vector<float>::iterator iter;
  
  symbVector new_symbols;
  for(iter = opt_values.begin(); iter != opt_values.end(); ++iter){
    term_holder.terminalsFromConstant(*iter, new_symbols);
    opt_val_symbols.push_back(new_symbols);
  }
  return numEpochsTrained;
}


///
/// Compresses the operator calculations for generating 
///
void NNSolutionCreator::compress_operator(vector<TerminalSymbol*> & postfix_stack,
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
      // when not an operator or constant take the argument from the stack and
      // convert it to a Constant and then put on the new_stack followed by the
      // original terminal
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
/// Returns vector with indexes of SNPs in the latest solution
/// @return vector of indexes
///
vector<int> NNSolutionCreator::getGeneIndexes(){
  vector<int> indexes;
  map<TerminalSymbol*, int>::iterator iter;
  GenotypeTerm* g;
  int i;
  for(iter = genos.begin(); iter != genos.end(); iter++){
    g = dynamic_cast<GenotypeTerm*>(iter->first);
    for(i=0; i<iter->second; i++){
      indexes.push_back(g->getIndex());
    }
  }
  return indexes;
}

///
/// Returns vector with indexes of Covariates in the latest solution
/// @return vector of indexes
///
inline vector<int> NNSolutionCreator::getCovarIndexes(){
  vector<int> indexes;
  map<TerminalSymbol*, int>::iterator iter;
  ContinVariable* c;
  int i;
  for(iter = covars.begin(); iter != covars.end(); iter++){
    c = dynamic_cast<ContinVariable*>(iter->first);
    for(i=0; i<iter->second; i++){
      indexes.push_back(c->getIndex());
    }
  }  
  return indexes;
}


///
/// Creates variable terminals in terminal holder.
/// @param
///
void NNSolutionCreator::set_variables(Dataset* set){
    unsigned int num_covars = (*set)[0]->num_covariates();
    unsigned int num_genos = (*set)[0]->num_genotypes();
    term_holder.create_terminals(num_genos, num_covars);
    terminals_set = true;
}


///
/// Evaluates single individual and returns value
/// @param ind Individual to evaluate
/// returns value
///
float NNSolutionCreator::evaluate_ind(Individual* ind){

    term_holder.set_ind(ind);
    
    vector<float> stack;
    deque<float> args;
    
    int num_args=0, k;
    unsigned int postfix_size = postfix_stack.size();
    
    for(unsigned int j=0; j<postfix_size; j++){
        args.clear();
        num_args = postfix_stack[j]->get_num_args();
        // evaluate and push any elements that don't take arguments
        if(!(num_args=postfix_stack[j]->get_num_args())){
            stack.push_back(postfix_stack[j]->evaluate(args));
        }
         // need to determine number of arguments
         // pop off arguments and then push back new result
         // when number of args for the element is < 0 
         // top number on stack should be the number of 
         // arguments for the function
        else{
            if(num_args<0){
                num_args = (unsigned int)stack.back();
                stack.pop_back();
            }
            // construct deque of arguments
            for(k=0; k<num_args; k++){
                args.push_front(stack.back());
                stack.pop_back();
            }
            // evaluate current operator and push onto stack
            stack.push_back(postfix_stack[j]->evaluate(args));
        }
    }
 
    // only value on stack should be final value
    return stack.back();
}


///
/// Checks to see if the individual has any missing data
/// @param ind
/// @param set
/// @return true if the individual should be used
///
bool NNSolutionCreator::use_ind(Individual* ind, Dataset* set){
    deque<float> args;
    
    map<TerminalSymbol*, int>::iterator iter;
 
    for(iter=covars.begin(); iter != covars.end(); iter++){
      if(iter->first->evaluate(args) == set->get_missing_covalue()){
        return false;
      }
    }
    
    for(iter=genos.begin(); iter != genos.end(); iter++){
      if(iter->first->evaluate(args) == set->get_missing_genotype()){
        return false;
      }
    }

    return true;
}

// #include <iostream>
// using namespace std;

///
/// Evaluates the neural network and returns the balanced accuracy
/// of the neural network.
/// @param set Dataset to use
///
float NNSolutionCreator::evaluate(Dataset* set){
   
    // to check for missing data use the list of variables gotten from the
    // creation of the stack and then check each ind to  make sure
    // there is a value in each
    Individual * ind;
    nIndsEvaluated = 0;

    calculator->reset();

    // when missing skip that ind
    for(unsigned int i=0; i < set->num_inds(); i++){
        ind = (*set)[i];
    
        ContinVariable::set_ind(ind);
        GenotypeTerm::set_ind(ind);
      
        if(!use_ind(ind, set))
            continue;

        nIndsEvaluated++;
        
        calculator->add_ind_score(evaluate_ind(ind), ind->get_status());
    }
    return calculator->get_score();
    
}

///
/// For Networks, the detailed logging determines the depth in nodes of the
/// model.
///
void NNSolutionCreator::detailed_logging(){
  ExpressionTree extree;
  extree.convert_postfix(postfix_stack);
  
  nn_depth = extree.get_max_depth();
}

///
/// Return the network depth
///
unsigned int NNSolutionCreator::get_detailed_log(){
    return nn_depth;
}

///
/// writes a dot compatible text file representing the network
/// @param os ostream to write to
///
void NNSolutionCreator::graphical_output(ostream& os, data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy, bool continmap_used){
  
  // network will be part of this tree
  ExpressionTree extree;
  
  extree.convert_postfix(postfix_stack);
  extree.output_dot(os,holder, map_used, ott_dummy, continmap_used);
  extree.clear_constants();
}


/// 
/// Performs evaluation on the dataset.  Each individual has their result passed
/// to the outptu stream provided.
/// @param set Dataset 
/// @param os ostream 
///
float NNSolutionCreator::evaluate_with_output(Dataset* set, ostream& os){

    // to check for missing data use the list of variables gotten from the
    // creation of the stack and then check each ind to  make sure
    // there is a value in each
    Individual * ind;

    calculator->reset();
    
    // when missing skip that ind
    for(unsigned int i=0; i < set->num_inds(); i++){
        ind = (*set)[i];
               
        ContinVariable::set_ind(ind);
        GenotypeTerm::set_ind(ind);       
        if(!use_ind(ind, set)){
            os << ind->get_id() << "\tMissing data\t" << ind->get_status() << endl;
            continue;
        }
        float score = evaluate_ind(ind);
        os << ind->get_id() << "\t" << score << "\t" << ind->get_status() << endl;
        calculator->add_ind_score(score, ind->get_status());
    }  
 
    return calculator->get_score();
}
