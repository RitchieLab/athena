#include "SymRegressSolutionCreator.h"
#include "TerminalSymbol.h"
#include "Terminals.h"

///
/// sets up a postfix stack for evaluation of the symbolic regression equation
/// @param pheno Phenotype from libGE that can be turned into a symbolic regression equation
/// @param set Dataset
///
void SymRegressSolutionCreator::establish_solution(vector<string>& symbols, Dataset* set){
    if(!terminals_set)
        set_variables(set);
    establish_solution(symbols);
}

///
/// Creates symbolic regression equation for evaluation
/// sets up a postfix stack for evaluation of the symbolic regression equation
/// @param symbols Phenotype from libGE that can be turned into a symbolic regression equation
///
void SymRegressSolutionCreator::establish_solution(vector<string>& symbols){
   
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
   
    postfix_stack.clear();
    
    vector<TerminalSymbol *> symbol_stack;
    
    vector<TerminalSymbol *> stack;
    
    deque<float> args, blank;
    TerminalSymbol* constant_terminal_ptr;

    unsigned int input_size = symbols.size();
    for(unsigned int input_count=0; input_count < input_size; input_count++){
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
          postfix_stack.push_back(constant_terminal_ptr);
          constants.push_back(constant_terminal_ptr);
          continue; // go to next string
        }
        
        if(!curr_term->get_priority()){
         if(curr_term == term_holder.left_paren()){
            stack.push_back(curr_term);
         }
         else if(curr_term == term_holder.right_paren()){
          while(stack.back() != term_holder.left_paren()){
            postfix_stack.push_back(stack.back());            
            stack.pop_back();
          }
          // pop off left paren
          stack.pop_back();
        }
        else if(curr_term == term_holder.comma()){
          while(stack.back() != term_holder.left_paren()){
            postfix_stack.push_back(stack.back());
            stack.pop_back();
          }
          // in this case do not pop off left paren as
          // it is being used to mark off arguments
          // to function calls
        }
        else{
          //  numbers and variables are pushed to output
          postfix_stack.push_back(curr_term);
        }
      }
      else{ // element is an operator
        while(stack.size() > 0 && 
              (curr_term->get_priority() <= stack.back()->get_priority())){
            postfix_stack.push_back(stack.back());
            stack.pop_back();
          }
        stack.push_back(curr_term);
      }
  }
 
 
  // after last terminal push anything on the stack to the postfix stack
  for(vector<TerminalSymbol *>::iterator stack_iter=stack.begin(); stack_iter != stack.end(); stack_iter++)
    postfix_stack.push_back(*stack_iter);

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
/// evaluates single individual and returns value for that individual
///
float SymRegressSolutionCreator::evaluate_ind(Individual* ind){
  return TerminalSymbol::ActivateSigmoid(NNSolutionCreator::evaluate_ind(ind));
}