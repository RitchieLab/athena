#include "Terminals.h"
#include <Stringmanip.h>
#include <sstream>

#include <iostream>

using namespace data_manage;


//////////////////////////////////////
// Constant class
//////////////////////////////////////

Constant::Constant(string symbol, int number_args)
  :TerminalSymbol(symbol, 0, TerminalSymbol::Constant){
    num_args = number_args;
    // must convert symbol into value
    stringstream ss(symbol);
    ss >> constant_value;
    label = symbol;
    shape = "box";
    style = "bold";
    type = "const";
}
  

Constant::Constant(float value)
  :TerminalSymbol("", 0, TerminalSymbol::Constant){
  num_args = 0;
  constant_value = value;
  stringstream ss;
  ss << value;
  ss >> label;
  name = label;
  priority = 0;
  shape = "box";
  style = "bold";   
  type = "const";
}

float Constant::evaluate(deque<float> & args){ 
  return constant_value;
} 

//////////////////////////////////////
// Genotype class
//////////////////////////////////////

Individual* GenotypeTerm::ind = NULL;


GenotypeTerm::GenotypeTerm(string gname, int var_index)
  :TerminalSymbol(gname, 0, TerminalSymbol::Genotype){
    num_args = 0;
    name = gname;
    // must convert symbol into value
    index_value = var_index-1;
    label = gname;
    style = "filled";
    shape = "box";
    type = label;
}

///
/// Returns value for the individual at the indicated location
/// @param args deque not used in this evaluation
/// @return value of the variable
/// @throws AthenaExcept when missing variable data at the indicated
/// locus for the indicated individual
///
float GenotypeTerm::evaluate(deque<float> & args){ 
  return ind->get_genotype(index_value);
} 


///
/// Set Individual pointer for class
/// @param i Individual*
///
void GenotypeTerm::set_ind(Individual* i){
    ind = i;
}

//////////////////////////////////////
// Continuous variable (covariate) class
//////////////////////////////////////

Individual* ContinVariable::ind = NULL;

ContinVariable::ContinVariable(string gname, int var_index)
  :TerminalSymbol(gname, 0, TerminalSymbol::Covariate){
    num_args = 0;
    name = gname;
    // must convert symbol into value
    label = gname;
    index_value = var_index-1;
    style = "filled";
    shape = "box";
    type = label;
    
}

///
/// Returns value for the individual at the indicated location
/// @param args deque not used in this evaluation
/// @return value of the variable
/// @throws AthenaExcept when missing variable data at the indicated
/// locus for the indicated individual
///
float ContinVariable::evaluate(deque<float> & args){ 
  return ind->get_covariate(index_value);
} 


///
/// Set Individual pointer for class
/// @param i Individual*
///
void ContinVariable::set_ind(Individual* i){
    ind = i;
}


//////////////////////////////////////
// Addition class
//////////////////////////////////////

Addition::Addition(string symbol, int number_args)
  :TerminalSymbol(symbol, 2, TerminalSymbol::Operator){
    num_args = number_args;
    shape = "diamond";
    style = "bold";
    label = "+";
    type = "Add";
}
    
float Addition::evaluate(deque<float> & args){
  float result = args[0];
  
  for(int i=1; i<num_args; i++){
    result += args[i];
  }
  return AdjustResult(result);
} 

//////////////////////////////////////
// Subtraction class
//////////////////////////////////////

Subtraction::Subtraction(string symbol, int number_args)
  :TerminalSymbol(symbol, 2, TerminalSymbol::Operator){
    num_args = number_args;
    shape = "diamond";
    style = "bold";
    label = "-";
    type = "Sub";
}
    
float Subtraction::evaluate(deque<float> & args){
  float result = args[0];
  
  for(int i=1; i<num_args; i++){
    result -= args[i];
  }
  return AdjustResult(result);
} 

//////////////////////////////////////
// Multiplication class
//////////////////////////////////////

Multiplication::Multiplication(string symbol, int number_args)
  :TerminalSymbol(symbol, 3, TerminalSymbol::Operator){
    num_args = number_args;
    shape = "diamond";
    style = "bold";
    label = "*";
    type = "Mult";
}
    
float Multiplication::evaluate(deque<float> & args){
  float result = args[0];
  
  for(int i=1; i<num_args; i++){
    if(args[i] == 0){
      result = 0;
      break;
    } 
    result *= args[i];
  }
  return AdjustResult(result);
} 

//////////////////////////////////////
// Division class
//////////////////////////////////////
Division::Division(string symbol, int number_args)
  :TerminalSymbol(symbol, 3, TerminalSymbol::Operator){
    num_args = number_args;
    shape = "diamond";
    style = "bold";
    label = "/";  
    type = "Div";
}
    
float Division::evaluate(deque<float> & args){
  float result = args[0];
  
  for(int i=1; i<num_args; i++){
    if(args[i] == 0)
      return 1.0;
    result /= args[i];
  }
  
  return AdjustResult(result);
} 


//////////////////////////////////////
// Power class
//////////////////////////////////////

Power::Power(string symbol, int number_args)
  :TerminalSymbol(symbol, 2, TerminalSymbol::Operator){
    num_args = number_args;
    shape = "diamond";
    style = "bold";
    label = "^";
    type = "Pow";
}
    
float Power::evaluate(deque<float> & args){
  return AdjustResult(pow(args[0], args[1]));
} 

//////////////////////////////////////
// pAdd class
//////////////////////////////////////
pAdd::pAdd(string symbol, int number_args)
  :TerminalSymbol(symbol, 4, TerminalSymbol::Neuron){
    num_args = number_args;
    var_args = true;
    shape = "doublecircle";
    style = "bold";
    label = "PADD";    
    type = label;
}

float pAdd::evaluate(deque<float> & args){
  float result = args[0];
  for(unsigned int i=1; i<args.size(); i++){
    result += args[i];
  }
  
  result = AdjustResult(result);
 
  return ActivateSigmoid(result);
} 

//////////////////////////////////////
// pSub class
//////////////////////////////////////
pSub::pSub(string symbol, int number_args)
  :TerminalSymbol(symbol, 4, TerminalSymbol::Neuron){
    num_args = number_args;
    var_args = true;
    shape = "doublecircle";
    style = "bold";
    label = "PSUB";    
    type = label;
}

float pSub::evaluate(deque<float> & args){
  float result = args[0]; 
  for(unsigned int i=1; i<args.size(); i++){
    result -= args[i];
  }
  result = AdjustResult(result);
  
  return ActivateSigmoid(result);
}

//////////////////////////////////////
// pMult class
//////////////////////////////////////
pMult::pMult(string symbol, int number_args)
  :TerminalSymbol(symbol, 4, TerminalSymbol::Neuron){
    num_args = number_args;
    var_args = true;
    shape = "doublecircle";
    style = "bold";
    label = "PMULT";
    type = label;
}

float pMult::evaluate(deque<float> & args){
  float result = args[0];
  
  for(unsigned int i=1; i<args.size(); i++){
    result *= args[i];
  }
  
  result = AdjustResult(result);
  return ActivateSigmoid(result);
}

//////////////////////////////////////
// pDiv class
//////////////////////////////////////
pDiv::pDiv(string symbol, int number_args)
  :TerminalSymbol(symbol, 4, TerminalSymbol::Neuron){
    num_args =number_args;
    var_args = true;
    shape = "doublecircle";
    style = "bold";
    label = "PDIV";    
    type = label;
}

float pDiv::evaluate(deque<float> & args){
  float result = args[0];
  
  for(unsigned int i=1; i<args.size(); i++){
    if(args[i] == 0){
//      result = 1.1;
      result = 10000000;
      break;
    }
    result /= args[i];
  }
  
  result = AdjustResult(result);
  
  return ActivateSigmoid(result);
}


//////////////////////////////////////
// pAnd class
//////////////////////////////////////

pAnd::pAnd(string symbol, int number_args)
  :TerminalSymbol(symbol, 4, TerminalSymbol::Neuron){
    num_args = number_args;
    var_args = true;
    shape = "doublecircle";
    style = "bold";
    label = "PAND";    
    type = label;
}

float pAnd::evaluate(deque<float> & args){  
  for(unsigned int i=0; i<args.size(); i++){
    if(args[i] <= 0){
      return 0.0;
    }
  }
  return 1.0;
  
}

//////////////////////////////////////
// pNand class
//////////////////////////////////////

pNand::pNand(string symbol, int number_args)
  :TerminalSymbol(symbol, 4, TerminalSymbol::Neuron){
    num_args = number_args;
    var_args = true;
    shape = "doublecircle";
    style = "bold";
    label = "PNAND";    
    type = label;
}

float pNand::evaluate(deque<float> & args){  
  for(unsigned int i=0; i<args.size(); i++){
    if(args[i] > 0){
      return 0.0;
    }
  }
  return 1.0;
}

//////////////////////////////////////
// pOr class
//////////////////////////////////////

pOr::pOr(string symbol, int number_args)
  :TerminalSymbol(symbol, 4, TerminalSymbol::Neuron){
    num_args = number_args;
    var_args = true;
    shape = "doublecircle";
    style = "bold";
    label = "POR";    
    type = label;
}

float pOr::evaluate(deque<float> & args){  
  for(unsigned int i=0; i<args.size(); i++){
    if(args[i] > 0){
      return 1.0;
    }
  }
  return 0.0;
}


//////////////////////////////////////
// pNor class
//////////////////////////////////////

pNor::pNor(string symbol, int number_args)
  :TerminalSymbol(symbol, 4, TerminalSymbol::Neuron){
    num_args = number_args;
    var_args = true;
    shape = "doublecircle";
    style = "bold";
    label = "PNOR";    
    type = label;
}

float pNor::evaluate(deque<float> & args){  
  for(unsigned int i=0; i<args.size(); i++){
    if(args[i] > 0){
      return 0.0;
    }
  }
  return 1.0;
}


//////////////////////////////////////
// pXor class
//////////////////////////////////////

pXor::pXor(string symbol, int number_args)
  :TerminalSymbol(symbol, 4, TerminalSymbol::Neuron){
    num_args = number_args;
    var_args = true;
    shape = "doublecircle";
    style = "bold";
    label = "PXOR";   
    type = label;
}

float pXor::evaluate(deque<float> & args){ 
  int value = (int)(args[0] > 0)?1:0;
  int nextval;
  
  for(unsigned int i=1; i<args.size(); i++){
    nextval = (int)(args[i] > 0)?1:0;
    if(nextval != value){
      value = 1;
    }
    else{
      value = 0;
    }
  }
  return (float) value;  
}


//////////////////////////////////////
// Weight class
//////////////////////////////////////
Weight::Weight(string symbol, int number_args)
  :TerminalSymbol(symbol, 4, TerminalSymbol::Weight){
    num_args = number_args;
    var_args = true;
    shape = "circle";
    style = "bold";
    label = "W";
    type = label;
}

float Weight::evaluate(deque<float> & args){
    float result = args[0];
  
  for(int i=1; i<num_args; i++){
    if(args[i] == 0){
      result = 0;
      break;
    } 
    result *= args[i];
  }
  
  return AdjustResult(result);
}

//////////////////////////////////////
// ConCat class
//////////////////////////////////////
ConCat::ConCat(string symbol, int number_args)
  :TerminalSymbol(symbol, 4){
    num_args = number_args;
    var_args = true;
}

// will be doubles for this implementation
// so add value of arg to 48
// '.' is 46 so that class has a value of -2
float ConCat::evaluate(deque<float> & args){
  float result;
  
  string newnumber;
  
  // check for negative value 
  if(args[0] < 0 && int(args[0]) != -2){
    newnumber += '-';
    // convert negative to positive
    args[0] = fabs(args[0]);
  }
  
  for(unsigned int i=0; i<args.size(); i++)
    newnumber += char(48+args[i]);

  // convert to float
  stringstream in(newnumber);
  in >> result;
  
  return result;
}


//////////////////////////////////////
// LogF class for log base 10 calculation
//////////////////////////////////////
LogF::LogF(string symbol, int number_args)
    :TerminalSymbol(symbol, 4, TerminalSymbol::PreOperator){
    num_args = number_args;
    shape = "diamond";
    style = "bold";
    label = "log";
    type = "log";    
  }
  
float LogF::evaluate(deque<float> & args){
    
  if(args[0] < 0){
    return 0;
  }
  else{
    return AdjustResult(log10(args[0]));
  }
}   


//////////////////////////////////////
// Sine class
//////////////////////////////////////
Sine::Sine(string symbol, int number_args)
    :TerminalSymbol(symbol, 4, TerminalSymbol::PreOperator){
    num_args = number_args;
    shape = "diamond";
    style = "bold";
    label = "sine";
    type = "sine";    
  }
  
float Sine::evaluate(deque<float> & args){  
  return AdjustResult(sin(args[0]));
} 


//////////////////////////////////////
// Sine class
//////////////////////////////////////
Cosine::Cosine(string symbol, int number_args)
    :TerminalSymbol(symbol, 4, TerminalSymbol::PreOperator){
    num_args = number_args;
    shape = "diamond";
    style = "bold";
    label = "cosine";
    type = "cosine";    
  }
  
float Cosine::evaluate(deque<float> & args){  
  return AdjustResult(cos(args[0]));
} 

//////////////////////////////////////
// Tangent class
//////////////////////////////////////
Tangent::Tangent(string symbol, int number_args)
    :TerminalSymbol(symbol, 4, TerminalSymbol::PreOperator){
    num_args = number_args;
    shape = "diamond";
    style = "bold";
    label = "cosine";
    type = "cosine";    
  }
  
float Tangent::evaluate(deque<float> & args){  
  return AdjustResult(tan(args[0]));
} 



//////////////////////////////////////
// Dot class
//////////////////////////////////////
Dot::Dot(string symbol, int number_args)
  :TerminalSymbol(symbol, 0){
    num_args = number_args;
    var_args = true;
}

// -2 is the offset to be used with the 
// ConCat operator to generate a '.'
float Dot::evaluate(deque<float> & args){
    return -2.0;
}

//////////////////////////////////////
// And class
//////////////////////////////////////

And::And(string symbol, int number_args)
  :TerminalSymbol(symbol, 2, TerminalSymbol::Operator){
    num_args = number_args;
    var_args = true;
    shape = "diamond";
    style = "bold";
    label = "and";  
    type = label;
}

float And::evaluate(deque<float> & args){  
  for(unsigned int i=0; i<args.size(); i++){
    if(args[i] <= 0){
      return 0.0;
    }
  }
  return 1.0;
  
}

//////////////////////////////////////
// Nand class
//////////////////////////////////////

Nand::Nand(string symbol, int number_args)
  :TerminalSymbol(symbol, 2, TerminalSymbol::Operator){
    num_args = number_args;
    var_args = true;
    shape = "diamond";
    style = "bold";
    label = "nand";  
    type = label;
}

float Nand::evaluate(deque<float> & args){  
  for(unsigned int i=0; i<args.size(); i++){
    if(args[i] > 0){
      return 0.0;
    }
  }
  return 1.0;
}

//////////////////////////////////////
// Or class
//////////////////////////////////////

Or::Or(string symbol, int number_args)
  :TerminalSymbol(symbol, 2,TerminalSymbol::Operator){
    num_args = number_args;
    var_args = true;
    shape = "diamond";
    style = "bold";
    label = "or";  
    type = label;
}

float Or::evaluate(deque<float> & args){  
  for(unsigned int i=0; i<args.size(); i++){
    if(args[i] > 0){
      return 1.0;
    }
  }
  return 0.0;
}


//////////////////////////////////////
// Nor class
//////////////////////////////////////

Nor::Nor(string symbol, int number_args)
  :TerminalSymbol(symbol, 2, TerminalSymbol::Operator){
    num_args = number_args;
    var_args = true;
    shape = "diamond";
    style = "bold";
    label = "nor";  
    type = label;
}

float Nor::evaluate(deque<float> & args){  
  for(unsigned int i=0; i<args.size(); i++){
    if(args[i] > 0){
      return 0.0;
    }
  }
  return 1.0;
}


//////////////////////////////////////
// Xor class
//////////////////////////////////////

Xor::Xor(string symbol, int number_args)
  :TerminalSymbol(symbol, 2, TerminalSymbol::Operator){
    num_args = number_args;
    var_args = true;
    shape = "diamond";
    style = "bold";
    label = "xor";  
    type = label;
}

float Xor::evaluate(deque<float> & args){ 
  int value = (int)(args[0] > 0)?1:0;
  int nextval;
  
  for(unsigned int i=1; i<args.size(); i++){
    nextval = (int)(args[i] > 0)?1:0;
    if(nextval != value){
      value = 1;
    }
    else{
      value = 0;
    }
  }
  return (float) value;  
}

