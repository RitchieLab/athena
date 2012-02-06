/* 
 * File:   TerminalSymbol.h
 * Author: dudeksm
 *
 * Created on November 14, 2008, 3:46 PM
 */

#ifndef _TERMINALSYMBOL_H
#define	_TERMINALSYMBOL_H

#include <math.h>
#include <deque>
#include <string>

using namespace std;

///
/// Base class for symbols used in neural networks
///

class TerminalSymbol{
    
public:
    
    enum TerminalType{
        Covariate,
        Genotype,
        Operator,
        Constant,
        NotVariable,
        Neuron,
        Weight,
        PreOperator
    };
    
    TerminalSymbol(){name = ""; priority = 0;}
    
    virtual ~TerminalSymbol(){}
    
    TerminalSymbol(std::string termname, int priority_level, TerminalType t_type=NotVariable)
      {name = termname; termtype = t_type; priority=priority_level;}
    
    /// Returns number of arguments needed by this element
    int get_num_args() const {return num_args;}
    
    /// Sets the number of arguments
    void set_num_args(int nargs){num_args = nargs;}
    
    /// Returns priority level for evaluating this element
    int get_priority() const {return priority;}
    
    /// Sets the priority level for evaluating the element
    void set_priority(int p){priority = p;}

    /// scale result using sigmoid function
    static float ActivateSigmoid(float x);    
    /// adjust results for infinite or nan
    static float AdjustResult(float x);
    
    /// nonterminals return 0 for evaluation
    virtual float evaluate(std::deque<float> & elements){return 0;}
       
    /// returns name of terminal symbol
    std::string get_name(){return name;}
    
    /// returns correct indicator for covariate or genotype or not
    TerminalType get_term_type(){return termtype;}
    
    /// returns label for dot file production
    virtual string get_label(){return label;}
    
    /// returns style for dot file production
    string get_style(){return style;}
    
    /// returns shape for dot file production
    string get_shape(){return shape;}
    
    /// returns type (used for dot)
    string get_type(){return type;}
    

    
  protected:
     std::string name, label, style, shape, type;
     bool var_args;
     int num_args, priority;    
     TerminalType termtype;
    
};



#endif	/* _TERMINALSYMBOL_H */


