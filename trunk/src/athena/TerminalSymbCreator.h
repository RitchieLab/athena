/* 
 * File:   TerminalSymbCreator.h
 * Author: dudeksm
 *
 * Created on November 19, 2008, 4:00 PM
 */

#ifndef _TERMINALSYMBCREATOR_H
#define	_TERMINALSYMBCREATOR_H

#include "TerminalSymbol.h"
#include <map>
#include <vector>
#include "HemannExcept.h"
#include <Individual.h>
#include "Structs.h"

using namespace data_manage;

class TerminalSymbCreator{
    
  public:
    TerminalSymbCreator();
    void create_terminals(int num_genotypes, int num_covariates);
    TerminalSymbol * create_constant(const std::string& symbol);
    
    void terminalsFromConstant(float value, symbVector& opt_symbols);
    
    void addGenotypeVariables(int num_variables);
    void addContinVariables(int num_variables);
    
    TerminalSymbol* get_term(std::string& symbol);
    
    inline TerminalSymbol* right_paren(){return rparen;}
    inline TerminalSymbol* left_paren(){return lparen;}
    inline TerminalSymbol* comma(){return commaptr;}
    inline TerminalSymbol* concaten(){return concat;}
    
    void set_ind(Individual* ind);
    
  private:
   
    void getStringFromNum(float value, symbVector& opt_symbols);
    std::map<std::string, TerminalSymbol*> TerminalMap;
    TerminalSymbol* rparen, *lparen, *commaptr, *concat;
    optSymbol left_paren_symb, right_paren_symb, concat_symb, period_symb;

};



#endif	/* _TERMINALSYMBCREATOR_H */

