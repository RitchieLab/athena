
#include "Terminals.h"

#include "GEGrammarSIAltTerm.h"

///
/// Replaces the Terminal symbols with the TerminalSymbol 
///
void GEGrammarSIAltTerm::replace_terminals(){
    
    unsigned int num_rules = (*this).size();
    unsigned int curr_prod, num_productions, curr_symbol, num_symbols;
    Symbol* old_symbol;
    
    for(unsigned int curr_rule=0;  curr_rule < num_rules; curr_rule++){
        num_productions = (*this)[curr_rule].size();
        for(curr_prod = 0; curr_prod < num_productions; curr_prod++){
            num_symbols = (*this)[curr_rule][curr_prod].size();
            for(curr_symbol=0; curr_symbol < num_symbols; curr_symbol++){
                old_symbol = (*this)[curr_rule][curr_prod][curr_symbol];
               if(old_symbol->getType() == TSymbol){             
                  (*this)[curr_rule][curr_prod][curr_symbol] = 
                     term_creator.create_terminal(*old_symbol);
                  delete old_symbol;
               }
            }
        }
   }
}
