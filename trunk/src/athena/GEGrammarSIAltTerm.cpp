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
