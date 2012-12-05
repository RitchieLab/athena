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
/* 
 * File:   GENNGrammerAdjuster.h
 * Author: dudeksm
 *
 * Created on December 5, 2008, 4:49 PM
 */

#ifndef _GENNGRAMMERADJUSTER_H
#define	_GENNGRAMMERADJUSTER_H


#include <GE/GEGrammarSI.h>
#include "AthenaGrammarSI.h"
#include <set>

///
/// Adjusts grammar for shorthand representation of continuous variables and 
/// for the use of dummy encoding for gentoypes
///

class GENNGrammarAdjuster{
    
public:

    /// Reads grammar file and applies any needed modifications to the grammar
    /// before proceeding
    void read_grammar_file(string filename);

    /// Doubles number of genotypes for ott dummy encoding
    void double_genotype_grammar();

    /// Includes all variables in set, overriding the grammar
    void include_all_vars(int nGenos, int nContin);
    
    /// passes the grammar string to the mapper
    void set_mapper(AthenaGrammarSI& mapper, int rank=0);
    
    ///  Expands shorthand for continuous variables and genotypes into useable format
    void expand_variables();
    
    /// Returns variables
    vector<string> get_variables();
   
    /// Clears variables held
    void clear_variables(){variables_included.clear();}
    
    /// Includes only variables in the variables_included set
    void add_variables(std::vector<std::string>& terminals);
    
    /// Creates new variable lines using only those in variables_included set
    void edit_only_var_included();
   
private:
    
    /// creates variable symbol for insertion in grammar
    string create_variable_symbol(string vPrefix, int num);
    
    vector<string>::iterator get_start_variables();
    
    int count_var_lines(vector<string>::iterator& start, vector<string>::iterator& last);
    
    /// contains the grammar
    vector<std::string> lines;
    std::set<std::string> variables_included;
};



#endif	/* _GENNGRAMMERADJUSTER_H */

