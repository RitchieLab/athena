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
 * File:   GEGrammarSIAltTerm.h
 * Author: dudeksm
 *
 * Created on November 17, 2008, 4:16 PM
 */

#ifndef _GEGRAMMARSIALTTERM_H
#define	_GEGRAMMARSIALTTERM_H

#include <GE/ge.h>
#include "Terminals.h"

///
/// Alternative GEGrammarSI class that replaces all basic terminals
/// with the TerminalSymbol class for use in evaluating the solutions
/// produced by the algorithm
///

class GEGrammarSIAltTerm: public GEGrammarSI{
    
public:
    
    void replace_terminals();

private:
    
    TerminalSymbCreator term_creator;
    
};



#endif	/* _GEGRAMMARSIALTTERM_H */

