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

