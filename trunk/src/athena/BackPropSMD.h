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
// BackPropSMD.h
#ifndef __BACKPROPSMD_H__
#define __BACKPROPSMD_H__

#include "Terminals.h"
#include <Dataset.h>
#include "BackPropTree.h"
#include "BackProp.h"

///
/// Runs back propagation on the usual GENN network representation.
/// Makes a few assumptions on the how back propagation will work in 
/// this case.  
///
class BackPropSMD:public BackProp{

  public:
  
    BackPropSMD();
    
    /// Receives a postfix stack representation of the neural network
    int runBackProp(vector<TerminalSymbol*> & postfix_stack, Dataset * set);
    
    /// Sets whether to use stochastic
    void useBatch(bool tf){runBatch = tf;}

  private:
    void initialize();
    float runBatchEpoch(Dataset* set);
    float runStochasticEpoch(Dataset* set);
    double calculateMSE(Dataset* set);

    bool runBatch;
    int num_epochs, stochastic_size;

};

#endif
