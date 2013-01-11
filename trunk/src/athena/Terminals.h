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
 * File:   Terminals.h
 * Author: dudeksm
 *
 * Created on November 14, 2008, 4:03 PM
 */

#ifndef _TERMINALS_H
#define	_TERMINALS_H

#include "TerminalSymbol.h"
#include <vector>
#include <Individual.h>

using namespace data_manage;

class Constant :public TerminalSymbol{
  
  public:
    Constant(std::string symbol, int number_args=0);
    Constant(float val);
    virtual ~Constant(){}
    
    virtual float evaluate(std::deque<float> & TerminalSymbols);

    private:
      float constant_value;
};

///////////////////////////////
// user-defined terminals below
///////////////////////////////


///
/// TerminalSymbol that acts as bias input for nodes
///
class BiasTerm :public TerminalSymbol{

	public:
		BiasTerm(std::string name, int number_args=0);
		
		virtual ~BiasTerm(){}
		
		virtual float evaluate(std::deque<float> & TerminalSymbols);

    private:
      float bias_value;		 

};



///
/// TerminalSymbol that contains an index for referencing 
/// genotypes in the dataset
///
class GenotypeTerm :public TerminalSymbol{
  
  public:
    GenotypeTerm(std::string name, int variable_index);
    
    virtual ~GenotypeTerm(){}
    
    virtual float evaluate(std::deque<float> & TerminalSymbols);
    
    inline int getIndex(){return index_value;}
    
    static void set_ind(Individual* ind);
    private:
      int index_value;
      
      static Individual* ind;
};


///
/// TerminalSymbol that contains an index for referencing 
/// genotypes in the dataset
///
class ContinVariable :public TerminalSymbol{
  
  public:
    ContinVariable(std::string name, int variable_index);
    
    virtual ~ContinVariable(){}
    
    virtual float evaluate(std::deque<float> & TerminalSymbols);
    static void set_ind(Individual* ind);
    
    inline int getIndex(){return index_value;}
    private:
      int index_value;
      static Individual* ind;
};


class Addition :public TerminalSymbol{
  
  public:
    Addition(std::string symbol, int number_args);
    virtual ~Addition(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Subtraction :public TerminalSymbol{
  
  public:
    Subtraction(std::string symbol, int number_args);
    virtual ~Subtraction(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Multiplication :public TerminalSymbol{
  
  public:
    Multiplication(std::string symbol, int number_args);
    virtual ~Multiplication(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Division :public TerminalSymbol{
  
  public:
    Division(std::string symbol, int number_args);
    virtual ~Division(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Power :public TerminalSymbol{
  
  public:
    Power(std::string symbol, int number_args);
    virtual ~Power(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class pAdd :public TerminalSymbol{
  
  public:
    pAdd(std::string symbol, int number_args);
    virtual ~pAdd(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class pSub :public TerminalSymbol{
  
  public:
    pSub(std::string symbol, int number_args);
    virtual ~pSub(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class pMult :public TerminalSymbol{
  
  public:
    pMult(std::string symbol, int number_args);
    virtual ~pMult(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class pDiv :public TerminalSymbol{
  
  public:
    pDiv(std::string symbol, int number_args);
    virtual ~pDiv(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Weight :public TerminalSymbol{
  
  public:
    Weight(std::string symbol, int number_args);
    virtual ~Weight(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};


class Dot :public TerminalSymbol{
  
  public:
    Dot(std::string symbol, int number_args=0);
    virtual ~Dot(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class ConCat :public TerminalSymbol{

  public:
    ConCat(std::string symbol, int number_args);
    virtual ~ConCat(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
    friend ostream & operator << (ostream & os, const TerminalSymbol & el);
};


class LogF  :public TerminalSymbol{
  public:
    LogF(std::string symbol, int number_args);
    virtual ~LogF(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Sine  :public TerminalSymbol{
  public:
    Sine(std::string symbol, int number_args);
    virtual ~Sine(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Cosine  :public TerminalSymbol{
  public:
    Cosine(std::string symbol, int number_args);
    virtual ~Cosine(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Tangent  :public TerminalSymbol{
  public:
    Tangent(std::string symbol, int number_args);
    virtual ~Tangent(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

   
// Following are the boolean operators for neural nets
class pAnd :public TerminalSymbol{
  
  public:
    pAnd(std::string symbol, int number_args);
    virtual ~pAnd(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class pNand :public TerminalSymbol{
  
  public:
    pNand(std::string symbol, int number_args);
    virtual ~pNand(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class pOr :public TerminalSymbol{
  
  public:
    pOr(std::string symbol, int number_args);
    virtual ~pOr(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class pNor :public TerminalSymbol{
  
  public:
    pNor(std::string symbol, int number_args);
    virtual ~pNor(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class pXor :public TerminalSymbol{
  
  public:
    pXor(std::string symbol, int number_args);
    virtual ~pXor(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

// Following are the boolean operators
class And :public TerminalSymbol{
  public:
    And(std::string symbol, int number_args);
    virtual ~And(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Nand :public TerminalSymbol{
  
  public:
    Nand(std::string symbol, int number_args);
    virtual ~Nand(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Or :public TerminalSymbol{
  
  public:
    Or(std::string symbol, int number_args);
    virtual ~Or(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Nor :public TerminalSymbol{
  
  public:
    Nor(std::string symbol, int number_args);
    virtual ~Nor(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};

class Xor :public TerminalSymbol{
  
  public:
    Xor(std::string symbol, int number_args);
    virtual ~Xor(){}
    virtual float evaluate(std::deque<float> & TerminalSymbols);
};


#endif	/* _TERMINALS_H */

