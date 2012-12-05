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
#include "NNSolution.h"
#include "TerminalSymbol.h"
#include "Terminals.h"
#include "CalculatorFactory.h"
#include <Stringmanip.h>
#include <sstream>

#include <iostream>

///
/// Clones current solution
///
Solution* NNSolution::clone(){
    NNSolution* newNNSol = new NNSolution;
 
    newNNSol->nn_depth=nn_depth;
    newNNSol->gram_depth=gram_depth;
    Solution* newSol = newNNSol;
    Solution* thiscopy = this;
    
    newSol->copy(thiscopy);

    return newSol;
}


///
/// Adjusts dummy encoding 
/// @parm genotype Index of genotype (1 is first)
/// @return 
///
int NNSolution::adjust_dummy_encoding(int genotype){
    return (genotype+1)/2;
}


///
/// Returns the genotypes in the solution
/// @param dummy_encoded when true will adjust genotype numbers
/// to reflect dummy encoding
/// @return vector of ints listing genotype index numbers
///
vector<int> NNSolution::get_genotypes(bool dummy_encoded){
    
    vector<int> genotypeIndexes;
    
    std::vector<std::string>::iterator symbolIter;
    for(symbolIter = symbols.begin(); symbolIter != symbols.end(); ++symbolIter){
      // look for letter 'G' followed by digits
      // when only one character long go to next symbol
      if((*symbolIter)[0] == 'G')            
          if((*symbolIter).size() >= 2){
             string num = (*symbolIter).substr(1, (*symbolIter).size()-1);
            if(Stringmanip::is_number(num)){
                int genoNum = Stringmanip::stoi(num);
                if(dummy_encoded)
                    genotypeIndexes.push_back(adjust_dummy_encoding(genoNum));
                else
                    genotypeIndexes.push_back(genoNum);
            }
          }
    }
    
    return genotypeIndexes;
}


///
/// Returns the covariates in the solution
/// @return vector of ints listing covariate numbers
///
vector<int> NNSolution::get_covariates(){
    
    vector<int> covariateIndexes;
    
    std::vector<std::string>::iterator symbolIter;
    for(symbolIter = symbols.begin(); symbolIter != symbols.end(); ++symbolIter){
      // look for letter 'C' followed by digits
      // when only one character long go to next symbol
      if((*symbolIter)[0] == 'C')            
          if((*symbolIter).size() >= 2){
            string num = (*symbolIter).substr(1, (*symbolIter).size()-1);
            if(Stringmanip::is_number(num)){
               covariateIndexes.push_back(Stringmanip::stoi(num));
            }
          }  
    }
    
    return covariateIndexes;
}


///
/// Cleans up output to remove extra rule information. Changes Concat operator 
/// into corresponding numbers for output
/// @param os ostream to write to
/// @param data
/// @param map_used
/// 
void NNSolution::output_clean(std::ostream& os, data_manage::Dataholder& data,
      bool map_used, bool ott_dummy, bool continmap_used){

  // concatenate numbers and output rest without spaces
  os << symbols[0];
  for(unsigned int symb=1; symb < symbols.size(); symb++){
    if(symbols[symb].compare("Concat")!=0){
      if(symbols[symb].compare("W")==0 || symbols[symb][0] == 'P'){
        os << " ";
        os << symbols[symb];
      }
      else if(map_used && symbols[symb][0] == 'G'){
        stringstream ss(symbols[symb].substr(1,symbols[symb].length()-1));
        int num;
        ss >> num;
        if(ott_dummy)
          num = (num-1)/2;
        else
          num -= 1;

        os << data.get_geno_name(num);
      }
      else if(continmap_used && symbols[symb][0] == 'C'){
        stringstream ss(symbols[symb].substr(1,symbols[symb].length()-1));
        int num;
        ss >> num;
        num -=1; // starts at 0 but labels start at 1
        os << data.get_covar_name(num);
      }
      else{
      	if(symbols[symb][0]=='G'){
      	  stringstream ss(symbols[symb].substr(1,symbols[symb].length()-1));
          int num;
          ss >> num;
          if(ott_dummy)
            num = (num-1)/2;
          else
            num -= 1;
      	  os << "G" << data.get_geno_name(num);
      	}
      	else      
	        os << symbols[symb];
      }
      
    }
    else{  // for Concat -- concatenate the numbers to create constant
      // skip next symbol -- must be '('
      // last number is not part of the concatenated number it is an indicator of the
      // number of arguments to pass to the concatentation operator so it is removed for
      // display purposes
      symb++;
      string num;
      while(symbols[++symb].compare(")")!=0){
        num += symbols[symb];
      }
      num = num.substr(0,num.size()-1);
      float realnum = Stringmanip::stodouble(num);
      os << realnum;
    }
  }


}

///
/// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
/// Works to alter mean squared error scores to r squared
/// @param tes
///
void NNSolution::adjust_score_out(Dataset* train_set, Dataset* test_set){

  float sstotal;
  int total_inds = calc_inds(train_set, sstotal);
  sol_fitness = alter_score(sol_fitness, total_inds, sstotal);

  total_inds = calc_inds(test_set, sstotal);

  test_score = alter_score(test_score, total_inds, sstotal);
}


///
/// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
/// Works to alter mean squared error scores to r squared
/// @param train_set Dataset
///
void NNSolution::adjust_score_out(Dataset* train_set){
  float sstotal;
  int total_inds = calc_inds(train_set, sstotal);

  sol_fitness = alter_score(sol_fitness, total_inds, sstotal);
}


///
/// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
/// Works to alter mean squared error scores to r squared
/// @param score
/// @param nIndsTested
/// @param sstotal
/// @return R-squared value
///
float NNSolution::adjust_score_out(float score, int nIndsTested, float sstotal){
  return alter_score(score, nIndsTested, sstotal);
}


///
/// Converts mean squared error 
/// @param set Dataset for this score
/// @param mse mean squared error
/// @param total_inds Total inds evaluated
/// @param sstotal
/// @return R squared score for the set
///
float NNSolution::alter_score(float mse, int total_inds, float sstotal){
   float rsquared = 1-((mse * total_inds) / sstotal);
  return rsquared;
}


///
/// Returns number of individuals from the set who are included in the calculation.
/// Individuals without data for one of the variables in the set are not included.
/// @param set Dataset
/// @param sstotal sets the sstotal
/// @return number of individuals
///
int NNSolution::calc_inds(Dataset* set, float& sstotal){

  // iterate through the symbols and look for anything with a 'G' or 'C'
  std::vector<std::string>::iterator iter;
  
  vector<int> genos, covars;
  int num;
  
  for(iter = symbols.begin(); iter != symbols.end(); iter++){
    string sym = *iter;
    if(sym[0] == 'G'){
      stringstream ss(sym.substr(1, sym.length()-1));
      ss >> num;
      genos.push_back(num-1);
    }
    else if(sym[0] == 'C' && sym.find_first_of("0123456789") != string::npos){
      stringstream ss(sym.substr(1, sym.length()-1));
      ss >> num;
      covars.push_back(num-1);
    }
  }

  int total_inds = 0;
  Individual* ind;
  bool any_missing;
  
  float diff=0.0, meanval, stat_total=0.0;
  vector<float> status_used;
  
  for(unsigned int currind=0; currind < set->num_inds(); currind++){
    ind = (*set)[currind];
    
    any_missing = false;
    for(unsigned int g=0; g < genos.size(); g++){
      if(ind->get_genotype(genos[g]) == set->get_missing_genotype()){
        any_missing = true;
        break;
      }
    }
    for(unsigned int c=0; c < covars.size(); c++){
      if(ind->get_covariate(covars[c]) == set->get_missing_covalue()){
        any_missing = true;
        break;
      }
    }
    
    if(!any_missing){
      total_inds++;
      stat_total += ind->get_status();
      status_used.push_back(ind->get_status());
    }
  }
  
  meanval = stat_total / total_inds;
  for(vector<float>::iterator iter=status_used.begin(); iter != status_used.end();
    ++iter){
    diff = diff + (*iter-meanval) * (*iter-meanval);
  }
  sstotal=diff;
  
  return total_inds;
  
}
