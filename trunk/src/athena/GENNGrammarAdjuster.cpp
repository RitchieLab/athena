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
#include "GENNGrammarAdjuster.h"
#include <sstream>
#include <map>


///
/// Reads in grammar file
/// @param filename
/// @param Mapper to pass grammar information to
///
void GENNGrammarAdjuster::read_grammar_file(string filename){
  ifstream file(filename.c_str(),ios::in|ios::binary);
	if(!file){
		cerr << "Could not open grammar file " << filename << ".\nExecution aborted.\n";
		exit(1);
	}
	char buffer[2048];
  lines.clear();
  string line;
	while(!file.eof()){
		file.getline(buffer, 2048);
		lines.push_back(buffer);
	}
	file.close();
}


///
/// Concatenates all lines and passes to the mapper for processing into grammar
/// @param mapper GEGrammarSI&
///
//void GENNGrammarAdjuster::set_mapper(GEGrammarSI& mapper){
void GENNGrammarAdjuster::set_mapper(AthenaGrammarSI& mapper, int rank){
  vector<string>::iterator iter;
  string program;
  
  for(iter = lines.begin(); iter != lines.end(); iter++){
    program += *iter + "\n";
  }
  program += "\n";

  mapper.readBNFString(program);
}



///
/// Reads the lines and selects the variables for storage and possible use later
/// returns the Genotype and Continuous variable strings
///
vector<string> GENNGrammarAdjuster::get_variables(){
  vector<string>::iterator iter;
  vector<string> vars;
  for(iter = lines.begin(); iter != lines.end(); ++iter){
    if(iter->find("<v>") == 0){
      break;
    }
  }
  
  std::size_t start, finish;
  
  // for each one look for G or C marking start of variables
  while((*iter).find_first_of("GC") != string::npos){
    string line = *iter;
    start = line.find_first_of("GC");
    finish = line.find_last_of("0123456789");
    vars.push_back(line.substr(start, finish));
    ++iter;
 }
  
  return vars;
}



///
/// Returns symbol name based on parameters 
/// @param vPrefix variable prefix 'G' for genotype and 'C' for continuous
/// @param num Number to include
/// @return string with symbol name
///
string GENNGrammarAdjuster::create_variable_symbol(string vPrefix, int num){
  stringstream ss;
  ss << num;
  return vPrefix + ss.str();
}


///
/// Returns iterator to start of variables line list in grammar
/// @return iterator to start
///
vector<string>::iterator GENNGrammarAdjuster::get_start_variables(){
  vector<string>::iterator start;
  for(start=lines.begin(); start != lines.end(); start++){
    if(start->find("<v>") == 0)
      break;
  }
  return start;
}


///
/// Returns number of lines associated with variables in 
/// current grammar
/// @param start iterator marking start of variable lines
/// @param last iterator containing last of the variable lines
/// @return number of lines in grammar
///
int GENNGrammarAdjuster::count_var_lines(vector<string>::iterator& start,
  vector<string>::iterator& last){

  int v_lines=1;
//  vector<string>::iterator iter;
  for(last=++start; last != lines.end(); last++){
    // look for ::= -- if there it means no more incrementing of lines
    if(last->find("::=") != string::npos){
      break;
    }
    if(last->find("|") != string::npos){
      v_lines++;
    }
    else{
      break;
    }
  }

  return v_lines;
}


///
/// Replaces any of the current genotypes or continuous variables specified in the grammar 
/// with rules placing every variable into the grammar
/// @param nGenos number of genotypes 
/// @param nContin number of continuous variables in set
///
void GENNGrammarAdjuster::include_all_vars(int nGenos, int nContin){

  vector<string>::iterator start = get_start_variables();
 
  vector<string>::iterator iter, last;
  int v_lines=count_var_lines(start,last);
 
  int c;
  start--; // move start back to correct row
  // first reset the first row
  if(nGenos >= 1){
    *start = "<v>      ::= G1";
    c = 1;
  }
  else{
    *start = "<v>      ::= C1";
    c = 2;
  }
  int v_used = 1;
  
  iter=start;
  
  string varType = "G", newSymbol, newLine;
  vector<string> holder;
  
  
  for(int g=2; g<=nGenos; g++){
    newSymbol = create_variable_symbol(varType, g);
    newLine = "           | " + newSymbol;
    if(v_used < v_lines){
      ++iter;
      *iter =  newLine;
      v_used++;
    }
    else{
      holder.push_back(newLine);
    }
  }
  
  varType = "C";
  for(; c<=nContin; c++){
    newSymbol = create_variable_symbol(varType, c);
    newLine = "           | " + newSymbol;
    if(v_used < v_lines){
      ++iter;
      *iter =  newLine;
      v_used++;
    }
    else{
      holder.push_back(newLine);
    }    
  }
 
  // if v_used is less than v_lines will need to erase the extra lines from the
  // grammar string vector
  if(v_used < v_lines){
    int num_to_delete = v_lines - v_used;
    vector<string>::iterator delete_start, delete_end;
    delete_start = start;
    int n;
    for(n=0; n < num_to_delete; n++){
      delete_start++;
    }
    delete_end = delete_start;
    for(;n < v_lines; n++){
      delete_end++;
    }
    lines.erase(delete_start, delete_end);
    }

  // insert in the holder strings if that vector contains anything
  if(!holder.empty()){
    // iter should hold location to insert vector after
    ++iter;
    lines.insert(iter, holder.begin(), holder.end());
    
  }
}


///
/// Expands shorthannd for continuous variables and genotypes into useable format.
/// The shorthand is for declaring multiple genotypes or continuous variables
/// at a time such as G1-10 or C5-7.
/// 
void GENNGrammarAdjuster::expand_variables(){
  
  // find line starting with "<v>" as start of genotypes and continuous variables
  vector<string>::iterator start = get_start_variables();
  vector<string>::iterator last;
  count_var_lines(start, last);

  start--; // move start back to first line

  // now start marks first line and last is last line in genotypes and continuous variables
  // need to process each line looking for '-' 
  map<string, bool> genos_included;
  map<string, bool> covars_included;
  std::size_t loc, dash, last_digit;
  string newSymbol;
  
  for(vector<string>::iterator iter=start; iter != last; ++iter){
    bool expand = false;
    last_digit = (*iter).find_last_of("0123456789");
    if(last_digit == string::npos)
      continue;
    
    if((dash=iter->find("-")) != string::npos){
      expand = true;
    }

    loc = iter->find_first_of("CG");
    string num;

    if(expand){
      stringstream ss;
      int first_num, last_num;
      for(std::size_t z=loc+1; z != dash; z++){
        num += (*iter)[z];
      }
      ss.str(num);
      ss >> first_num;
      num = "";
      for(std::size_t z=dash+1; z <= last_digit; z++){
        num += (*iter)[z];
      }
      stringstream lss;
      lss.str(num);
      lss >> last_num;
      string vPrefix;
      vPrefix += (*iter)[loc];
      for(int n = first_num; n <= last_num; n++){
        newSymbol = create_variable_symbol(vPrefix, n);
        if((*iter)[loc]=='G'){
          genos_included[newSymbol] = true;
        }
        else 
          covars_included[newSymbol] = true;
      }
    }
    else{
      for(std::size_t z=loc+1; z <= last_digit; z++){
        num+=(*iter)[z];
      }
      newSymbol = (*iter)[loc] + num;

      if((*iter)[loc]=='G'){
        genos_included[newSymbol] = true;
      }
      else{ 
        covars_included[newSymbol] = true;
      }
    }
  }
  
  // delete all current genotype and continuous lines
  // then insert new lines matching genotypes and continuous variables to include
  last--;  // move last to point back to final element
  vector<string>::iterator start_delete = start++; // move to first past start
  start--;
  if(start_delete != lines.end())
    lines.erase(start_delete, last);
  
  map<string, bool>::iterator genoiter = genos_included.begin();
  map<string, bool>::iterator covariter = covars_included.begin();
  
  // need to reset the first line
  if(!genos_included.empty()){
    *start = "<v>      ::= " + genoiter->first;
    genoiter++;
  }
  else{
    *start = "<v>      ::= " + covariter->first;
    covariter++;
  }
  
  start++;
  
  // create temp vector
  vector<string> holder;
  for(; genoiter != genos_included.end(); genoiter++){
    holder.push_back("           | " + genoiter->first);
  }
  
  for(; covariter != covars_included.end(); covariter++){
    holder.push_back("           | " + covariter->first);
  }
  
  if(!holder.empty()){
    lines.insert(start, holder.begin(), holder.end());
  }
  
}

///
/// Doubles number of genotypes for ott dummy encoding in the grammar
/// @param mapper GEGrammarSI
///
void GENNGrammarAdjuster::double_genotype_grammar(){
  vector<string>::iterator start = get_start_variables();

  vector<string>::iterator last;
  count_var_lines(start, last);
  
  // need to make a list of all the genotypes and covariates in the grammar
  vector<int> genos;
  vector<string> covars;  // no need to manipulate these
  
  std::size_t loc, last_digit;
  string newSymbol;

  start--;
  
  for(vector<string>::iterator iter=start; iter != last; ++iter){   
    loc = iter->find_first_of("CG");
    string num;
    last_digit = iter->find_last_of("0123456789");
    
    for(std::size_t z=loc+1; z <= last_digit; z++){
      num+=(*iter)[z];
    }
    
    if((*iter)[loc]=='G'){
      stringstream ss(num);
      int n;
      ss >> n;
      genos.push_back(n);
    }
    else {
      covars.push_back("C" + num);
    }
 
  }
  
  // delete all current genotype and continuous lines
  // then insert new lines matching genotypes and continuous variables to include
  last--;  // move last to point back to final element
  vector<string>::iterator start_delete = start++; // move to first past start
  start--;
  if(start_delete != lines.end())
    lines.erase(start_delete, last);

  vector<string> holder;
  // now begin adding back by setting start position
  // need to reset the first line
  unsigned int covarstart = 0;
  if(!genos.empty()){
    *start = "<v>      ::= " + create_variable_symbol("G", genos[0]*2-1);
    holder.push_back("           | " +create_variable_symbol("G", genos[0]*2));
  }
  else{
    *start = "<v>      ::= " + covars[0];
    covarstart++;
  }

  start++;

  
  // now begin adding new ones -- adjusting the genotypes 
  for(unsigned int i=1; i < genos.size(); i++){
    holder.push_back("           | " +create_variable_symbol("G", genos[i]*2-1));
    holder.push_back("           | " +create_variable_symbol("G", genos[i]*2));
  }
  // add covariates
  for(; covarstart < covars.size(); covarstart++){
    holder.push_back("           | " + covars[covarstart]);
  }
  
  // now insert all new lines
  if(!holder.empty()){
    lines.insert(start, holder.begin(), holder.end());
  }  
  
}



///
/// Edits grammar to include only the variables contained in the
/// variables_included set
///
void GENNGrammarAdjuster::edit_only_var_included(){
 
  vector<string>::iterator start = get_start_variables();
  vector<string>::iterator last;
  int v_lines=count_var_lines(start, last);
  
  // move back to start row
  start--;
 
  // get first variable from set
  set<string>::iterator varIter = variables_included.begin();

  *start = "<v>      ::= " + *varIter;
  varIter++;
  start++;
  int v_used = 1;
  
  for(;varIter != variables_included.end() && start != lines.end(); varIter++){
    *start = "           | " + *varIter;
    start++;
    v_used++;
  }

  // need to add more lines
  if(varIter != variables_included.end()){
    // create temp vector and then insert
    vector<string> temp;
    for(;varIter != variables_included.end(); varIter++){
      temp.push_back("           | " + *varIter);
    }
    lines.insert(start, temp.begin(), temp.end());
  }
  // delete extra variables from the grammar
  else if(v_used < v_lines){
    int num_to_delete = v_lines - v_used;
    vector<string>::iterator delete_start, delete_end;
    delete_start = delete_end = start;
    int n;
    for(n=0; n < num_to_delete; n++){
      delete_end++;
    }
    lines.erase(delete_start, delete_end);
  } 
  
}


///
/// Extracts variables from string vector and keeps set of
/// variables.  This set can be used to construct a new grammar
/// @param terminals string vector containing terminals
///
void GENNGrammarAdjuster::add_variables(vector<string>& terminals){
  
  vector<string>::iterator iter;
  
  std::size_t loc;
  
  for(iter = terminals.begin(); iter != terminals.end(); iter++){
    loc = iter->find_first_of("CG");
    // must start with C or G to be a variable
    if(loc == 0){
      // rest of string must only be digits
      loc =  iter->find_first_not_of("1234567890", 1);
      if(loc == string::npos){
        variables_included.insert(*iter);
      }
    }
  } 
}




