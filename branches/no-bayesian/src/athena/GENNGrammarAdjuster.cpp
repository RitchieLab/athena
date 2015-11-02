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
#include <Stringmanip.h>


///
/// Reads in grammar file
/// @param filename
/// @param Mapper to pass grammar information to
///
void GENNGrammarAdjuster::readGrammarFile(string filename){
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
		string newline(buffer);
// 		if(newline.find_first_of("01234567890abcdedfghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ<>,")!=string::npos)
			lines.push_back(buffer);
	}
	file.close();
}



///
/// Concatenates all lines and passes to the mapper for processing into grammar
/// @param mapper AthenaGrammarSI&
///
void GENNGrammarAdjuster::setMapper(AthenaGrammarSI& mapper, int rank){
	vector<string>::iterator iter;
	string program;
	
	for(iter = lines.begin(); iter != lines.end(); iter++){
		program += *iter + "\n";
	}
	program += "\n";
// cout << program << endl;
// exit(1);
	mapper.readBNFString(program);
}



///
/// Reads the lines and selects the variables for storage and possible use later
/// returns the Genotype and Continuous variable strings
///
vector<string> GENNGrammarAdjuster::getVariables(){
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
string GENNGrammarAdjuster::createVariableSymbol(string vPrefix, int num){
	stringstream ss;
	ss << num;
	return vPrefix + ss.str();
}


///
/// Returns iterator to start of variables line list in grammar
/// @return iterator to start
///
vector<string>::iterator GENNGrammarAdjuster::getStartVariables(){
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
int GENNGrammarAdjuster::countVarLines(vector<string>::iterator& start,
	vector<string>::iterator& last){

	int vLines=1;
	for(last=++start; last != lines.end(); last++){
		// look for ::= -- if there it means no more incrementing of lines
		if(last->find("::=") != string::npos){
			break;
		}
		if(last->find("|") != string::npos){
			vLines++;
		}
		else{
			break;
		}
	}

	return vLines;
}


///
/// Replaces any of the current genotypes or continuous variables specified in the grammar 
/// with rules placing every variable into the grammar
/// @param nGenos number of genotypes 
/// @param nContin number of continuous variables in set
///
void GENNGrammarAdjuster::includeAllVars(int nGenos, int nContin){

	vector<string>::iterator start = getStartVariables();
 
	vector<string>::iterator iter, last;
	int vLines=countVarLines(start,last);
 
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
	int vUsed = 1;
	
	iter=start;
	
	string varType = "G", newSymbol, newLine;
	vector<string> holder;
	
	
	for(int g=2; g<=nGenos; g++){
		newSymbol = createVariableSymbol(varType, g);
		newLine = "           | " + newSymbol;
		if(vUsed < vLines){
			++iter;
			*iter =  newLine;
			vUsed++;
		}
		else{
			holder.push_back(newLine);
		}
	}
	
	varType = "C";
	for(; c<=nContin; c++){
		newSymbol = createVariableSymbol(varType, c);
		newLine = "           | " + newSymbol;
		if(vUsed < vLines){
			++iter;
			*iter =  newLine;
			vUsed++;
		}
		else{
			holder.push_back(newLine);
		}    
	}
 
	// if vUsed is less than vLines will need to erase the extra lines from the
	// grammar string vector
	if(vUsed < vLines){
		int numToDelete = vLines - vUsed;
		vector<string>::iterator deleteStart, deleteEnd;
		deleteStart = start;
		int n;
		for(n=0; n < numToDelete; n++){
			deleteStart++;
		}
		deleteEnd = deleteStart;
		for(;n < vLines; n++){
			deleteEnd++;
		}
		lines.erase(deleteStart, deleteEnd);
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
void GENNGrammarAdjuster::expandVariables(){
	
	// find line starting with "<v>" as start of genotypes and continuous variables
	vector<string>::iterator start = getStartVariables();
	vector<string>::iterator last;
	countVarLines(start, last);

	start--; // move start back to first line

	// now start marks first line and last is last line in genotypes and continuous variables
	// need to process each line looking for '-' 
	map<string, bool> genosIncluded;
	map<string, bool> covarsIncluded;
	std::size_t loc, dash, lastDigit;
	string newSymbol;
	
	for(vector<string>::iterator iter=start; iter != last; ++iter){
		bool expand = false;
		lastDigit = (*iter).find_last_of("0123456789");
		if(lastDigit == string::npos)
			continue;
		
		if((dash=iter->find("-")) != string::npos){
			expand = true;
		}

		loc = iter->find_first_of("CG");
		string num;

		if(expand){
			stringstream ss;
			int firstNum, lastNum;
			for(std::size_t z=loc+1; z != dash; z++){
				num += (*iter)[z];
			}
			ss.str(num);
			ss >> firstNum;
			num = "";
			for(std::size_t z=dash+1; z <= lastDigit; z++){
				num += (*iter)[z];
			}
			stringstream lss;
			lss.str(num);
			lss >> lastNum;
			string vPrefix;
			vPrefix += (*iter)[loc];
			for(int n = firstNum; n <= lastNum; n++){
				newSymbol = createVariableSymbol(vPrefix, n);
				if((*iter)[loc]=='G'){
					genosIncluded[newSymbol] = true;
				}
				else 
					covarsIncluded[newSymbol] = true;
			}
		}
		else{
			for(std::size_t z=loc+1; z <= lastDigit; z++){
				num+=(*iter)[z];
			}
			newSymbol = (*iter)[loc] + num;

			if((*iter)[loc]=='G'){
				genosIncluded[newSymbol] = true;
			}
			else{ 
				covarsIncluded[newSymbol] = true;
			}
		}
	}
	
	// delete all current genotype and continuous lines
	// then insert new lines matching genotypes and continuous variables to include
	last--;  // move last to point back to final element
	vector<string>::iterator startDelete = start++; // move to first past start
	start--;
	if(startDelete != lines.end())
		lines.erase(startDelete, last);
	
	map<string, bool>::iterator genoiter = genosIncluded.begin();
	map<string, bool>::iterator covariter = covarsIncluded.begin();
	
	// need to reset the first line
	if(!genosIncluded.empty()){
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
	for(; genoiter != genosIncluded.end(); genoiter++){
		holder.push_back("           | " + genoiter->first);
	}
	
	for(; covariter != covarsIncluded.end(); covariter++){
		holder.push_back("           | " + covariter->first);
	}
	
	if(!holder.empty()){
		lines.insert(start, holder.begin(), holder.end());
	}
}


///
/// Excludes indicated variables from grammar
///
void GENNGrammarAdjuster::excludeVariables(std::set<std::string>& varSet){

// vector<string>::iterator tmpIter;
// for(tmpIter=lines.begin(); tmpIter != lines.end(); ++tmpIter){
// cout << *tmpIter << endl;
// }

	// find line starting with "<v>" as start of genotypes and continuous variables
	vector<string>::iterator start = getStartVariables();
	vector<string>::iterator last;
	countVarLines(start, last);
// cout << "last=" << *last << endl;
	vector<string> tmp;
	
	vector<string>::iterator iter = lines.begin();
	// copy up to the start of the variables
	for(;iter != start-1; ++iter){
		tmp.push_back(*iter);
	}
	
// 	vector<string>::iterator final = last+1;
	// check start and see if it is to be excluded
	string varname = getVarname(*iter);
// cout << "varname=" << varname << endl;
	if(varSet.find(varname) == varSet.end()){
		tmp.push_back(*iter);
		++iter;
	}
	else{
		while(varSet.find(varname) != varSet.end() && iter != last){
			++iter;
			varname = getVarname(*iter);
		}
		if(iter != last){
			tmp.push_back("<v>      ::= " + varname);
			++iter;
		}
	}
	
	while(iter != last){
		varname = getVarname(*iter);
		if(varSet.find(varname) == varSet.end()){
			tmp.push_back(*iter);
		}
		++iter;
	}
	
	for(;iter != lines.end(); ++iter){
		tmp.push_back(*iter);
	}
	
	lines = tmp;
// cout << "--------------------------------" << endl;
// for(tmpIter=lines.begin(); tmpIter != lines.end(); ++tmpIter){
// cout << *tmpIter << endl;
// }
// exit(1);
}

///
/// extracts variable name from a variable line
/// @returns name of variable from grammar
///
string GENNGrammarAdjuster::getVarname(string varLine){
	size_t nameStart, nameEnd;
	nameStart=varLine.find_first_of("CG");
	nameEnd = varLine.find_last_of("0123456789");
	return varLine.substr(nameStart, nameEnd-nameStart+1);
}

///
/// Doubles number of genotypes for ott dummy encoding in the grammar
///
void GENNGrammarAdjuster::doubleGenotypeGrammar(){
	vector<string>::iterator start = getStartVariables();

	vector<string>::iterator last;
	countVarLines(start, last);
	
	// need to make a list of all the genotypes and covariates in the grammar
	vector<int> genos;
	vector<string> covars;  // no need to manipulate these
	
	std::size_t loc, lastDigit;
	string newSymbol;

	start--;
	
	for(vector<string>::iterator iter=start; iter != last; ++iter){   
		loc = iter->find_first_of("CG");
		string num;
		lastDigit = iter->find_last_of("0123456789");
		
		for(std::size_t z=loc+1; z <= lastDigit; z++){
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
	vector<string>::iterator startDelete = start++; // move to first past start
	start--;
	if(startDelete != lines.end())
		lines.erase(startDelete, last);

	vector<string> holder;
	// now begin adding back by setting start position
	// need to reset the first line
	unsigned int covarStart = 0;
	if(!genos.empty()){
		*start = "<v>      ::= " + createVariableSymbol("G", genos[0]*2-1);
		holder.push_back("           | " +createVariableSymbol("G", genos[0]*2));
	}
	else{
		*start = "<v>      ::= " + covars[0];
		covarStart++;
	}

	start++;

	
	// now begin adding new ones -- adjusting the genotypes 
	for(unsigned int i=1; i < genos.size(); i++){
		holder.push_back("           | " +createVariableSymbol("G", genos[i]*2-1));
		holder.push_back("           | " +createVariableSymbol("G", genos[i]*2));
	}
	// add covariates
	for(; covarStart < covars.size(); covarStart++){
		holder.push_back("           | " + covars[covarStart]);
	}
	
	// now insert all new lines
	if(!holder.empty()){
		lines.insert(start, holder.begin(), holder.end());
	}  
	
}



///
/// Edits grammar to include only the variables contained in the
/// variablesIncluded set
///
void GENNGrammarAdjuster::editOnlyVarIncluded(){
 
	vector<string>::iterator start = getStartVariables();
	vector<string>::iterator last;
	int vLines=countVarLines(start, last);
	
	// move back to start row
	start--;
 
	// get first variable from set
	set<string>::iterator varIter = variablesIncluded.begin();

	*start = "<v>      ::= " + *varIter;
	varIter++;
	start++;
	int vUsed = 1;
	
	for(;varIter != variablesIncluded.end() && start != lines.end(); varIter++){
		*start = "           | " + *varIter;
		start++;
		vUsed++;
	}

	// need to add more lines
	if(varIter != variablesIncluded.end()){
		// create temp vector and then insert
		vector<string> temp;
		for(;varIter != variablesIncluded.end(); varIter++){
			temp.push_back("           | " + *varIter);
		}
		lines.insert(start, temp.begin(), temp.end());
	}
	// delete extra variables from the grammar
	else if(vUsed < vLines){
		int numToDelete = vLines - vUsed;
		vector<string>::iterator deleteStart, deleteEnd;
		deleteStart = deleteEnd = start;
		int n;
		for(n=0; n < numToDelete; n++){
			deleteEnd++;
		}
		lines.erase(deleteStart, deleteEnd);
	} 
}


///
/// Sets parent and child grammar to restrict allowed sizes
/// @param minP minimum number of parents in Bayesian network for a node
/// @param maxP maximum number of parents in Bayesian network for a node
/// @param minC minimum number of children in Bayesian network for a node
/// @param maxC maximum number of children in Bayesian network for a node
///
void GENNGrammarAdjuster::setBayesianSize(unsigned int minP, unsigned int maxP, unsigned int minC,
	unsigned int maxC){
	
	vector<string>::iterator lineIter = lines.begin(), startIter, endIter;
	for(;lineIter != lines.end(); ++lineIter){
		if(lineIter->find("<parents>")==0){
			startIter=lineIter;
			break;
		}
	}
	for(endIter=startIter+1; endIter != lines.end(); ++endIter){
		if(endIter->find("<")==0){
			break;
		}
	}
	
	// delete extra lines
	if(endIter != startIter){
		endIter = lines.erase(startIter, endIter);
	}
	else{
		endIter = lines.erase(startIter);
	}
	
	vector<string> newLines;
	unsigned int i = minP;
	newLines.push_back("<parents>  ::= ");
	for(unsigned int j=1; j<=i; j++){
		newLines.back() += "<v>";
	}
	for(i=minP+1; i<=maxP; i++){
		newLines.push_back("           | ");
		for(unsigned int j=1; j<=i; j++){
			newLines.back() += "<v>";
		}
	}
	
	// insert new lines before the endIter
	lines.insert(endIter, newLines.begin(), newLines.end());
	
	for(lineIter = lines.begin();lineIter != lines.end(); ++lineIter){
		if(lineIter->find("<children>")==0){
			startIter=lineIter;
			break;
		}
	}
	for(endIter=startIter+1; endIter != lines.end(); ++endIter){
		if(endIter->find("<")==0){
			break;
		}
	}
	
	// delete extra lines
	if(endIter != startIter){
		endIter = lines.erase(startIter, endIter);
	}
	else{
		endIter = lines.erase(startIter);
	}
	
	newLines.clear();
	i=minC;
	newLines.push_back("<children>  ::= ");
	for(unsigned int j=1; j<=i; j++){
		newLines.back() += "<child>";
	}
	for(i=minC+1; i<=maxC; i++){
		newLines.push_back("           | ");
		for(unsigned int j=1; j<=i; j++){
			newLines.back() += "<child>";
		}
	}	
	
	// insert new lines before the endIter
	lines.insert(endIter, newLines.begin(), newLines.end());
	
}


///
/// Extracts variables from string vector and keeps set of
/// variables.  This set can be used to construct a new grammar
/// @param terminals string vector containing terminals
///
void GENNGrammarAdjuster::addVariables(vector<string>& terminals){
	
	vector<string>::iterator iter;
	std::size_t loc;
	
	for(iter = terminals.begin(); iter != terminals.end(); iter++){
		loc = iter->find_first_of("CG");
		// must start with C or G to be a variable
		if(loc == 0){
			// rest of string must only be digits
			loc =  iter->find_first_not_of("1234567890", 1);
			if(loc == string::npos){
				variablesIncluded.insert(*iter);
			}
		}
	} 
}


///
/// Create and include constants in the grammar.  Append them to the existing grammar
/// and remove the original Constant lines.
/// @param min constant value
/// @param max constant value
/// @param interval between each constant value
///
void GENNGrammarAdjuster::setConstants(float min, float max, float interval){

	float center = (max+min)/2.0;

	vector<string> constantLines;
	lines.push_back("<Constant> ::= " + data_manage::Stringmanip::numberToString(center));

	float increasingValue=center, decreasingValue=center;
	increasingValue += interval;
	decreasingValue -= interval;
	char buffer[20];
	constantsIncluded.clear();
	
	while(increasingValue <= max){
// 		constantLines.push_back("           | " + data_manage::Stringmanip::numberToString(increasingValue));
// 		constantLines.push_back("           | " + data_manage::Stringmanip::numberToString(decreasingValue));
		sprintf(buffer, "%.3f",increasingValue);
		constantLines.push_back("           | " + string(buffer));
		constantsIncluded.push_back(buffer);
		sprintf(buffer, "%.3f",decreasingValue);
		constantLines.push_back("           | " + string(buffer));
		constantsIncluded.push_back(buffer);
		increasingValue += interval;
		decreasingValue -= interval;
	}
	
	vector<string>::iterator start, last, iter;
	// delete the original constant lines
	for(iter=lines.begin(); iter != lines.end(); ++iter){
		if(iter->find("<Constant>")==0){
			start=iter;
			break;
		}
	}
	for(last=start+1;last!= lines.end(); ++last){
		if(last->find("::=") != string::npos){
			break;
		}
	}
	
	lines.erase(start, last);
	
	// append the new ones
	lines.insert(lines.end(), constantLines.begin(), constantLines.end());
	
}
