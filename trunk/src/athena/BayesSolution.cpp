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
#include "BayesSolution.h"
#include "TerminalSymbol.h"
#include "Terminals.h"
#include "CalculatorFactory.h"
#include <Stringmanip.h>
#include <sstream>

#include <iostream>

///
/// Clones current solution
///
Solution* BayesSolution::clone(){
		BayesSolution* newNNSol = new BayesSolution;
 
// 		newNNSol->nnDepth=nnDepth;
// 		newNNSol->gramDepth=gramDepth;
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
int BayesSolution::adjustDummyEncoding(int genotype){
		return (genotype+1)/2;
}


///
/// Returns the genotypes in the solution
/// @param dummyEncoded when true will adjust genotype numbers
/// to reflect dummy encoding
/// @return vector of ints listing genotype index numbers
///
vector<int> BayesSolution::getGenotypes(bool dummyEncoded){
		
		vector<int> genotypeIndexes;
		
		std::vector<std::string>::iterator symbolIter;
		for(symbolIter = symbols.begin(); symbolIter != symbols.end(); ++symbolIter){
			// look for letter 'G' followed by digits
			// when only one character long go to next symbol
			if((*symbolIter)[0] == 'G')            
					if((*symbolIter).size() >= 2){
						 string num = (*symbolIter).substr(1, (*symbolIter).size()-1);
						if(Stringmanip::is_number(num)){
								int genoNum = Stringmanip::stringToNumber<int>(num);
								if(dummyEncoded)
										genotypeIndexes.push_back(adjustDummyEncoding(genoNum));
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
vector<int> BayesSolution::getCovariates(){
		
		vector<int> covariateIndexes;
		
		std::vector<std::string>::iterator symbolIter;
		for(symbolIter = symbols.begin(); symbolIter != symbols.end(); ++symbolIter){
			// look for letter 'C' followed by digits
			// when only one character long go to next symbol
			if((*symbolIter)[0] == 'C')            
					if((*symbolIter).size() >= 2){
						string num = (*symbolIter).substr(1, (*symbolIter).size()-1);
						if(Stringmanip::is_number(num)){
							 covariateIndexes.push_back(Stringmanip::stringToNumber<int>(num));
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
/// @param mapUsed
/// 
void BayesSolution::outputClean(std::ostream& os, data_manage::Dataholder& data,
			bool mapUsed, bool ottDummy, bool continMapUsed){
	// concatenate numbers and output rest without space
// 	os << symbols[0];
	for(unsigned int symb=0; symb < symbols.size(); symb++){
		if(symbols[symb].compare("Concat")!=0){
			if(symbols[symb].compare("W")==0 || symbols[symb][0] == 'P'){
				os << " ";
				os << symbols[symb];
			}
			else if(mapUsed && symbols[symb][0] == 'G'){
				stringstream ss(symbols[symb].substr(1,symbols[symb].length()-1));
				int num;
				ss >> num;
				if(ottDummy)
					num = (num-1)/2;
				else
					num -= 1;
				os << data.getGenoName(num);
			}
			else if(continMapUsed && symbols[symb][0] == 'C'){
				stringstream ss(symbols[symb].substr(1,symbols[symb].length()-1));
				int num;
				ss >> num;
				num -=1; // starts at 0 but labels start at 1
				os << data.getCovarName(num);
			}
			else{
				if(symbols[symb][0]=='G'){
				  stringstream ss(symbols[symb].substr(1,symbols[symb].length()-1));
					int num;
					ss >> num;
					if(ottDummy)
						num = (num-1)/2;
					else
						num -= 1;
				  os << "G" << data.getGenoName(num);
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
			float realnum = Stringmanip::stringToNumber<double>(num);
			os << realnum;
		}
	}
}



///
/// Adjusts output of the scores when needed (e.g. meansquared to rSquared)
/// Works to alter mean squared error scores to r squared
/// @param tes
///
void BayesSolution::adjustScoreOut(Dataset* trainSet, Dataset* testSet){
	solFitness = alterScore(solFitness, trainSet->getConstant());
	testScore = alterScore(testScore, testSet->getConstant());
}


///
/// Adjusts output of the scores when needed (e.g. meansquared to rSquared)
/// Works to alter mean squared error scores to r squared
/// @param trainSet Dataset
///
void BayesSolution::adjustScoreOut(Dataset* trainSet){
	solFitness = alterScore(solFitness, trainSet->getConstant());
}


///
/// Adjusts score passed and returns value
///
float BayesSolution::adjustScoreOut(float score, int nIndsTested, 
	float constant, std::string calcName){
	return alterScore(score, constant);
}

///
/// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
///
void BayesSolution::adjustScoreOut(Dataset* trainSet, Dataset* testSet, 
	std::string calcName){
	solFitness = alterScore(solFitness, trainSet->getConstant());
	testScore = alterScore(testScore, testSet->getConstant());
}
		
///		
/// Adjusts output of the scores when needed
///
void BayesSolution::adjustScoreOut(Dataset* trainSet, std::string calcName){
	solFitness = alterScore(solFitness, trainSet->getConstant());
}


///
/// Adjusts output of the scores when needed (e.g. meansquared to rSquared)
/// Works to alter mean squared error scores to r squared
/// @param score
/// @param nIndsTested
/// @param ssTotal
/// @return R-squared value
///
float BayesSolution::adjustScoreOut(float score, int nIndsTested, float ssTotal){
	return alterScore(score, ssTotal);
}


///
/// Converts difference in network to total score
/// @param score
/// @param c value of unconnected network
/// @return R squared score for the set
///
float BayesSolution::alterScore(float score, double c){
// 	if(score > 0.0){
		score = c + score + c;
// 	}
	return score;
}


///
/// Returns number of individuals from the set who are included in the calculation.
/// Individuals without data for one of the variables in the set are not included.
/// @param set Dataset
/// @param ssTotal sets the ssTotal
/// @return number of individuals
///
// int BayesSolution::calcInds(Dataset* set, float& ssTotal){
// 
// 	// iterate through the symbols and look for anything with a 'G' or 'C'
// 	std::vector<std::string>::iterator iter;
// 	
// 	vector<int> genos, covars;
// 	int num;
// 	
// 	for(iter = symbols.begin(); iter != symbols.end(); iter++){
// 		string sym = *iter;
// 		if(sym[0] == 'G'){
// 			stringstream ss(sym.substr(1, sym.length()-1));
// 			ss >> num;
// 			genos.push_back(num-1);
// 		}
// 		else if(sym[0] == 'C' && sym.find_first_of("0123456789") != string::npos){
// 			stringstream ss(sym.substr(1, sym.length()-1));
// 			ss >> num;
// 			covars.push_back(num-1);
// 		}
// 	}
// 
// 	int totalInds = 0;
// 	Individual* ind;
// 	bool anyMissing;
// 	
// 	float diff=0.0, meanVal, statTotal=0.0;
// 	vector<float> statusUsed;
// 	
// 	for(unsigned int currind=0; currind < set->numInds(); currind++){
// 		ind = (*set)[currind];
// 		
// 		anyMissing = false;
// 		for(unsigned int g=0; g < genos.size(); g++){
// 			if(ind->getGenotype(genos[g]) == set->getMissingGenotype()){
// 				anyMissing = true;
// 				break;
// 			}
// 		}
// 		for(unsigned int c=0; c < covars.size(); c++){
// 			if(ind->getCovariate(covars[c]) == set->getMissingCoValue()){
// 				anyMissing = true;
// 				break;
// 			}
// 		}
// 		
// 		if(!anyMissing){
// 			totalInds++;
// 			statTotal += ind->getStatus();
// 			statusUsed.push_back(ind->getStatus());
// 		}
// 	}
// 	
// 	meanVal = statTotal / totalInds;
// 	for(vector<float>::iterator iter=statusUsed.begin(); iter != statusUsed.end();
// 		++iter){
// 		diff = diff + (*iter-meanVal) * (*iter-meanVal);
// 	}
// 	ssTotal=diff;
// 	
// 	return totalInds;
// }
