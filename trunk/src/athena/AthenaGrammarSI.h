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
//AthenaGrammarSI.h

#ifndef _AthenaGrammarSI_H
#define	_AthenaGrammarSI_H

#include <GE/GEGrammarSI.h>
#include "Structs.h"
#include "SolutionCreator.h"
#include <map>
#include <set>


///
/// Adds functionality for incorporating biofilter models into
/// initialized models.
///
class AthenaGrammarSI: public GEGrammarSI{

	public:
		void addModel(std::vector<int> indexes);
		
		bool growTree(DerivationTree &tree,const bool &growMethod,const unsigned int &maximumDepth);
		bool init(const unsigned int index=UINT_MAX);
		
		void setRestrictRule(std::string ruleString);
		
		/// Sets any codon values in the Grammar Models to reflect the replacement that should be made
		void setModelCodons(bool dummyEncoded);
		
		/// Resets iterator for models to point to first one in list
		void resetGrammarModels(){currGramModel = gramModels.begin();}
		
		/// Clears biofilter models
		void clearModels(){gramModels.clear();}
		
		/// Returns the codon value for the variable requested
		int getCodonVarValue(std::string variable, int maxCodonValue);
		
		/// Establishes variable_codon_map 
		void setVariableCodonMap();
		
		/// Takes a genome translates it using current rules and then alters the variables for new mapper rules
		void convertGenomeVariables(AthenaGrammarSI& newMapper, const GA1DArrayGenome<int> &genome);
		
		/// Determines size of genome block from starting codon passed
		int determineBlockLength(int startCodon);
		
		/// Returns codon matching
		int getMatchingCodon(std::string rule);
		
		/// Establishes current genome 
		void establishCodons(const GA1DArrayGenome<int> &genome);
		
		/// Returns string held in symbol for the codon vector
		std::string getRuleString(int codon){return codonVector[codon];}
		
		/// Return codons
		std::vector<std::string> getCodons(){return codonVector;}
		
		struct codonBlocks{
			int start,end;
		};
		
		/// Sets start symbol
		inline void setOptStartSymbol(std::string stSymb){
			optStartSymbol = stSymb;
			Symbol searchSymbol(optStartSymbol, NTSymbol);
			startOptRulePtr = findRule(searchSymbol);
		}
		
		/// Sets symbol for variable search
		inline void setVariableSymbol(std::string varSymb){
		  Symbol searchSymbol(varSymb, NTSymbol);
		  varRulePtr = findRule(searchSymbol);
		}
		
		/// Sets the set of symbols that define the optimization parameters in grammar
		inline void setOptSymbolSet(std::set<string> optSet){optSymbols = optSet;}
		
		/// Sets genotype for use in optimization
		vector<codonBlocks> setGenotypeOpt(const GA1DArrayGenome<int> &genome,
			std::set<string> compressedSet, string startSymbol, bool singelOpt=false);
		
		/// Sets genotype for use in optimization
		vector<codonBlocks> setGenotypeOpt(const GA1DArrayGenome<int> &genome);
		
		/// Returns vector of ints corresponding to block 
		vector<int> translateOptValue(symbVector& opt_symbol);
		
		/// Sets left bound on optimized values usually '(' 
		inline void setLeftOptBound(char l){leftOptBound = l;}
		
		/// Sets right bound on optimized values usually ')'
		inline void setRightOptBound(char r){rightOptBound = r;}
		
		/// Sets symbols that are used as an argument to other symbols (such as num for Concat)
		inline void setArgSymbols(std::set<string> symbSet){isArg=symbSet;}
		
		/// Constructs a reverse grammer that will give the left hand side when passed the right side
		void constructReverseGrammar();
		
		/// Builds the derivation tree and returns max depth of the tree
		int buildDerivationTree();
		
		
		/// Return a new variable
		std::string getNewVariable(int& newCodon);
		
		/// change genome variables as specified in variableMap
    void changeVariables(GA1DArrayGenome<int>& genome, map<int, SolutionCreator::TerminalInfo> variableMap);
		
	protected:
		struct GrammarModel{
			vector<int> datasetIndexes; //location in original array (corresponds to V1,V2, etc. in grammar)
			vector<int> codonValues; // contains codon values to replace the existing ones in genome with
		};
	
    int genotype2PhenotypeStepReplaceVar(stack<const Symbol*> &nonterminals,
		  int& genoIt, GA1DArrayGenome<int>& genome, map<int, SolutionCreator::TerminalInfo>& variableMap);
	
		int getMax(DerivationTree& tree);
	
		bool genotype2PhenotypeConvert(AthenaGrammarSI& newMapper, const bool buildDerivationTree=false);
	
		int genotype2PhenotypeStepConvert(AthenaGrammarSI& newMapper, stack<const Symbol*> &nonterminals,
			Genotype::iterator genoIt, bool buildDerivationTree);
		
		bool genotype2PhenotypeOpt(vector<codonBlocks>& blocks);
		
		int genotype2PhenotypeStepOpt(stack<const Symbol*> &nonterminals,
			Genotype::iterator genoIt);
		
		bool getRuleFromProduction(string& production_string, string& rule_left_side, int& codon);
		
		int fillCodons();
		
		vector<GrammarModel> gramModels;
		vector<GrammarModel>::iterator currGramModel;
		int currModelIndex, currVarIndex;
		unsigned int currModelCodonIndex;
		Rule *restrictRulePtr, *startOptRulePtr, *varRulePtr;
		char leftOptBound, rightOptBound;
		
		struct reverseRule{
			reverseRule(string rStr, int c, int nCodons){
				rule = rStr;
				codon = c;
				numCodons = nCodons;
			}
			
			reverseRule(){
				rule ="";
				codon = -1;
				numCodons = 0;
			}
			string rule;
			int codon, numCodons;
		};
		
		struct GramElement{
			string symbol;
			vector<int> codons;
			bool noNT;
		};
		
		map<string, reverseRule> reverseRules;
	 
		std::set<string> optSymbols, isArg;
		std::string optStartSymbol;
		bool optIsSingle;

		// map has key as variable name and value as codon value appropriate for it
		std::map<std::string, int> variableCodonMap;
		
		std::multimap<std::string, int>  codonMap;
		std::vector<std::string> codonVector;
	

};

#endif

