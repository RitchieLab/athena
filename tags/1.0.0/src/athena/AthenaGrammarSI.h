//AthenaGrammarSI.h

#ifndef _AthenaGrammarSI_H
#define	_AthenaGrammarSI_H

#include <GE/GEGrammarSI.h>
#include "Structs.h"
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
    
    void setRestrictRule(std::string rule_string);
    
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
    std::string getRuleString(int codon){return codon_vector[codon];}
    
    /// Return codons
    std::vector<std::string> getCodons(){return codon_vector;}
    
    struct codonBlocks{
      int start,end;
    };
    
    /// Sets start symbol
    inline void setOptStartSymbol(std::string st_symb){
      opt_start_symbol = st_symb;
      Symbol searchSymbol(opt_start_symbol, NTSymbol);
      startOptRulePtr = findRule(searchSymbol);
    }
    
    /// Sets the set of symbols that define the optimization parameters in grammar
    inline void setOptSymbolSet(std::set<string> opt_set){opt_symbols = opt_set;}
    
    /// Sets genotype for use in optimization
    vector<codonBlocks> setGenotypeOpt(const GA1DArrayGenome<int> &genome,
      std::set<string> compressed_set, string startSymbol);
    
    /// Sets genotype for use in optimization
    vector<codonBlocks> setGenotypeOpt(const GA1DArrayGenome<int> &genome);
    
    /// Returns vector of ints corresponding to block 
    vector<int> translateOptValue(symbVector& opt_symbol);
    
    /// Sets left bound on optimized values usually '(' 
    inline void setLeftOptBound(char l){left_opt_bound = l;}
    
    /// Sets right bound on optimized values usually ')'
    inline void setRightOptBound(char r){right_opt_bound = r;}
    
    /// Sets symbols that are used as an argument to other symbols (such as num for Concat)
    inline void setArgSymbols(std::set<string> symb_set){is_arg=symb_set;}
    
    /// Constructs a reverse grammer that will give the left hand side when passed the right side
    void constructReverseGrammar();
    
  private:
    struct GrammarModel{
      vector<int> dataset_indexes; //location in original array (corresponds to V1,V2, etc. in grammar)
      vector<int> codon_values; // contains codon values to replace the existing ones in genome with
    };
  
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
    Rule *restrictRulePtr, *startOptRulePtr;
    char left_opt_bound, right_opt_bound;
    
    struct reverseRule{
      reverseRule(string rStr, int c, int ncodons){
        rule = rStr;
        codon = c;
        numcodons = ncodons;
      }
      
      reverseRule(){
        rule ="";
        codon = -1;
        numcodons = 0;
      }
      string rule;
      int codon, numcodons;
    };
    
    struct gramelement{
      string symbol;
      vector<int> codons;
      bool noNT;
    };
    
    map<string, reverseRule> reverseRules;
   
    std::set<string> opt_symbols, is_arg;
    std::string opt_start_symbol;

    // map has key as variable name and value as codon value appropriate for it
    std::map<std::string, int> variable_codon_map;
    
    std::multimap<std::string, int>  codon_map;
    std::vector<std::string> codon_vector;

};

#endif

