// santafe-gcc.cpp -*- C++ -*-
#ifndef _SANTAFE_GCC_CPP_
#define _SANTAFE_GCC_CPP_

#include<cstdio>
#include<cstdlib>
#include<cstring>

#include<ga/ga.h>
#include<GE/ge.h>
#include"GEListGenome.h"

#include <map>

using namespace std;

// GE mapper, declared in main.cpp
extern GEGrammarSI mapper;

void expand_variables(Rule& varRule);
string create_variable_symbol(string vPrefix, int num);
void double_genotype_grammar(Rule& varRule);

// Numer of Objective function calls, declared in main.cpp
extern int obj_calls;

//Buffers for start and end code of the individual.c file
string SFstart;
string SFend;

//Loads the start and end code required to compile phenotypes,
//and initialises the mapper
void app_init(unsigned int wrappingEvents,string grammarFile){
	char buffer[1000];
	FILE *c_code;

	/* Read buffer code files */
	SFstart="";
	if(!(c_code=fopen("santafestart.c","r"))){
		cerr << "Could not read santafestart.c\n";
		cerr << "Execution aborted.\n";
		exit(0);
	}
	while(fgets(buffer,1000,c_code))
		SFstart+=buffer;
	fclose(c_code);
	SFend="";
	if(!(c_code=fopen("santafeend.c","r"))){
		cerr << "Could not read santafeend.c\n";
		cerr << "Execution aborted.\n";
		exit(0);
	}
	while(fgets(buffer,1000,c_code))
		SFend+=buffer;
	fclose(c_code);

	/* Init GE mapper */
	//Set maximum number of wrapping events per mapping
	mapper.setMaxWraps(wrappingEvents);
	//Load grammar
	if (!mapper.readBNFFile(grammarFile)){
                cerr << "Could not read " << grammarFile << "\n";
		cerr << "Execution aborted.\n";
		exit(0);
	}
	
	cout << "read map file " << grammarFile << endl;
	cout << "Total # of rules = " << mapper.size() << endl;
	
	
// 	for(uint i=0; i<mapper.size(); i++){
// 	  cout << "rule #" << i << " " << *(mapper[i].lhs.back()) << endl;
// 	  
// 	  // show productions of the rules
// 	  for(uint j=0; j<mapper[i].size(); j++){
// 	    cout << "j=" << j;
// 	    for(uint k=0; k<mapper[i][j].size(); k++){
// 	      cout << " -- " << *(mapper[i][j][k]) << " ";
// 	    }
// 	    cout << endl;
// 	  }
// 	  
// 	  
// 	}

/// /   string v = "<v>";
//   Symbol var_symbol = v;	
// 	Rule * variableRule = mapper.findRule(var_symbol);
// 	
// 	cout << "variableRule=" << variableRule << " NULL=" << NULL << endl;
// 	exit(1);
	
	string v = "<v>";
	// find rule with <v>
	uint rule_num;
	for(uint i=0; i<mapper.size(); i++){
	  if(v.compare(*(mapper[i].lhs.back())) == 0){
	    rule_num = i;
	  }
	}
	cout << "rule_num=" << rule_num << endl;
cout << "before expanded vars" << endl;
	expand_variables(mapper[rule_num]);

cout << "after expanded vars" << endl;

	  // show productions of the rules
	  for(uint j=0; j<mapper[rule_num].size(); j++){
	    cout << "j=" << j;
	    for(uint k=0; k<mapper[rule_num][j].size(); k++){
	      cout << " -- " << *(mapper[rule_num][j][k]) << " ";
	    }
	    cout << endl;
	  }

  double_genotype_grammar(mapper[rule_num]);
 cout << "after doubled genotypes" << endl;
	  // show productions of the rules
	  for(uint j=0; j<mapper[rule_num].size(); j++){
	    cout << "j=" << j;
	    for(uint k=0; k<mapper[rule_num][j].size(); k++){
	      cout << " -- " << *(mapper[rule_num][j][k]) << " ";
	    }
	    cout << endl;
	  }  
  
  
// 	for(uint i=0; i<mapper.size(); i++){
// 	  cout << "rule #" << i << " " << *(mapper[i].lhs.back()) << endl;
// 	  
// 	  // show productions of the rules
// 	  for(uint j=0; j<mapper[i].size(); j++){
// 	    cout << "j=" << j;
// 	    for(uint k=0; k<mapper[i][j].size(); k++){
// 	      cout << " -- " << *(mapper[i][j][k]) << " ";
// 	    }
// 	    cout << endl;
// 	  }
// 	  
// 	  
// 	}
	
	// get rule for Variables
// 	Symbol searchSymbol("V");
// 	Rule* var_rule = mapper.findRule(searchSymbol);
// 	cout << "rule productions" << endl;
// 	for(uint i=0; i<var_rule->size(); i++){
// 	  cout << "i=" << i << " prod=";
// 	  for(uint j=0; j<(*var_rule)[i].size(); j++){
// 	    cout << (*var_rule)[i][j] << " ";
// 	  }
// 	  cout << endl;
// 	}
	
  }
  
  // should be in format V1-20 or so
//   string start, end;
//   std::size_t dash_pos = symbol.find("-");
//   start = symbol.substr(1, dash_pos-1);
//   end = symbol.substr(dash_pos+1, symbol.size() - dash_pos);
//   
//   cout << "variable splitting method start=>" << start << " end=>" << end << "<== " << endl;
//   stringstream ss;
//   int begin, final;
//   ss.str(start);
//   ss >> begin;
// //   ss.str("");
// // 
// //   ss.str(end);
//   stringstream ssend(end);
// cout << "setting ss to be ==>" << ss.str() << "<== " << endl;
//   ssend >> final;
//   string newSymbol;
// cout << "begin=" << begin << " final=" << final << endl;
//   for(int i=begin; i<=final; i++){
// //     ss.str("");
//     Production prod;
//     stringstream stream;
//     stream << i;
//     newSymbol = "V" + stream.str();
// cout << "newSymbol in multvar=" << newSymbol << endl;
//     prod.push_back(new Symbol(newSymbol));
//     currentRule->push_back(prod);
//   }

// Expands the variables from shorthand in the rules 
void expand_variables(Rule& varRule){

  string start, end, varType, newSymbol;
  int begin, final;
  
  // loop through and find any rules 
  for(uint p=0; p < varRule.size(); p++){
  
    string symbol = *(varRule[p][0]);
    std::size_t dash_pos = symbol.find("-");
    if(dash_pos != string::npos){
cout << "found " << symbol << endl;

      stringstream ss(symbol);
      varType = symbol[0];
      start = symbol.substr(1, dash_pos-1);
      end = symbol.substr(dash_pos+1, symbol.size() - dash_pos);
cout << "start=" << start << " end=" << end << endl;     
      stringstream ss_start(start);
      stringstream ss_end(end);
      ss_start >> begin;
      ss_end >> final;
cout << "begin=" << begin << endl;
cout << "final=" << final << endl;

      newSymbol = create_variable_symbol(varType, begin);
//       Production product;
//       product.push_back(new Symbol(newSymbol));
cout << "p=" << p << endl;
cout << *(varRule[p][0]) << endl;
//       varRule[p] = product;
      *(varRule[p][0]) = newSymbol;
cout << *(varRule[p][0]) << endl;    
      for(uint i=begin+1; i<=final; i++){
        newSymbol = create_variable_symbol(varType, i);
        Production prod;
        prod.push_back(new Symbol(newSymbol));
        varRule.push_back(prod);
      }
      
    }

  }
 
 for(uint j=0; j<varRule.size(); j++){
   cout << "j=" << j;
  for(uint k=0; k<varRule[j].size(); k++){
    cout << " -- " << *(varRule[j][k]) << " ";
  }
  cout << endl;
} 
 
 
}

// doubles all genotype variables for use with ott dummy
// only ones included in the list are doubled so need to check for 
// duplicates
void double_genotype_grammar(Rule& varRule){
  
  // create first pass list
  vector<int> genos_included;
  map<int, bool> geno_map;
  
  for(uint i=0; i<varRule.size(); i++){
    if((*varRule[i][0])[0] == 'V'){  // change to G for new version
      string num  = (*varRule[i][0]).substr(1, (*varRule[i][0]).size()-1);
      stringstream ss(num);
      int geno;
      ss >> geno;
      genos_included.push_back(geno);
      geno_map[geno] = true;
    }
  }
  
  string symbol;
  vector<int> to_include(2,0);
  
  // now need to add the missing genotypes to list
  for(uint i=0; i<genos_included.size(); i++){
    to_include[0] = genos_included[i] *2 -1;
    to_include[1] = genos_included[i] *2; 
    
    for(uint j=0; j<2; j++){
      // check to see if the first geno is included
      if(geno_map.find(to_include[j]) == geno_map.end()){
        stringstream ss;
        ss << to_include[j];
        symbol = 'V' + ss.str();
        Production newprod;
        newprod.push_back(new Symbol(symbol));
        varRule.push_back(newprod);
        geno_map[to_include[j]] = true;
      }  
    }
    
  }
  
  
}

// creates symbol for storing in rule
string create_variable_symbol(string vPrefix, int num){
  stringstream ss;
  ss << num;
  return vPrefix + ss.str();
}
  

//Attributes a fitness score to a genome
float objfunc(GAGenome &g){
	obj_calls++;
	GEListGenome &genome = static_cast<GEListGenome&>(g);
	//Assign genotype to mapper
	mapper.setGenotype(genome);
	//Grab phenotype
        Phenotype const *phenotype=mapper.getPhenotype();
	if(phenotype->getValid()){
cout << phenotype->getString() << endl;
for(uint i=0; i<phenotype->size(); i++){
  cout  << *(*phenotype)[i] << " -- " << endl;
}
exit(1);
		FILE *file;
		int fitness;
		//Create output file
		if(!(file=fopen("individual.c","w"))){
			cerr << "Could not open individual.c.\n";
			cerr << "Execution aborted.\n";
			exit(0);
		}
		//Write start buffer to file
		fprintf(file,"%s",SFstart.c_str());
		//Write phenotype code to file
		fprintf(file,"%s",phenotype->getString().c_str());
		//Write end buffer to file
		fprintf(file,"%s",SFend.c_str());
		fclose(file);

#if defined GECART_USE_TCC_COMP

		//Compile and execute file with TCC
		if(system("tcc individual.c tcc_GEant.o\
				tcc_GEtrail.o && ./a.out > result")==-1){
			cerr << "Compilation or execution failed.\n";
			cerr << "Execution aborted.\n";
			exit(0);
		}
#else
		//Compile and execute file with GCC
		if(system("gcc -pipe individual.c gcc_GEant.o\
				gcc_GEtrail.o && ./a.out > result")==-1){
			cerr << "Compilation or execution failed.\n";
			cerr << "Execution aborted.\n";
			exit(0);
		}
#endif

		//Open result file, containing fitness score
		if(!(file=fopen("result","r"))){
			cerr << "Could not open result file.\n";
			cerr << "Execution aborted.\n";
			exit(0);
		}
		fscanf(file,"%d",&fitness);
		fclose(file);
		// Set effective size of genome
		genome.setEffectiveSize(mapper.getGenotype()->getEffectiveSize());
		return fitness;
	}
	else
		return 0;
  }

//Print an individual to stdout
void print_individual(const GAGenome &g){
	GAListGenome<unsigned char> &genome =
		(GAListGenome<unsigned char> &) g;
	//Assign genotype to mapper
	mapper.setGenotype(genome);
	//Print phenotype
	cout << *(mapper.getPhenotype());
	cout << endl;
	cout << "Genotype = " << *mapper.getGenotype() << "\n";
	cout << "Total length     = "
		<< mapper.getGenotype()->size() << "\n";
	cout << "Effective length = "
		<< mapper.getGenotype()->getEffectiveSize() << "\n";
	cout << "Wrapping events = "
		<< mapper.getGenotype()->getWraps() << "\n";
}

#endif

