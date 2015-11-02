//NNSolutionCreatorTest.cpp

#include "../NNSolutionCreatorIncludeAll.h"
#include "../NNSolutionCreatorIncludeOnce.h"
#include <MDRFileHandler.h>
#include <OttDummyConvert.h>
#include <ScaleContinuous.h>
#include <gtest/gtest.h>

#include <iostream>
using namespace std;
using namespace data_manage;


void FillDataSet(Dataholder* holder, Dataset* set, string filename);

//{"PA","(","W","(","Concat","(","2",".","7","3",")",",","G2",")",",","W","(","Concat","(","5",".","3","6","4",")",",","G1",")",",","2",")"}

//rank=0 Old genome score is 0.217632
//{"PA","(","W","(","(","Concat","(","7",".","2","1","4",")","/","(","Concat","(","0","1",")","-","Concat","(","0","1",")",")",")",",","G13",")",",","W","(","Concat","(","1","9",".","6","7","5",")",",","G10",")",",","2",")"}
// rank=0 New genome score is 0.162289

class NNSolutionCreatorTest:public testing::Test{
  protected:
  virtual void SetUp(){
    // create symbol list describing network to test on
    string symbs[49] = {"PA","(","W","(","(","Concat","(","7",".","2","1","4",")",
      "/","(","Concat","(","0","1",")","-","Concat","(","0","1",")",")",")",",",
      "G10",")",",","W","(","Concat","(","1","9",".","6","7","5",")",",","G19",
      ")",",","2",")"};

    for(int i=0; i<49; i++){
      symbols.push_back(symbs[i]);
    }
    
    // create dataset to work on -- first need dataholder filled
    FillDataSet(&holder, &set, "sample_data/config.sim.0.1.1-DATASET.mdr");
  }  
  vector<string> symbols;
  Dataholder holder;
  Dataset set;
};

class NNSolutionCreatorIncludeAllTest:public testing::Test{
  protected:
  virtual void SetUp(){
    // create symbol list describing network to test on
    string symbs[49] = {"PA","(","W","(","(","Concat","(","7",".","2","1","4",")",
      "/","(","Concat","(","0","1",")","-","Concat","(","0","1",")",")",")",",",
      "G10",")",",","W","(","Concat","(","1","9",".","6","7","5",")",",","G19",
      ")",",","2",")"};

    for(int i=0; i<49; i++){
      symbols.push_back(symbs[i]);
    }
    
    // create dataset to work on -- first need dataholder filled
    FillDataSet(&holder, &set, "sample_data/config.sim.0.1.1-DATASET.mdr");
  }  
  vector<string> symbols;
  Dataholder holder;
  Dataset set;
};


class NNSolutionCreatorIncludeOnceTest:public testing::Test{
  protected:
  virtual void SetUp(){
    // create symbol list describing network to test on
    string symbs[49] = {"PA","(","W","(","(","Concat","(","7",".","2","1","4",")",
      "/","(","Concat","(","0","1",")","-","Concat","(","0","1",")",")",")",",",
      "G10",")",",","W","(","Concat","(","1","9",".","6","7","5",")",",","G19",
      ")",",","2",")"};

    for(int i=0; i<49; i++){
      symbols.push_back(symbs[i]);
    }
    
    // create dataset to work on -- first need dataholder filled
    FillDataSet(&holder, &set, "sample_data/config.sim.0.1.1-DATASET.mdr");
  }  
  vector<string> symbols;
  Dataholder holder;
  Dataset set;
};

// fills dataholder and dataset for use in tests
void FillDataSet(Dataholder* holder, Dataset* set, string filename){
  MDRFileHandler mdrreader;
  mdrreader.parse_file(filename, holder, -1, -9999, false);
  OttDummyConvert ott;
  ott.convert_genotypes(holder);
  
  // normalize status values
  ScaleContinuous scaler;
  scaler.adjust_status(holder);

  for(unsigned int i=0; i < holder->num_inds(); i++){
    set->add_ind(holder->get_ind(i));
  }
  
  set->set_missing_genotype(holder->get_missing_genotype());
  set->set_missing_covalue(holder->get_missing_covalue());  
  
// output set for script testing
// for(size_t i=0; i < holder->num_inds(); i++){
//   Individual * ind = holder->get_ind(i);
//   cout << ind->get_status();
//   for(unsigned int j=0; j<ind->num_genotypes(); j++){
//     cout << " " << ind->get_genotype(j);
//   }
//   cout << endl;
// }
// exit(1);
  
}

// check establish solution
TEST_F(NNSolutionCreatorTest, EstablishSolution){
  NNSolutionCreator creator;
  creator.establish_solution(symbols, &set);
  creator.set_calculator("RSQUARED");
  EXPECT_EQ(2, int(creator.get_num_genes()));
// cout << "first gene is " << creator.getGeneIndexes().at(0) << endl;
// cout << "second snp is " << creator.getGeneIndexes().at(1) << endl;
  EXPECT_EQ(9, int(creator.getGeneIndexes().at(0))) << "First snp in network should have index 9";
  EXPECT_EQ(18, int(creator.getGeneIndexes().at(1))) << "Second snp in network should have index of 12";
  EXPECT_EQ(0, int(creator.get_num_covars()));
}

// check evaluate
TEST_F(NNSolutionCreatorTest, Evaluate){
  NNSolutionCreator creator;
  creator.set_calculator("RSQUARED");
  creator.establish_solution(symbols, &set);
  float score = creator.evaluate(&set);
//  cout << creator.evaluate(&set) << endl;
// cout << "original score " << score << endl;
  EXPECT_GT(0.0001, fabs(0.117177 - score));
// cout << "num inds eval =" << creator.getNumIndsEvaluated() << endl;
  EXPECT_EQ(1000, creator.getNumIndsEvaluated());
}

// check backpropagation procedure for the optimization
TEST_F(NNSolutionCreatorTest, OptimizeSolution){
  NNSolutionCreator creator;
  creator.set_calculator("RSQUARED");
  int epochs = creator.optimizeSolution(symbols, &set);
  EXPECT_EQ(2, epochs) << "Should only run 2 epochs on this network";
cout << "Optimized score=" << creator.getOptimizedScore() << endl;
  EXPECT_GT(0.000001, fabs(0.0938615 - creator.getOptimizedScore()));
cout << "Optimized score=" << creator.getOptimizedScore() << endl;
cout << "Number of epochs=" << epochs << endl;
  // retrieve new weights and check
  vector<float> optimizedWeights = creator.getOptimizedValues();
  EXPECT_EQ(2, int(optimizedWeights.size()));
cout << optimizedWeights[0] << " " << optimizedWeights[1] << endl;
  float weights[] = {0.165306, 19.67};
  for(int i=0; i<2; i++){
    EXPECT_GT(0.00001, fabs(weights[i] - optimizedWeights[i]));
  }
}


// check optimization constants
TEST_F(NNSolutionCreatorTest, OptimizationConstants){
  NNSolutionCreator creator;
  EXPECT_EQ("<cop>", creator.getStartOptSymbol());
  char left = creator.getLeftOptBound();
  EXPECT_EQ('(', left);
  char right = creator.getRightOptBound();
  EXPECT_EQ(')', right);
  std::set<string> optargsymbols = creator.getOptArgSymbols();
  EXPECT_EQ("(<num>)", *(optargsymbols.begin()));
  std::set<string> optSymbolsIncluded = creator.getOptIncluded();
  EXPECT_EQ(22, int(optSymbolsIncluded.size()));
  
  string optsymbols[] = {"<cop>","<op>","<Concat>","(",")", "<num>","Concat","1","2","3","4",
    "5","6","7","8","9","0","<dig>","+","-","*","/"};
  
  for(int i=0; i<22; i++){
    std::set<string>::iterator iter = optSymbolsIncluded.find(optsymbols[i]);
//    EXPECT_NE(optSymbolsIncluded.end(), iter);// << "optsymbol " << optsymbols[i]
//      << " not in set";
    EXPECT_FALSE(optSymbolsIncluded.end() == optSymbolsIncluded.find(optsymbols[i])) << "optsymbol " << optsymbols[i]
      << " not in set or argument symbols";
  } 
}

// check the Include All type of network 
TEST_F(NNSolutionCreatorIncludeAllTest, Evaluate){
  NNSolutionCreatorIncludeAll creator;
  creator.set_calculator("RSQUARED");
  vector<string> required;
  required.push_back("G10");
  required.push_back("G19");
  creator.restrict(required);
  creator.establish_solution(symbols, &set);
  EXPECT_EQ(2, int(creator.get_num_genes()));
// cout << "first gene is " << creator.getGeneIndexes().at(0) << endl;
// cout << "second snp is " << creator.getGeneIndexes().at(1) << endl;
  EXPECT_EQ(9, int(creator.getGeneIndexes().at(0))) << "First snp in network should have index 9";
  EXPECT_EQ(18, int(creator.getGeneIndexes().at(1))) << "Second snp in network should have index of 12";
  EXPECT_EQ(0, int(creator.get_num_covars()));  
  float score = creator.evaluate(&set);
  EXPECT_GT(0.0001, fabs(0.117177 - score));
  
  // create with different list of required variables
  NNSolutionCreatorIncludeAll creator2;
  creator2.set_calculator("RSQUARED");
  vector<string> required2;
  required2.push_back("G10");
  required2.push_back("G18");
  creator2.restrict(required2);
  creator2.establish_solution(symbols, &set);
  EXPECT_EQ(2, int(creator2.get_num_genes()));
// cout << "first gene is " << creator.getGeneIndexes().at(0) << endl;
// cout << "second snp is " << creator.getGeneIndexes().at(1) << endl;
  EXPECT_EQ(9, int(creator2.getGeneIndexes().at(0))) << "First snp in network should have index 9";
  EXPECT_EQ(18, int(creator2.getGeneIndexes().at(1))) << "Second snp in network should have index of 12";
  EXPECT_EQ(0, int(creator2.get_num_covars()));  
  score = creator2.evaluate(&set);
  EXPECT_GT(0.0001, fabs(1e+06 - score)) << "Network should give worst score because it is missing G19 in network";
}



// check the Include Once type of network 
// where each variable should appear once and only once in the network
TEST_F(NNSolutionCreatorIncludeOnceTest, Evaluate){
  NNSolutionCreatorIncludeOnce creator;
  creator.set_calculator("RSQUARED");
  vector<string> required;
  required.push_back("G10");
  required.push_back("G19");
  creator.restrict(required);
  creator.establish_solution(symbols, &set);
  EXPECT_EQ(2, int(creator.get_num_genes()));
// cout << "first snp is " << creator.getGeneIndexes().at(0) << endl;
// cout << "second snp is " << creator.getGeneIndexes().at(1) << endl;
  EXPECT_EQ(9, int(creator.getGeneIndexes().at(0))) << "First snp in network should have index 9";
  EXPECT_EQ(18, int(creator.getGeneIndexes().at(1))) << "Second snp in network should have index of 12";
  EXPECT_EQ(0, int(creator.get_num_covars()));  
  float score = creator.evaluate(&set);
  EXPECT_GT(0.0001, fabs(0.117177 - score));

  // create symbol list describing network to test on for failure
  string symbs2[49] = {"PA","(","W","(","(","Concat","(","7",".","2","1","4",")",
    "/","(","Concat","(","0","1",")","-","Concat","(","0","1",")",")",")",",",
    "G19",")",",","W","(","Concat","(","1","9",".","6","7","5",")",",","G19",
    ")",",","2",")"};

  vector<string> symbols2;
  for(int i=0; i<49; i++){
    symbols2.push_back(symbs2[i]);
  }
  
  // create with different list of required variables
  NNSolutionCreatorIncludeOnce creator2;
  creator2.set_calculator("RSQUARED");
  vector<string> required2;
  required2.push_back("G19");
  required2.push_back("G18");
  creator2.restrict(required2);
  creator2.establish_solution(symbols2, &set);
  EXPECT_EQ(2, int(creator2.get_num_genes()));
// cout << "first gene is " << creator.getGeneIndexes().at(0) << endl;
// cout << "second snp is " << creator.getGeneIndexes().at(1) << endl;
  EXPECT_EQ(18, int(creator2.getGeneIndexes().at(0))) << "First snp in network should have index 9";
  EXPECT_EQ(18, int(creator2.getGeneIndexes().at(1))) << "Second snp in network should have index of 12";
  EXPECT_EQ(0, int(creator2.get_num_covars()));  
  score = creator2.evaluate(&set);
  EXPECT_GT(0.0001, fabs(1e+06 - score)) << "Network should give worst score because it is missing G19 in network";
}

