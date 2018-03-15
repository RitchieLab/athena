//TerminalsTest.cpp

#include "../Terminals.h"
#include <gtest/gtest.h>
#include <iostream>
#include <deque>
using namespace std;

// Test constructors
TEST(ConstantTerminalTest, ConstructorTest){
  Constant const1("8.35", 0);
  deque<float> args;
  EXPECT_GT(0.00001, fabs(8.35 - const1.evaluate(args))) << "Standard constructor for Constant should have value of 8.35";
  EXPECT_EQ(0, const1.get_priority()) << "Constant priority should be 0";
  Constant const2(-0.5);
  EXPECT_EQ(-0.5, const2.evaluate(args)) << "Constant constructor from float";
}

// Establish individuals for use in testing GenotypeTerms and ContinVariable terminals
class VariableTest:public testing::Test{
  protected:
  virtual void SetUp(){
    ind1.set_status(1);
    ind1.add_genotype(0);
    ind1.add_genotype(1);
    ind1.add_genotype(2);
    ind2.set_status(0.5);    
    ind2.add_genotype(2);
    ind2.add_genotype(1);
    ind2.add_genotype(0);
    ind3 = ind2;
    ind3.add_covariate(0.2);
    ind3.add_covariate(0.1);
    ind3.add_covariate(-1.3);
    ind2.add_covariate(1.8);
    ind2.add_covariate(99.9);
    ind2.add_covariate(0);
    
  }
  Individual ind1, ind2, ind3;
};

// Test genotype term
TEST_F(VariableTest, GenotypeTermTest){
  GenotypeTerm::set_ind(&ind1);
  GenotypeTerm g0("G1", 1);
  GenotypeTerm g1("G2", 2);
  GenotypeTerm g2("G3", 3);
  deque<float> args;
  EXPECT_EQ(0, g0.evaluate(args)) << "First genotype for first ind";
  EXPECT_EQ(1, g1.evaluate(args)) << "Second genotype for first ind";
  EXPECT_EQ(2, g2.evaluate(args)) << "Third genotype for first ind";
  
  GenotypeTerm::set_ind(&ind2);
  EXPECT_EQ(2, g0.evaluate(args)) << "First genotype for first ind";
  EXPECT_EQ(1, g1.evaluate(args)) << "Second genotype for first ind";
  EXPECT_EQ(0, g2.evaluate(args)) << "Third genotype for first ind";
}

// Test genotype term
TEST_F(VariableTest, ContinVariableTest){
  ContinVariable::set_ind(&ind2);
  ContinVariable c0("C1", 1);
  ContinVariable c1("C2", 2);
  ContinVariable c2("C3", 3);
  deque<float> args;
  EXPECT_FLOAT_EQ(1.8,c0.evaluate(args)) << "First covar for first ind";
  EXPECT_FLOAT_EQ(99.9,c1.evaluate(args)) << "Second covar for first ind";
  EXPECT_FLOAT_EQ(0,c2.evaluate(args)) << "Third covar for first ind";
  
  ContinVariable::set_ind(&ind3);
  EXPECT_FLOAT_EQ(0.2,c0.evaluate(args)) << "First covar for second ind";
  EXPECT_FLOAT_EQ(0.1,c1.evaluate(args)) << "Second covar for second ind";
  EXPECT_FLOAT_EQ(-1.3,c2.evaluate(args)) << "Third covar for second ind";  
}

// Test basic addition
TEST(AdditionTerminal, Evaluation){
  Addition add("+", 2);
  deque<float> args;
  args.push_back(2);
  args.push_back(1);

  EXPECT_EQ(3, add.evaluate(args));
  args[1] = -1;
  EXPECT_EQ(1, add.evaluate(args));
  
  args.clear();
  args.push_back(-9.2);
  args.push_back(0);
  
//   EXPECT_EQ(0.00001, fabs(-9.2-add.evaluate(args)));
  EXPECT_FLOAT_EQ(-9.2, add.evaluate(args));
  
}

// Test subtraction
TEST(SubtractionTerminal, Evaluation){
  Subtraction sub("-", 2);
  deque<float> args;
  args.push_back(2);
  args.push_back(1);
  EXPECT_FLOAT_EQ(1,sub.evaluate(args));
  
  args[0] = -1;
  args[1] = 400;
  EXPECT_FLOAT_EQ(-401, sub.evaluate(args));
  
  args[0] = 0.22;
  args[1] = -0.1;
  EXPECT_FLOAT_EQ(0.32, sub.evaluate(args));
  
}


// Test subtraction
TEST(MultiplicationTerminal, Evaluation){
  Multiplication mult("*", 2);
  deque<float> args;
  args.push_back(2);
  args.push_back(1);
  EXPECT_FLOAT_EQ(2,mult.evaluate(args));
  
  args[0] = -1;
  args[1] = 400;
  EXPECT_FLOAT_EQ(-400, mult.evaluate(args));
  
  args[0] = 0.22;
  args[1] = -2;
  EXPECT_FLOAT_EQ(-0.44, mult.evaluate(args));
  
}

// Test Division
TEST(DivisionTerminal, Evaluation){
  Division div("/", 2);
  deque<float> args;
  args.push_back(20);
  args.push_back(10);
  EXPECT_FLOAT_EQ(2,div.evaluate(args));
  
  // when trying to divide by zero should return 1.0
  args[0] = 50;
  args[1] = 0;
  EXPECT_FLOAT_EQ(1.0, div.evaluate(args));
  
  args[0] = 0.22;
  args[1] = -2;
  EXPECT_FLOAT_EQ(-0.11, div.evaluate(args)); 
}

// Test node Addition
TEST(pAddTerminal, Evaluation){
  pAdd padd("PADD", -1);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(0.99004817, padd.evaluate(args));
  args.push_back(-5);
  EXPECT_FLOAT_EQ(0.40131232, padd.evaluate(args));
  args.clear();
  args.push_back(0);
  args.push_back(0);
  EXPECT_FLOAT_EQ( 0.5, padd.evaluate(args));
}


// Test node subtraction
TEST(pSubTerminal, Evaluation){
  pSub psub("PSUB", -1);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(0.083172686, psub.evaluate(args));
  args.push_back(-5);
  EXPECT_FLOAT_EQ( 0.93086159, psub.evaluate(args));
  args.clear();
  args.push_back(0);
  args.push_back(0);
  EXPECT_FLOAT_EQ( 0.5, psub.evaluate(args));
}

// Test node Multiplication
TEST(pMultTerminal, Evaluation){
  pMult pmult("PMULT", -1);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ( 0.97916365, pmult.evaluate(args));
  args.push_back(-5);
  EXPECT_FLOAT_EQ(4.363462e-09, pmult.evaluate(args));
  args.clear();
  args.push_back(1);
  args.push_back(0.245);
  args.push_back(0);
  EXPECT_FLOAT_EQ( 0.5, pmult.evaluate(args));
}

// Test node Division
TEST(pDivTerminal, Evaluation){
  pDiv pdiv("PDIV", -1);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ( 0.57793099, pdiv.evaluate(args));
  args.push_back(-5);
  EXPECT_FLOAT_EQ(0.4842909, pdiv.evaluate(args));
  args.clear();
  // check division by zero
  args.push_back(1);
  args.push_back(0.245);
  args.push_back(0);
  EXPECT_FLOAT_EQ(1, pdiv.evaluate(args));
}


// Test boolean And operator
TEST(pAndTerminal, Evaluation){
  pAnd pand("PAND", -1);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(1, pand.evaluate(args));
  args.push_back(-5);
  EXPECT_FLOAT_EQ(0, pand.evaluate(args));
  args.clear();
  // check with zero
  args.push_back(1);
  args.push_back(0.245);
  args.push_back(0);
  EXPECT_FLOAT_EQ(0, pand.evaluate(args));
}


// Test boolean pNand operator
TEST(pNandTerminalTest, Evaluation){
  pNand pnand("PNAND", -1);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(0, pnand.evaluate(args));
  args.push_back(-5);
  EXPECT_FLOAT_EQ(0, pnand.evaluate(args));
  args.clear();
  // check with zero
  args.push_back(-1);
  args.push_back(-0.245);
  args.push_back(0);
  EXPECT_FLOAT_EQ(1, pnand.evaluate(args));
}


// Test boolean pOr operator
TEST(pOrTerminalTest, Evaluation){
  pOr por("POR", -1);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(1, por.evaluate(args));
  args.push_back(-5);
  EXPECT_FLOAT_EQ(1, por.evaluate(args));
  args.clear();
  // check with zero
  args.push_back(-1);
  args.push_back(-0.245);
  args.push_back(0);
  EXPECT_FLOAT_EQ(0, por.evaluate(args));
}

// Test boolean pNor operator
TEST(pNorTerminalTest, Evaluation){
  pNor pnor("PNOR", -1);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(0, pnor.evaluate(args));
  args.push_back(-5);
  EXPECT_FLOAT_EQ(0, pnor.evaluate(args));
  args.clear();
  // check with zero
  args.push_back(-1);
  args.push_back(-0.245);
  args.push_back(0);
  EXPECT_FLOAT_EQ(1, pnor.evaluate(args));
}

// Test boolean pXor operator
TEST(pXorTerminalTest, Evaluation){
  pXor pxor("PXOR", -1);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(0, pxor.evaluate(args));
  args.push_back(-5);
  EXPECT_FLOAT_EQ(0, pxor.evaluate(args));
  args.clear();
  // check with zero
  args.push_back(0);
  args.push_back(1);
  EXPECT_FLOAT_EQ(1, pxor.evaluate(args));
}

// Test weight class -- basically same as multiplication
TEST(WeightTerminal, Evaluation){
  Weight wt("W", 2);

  deque<float> args;
  args.push_back(2);
  args.push_back(1);
  EXPECT_FLOAT_EQ(2,wt.evaluate(args));
  
  args[0] = -1;
  args[1] = 400;
  EXPECT_FLOAT_EQ(-400, wt.evaluate(args));
  
  args[0] = 0.22;
  args[1] = -2;
  EXPECT_FLOAT_EQ(-0.44, wt.evaluate(args));
}

// Test concatenation class
TEST(ConCatTerminalTest, Evaluation){
  ConCat concat("Concat", -1);
  
  deque<float> args;
  args.push_back(5);
  args.push_back(-2);
  args.push_back(3);
  
  EXPECT_FLOAT_EQ(5.3, concat.evaluate(args));
  
  args.clear();
  args.push_back(-2);
  args.push_back(1);
  args.push_back(7);
  EXPECT_FLOAT_EQ(0.17, concat.evaluate(args));
  
  args.clear();
  args.push_back(-1);
  EXPECT_FLOAT_EQ(-1, concat.evaluate(args));
  
}

// Test xor class
TEST(XorTerminalTest, Evaluation){
  Xor XOR("XOR", 2);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(0, XOR.evaluate(args));
  
  args[1] = -1;
  EXPECT_FLOAT_EQ(1, XOR.evaluate(args));
  args.clear();
  // check with zero
  args.push_back(0);
  args.push_back(1);
  EXPECT_FLOAT_EQ(1, XOR.evaluate(args));
}


// Test or class
TEST(OrTerminalTest, Evaluation){
  Or OR("OR", 2);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(1, OR.evaluate(args));
  
  args[1] = -1;
  EXPECT_FLOAT_EQ(1, OR.evaluate(args));
  args.clear();
  // check with zero
  args.push_back(0);
  args.push_back(1);
  EXPECT_FLOAT_EQ(1, OR.evaluate(args));
}


// Test Nor class
TEST(NorTerminalTest, Evaluation){
  Nor nor("NOR", 2);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(0, nor.evaluate(args));
  
  args[1] = -1;
  EXPECT_FLOAT_EQ(0, nor.evaluate(args));
  args.clear();
  // check with zero
  args.push_back(0);
  args.push_back(-1);
  EXPECT_FLOAT_EQ(1, nor.evaluate(args));
}

// Test And class
TEST(AndTerminalTest, Evaluation){
  And AND("NOR", 2);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(1, AND.evaluate(args));
  
  args[1] = -1;
  EXPECT_FLOAT_EQ(0, AND.evaluate(args));
  args.clear();
  // check with zero
  args.push_back(0);
  args.push_back(-1);
  EXPECT_FLOAT_EQ(0, AND.evaluate(args));
}

// Test And class
TEST(NandTerminalTest, Evaluation){
  Nand NAND("NAND", 2);
  // first argument is number of terms to evaluate
  deque<float> args;
  args.push_back(1.1);
  args.push_back(3.5);
 
  EXPECT_FLOAT_EQ(0, NAND.evaluate(args));
  
  args[1] = -1;
  EXPECT_FLOAT_EQ(0, NAND.evaluate(args));
  args.clear();
  // check with zero
  args.push_back(0);
  args.push_back(-1);
  EXPECT_FLOAT_EQ(1, NAND.evaluate(args));
}


