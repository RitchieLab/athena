# Location of original source code
USER_DIR = ..

# Location of stat calculator header
CALC_DIR = ../../stat_calc

# Location/name of calc library
LIB_DIR = ../../lib

GADIR = ../../galib247
GEDIR = ../../libGE-0.26

GALIBLIB = ../../galib247/ga
GELIBLIB = ../../libGE-0.26

# Point to googletest root directory
GTEST_DIR = ../../gtest-1.5.0

# Point to data_manage headers
DATA_MANAGE_DIR = ../../data_manage

# Point to annie backprop library
ANNIE_DIR = ../../annie

# Flags for preprocessor
CPPFLAGS += -I$(GTEST_DIR) -I$(GTEST_DIR)/include -I$(DATA_MANAGE_DIR) \
 -I$(ANNIE_DIR) -I$(GADIR) -I$(GEDIR)

# Flags for the compiler
CXXFLAGS += -g -Wall 

# Flags for linking
LIBS = -L$(LIB_DIR) -L$(GALIBLIB) -L$(GELIBLIB) -lm -lga -lGE -ldata_manage -lstat_calc -lannie

# compiler
CXX = g++

# All tests produced by this makefile. Remember to add new tests
# you created to this list
TESTS =  BalAccCalculatorTest MeanSquaredErrCalculatorTest NNSolutionCreatorTest TerminalsTest

# All Google Test headers.  Usually you shouldn't change this definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
	$(GTEST_DIR)/include/gtest/internal/*.h

all : $(TESTS)

clean : 
	rm -rf $(TESTS) gtest.a gtest_main.a *.o

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Build tests.  Test should link with gtest.a if it defines own main()
# function or gtest_main.a if it doesn't.
BalAccCalculator.o: $(USER_DIR)/BalAccCalculator.cpp $(USER_DIR)/BalAccCalculator.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/BalAccCalculator.cpp

BalAccCalculatorTest.o: BalAccCalculatorTest.cpp $(USER_DIR)/BalAccCalculator.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c BalAccCalculatorTest.cpp

BalAccCalculatorTest: BalAccCalculatorTest.o BalAccCalculator.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# MeanSquaredErrCalculatorTest
MeanSquaredErrCalculator.o: $(USER_DIR)/MeanSquaredErrCalculator.cpp $(USER_DIR)/MeanSquaredErrCalculator.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/MeanSquaredErrCalculator.cpp

MeanSquaredErrCalculatorTest.o: MeanSquaredErrCalculatorTest.cpp $(USER_DIR)/MeanSquaredErrCalculator.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c MeanSquaredErrCalculatorTest.cpp

MeanSquaredErrCalculatorTest: MeanSquaredErrCalculatorTest.o MeanSquaredErrCalculator.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# NNSolutionCreatorTest
NNSOLUTIONCREATORTESTOBJS = NNSolutionCreator.o NNLog.o NNSolution.o MeanSquaredErrCalculator.o \
	BalAccCalculator.o Terminals.o BackPropAnnie.o ExpressionTree.o TerminalSymbol.o \
	BackPropAnnieTree.o BackProp.o BackPropTree.o TerminalSymbCreator.o HemannGrammarSI.o \
	Solution.o CalculatorFactory.o RSquaredCalculator.o NNSolutionCreatorIncludeAll.o \
	NNSolutionCreatorIncludeOnce.o

NNSOLUTIONCREATORTESTHEADERS = $(USER_DIR)/NNSolutionCreator.h $(USER_DIR)/SolutionCreator.h \
	$(USER_DIR)/NNSolution.h $(USER_DIR)/Solution.h $(USER_DIR)/MeanSquaredErrCalculator.h \
	$(USER_DIR)/BalAccCalculator.h $(USER_DIR)/Terminals.h $(USER_DIR)/BackPropAnnie.h \
	$(USER_DIR)/ExpressionTree.h $(USER_DIR)/TerminalSymbol.h $(USER_DIR)/BackPropAnnieTree.h \
	$(USER_DIR)/BackProp.h $(USER_DIR)/BackPropTree.h $(USER_DIR)/HemannGrammarSI.h \
	$(USER_DIR)/Solution.h $(USER_DIR)/CalculatorFactory.h $(USER_DIR)/RSquaredCalculator.h \
	$(USER_DIR)/NNSolutionCreatorIncludeAll.h $(USER_DIR)/NNSolutionCreatorIncludeOnce.h

NNSolutionCreator.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/NNSolutionCreator.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/NNSolutionCreator.cpp

NNSolution.o:  $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/NNSolution.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/NNSolution.cpp	

NNLog.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/NNLog.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/NNLog.cpp

Terminals.o: $(NNSOLUTIONCREATORTESTHEADERS)  $(GTEST_HEADERS) $(USER_DIR)/Terminals.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/Terminals.cpp

BackPropAnnie.o: $(NNSOLUTIONCREATORTESTHEADERS)  $(GTEST_HEADERS) $(USER_DIR)/BackPropAnnie.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/BackPropAnnie.cpp

ExpressionTree.o:  $(NNSOLUTIONCREATORTESTHEADERS)  $(GTEST_HEADERS) $(USER_DIR)/ExpressionTree.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/ExpressionTree.cpp

TerminalSymbol.o:  $(NNSOLUTIONCREATORTESTHEADERS)  $(GTEST_HEADERS) $(USER_DIR)/TerminalSymbol.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/TerminalSymbol.cpp

BackPropAnnieTree.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/BackPropAnnieTree.cpp
	 $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/BackPropAnnieTree.cpp

BackProp.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/BackProp.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/BackProp.cpp

BackPropTree.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/BackPropTree.cpp
	 $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/BackPropTree.cpp

TerminalSymbCreator.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/TerminalSymbCreator.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/TerminalSymbCreator.cpp

RSquaredCalculator.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/RSquaredCalculator.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/RSquaredCalculator.cpp

HemannGrammarSI.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/HemannGrammarSI.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/HemannGrammarSI.cpp

CalculatorFactory.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/CalculatorFactory.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/CalculatorFactory.cpp

Solution.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/Solution.cpp
	 $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/Solution.cpp

NNSolutionCreatorIncludeAll.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/NNSolutionCreatorIncludeAll.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/NNSolutionCreatorIncludeAll.cpp

NNSolutionCreatorIncludeOnce.o:  $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/NNSolutionCreatorIncludeOnce.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/NNSolutionCreatorIncludeOnce.cpp

NNSolutionCreatorTest.o: $(NNSOLUTIONCREATORTESTHEADERS) $(GTEST_HEADERS) NNSolutionCreatorTest.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c NNSolutionCreatorTest.cpp

NNSolutionCreatorTest: $(NNSOLUTIONCREATORTESTOBJS) gtest_main.a NNSolutionCreatorTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@ $(LIBS)

TERMINALSTESTOBJS = Terminals.o TerminalSymbol.o 

TERMINALSTESTHEADERS = $(USER_DIR)/Terminals.h $(USER_DIR)/TerminalSymbol.h $(DATA_MANAGE_DIR)/Individual.h

TerminalsTest.o: $(TERMINALSTESTHEADERS) $(GTEST_HEADERS) TerminalsTest.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c TerminalsTest.cpp

TerminalsTest: $(TERMINALSTESTOBJS) gtest_main.a TerminalsTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@ $(LIBS)

