# -----------------------------------------------------------------------------
# Make ATHENA
# -----------------------------------------------------------------------------

# compiler flags


STATDIR = ../stat_calc
DATADIR = ../data_manage
GADIR = ../galib247
GEDIR = ../libGE-0.26
ANNIEDIR = ../annie

INCLUDEDIR = -I. -I$(STATDIR) -I$(DATADIR) -I$(GADIR) -I$(GEDIR) -I$(ANNIEDIR)
DFLAGS = -DPARALLEL
#DFLAGS = 
CCFLAGS = -Wall $(DFLAGS) -O2 $(INCLUDEDIR)
#CCFLAGS = -Wall $(DFLAGS) -g $(INCLUDEDIR)

GALIBLIB = ../galib247/ga
GELIBLIB = ../libGE-0.26

LIBDIR = -L$(GALIBLIB) -L$(GELIBLIB) -L../lib 
LIBS = -lm -lga -lGE -ldata_manage -lstat_calc -lannie
PROG_NAME = athena
EXECUTABLE = athena
TEST_NAME = athena_test

SYSOBJS = AlgorithmFactory.o BalAccCalculator.o CalculatorFactory.o Config.o \
	ConfigFileReader.o GE1DArrayGenome.o GENNGrammarAdjuster.o GEObjective.o \
	athena.o InitGEgenome.o MeanSquaredErrCalculator.o NNSolution.o \
	NNSolutionCreator.o OutputManager.o Population.o Solution.o \
	SolutionFactory.o Terminals.o TerminalSymbCreator.o TerminalSymbol.o \
	GENNAlg.o NNSolutionCreatorIncludeAll.o NNSolutionCreatorIncludeOnce.o \
	ExpressionTree.o NNLog.o Algorithm.o RSquaredCalculator.o BioFilterReader.o \
	HemannGrammarSI.o BioFilterModelCollection.o BackProp.o BackPropTree.o \
	BioFilterBinReader.o genegenemodel.o BackPropSMD.o BackPropAnnie.o \
	BackPropAnnieTree.o SymRegressSolutionCreator.o GESymbReg.o

SYSHEADERS = AlgorithmFactory.h Algorithm.h AlgorithmList.h BalAccCalculator.h \
	CalculatorFactory.h CalculatorList.h ConfigFileReader.h Config.h \
	GE1DArrayGenome.h GENNAlg.h GENNGrammarAdjuster.h GEObjective.h \
	HemannExcept.h InitGEgenome.h MeanSquaredErrCalculator.h NNSolution.h \
	NNSolutionCreator.h OutputManager.h Population.h rbtree.h \
	SolutionCalculator.h SolutionCreator.h SolutionFactory.h Solution.h \
	Terminals.h TerminalSymbCreator.h TerminalSymbol.h RSquaredCalculator.h \
	NNSolutionCreatorIncludeAll.h NNSolutionCreatorIncludeOnce.h \
	ExpressionTree.h Tree.hh AlgorithmLog.h NNLog.h BioFilterReader.h \
	HemannGrammarSI.h BioFilterModelCollection.h BackProp.h \
	BackPropTree.h Structs.h BioFilterBinReader.h BioReader.h \
	genegenemodel.h genegenemodelreader.h region.h snpholder.h \
	snpsnpmodel.h biomodel.h BackPropSMD.h BackPropAnnie.h BackPropAnnieTree.h \
	SymRegressSolutionCreator.h GESymbReg.h

ifeq ($(DFLAGS),-DPARALLEL)
CC = mpiCC
SYSOBJS += TransferData.o
SYSHEADERS += TransferData.h
EXECUTABLE = athenap
TEST_NAME = athenap_test
else
CC = g++
endif

all: $(TEST_NAME)

$(PROG_NAME): $(SYSOBJS) $(SYSHEADERS)
	$(CC) -o ../../bin/$(EXECUTABLE) $(SYSOBJS) $(LIBDIR) $(LIBS)

%.o: %.cpp $(SYSHEADERS) 
	$(CC) -c $(CCFLAGS) $< -o $@

final: $(PROG_NAME)

$(TEST_NAME): $(SYSOBJS) $(SYSHEADERS)
	$(CC) -o ../../bin/$(TEST_NAME) $(SYSOBJS) $(LIBDIR) $(LIBS)

tarball:
	tar cvfz hemann_source_`date '+%Y%m%d'`.tar.gz \
	*.cpp *.h* Makefile* 

clean:
	rm -f *.o *~ core

