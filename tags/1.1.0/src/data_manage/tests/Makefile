# Location of original source code
USER_DIR = ..

# Location of stat calculator header
CALC_DIR = ../../stat_calc

# Location/name of calc library
CALC_LIB = ../../lib

# Point to googletest root directory
GTEST_DIR = ../../gtest-1.5.0

# Flags for preprocessor
CPPFLAGS += -I$(GTEST_DIR) -I$(GTEST_DIR)/include

# Flags for the compiler
CXXFLAGS += -g -Wall -Wextra

# compiler
CXX = g++

# All tests produced by this makefile. Remember to add new tests
# you created to this list
TESTS = IndividualTest FileHandlerTest ContinFileReaderTest DatasetTest \
	CVIntervalTest CVSetTest CrossValidatorTest DataholderTest \
	MapFileReaderTest NormalizeContinuousTest ScaleContinuousTest \
	OttDummyConvertTest StephenDummyConvertTest StringmanipTest

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

Individual.o: $(USER_DIR)/Individual.cpp $(USER_DIR)/Individual.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/Individual.cpp

IndividualTest.o: IndividualTest.cpp $(USER_DIR)/Individual.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c IndividualTest.cpp

IndividualTest: IndividualTest.o Individual.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

FILEHANDLERHEADERS = $(USER_DIR)/Dataholder.h $(USER_DIR)/FileHandler.h $(USER_DIR)/MDRFileHandler.h \
	$(USER_DIR)/Stringmanip.h $(USER_DIR)/DataExcept.h

FILEHANDLERTESTOBJS = Dataholder.o FileHandler.o MDRFileHandler.o Individual.o Stringmanip.o \
	DataExcept.o

Dataholder.o:  $(USER_DIR)/Dataholder.cpp  $(USER_DIR)/Dataholder.h $(GTEST_HEADERS) \
	$(USER_DIR)/Individual.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/Dataholder.cpp

FileHandler.o: $(USER_DIR)/FileHandler.cpp $(USER_DIR)/FileHandler.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/FileHandler.cpp

MDRFileHandler.o: $(USER_DIR)/MDRFileHandler.cpp $(USER_DIR)/MDRFileHandler.h \
	$(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/MDRFileHandler.cpp

Stringmanip.o: $(USER_DIR)/Stringmanip.cpp $(USER_DIR)/Stringmanip.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/Stringmanip.cpp

DataExcept.o: $(USER_DIR)/DataExcept.cpp $(USER_DIR)/DataExcept.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/DataExcept.cpp

FileHandlerTest.o: FileHandlerTest.cpp $(FILEHANDLERHEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c FileHandlerTest.cpp

FileHandlerTest: $(FILEHANDLERTESTOBJS) gtest_main.a FileHandlerTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@
	
ContinFileReader.o: $(USER_DIR)/ContinFileReader.cpp $(USER_DIR)/ContinFileReader.h \
	$(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/ContinFileReader.cpp

ContinFileReaderTest.o: ContinFileReaderTest.cpp $(USER_DIR)/ContinFileReader.h 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c ContinFileReaderTest.cpp

ContinFileReaderTest: ContinFileReader.o Individual.o Dataholder.o ContinFileReaderTest.o \
	gtest_main.a DataExcept.o Stringmanip.o
	 $(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# Dataset Test
DATASETTESTOBJS = Dataset.o Individual.o Stringmanip.o

DATASETTESTSOURCE = $(USER_DIR)/Dataset.cpp $(USER_DIR)/Individual.cpp \
	$(USER_DIR)/Stringmanip.cpp	

DATASETTESTHEADERS = $(USER_DIR)/Dataset.h $(USER_DIR)/Individual.h \
	$(USER_DIR)/Stringmanip.h

#$(DATASETTESTOBJS): $(DATASETHEADERS) $(DATASETTESTSOURCE) $(GTEST_HEADERS)
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

#Dataset.o: $(USER_DIR)/Dataset.cpp $(DATASETHEADERS) $(GTEST_HEADERS)
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/Dataset.cpp

Dataset.o: $(USER_DIR)/Dataset.cpp $(USER_DIR)/Dataset.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/Dataset.cpp

DatasetTest.o: DatasetTest.cpp $(DATASETTESTHEADERS) $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c DatasetTest.cpp

DatasetTest: $(DATASETTESTOBJS) gtest_main.a DatasetTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# CVInterval Test

CVINTERVALTESTOBJS = Dataset.o Individual.o Stringmanip.o CVInterval.o

CVINTERVALTESTHEADERS = $(USER_DIR)/Dataset.h $(USER_DIR)/Individual.h \
	$(USER_DIR)/Stringmanip.h $(USER_DIR)/CVInterval.h

CVInterval.o: $(CVINTERVALTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/CVInterval.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/CVInterval.cpp

CVIntevalTest.o: $(CVINTERVALTESTHEADERS) $(GTEST_HEADERS) CVIntervalTest.cpp
	 $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c CVIntervalTest.cpp

CVIntervalTest: $(CVINTERVALTESTOBJS) gtest_main.a CVIntervalTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@


# CVSet Test
CVSETTESTOBJS = CVSet.o CVInterval.o Individual.o Dataset.o

CVSETTESTHEADERS = $(USER_DIR)/CVSet.h $(USER_DIR)/CVInterval.h $(USER_DIR)/Individual.h \
	$(USER_DIR)/Dataset.h

CVSet.o: $(CVSETTESTHEADERS) $(USER_DIR)/CVSet.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/CVSet.cpp

CVSetTest.o: $(CVSETTESTHEADERS) $(GTEST_HEADERS) CVSetTest.cpp
	$(CXX)  $(CPPFLAGS) $(CXXFLAGS) -c CVSetTest.cpp

CVSetTest: $(CVSETTESTOBJS) gtest_main.a CVSetTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# CrossValidatorTest
CROSSVALIDATORTESTOBJS = $(CVSETTESTOBJS) MDRFileHandler.o Dataholder.o CrossValidator.o \
	DataExcept.o FileHandler.o Stringmanip.o

CROSSVALIDATORTESTHEADERS = $(CVSETTESTHEADERS) $(USER_DIR)/MDRFileHandler.h \
	$(USER_DIR)/Dataholder.h $(USER_DIR)/CrossValidator.h $(USER_DIR)/DataExcept.h \
	$(USER_DIR)/FileHandler.h $(USER_DIR)/Stringmanip.h

CrossValidator.o: $(CROSSVALIDATORHEADERS) $(GTEST_HEADERS) $(USER_DIR)/CrossValidator.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/CrossValidator.cpp

CrossValidatorTest: $(CROSSVALIDATORTESTOBJS) gtest_main.a CrossValidatorTest.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# DataholderTest
DATAHOLDERTESTOBJS = Dataholder.o Individual.o DataExcept.o Stringmanip.o

DATAHOLDERTESTHEADERS = $(USER_DIR)/Dataholder.h $(USER_DIR)/Individual.h \
	$(USER_DIR)/DataExcept.h $(USER_DIR)/Stringmanip.h

DataholderTest.o: $(DATAHOLDERTESTHEADERS) $(GTEST_HEADERS) DataholderTest.cpp
	$(CXX)  $(CPPFLAGS) $(CXXFLAGS) -c DataholderTest.cpp

DataholderTest: $(DATAHOLDERTESTOBJS) gtest_main.a DataholderTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# MapFileReaderTest
MAPFILEREADERTESTOBJS = MapFileReader.o DataExcept.o Dataholder.o Individual.o \
	Stringmanip.o

MAPFILEREADERTESTHEADERS = $(USER_DIR)/Dataholder.h $(USER_DIR)/Individual.h \
	$(USER_DIR)/DataExcept.h $(USER_DIR)/MapFileReader.h $(USER_DIR)/Stringmanip.h

MapFileReader.o: $(MAPFILEREADERTESTHEADERS) $(GTEST_HEADERS) $(USER_DIR)/MapFileReader.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/MapFileReader.cpp

MapFileReaderTest.o: $(MAPFILEREADERTESTHEADERS) $(GTEST_HEADERS) MapFileReaderTest.cpp
	 $(CXX)  $(CPPFLAGS) $(CXXFLAGS) -c MapFileReaderTest.cpp

MapFileReaderTest: $(MAPFILEREADERTESTOBJS) gtest_main.a MapFileReaderTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# NormalizeContinuousTest
NORMALIZECONTINUOUSTESTOBJS = NormalizeContinuous.o ScaleData.o Dataholder.o Individual.o \
	Stringmanip.o

NORMALIZECONTINUOUSTESTHEADERS = $(USER_DIR)/NormalizeContinuous.h $(USER_DIR)/ScaleData.h \
	$(USER_DIR)/Dataholder.h $(USER_DIR)/Individual.h $(CALC_DIR)/Descriptive.h \
	$(USER_DIR)/Stringmanip.h

NormalizeContinuous.o: $(NORMALIZECONTINUOUSTESTHEADERS) $(GTEST_HEADERS) \
	$(USER_DIR)/NormalizeContinuous.cpp $(CALC_DIR)/Descriptive.h
	 $(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(CALC_DIR) -c $(USER_DIR)/NormalizeContinuous.cpp

ScaleData.o: $(NORMALIZECONTINUOUSTESTHEADERS) $(GTEST_HEADERS) \
	$(USER_DIR)/ScaleData.cpp 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/ScaleData.cpp

NormalizeContinuousTest.o: $(NORMALIZECONTINUOUSTESTHEADERS) $(GTEST_HEADERS) \
	NormalizeContinuousTest.cpp
	$(CXX)  $(CPPFLAGS) $(CXXFLAGS) -c NormalizeContinuousTest.cpp

NormalizeContinuousTest: $(NORMALIZECONTINUOUSTESTOBJS) gtest_main.a $(CALC_LIB)/libstat_calc.a \
	NormalizeContinuousTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@
	
	
# ScaleContinuousTest
SCALECONTINUOUSTESTOBJS = ScaleContinuous.o ScaleData.o Dataholder.o Individual.o \
        Stringmanip.o

SCALECONTINUOUSTESTHEADERS = $(USER_DIR)/ScaleContinuous.h $(USER_DIR)/ScaleData.h \
        $(USER_DIR)/Dataholder.h $(USER_DIR)/Individual.h $(CALC_DIR)/Descriptive.h \
	$(USER_DIR)/Stringmanip.h

ScaleContinuous.o: $(SCALECONTINUOUSTESTHEADERS) $(GTEST_HEADERS) \
        $(USER_DIR)/ScaleContinuous.cpp $(CALC_DIR)/Descriptive.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(CALC_DIR) -c $(USER_DIR)/ScaleContinuous.cpp

ScaleContinuousTest.o: $(SCALECONTINUOUSTESTHEADERS) $(GTEST_HEADERS) \
        ScaleContinuousTest.cpp
	$(CXX)  $(CPPFLAGS) $(CXXFLAGS) -c ScaleContinuousTest.cpp

ScaleContinuousTest: $(SCALECONTINUOUSTESTOBJS) gtest_main.a $(CALC_LIB)/libstat_calc.a \
        ScaleContinuousTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# OttDummyConvertTest
OTTDUMMYCONVERTTESTOBJS = Dataholder.o Individual.o Stringmanip.o OttDummyConvert.o

OTTDUMMYCONVERTTESTHEADERS = $(USER_DIR)/Dataholder.h $(USER_DIR)/Individual.h \
	$(USER_DIR)/Stringmanip.h $(USER_DIR)/OttDummyConvert.h

OttDummyConvert.o: $(OTTDUMMYCONVERTTESTHEADERS) $(GTEST_HEADERS) \
	$(USER_DIR)/OttDummyConvert.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(CALC_DIR) -c $(USER_DIR)/OttDummyConvert.cpp

OttDummyConvertTest.o: $(OTTDUMMYCONVERTTESTHEADERS) $(GTEST_HEADERS) \
	OttDummyConvertTest.cpp
	$(CXX)  $(CPPFLAGS) $(CXXFLAGS) -c OttDummyConvertTest.cpp

OttDummyConvertTest: $(OTTDUMMYCONVERTTESTOBJS)  gtest_main.a OttDummyConvertTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# StephenDummyConvertTest
STEPHENDUMMYCONVERTTESTOBJS = Dataholder.o Individual.o Stringmanip.o StephenDummyConvert.o

STEPHENDUMMYCONVERTTESTHEADERS = $(USER_DIR)/Dataholder.h $(USER_DIR)/Individual.h \
	$(USER_DIR)/Stringmanip.h $(USER_DIR)/StephenDummyConvert.h

StephenDummyConvert.o: $(STEPHENDUMMYCONVERTTESTHEADERS) $(GTEST_HEADERS) \
	$(USER_DIR)/StephenDummyConvert.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(CALC_DIR) -c $(USER_DIR)/StephenDummyConvert.cpp

StephenDummyConvertTest.o: $(STEPHENDUMMYCONVERTTESTHEADERS) $(GTEST_HEADERS) \
	StephenDummyConvertTest.cpp
	$(CXX)  $(CPPFLAGS) $(CXXFLAGS) -c StephenDummyConvertTest.cpp

StephenDummyConvertTest: $(STEPHENDUMMYCONVERTTESTOBJS)  gtest_main.a StephenDummyConvertTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

#StringmanipTest
STRINGMANIPTESTOBJS = Stringmanip.o

STRINGMANIPTESTHEADERS = $(USER_DIR)/Stringmanip.h

StringmanipTest.o: $(STRINGMANIPTESTHEADERS) $(GTEST_HEADERS) \
        StringmanipTest.cpp
	$(CXX)  $(CPPFLAGS) $(CXXFLAGS) -c StringmanipTest.cpp

StringmanipTest: $(STRINGMANIPTESTOBJS)  gtest_main.a StringmanipTest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

