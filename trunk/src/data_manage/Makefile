
LIBDIR = ../lib/
LIB= $(LIBDIR)libdata_manage.a
CFLAGS=  -Wall -Wno-deprecated

OBJS = ContinFileReader.o CrossValidator.o CVInterval.o CVSet.o DataExcept.o \
	Dataholder.o Dataset.o FileHandler.o Individual.o MapFileReader.o \
	MDRFileHandler.o NormalizeContinuous.o OttDummyConvert.o \
	ScaleContinuous.o ScaleData.o Stringmanip.o OutputSet.o \
	ScaledDataFactory.o RandomNoReplace.o StephenDummyConvert.o \
	ContinMapFileReader.o ScaleGroupContinuous.o CustomDummyConvert.o \
	EncodingFactory.o
	
HEADERS = ContinFileReader.h CrossValidator.h CVInterval.h CVSet.h DataExcept.h \
	Dataholder.h Dataset.h FileHandler.h Individual.h MapFileReader.h \
	MDRFileHandler.h NormalizeContinuous.h OttDummyConvert.h \
	ScaleContinuous.h ScaleData.h Stringmanip.h OutputSet.h \
	ScaledDataFactory.h RandomNoReplace.h StephenDummyConvert.h \
	ContinMapFileReader.h ScaleGroupContinuous.h DummyConvert.h \
	CustomDummyConvert.h EncodingFactory.h Utility.h

CC          = g++ -O2
#CC          = g++ -g
PRELINK     =
LD          = g++ -w
AR          = ar rv
INSTALL     = install -c
RM          = rm
INC_DIRS    = -I. -I../stat_calc


%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INC_DIRS) -c $< -o $@
$(LIB): $(OBJS)
	$(PRELINK)
	$(AR) $(LIB) $?
	@echo $(LIB) is now up-to-date
clean:
	rm -f *.o *~ core 
