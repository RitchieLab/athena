
LIBDIR = ../lib/
LIB= $(LIBDIR)libstat_calc.a
CFLAGS=  -Wall -Wno-deprecated

OBJS = Descriptive.o AUCcalc.o

HEADERS = Descriptive.h AUCcalc.cpp

CC          = g++ -O2
CC          = g++ -g
PRELINK     =
LD          = g++ -w
AR          = ar rv
INSTALL     = install -c
RM          = rm
INC_DIRS    = -I.


%.o: %.cpp $(INTERFACES)
	$(CC) $(CFLAGS) $(INC_DIRS) -c $< -o $@
$(LIB): $(OBJS)
	$(PRELINK)
	$(AR) $(LIB) $?
	@echo $(LIB) is now up-to-date
clean:
	rm -f *.o *~ core 
