#
# Makefile for Simple Genetic Algorithm C code
#
#####################################
# Define implicit compilation rules #
#####################################

# uncomment following lines for UNIX computers
# nominally tested on Sun, Cray UNICOS, and Vax Ultrix systems
CC=ccc
CFLAGS  = -I/opt/builds/libGE-cvs-0.26alpha
LDFLAGS = -L/opt/builds/libGE-cvs-0.26alpha
LDLIBS  = -lm -lGE

# uncomment two of the following lines to determine selection routine
# roulette wheel, stochastic remainder, and tournament selection
# are distributed with SGA-C
SELECTCODE = rselect.c
SELECTOBJ = rselect.o
#SELECTCODE = srselect.c
#SELECTOBJ = srselect.o
#SELECTCODE = tselect.c
#SELECTOBJ = tselect.o


.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<


SRCCODE=main.c app.c generate.c initial.c memory.c\
	operators.c random.c report.c $(SELECTCODE)\
	statistic.c utility.c

OBJECTS=main.o app.o generate.o initial.o\
	memory.o operators.o random.o report.o\
	$(SELECTOBJ) statistic.o utility.o

sga:    $(OBJECTS)
	$(CC) -c GEant.c GEtrail.c
	g++ $(CFLAGS) -o $@ $(OBJECTS) $(LDLIBS)
	chmod 755 $@

$(OBJECTS): sga.h external.h

listing:
	cat *.h > sga.list
	cat $(SRCCODE) >> sga.list

clean:
	rm -f sga 
	rm -f *.on sga.sym
	rm -f *.o
