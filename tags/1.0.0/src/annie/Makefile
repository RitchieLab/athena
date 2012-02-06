LIBDIR = ../lib/
LIB= $(LIBDIR)libannie.a
CFLAGS=  -Wall -Wno-deprecated

OBJS = AbstractNeuron.o File.o Layer.o Network.o SimpleNeuron.o \
	CenterNeuron.o InputNeuron.o MultiLayerNetwork.o Neuron.o \
	Exception.o Link.o Control.o Vector.o random.o auxf.o \
	InputLayer.o TrainingSet.o examples.o TwoLayerNetwork.o \
	args.o listeners.o

HEADERS = AbstractNeuron.h File.h Layer.h Network.h SimpleNeuron.h \
	CenterNeuron.h InputNeuron.h MultiLayerNetwork.h Neuron.h \
	Exception.h Link.h Control.h Vector.h random.h auxf.h \
	InputLayer.h TrainingSet.h examples.h TwoLayerNetwork.h \
	args.h listeners.h defines.h config.h

CC          = g++ -O2
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
