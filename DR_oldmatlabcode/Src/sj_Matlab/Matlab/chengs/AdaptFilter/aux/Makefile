
TARGETS=  mexversion.o libAFaux.a 
SRCS=  TData.cc hippoIO.cc numerics.cc Grad.cc 
OBJS=   $(SRCS:.cc=.o)
AUX= Makefile

### TARGETS ###
all: $(TARGETS)

include ../makeconfig
include ../makerules

libAFaux.a: $(OBJS) $(AUX)
	ar rvs libAFaux.a $(OBJS)

mexversion.o:
	g++ $(CFLAGS) -c -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE -pthread -DNDEBUG \
	-I$(MDIR)/extern/include -I$(MDIR)/simulink/include \
	$(MDIR)/extern/src/mexversion.c
