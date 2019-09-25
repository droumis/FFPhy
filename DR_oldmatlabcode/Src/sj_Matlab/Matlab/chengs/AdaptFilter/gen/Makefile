TARGETS=  libAFgen.a genSpikes.mexglx 

SRCS= genSpikes.cc genFactory.cc PoissonGen.cc RescalingGen.cc VGenerator.cc 
OBJS=   $(SRCS:.cc=.o)
AUX= Makefile

GENSPIKES_MEX_OBJS=  genSpikes.o libAFgen.a  ../model/libAFmodel.a ../aux/libAFaux.a 

### TARGETS ###
all: $(TARGETS) $(AUX)

genSpikes.mexglx: $(GENSPIKES_MEX_OBJS) $(MEX_AUX) 
	$(CC) -shared $(MEX_LDFLAGS) $(GENSPIKES_MEX_OBJS) -o genSpikes.mexglx
	cp genSpikes.mexglx ../main  -u

include ../makeconfig
include ../makerules

libAFgen.a: $(OBJS)
	ar rvs libAFgen.a $(OBJS)


