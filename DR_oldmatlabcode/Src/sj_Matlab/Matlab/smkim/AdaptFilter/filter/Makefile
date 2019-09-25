TARGETS=  libAFfilter.a adaptFilter.mexglx 
# SRCS= adaptFilter.cc AscentFilter.cc  AFSmooth.cc EM.cc KalmanSmooth.cc filterFactory.cc 
SRCS= adaptFilter.cc AscentFilter.cc  filterFactory.cc 
OBJS=   $(SRCS:.cc=.o)
AUX= Makefile

ADAPTFILTER_MEX_OBJS=  adaptFilter.o ../aux/mexversion.o ../aux/libAFaux.a ../model/libAFmodel.a libAFfilter.a  

all: $(TARGETS)

adaptFilter.mexglx: $(ADAPTFILTER_MEX_OBJS) $(MEX_AUX) 
	$(CC) -shared $(MEX_LDFLAGS) $(ADAPTFILTER_MEX_OBJS) -o adaptFilter.mexglx
	cp adaptFilter.mexglx ../main  -u

include ../makeconfig
include ../makerules

libAFfilter.a: $(OBJS) $(AUX)
	ar rvs libAFfilter.a $(OBJS)

