TARGETS=  libAFstats.a calcStatistics.mexglx 
SRCS= calcStatistics.cc statsFactory.cc VStat.cc VStat1d.cc VStat2d.cc VStatAM.cc StatLimits.cc MutualInfo.cc RescaledSpikes.cc LinCorr.cc LinCorrOffset.cc LinCorrShift.cc Integral2d.cc Moments2d.cc CircMeasures.cc Moments.cc 
OBJS=   $(SRCS:.cc=.o)
AUX= Makefile

CALCSTATISTICS_MEX_OBJ=  calcStatistics.o ../aux/mexversion.o ../aux/libAFaux.a ../model/libAFmodel.a libAFstats.a


### TARGETS ###
all: $(TARGETS)

calcStatistics.mexglx: $(CALCSTATISTICS_MEX_OBJ) $(MEX_AUX) 
	$(CC) -shared $(MEX_LDFLAGS) $(CALCSTATISTICS_MEX_OBJ) -o calcStatistics.mexglx
	cp calcStatistics.mexglx ../main  -u

include ../makeconfig
include ../makerules

libAFstats.a: $(OBJS) $(AUX)
	ar rvs libAFstats.a $(OBJS)


