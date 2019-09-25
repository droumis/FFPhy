TARGETS=  # adaptFilterIni.mexglx 

# COMMON_OBJ= $(AFLIB) mexversion.o 
# COMMON_OBJ= mexversion.o ../gen/libAFgen.a ../filter/libAFfilter.a  ../stats/libAFstats.a ../model/libAFmodel.a ../aux/libAFaux.a
COMMON_OBJ= ../gen/libAFgen.a ../filter/libAFfilter.a  ../stats/libAFstats.a ../model/libAFmodel.a ../aux/libAFaux.a
MEX_AUX= Makefile

# SRCS= calcStatistics.cc genSpikes.cc eval2DModel.cc adaptFilterIni.cc  

### TARGETS ###
all: $(TARGETS)

include ../makeconfig
include ../makerules

try:  tryg tryf trys

tryf: $(ADAPTFILTER_MEX_OBJS) $(MEX_AUX) $(COMMON_OBJ)
	$(CC) $(EXEBIT) -g -Wl,--rpath,$(MDIR)/sys/os/glnx86:$(MDIR)/bin/glnx86 -L$(MDIR)/sys/os/glnx86 -L$(MDIR)/bin/glnx86 -lmx -lmex -lm $(ADAPTFILTER_MEX_OBJS)  $(COMMON_OBJ) try.cc -o tryf

tryg: $(GENSPIKES_MEX_OBJS) $(MEX_AUX) $(COMMON_OBJ)
	$(CC) $(EXEBIT) -g  -L$(MDIR)/sys/os/glnx86 -L$(MDIR)/bin/glnx86 -lmx -lmex -lm $(GENSPIKES_MEX_OBJS)  $(COMMON_OBJ) try.cc -o tryg

trys: $(MEX_AUX) $(COMMON_OBJ) calcStatistics.o
	$(CC) $(EXEBIT) -g  -Wl,--rpath,$(MDIR)/sys/os/glnx86:$(MDIR)/bin/glnx86 -L$(MDIR)/sys/os/glnx86 -L$(MDIR)/bin/glnx86 -lmx -lmex -lm calcStatistics.o $(COMMON_OBJ) try.cc -o trys

trye: eval2DModel.o $(MEX_AUX) $(COMMON_OBJ)
	$(CC) $(EXEBIT) -g -Wl,--rpath,$(MDIR)/sys/os/glnx86:$(MDIR)/bin/glnx86 -L$(MDIR)/sys/os/glnx86 -L$(MDIR)/bin/glnx86 -lmx -lmex -lm eval2DModel.o  $(COMMON_OBJ) try.cc -o trye

tryv: $(MEX_AUX) $(COMMON_OBJ) ../newvis/visModel.o ../newvis/visModel.o
	$(CC) $(EXEBIT) -g  -Wl,--rpath,$(MDIR)/sys/os/glnx86:$(MDIR)/bin/glnx86 -L$(MDIR)/sys/os/glnx86 -L$(MDIR)/bin/glnx86 -lmx -lmex -lm ../newvis/visModel.o $(COMMON_OBJ) try.cc -o tryv

trye: eval2DModel.o $(MEX_AUX) $(COMMON_OBJ)

# mexversion.o:
#     g++ $(EXEBIT) -c -g -DMATLAB_MEX_FILE -fPIC -ansi -D_GNU_SOURCE -pthread -O -DNDEBUG \
#     -I$(MDIR)/extern/include -I$(MDIR)/simulink/include -Wall \
#     $(MDIR)/extern/src/mexversion.c

