

############ Matlab ############

############ Compiler and linker flags ############


############ Sources and Objects ############
FLSRCS		= GUIMain.fl 

CCSRCS_VO 	= Geometry.cc Plot1d.cc Plot2d.cc
CCSRCS_APP	= camera.cc ModelVis.cc EventHandler.cc Display.cc MovieControl.cc DataCache.cc Activity.cc
CCSRCS_UT   = glUtils.cc axes.cc colormapTexture.cc visModel.cc
CCSRCS 		= $(CCSRCS_APP) $(CCSRCS_VO) $(CCSRCS_UT)
SRCS		= $(FLSRCS:.fl=.cxx) $(CCSRCS) 
OBJS		= $(FLSRCS:.fl=.o) $(CCSRCS:.cc=.o)
AUX		=Makefile
MEX_OBJ=  ../aux/mexversion.o ../aux/libAFaux.a ../model/libAFmodel.a

############ Targets ############

TARGETS			= $(FLSRCS:.fl=.cxx) visModel.mexglx

all: $(TARGETS)

include ../makeconfig
include ../makerules


visModel.mexglx: $(AUX) $(MEX_OBJ) $(OBJS) 
	$(CC) -shared $(MEX_LDFLAGS) $(OBJS) $(MEX_OBJ) $(UI_LIBS) -o visModel.mexglx
	cp -f visModel.mexglx ../main -u
# visModel.mexglx: $(AUX) $(MEX_OBJ) $(OBJS) 
#     $(CC) -shared $(MEX_LDFLAGS) -L/usr/local/lib -lfltk -ldl -lm -lXext -lX11 -lGL $(OBJS) $(MEX_OBJ) -o visModel.mexglx
#     cp -f visModel.mexglx ../main -u

tryv: $(OBJS) $(AUX) $(MEX_OBJ)
	$(CC) $(OTIM) $(EXEBIT) -Wl,--rpath,$(MDIR)/sys/os/glnx86:$(MDIR)/bin/glnx86 -L$(MDIR)/sys/os/glnx86 -L$(MDIR)/bin/glnx86 -lmx -lmex -lm  $(UI_LIBS) $(MEX_OBJ) $(OBJS)  ../main/try.cc -o tryv



# mexversion.o:
#     g++ $(EXEBIT) -c -g -DMATLAB_MEX_FILE -fPIC -ansi -D_GNU_SOURCE -pthread -O -DNDEBUG \
#     -I$(MDIR)/extern/include -I$(MDIR)/simulink/include -Wall \
#     $(MDIR)/extern/src/mexversion.c

############ Suffix rules ############

.SUFFIXES: .fl .cxx

.fl.cxx .fl.h:
	$(FLC) $<

.cxx.o:
	$(CC) -c $(CFLAGS) $(INCL) -o $@ $<


############ include makefiles ############
