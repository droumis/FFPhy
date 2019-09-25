############ include makefiles ############
include ../makeinclude

############ FLTK and other graphics ############
FLC		=fluid -c

X_LIBS		=-L/usr/lib -L/usr/X11R6/lib -Wl,--rpath,/usr/lib:/usr/X11R6/lib  -lXext -lX11 -lm #-lsupc++
FLTK_LIBS	=-L/usr/lib -Wl,--rpath,/usr/lib -lfltk -lfltk_gl -lfltk_forms
GL_LIBS		=-L/usr/lib -L/usr/X11R6/lib -Wl,--rpath,/usr/lib:/usr/X11R6/lib -lGL -lGLU #-lsupc++
UI_LIBS		= $(X_LIBS) $(FLTK_LIBS) $(GL_LIBS) 

############ Matlab ############
MEX_AUX		=Makefile

############ Compiler and linker flags ############


############ Sources and Objects ############
FLSRCS		= afViewUI.fl afMeshViewUI.fl 

CCSRCS_VO 	= afTrajMesh.cc afTrajMovie.cc afIsiPlot.cc afGeomPlot.cc 
CCSRCS_APP	= afDataModel.cc trajInfo.cc viewWindow.cc meshWindow.cc movieWindow.cc camera.cc movieContext.cc
CCSRCS_UT   = glUtils.cc axes.cc colormapTexture.cc mexMain.cc
CCSRCS 		= $(CCSRCS_VO) $(CCSRCS_APP) $(CCSRCS_UT)
SRCS		= $(FLSRCS:.fl=.cxx) $(CCSRCS)
OBJS		= $(FLSRCS:.fl=.o) $(CCSRCS:.cc=.o)

############ Targets ############

#TARGETS		= test vis.mexglx 

TARGETS			= vis.mexglx

all: $(TARGETS)

#test: main.o $(OBJS) 
#	$(CC) $(EXEBIT) -g -shared $(OBJS) $(UI_LIBS) main.o -o test	

vis.mexglx: $(OBJS) mexversion.o $(MEX_AUX) \
		../main/libAFaux.a ../main/libAFgen.a ../main/libAFmodel.a ../main/libAFfilter.a
	$(CC) -shared $(MEX_LDFLAGS) $(OBJS) mexversion.o $(UI_LIBS) \
		../main/libAFaux.a ../main/libAFgen.a ../main/libAFmodel.a ../main/libAFfilter.a  -o vis.mexglx
	cp -f vis.mexglx ../main

mexversion.o:
	g++ $(EXEBIT) -c -g -DMATLAB_MEX_FILE -fPIC -ansi -D_GNU_SOURCE -pthread -O -DNDEBUG \
	-I$(MDIR)/extern/include -I$(MDIR)/simulink/include -Wall \
	$(MDIR)/extern/src/mexversion.c

############ Suffix rules ############

.SUFFIXES: .fl .cxx

.fl.cxx .fl.h:
	$(FLC) $<

.cxx.o:
	$(CC) -c $(CFLAGS) -o $@ $<

.cc.o:
	$(CC) -c $(CFLAGS) -o $@ $<

############ Utilities ############

clean:
	rm -f $(TARGETS) *.o *~ $(FLSRCS:.fl=.cxx) $(FLSRCS:.fl=.h) makedepend*

############ Dependencies ############

depend: $(FLSRCS:.fl=.cxx)
	rm -f makedepend
	echo "" > makedepend
	makedepend -- $(CFLAGS) -- $(SRCS) -f makedepend

makedepend:
	echo "" > makedepend
	makedepend -- $(CFLAGS) -- $(SRCS) -f makedepend

include makedepend

