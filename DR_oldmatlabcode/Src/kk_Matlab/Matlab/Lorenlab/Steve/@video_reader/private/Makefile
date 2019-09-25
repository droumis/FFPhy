MEXSUFFIX = mexglx

CXX       = g++
CFLAGS    = -ggdb -c -fpic

MEX       = mex
MEXFLAGS  = -cxx -g

LIBS      = -lmpeg2 -lmpeg2convert
INCLUDE   = -I/usr/include/mpeg2dec

all: video_reader_mex.$(MEXSUFFIX)

video_reader_mex.$(MEXSUFFIX): VideoReader.o
	$(MEX) $(MEXFLAGS) video_reader_mex.cpp VideoReader.o $(LIBS) $(INCLUDE)

VideoReader.o:
	$(CXX) $(CFLAGS) VideoReader.cpp $(LIBS) $(INCLUDE)

clean:
	rm video_reader_mex.$(MEXSUFFIX) VideoReader.o

