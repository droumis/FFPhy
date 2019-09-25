# Makefile for extract_frames
# last edited by SMK, 2 April 2008

CC = g++
CFLAGS = -ggdb -O2

LIBS = -lpng -lmpeg2 -lmpeg2convert
INCDIR = /usr/include/mpeg2dec/

extract_frames: extract_frames.cpp
	$(CC) $(CFLAGS) $(LIBS) -I$(INCDIR) -o extract_frames extract_frames.cpp



