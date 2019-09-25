GL_LIB = libMesaGL.a
GLU_LIB = libMesaGLU.a
GLUT_LIB = libglut.a
CC = gcc
#CFLAGS = -O2 -funroll-loops -ansi -Wall -pedantic -ffast-math -D_SVID_SOURCE -D_BSD_SOURCE -DSHM -DUSE_X86_ASM -I/usr/X11R6/include
#CFLAGS = -g -funroll-loops -ansi -Wall -pedantic -ffast-math -D_SVID_SOURCE -D_BSD_SOURCE -DSHM -DUSE_X86_ASM -I/usr/X11R6/include
CFLAGS = -g -funroll-loops -Wall -pedantic -ffast-math -D_SVID_SOURCE -D_BSD_SOURCE -DSHM -DUSE_X86_ASM -I/usr/X11R6/include
MAKELIB = mklib.ar-ruv
XLIBS = -L/usr/X11/lib -L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE
ASM_SOURCES = asm_386.S
TARGETS = poisadapt

INCDIR = /usr/local/include
LIBDIR = /usr/local/lib

GL_LIBS = -L$(LIBDIR) -lglut -lMesaGLU -lMesaGL -lm $(XLIBS)		
LIB_DEP = $(LIBDIR)/$(GL_LIB) $(LIBDIR)/$(GLU_LIB) $(LIBDIR)/$(GLUT_LIB)

PROGS = poisadapt

viewadaptsuzuki: viewadaptsuzuki.o adapt.h callbacks.o readdata.o cardinal.o $(LIB_DEP)
	gcc -I$(INCDIR) $(CFLAGS) -g -o viewadaptsuzuki viewadaptsuzuki.o callbacks.o readdata.o cardinal.o $(GL_LIBS)

callbacks.o:  callbacks.c adapt.h $(LIB_DEP)
	gcc -I$(INCDIR) $(CFLAGS) -g -c -o callbacks.o callbacks.c

viewadaptsuzuki.o: viewadaptsuzuki.c adapt.h $(LIB_DEP)
	gcc -I$(INCDIR) $(CFLAGS) -g -c -o viewadaptsuzuki.o viewadaptsuzuki.c

readdata.o: readdata.c adapt.h $(LIB_DEP)
	gcc -I$(INCDIR) $(CFLAGS) -g -c -o readdata.o readdata.c

cardinal.o: cardinal.c adapt.h $(LIB_DEP)
	gcc -I$(INCDIR) $(CFLAGS) -g -c -o cardinal.o cardinal.c


