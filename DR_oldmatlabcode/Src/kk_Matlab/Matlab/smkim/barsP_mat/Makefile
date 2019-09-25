CC=gcc
CFLAGS=-ansi

com.o: com.c
	$(CC) -c com.c

ranlibsub.o: ranlibsub.c ranlib.h
	$(CC) -c ranlibsub.c

logspline.o: logspline.c logspline.h
	$(CC) -c logspline.c

barsP_utils.o: barsP_utils.c barsP_utils.h
	$(CC) -c barsP_utils.c

barsP_funcs.o: barsP_funcs.c barsP_funcs1.h barsP_funcs.h
	$(CC) -c barsP_funcs.c

barsP.out: barsP.c barsP_funcs.h barsP_funcs.o barsP_utils.o logspline.o ranlibsub.o com.o
	$(CC) -o barsP.out $(CFLAGS) barsP.c barsP_funcs.o barsP_utils.o logspline.o ranlibsub.o com.o -llapack -lblas -lg2c


