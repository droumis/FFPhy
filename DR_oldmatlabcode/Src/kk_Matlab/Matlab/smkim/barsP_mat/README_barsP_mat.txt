BarsP Matlab wrappers 1.0
Ryan Kelly 4/11/08


This package contains a matlab wrapper for BARS (Bayesian Adaptive Regression Splines).  More details on the software are in the subdirectories here, or on the BARS web page: 
http://lib.stat.cmu.edu/~kass/bars/

I've compiled the BARS source code and added a wrapper function, barsP_mat.c, which may be compiled with the Matlab compiler, mex.  It has compiled successfully on at least one OS X, other unix, and windows machine.


The changes made to the original source code are:
- in barsP_funcs1.h, added a preprocessor directive to add trailing underscores for calling external Fortran programs
- in logspline.c, added a preprocessor directive to include malloc.h only for systems which have it, otherwise stdlib.h is included.
- in logspline.c, in function rpqlsd(), defined a float "zero" to get around a windows compiler issue.
- added testBars.m, barsP_mat.c, and barsP_mat.m to the directory
- All memory allocations and frees are replaced with their mx* versions.

Matlab Compilation
In Matlab, run the command
>> mex barsP_mat.c barsP_utils.c barsP_funcs.c com.c logspline.c ranlibsub.c
This will produce a binary barsP_mat which may be called directly from Matlab.  
IMPORTANT: In Windows, and some unix setups, you must include the lapack Fortran library explicitly when compiling. 
e.g.
Windows (depending on the version of matlab, it's one of these.  libmwblas.lib may not be necessary.)
>> mex barsP_mat.c barsP_utils.c barsP_funcs.c com.c logspline.c ranlibsub.c libmwlapack.lib libmwblas.lib
>> mex barsP_mat.c barsP_utils.c barsP_funcs.c com.c logspline.c ranlibsub.c <your_matlab_directory>/extern/lib/win32/lcc/libmwlapack.lib 
64 bit Unix
>> mex barsP_mat.c barsP_utils.c barsP_funcs.c com.c logspline.c ranlibsub.c <your_matlab_directory>/bin/glnxa64/libmwlapack.so
Or possibly
>> mex barsP_mat.c barsP_utils.c barsP_funcs.c com.c logspline.c ranlibsub.c <your_matlab_directory>/bin/glnxa64/libmwlapack.so <your_matlab_directory>/bin/glnxa64/libmwblas.so

Starting with R2007, the BLAS and LAPACK libraries are separate .lib files.  They both need to be included. There is a problem with the windows LCC compiler for R2007b which will prevent you from compiling.  R2008a seems to work fine.

Execution
Once you have the binary, the following matlab commands will get you started.
>> load example.data
>> testBars(example(2:end,:),60)

For a little more information, type 
>> help barsP_mat

And also look inside testBars.m for an example of plotting the results.