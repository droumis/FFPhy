/* $Id: hippoIO.h,v 1.7 2008/08/20 21:44:11 chengs Exp $
   
   Sen Cheng, 2004/05/12
   Header file for io functions.
   
*/

#ifndef HIPPOIO_H
#define HIPPOIO_H

#include <iostream>
#include <stdexcept>
using namespace std;

/* read vector from binary file that only contains one vector of floats into
   data 
   returns: the number of elements read in
*/
int readdata(char *filename, float *data, int nsamples);

/* returns the number of elements of the specificied size in the binary file */
int getnelements(char *filename, int esize);

#define hippo_Virtual							\
    cerr << "Error: called virtual function '" << __FUNCTION__  << "'";	\
    cerr << "\n\tin file " __FILE__ ", l. " << __LINE__  << endl; \
    throw 5;
#define hippo_Empty							\
    cerr << "Error:  function not defined, yet.";			\
    cerr << ", " __FILE__ ", l. " << __LINE__  << endl;
#define hippo_Warning(str)\
    cerr <<  str << ", " __FILE__ ", l. " << __LINE__ << endl;
#define hippo_Print(var)						\
    cerr <<  #var "= " << var;						\
    cerr <<  ", " __FILE__ ", l. " << __LINE__ << endl;
#define hippo_Error(msg,arg)						\
    do{									\
        cerr << "Error: " << msg << arg;						\
        cerr << " " __FILE__ ", l. " << __LINE__  << "." <<endl;	\
        throw 1; \
    } while(0);
#define hippo_ErrorMsg(msg) \
    do{									\
        cerr << "Error: " << msg ;						\
        cerr << " " __FILE__ ", l. " << __LINE__  << "." << endl;	\
        throw 1; \
    } while(0);

#define hippo_Exit \
    do{									\
        cerr << "Exit  in " __FILE__ ", l. " << __LINE__  << endl;	\
        throw 1; \
    } while(0);

#define hippo_Msg(msg) \
    do{									\
        cerr << "    " __FILE__ ", l. " << __LINE__ << ": " << msg << endl; \
    } while(0);

#ifdef DEBUG
#define hippo_Message(string)				\
    cerr << "Message: "  string;				\
    cerr << ", " __FILE__ ", l. " << __LINE__  << endl;

#define hippo_Assert(cond, string)					\
    if (!(cond)) {							\
        cerr << "assertion failed: " << #cond  \
            << "\n  " << string		\
            << "\n  in " __FILE__ ", l. " << __LINE__  << endl; \
        throw 1; \
    }
#define hippo_Mark						\
    cerr <<  "mark:  " __FILE__ ", l. " << __LINE__ << endl;

#else // ifdef DEBUG

#define hippo_Message(string) 
#define hippo_Assert(cond, string) 
#define hippo_Mark 

#endif // ifdef DEBUG


#define hippo_CAlloc(ptr,N,TYPE) do {					\
	(ptr)= (TYPE *) calloc((N), sizeof(TYPE));			\
	hippo_Assert(ptr, "memory allocation failed");			\
    } while(0)
					  

// return determinant of 2x2 matrix A{11,12,21,22}
double calcMatrixDet(double A[]);


// return inverse of 2x2 matrix A in Ainv. A is vector with A{11,12,21,22}
void calcMatrixInv(double A[], double Ainv[]);

#endif  // HIPPOIO_H
