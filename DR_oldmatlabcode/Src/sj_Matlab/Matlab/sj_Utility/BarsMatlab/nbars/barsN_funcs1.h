/* barsN_funcs1.h */
#ifndef _intbarsNfuncs
#define _intbarsNfuncs

/*
#define dgeqrf dgeqrf_
#define dormqr dormqr_
#define dpotrf dpotrf_
#define dpotrs dpotrs_
#define dgemm dgemm_
#define dpocon dpocon_
#define dpotri dpotri_
#define dtrtri dtrtri_
#define dtrtrs dtrtrs_
#define dtrmv dtrmv_


void dgeqrf( int *, int *, double *, int *, double *, double *, int *, int * );
void dormqr( char *, char *, int *, int *, int *, double *, int *, double *,    double *, int *, double *, int *, int * );
void dpotrf( char *, int *, double *, int *, int * );
void dpotrs( char *, int *, int *, double *, int *, double *, int *, int * );
void dgemm( char*, char*, int*, int*, int*, double*, double*, int*, double*,  int*, double*, double*, int* );
void dpocon( char *, int *, double *, int *, double *, double *, double *, int *, int * );
void dpotri( char *, int *, double *, int *, int *);
void dtrtri( char *, char *, int *, double *, int *, int *);
void dtrtrs( char *, char *, char *, int *, int *, double *, int *,
        double *, int *, int *);
void dtrmv(char *, char *, char *,int *, double *,int *,
	   double *, int *);
*/
/* Some constants for lapack and blas */
static  char Left = 'L';
static  char Right = 'R';
static  char Trans = 'T';
static  char NoTrans = 'N';
static  char Upper = 'U';
static  char Lower = 'L';
static  char Unit = 'U';
static  char NoUnit = 'N';
static  char No = 'N';
static int iZero = 0;
static int iOne = 1;
static int iTwo = 2;
static int iFour = 4;
static double Zero = 0.0;
static double One = 1.0;

#endif
