/* barsP_utils.h */
#ifndef _barsPutils
#define _barsPutils

/* String utility functions */
int streq(char *,char *);
int strbool(char *);
int isNotNone(char *);
void alert(char *);

/* Memory allocation */
int *tssivec(int);
int *tssivec0(int);
char *tsscvec(int);
char **tsscmat(int,int);
void free_cmat(char **, int);
double *tssdvec(int);
double *tssdvec0(int);
double **tssdmat(int,int);
double *tssdmat2(int,int);

/* Bounding functions */
void asSafePositive(double *);
double min(double,double);
double max(double,double);
int imin(int, int);
int imax(int, int);

/* Random Number Generation functions */
double rpoisson(double);
void gen_poisson(int, double*,double*,int, int);
void setseed(time_t *);
void setseed_mpi(time_t *, int);
void rmvnorm(int, double *);
double gammln(double);
double betaln(double,double);
double dbeta(double,double,double);

/* Knot utilities */
int check_cand(double, int, double *);
void set_default_knots( int, double *);
void boundaries(int,int,double*,double*);
void get_logspline_knots(int,double *,double *,double *,int *,double);

/* Generic utilities */
void swap (double **,double **);
int binary_search(double ,double *,int);

/* X Values utilities */
void normalize(int, double *, double *, double *, double *);
void normalize_bin(int, double *, double *, double *, double *, double *);
void lin_interp(double *, double *, int, int *, double *);
void mymult(double *, double *, int, int, int,int *, double *);

/* Basic Math */
void scale(int, double, double *);
void affine(int, double, double *, double *);
void transpose_matrix (int, int, double *, double *);
double squared_norm (double *,int);

/* IO functions */
void write_matrix(FILE *, double *, int, int);
void write_vector(FILE *, double *, int);

/* Poisson Regression functions */
double glmp_fitloglik(int, double *, double *, int);
void exb(int, int, double *, double *, double *);
double glmp_getsig2(int, int, double *, double *);
double glmp_getsigma(int, int, double *, double *);
void glmp_getz(int, double *, double *, double *);
void glmp_geth(int, int, double *, double *, double *);
void glmp_getmu(int, double *, double *);
double glmp_loglik(int, double *,double *, double *);
double glmp_loglik_3(int, double *, double *);

/* Parameters functions */
void cumsum(int, double *,double *);
double get_sum(int, double *);
double get_max(int, double *);
double get_min(int, double *);
double get_mode(int, double *, double *, double, double);
void get_mode_stats(int , double *, double *, double, 
         double, double *);
void get_mode_stats0(int, double *, double *, 
         double *);
void get_mustats(int, double *, double *, double, double, double *);
void get_mustats2(int, double *, double *, double, double, double *, int, int);
double getNorm(double *,int);
void conf_int(double *, int, double, double *);
void conf_int95(double *, int, double *);
void sum_stats(int, double *,double *);
void max_spline(int, double *, double *, double *,
           double *, double *, double *, double *);

/* Basis functions */
void outer_plus(int, int*,int, int*);
void basis_cw_trans(double *, double *, int *, int, int);
void spline_basis(double *,int *, int *, double *, int *, int *,
          double *,int *);
          
#endif

