#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "logspline.h"
#include "ranlib.h"
#include "barsP_utils.h"
#include "mex.h"

/*
barsP_utils.c

utility functions for barsP: Bars for Poisson Count Data

*/


/* BEGIN String utility functions */

int streq(char *a, char *b){
  return (strcmp(a,b) == 0);
}

int strbool(char *a){
  return ((streq(a,"true")) || (streq(a,"TRUE")) || (streq(a,"True")) ||
    (streq(a,"T")) || (streq(a,"1")) || (streq(a,"t")) || (streq(a,"yes")) ||
    (streq(a,"YES")) || (streq(a,"Y")) || (streq(a,"y")) || (streq(a,"Yes")));
}

int isNotNone(char *c){
    return (!((streq(c,"none")) || (streq(c,"NONE")) || (streq(c,"None"))));
  /* return ((strcmp(c,"none") == 0) ? 0 : 1);*/
}

void alert(char *c){
  printf(c);
  fflush(stdout);
}

/* END String utility functions */


/* BEGIN Memory allocation */

int *tssivec(nh)
int nh;
{
   return (int *)calloc((long) (nh),sizeof(int));
}

int *tssivec0(nh)
int nh;
{
  int *u;
  int i;
  u = tssivec(nh);
  for(i=0;i<nh;i++){
    u[i] = 0;
  }
  return u;
}

double *tssdvec(nh)
int nh;
{
   return (double *)calloc((long) (nh),sizeof(double));
}

double *tssdvec0(nh)
int nh;
{
  double *u;
  int i;
  u = tssdvec(nh);
  for(i=0;i<nh;i++){
    u[i] = 0.0;
  }
  return u;
}


char *tsscvec(nh)
int nh;
{
  return (char *)calloc((long) nh,sizeof(char));
}

char **tsscmat(nrh, nch)
int nrh,nch;
{
  int i;
  char **c;
  c = (char **) calloc((long) nrh, sizeof(char*));
  for(i=0;i<nrh;i++){
    c[i] = (char *) calloc((long) nch,sizeof(char));
  }
  return c;
}

void free_cmat(char **c, int nrh){
  int i;
  for(i=0;i<nrh;i++){
    free(c[i]);
  }
  free(c);
}

double **tssdmat(nrh,nch)
int nrh,nch;
{
   int i;
   double **m;
   m=(double **) calloc((long) nrh,sizeof(double*));
   for(i=0;i<nrh;i++)
      m[i]=(double *) calloc((long) nch,sizeof(double));
   return m;
}

double *tssdmat2(nrh,nch)
int nrh,nch;
{
return( (double *) calloc((long) (nrh*nch),sizeof(double)));
}

/* END Memory Allocation */

/* BEGIN Bounding functions */

void asSafePositive(double *d){
    double sd = sqrt(DBL_EPSILON);
    if ((*d) < sd) *d = sd;
    if ((*d) > (1.0/sd)) *d = (1.0/sd);
}

double min(double a, double b){
  return ((a < b) ? a : b);
}
double max(double a, double b){
  return ((a > b) ? a : b);
}
int imin(int a, int b){
  return ((a < b) ? a : b);
}
int imax(int a, int b){
  return ((a > b) ? a : b);
}

/* END Bounding functions */

/* BEGIN Random Number Generation functions */

double rpoisson(double lam){
  float lam2;
  long int rp2;
  lam2 = (float) lam;
  rp2 = ignpoi(lam2);
  return ((double)rp2);
}

void gen_poisson_data(int n, double *z, double *y, int bps, int trials){
  double mfactor;
  mfactor = (double)trials / (double)bps;
  while(n--){
    *y++ = rpoisson((*z++) * mfactor);
  }
}

void setseed(time_t *t1){
  struct tm *mytm;
  mytm = gmtime(t1);
  setall((long)((mytm->tm_sec + 1) * (mytm->tm_hour + 1) *
        (mytm->tm_mon + 1) *
        (mytm->tm_wday + 1) * (mytm->tm_mday)),
     (long)((mytm->tm_sec + 1) * (mytm->tm_min + 1) *
        (mytm->tm_hour + 1) * (mytm->tm_yday + 1)));
}

void setseed_mpi(time_t *t1,int prnum){
  int pnum = prnum + 1;
  struct tm *mytm;
  mytm = gmtime(t1);

  setall((long)((mytm->tm_sec + pnum) * (mytm->tm_hour + 1) *
        (mytm->tm_mon + 1) *
        (mytm->tm_wday + pnum) * (mytm->tm_mday)),
     (long)((mytm->tm_sec + 1) * (mytm->tm_min + pnum) *
        (mytm->tm_hour + 1) * (mytm->tm_yday + pnum)));
}

void rmvnorm(int n, double *z){
  while(n--){
    *z++ = snorm();
  }
}
double gammln(double a){
  double gx,gy,tmp,ser;
  static double cof[6] = {76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
                0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  gy = gx=a;
  tmp = gx + 5.5;
  tmp -= (gx+0.5)*log(tmp);
  ser = 1.00000000019015;
  for(j=0;j<6;j++) ser += cof[j]/++gy;
  return -tmp+log(2.5066282746310005*ser/gx);
}

double betaln(double a,double b){
  return (gammln(a) + gammln(b) - gammln((a+b)));
}

double dbeta(double x,double a,double b){
  return ((x <= 0.0) || (x >= 1.0)) ? 0.0 :
    exp((a-1.0)*log(x) + (b-1.0)*log(1.0-x) - betaln(a,b));
}

/* END Random Number Generation functions */

/* BEGIN Knot utilities */

int check_cand1(double cand, double used){
  return (fabs(cand - used) < DBL_EPSILON) ? 1 : 0;
}

int check_cand(double cand, int k, double *knots){
  int i,iout = 0;
  iout += check_cand1(cand,0.0);
  iout += check_cand1(cand,1.0);
  for(i=0;i<k;i++){
    iout += check_cand1(cand,knots[i]);
  }
  return iout;
}

void set_default_knots(int k, double *knots){
  int i;
  for(i=0;i<k;i++)
    *knots++ = ((double)(i+1))/((double)(k+1));
}

void boundaries (int k, int order, double *ki, double *ka){
  int od=order;
  while(od--){
    *ka++ = 0.0;
  }
  while(k--){
    *ka++ = *ki++;
  }
  while(order--){
    *ka++ = 1.0;
  }
}

double get_sum(int n, double *y){
  double dout = 0.0;
  while(n--){
    dout += *y++;
  }
  return dout;
}

void get_logspline_knots(int n,double *x,double *y,double *st_knots,int *k,double binwidth){

  int i,j,nn,*ls_ip,*ls_ad,m,ww;
  double *ls_coef,*ls_dp,*ls_logl,*ls_kts,w,*ls_y;
  double *my_kts,xi_sub,x_min,x_max;
/*  double *xxxx,*yyyy;
  FILE *out_file;
*/
  nn = (int) get_sum(n,y);

  ls_y = tssdvec0(n);
  ls_ip = tssivec0(7);
  ls_ad = tssivec0(1000);
  ls_coef = tssdvec0(nn);
  ls_dp = tssdvec0(3);
  ls_logl = tssdvec0(1000);
  ls_kts = tssdvec0(1000);

/*
   xxxx = tssdvec0(nn);
   yyyy = tssdvec0(nn);
*/
  my_kts = tssdvec0(3);
  ls_ip[0] = nn;
  ls_ip[1] = 0;
  ls_ip[2] = 3;
  ls_ip[3] = 0;
  ls_ip[4] = 1;
  ls_ip[5] = 1;
  ls_ip[6] = -1;
  m = 0;

  x_min = 100.0;
  x_max = -100.0;

  for(i=0;i<n;i++){
    ww = (int)(y[i] + 0.5);
    if ((x[i] > x_max) && (y[i] > 0.5)) x_max = x[i];
    if ((x[i] < x_min) && (y[i] > 0.5)) x_min = x[i];
    for(j=0;j<ww;j++){
      ls_coef[m] = ((double)j/(double)(ww - 1)) - 0.5;
      ls_coef[m] *= binwidth;
      ls_coef[m] += x[i];
      if ((ls_coef[m] > DBL_EPSILON) && (1.0 - ls_coef[m] > DBL_EPSILON))
  m++;
    }
  }

/*  printf("xmin = %le, xmax = %le\n",x_min,x_max); */

/*  memcpy(xxxx,ls_coef,sizeof(double)*m); */


  ls_dp[0] = -1.0;
  ls_dp[1] = 0.0;
  ls_dp[2] = 0.0;
  ls_kts[0] = ls_coef[0];
  ls_kts[2] = ls_coef[(m-1)];
  ls_kts[1] = 0.5 * (ls_kts[0] + ls_kts[2]);
  for(i=0;i<3;i++){
    my_kts[i] = ls_kts[i];
  }
  ls_ip[0] = m;
  ls_ip[2] = 0;

  nlogcensor(ls_ip,ls_coef,ls_dp,ls_logl,ls_ad,ls_kts);
/*
  out_file = fopen("abc","w");
  for(i=0;i<m;i++){
    yyyy[i] = ls_coef[0] + (xxxx[i] * ls_coef[1]);
  for(j=0;j<ls_ip[1];j++){
    if (xxxx[i] > ls_kts[j]) yyyy[i] += (ls_coef[(j+2)] * pow((xxxx[i] - ls_kts[j]),3));
  }
    yyyy[i] = exp(yyyy[i]) * (m / 60.0);
  fprintf(out_file,"%lf %lf\n",xxxx[i],yyyy[i]);
}
   fclose(out_file);
*/

  if ((ls_ip[0] > 0) && (ls_ip[0] < 100)){
    j = 1;
    st_knots[0] = my_kts[1];
  } else {
    x_min += DBL_EPSILON;
    x_max -= DBL_EPSILON;
    j = 0;
    for(i=0;i<ls_ip[1];i++){
      if ((ls_kts[i] > x_min) && (ls_kts[i] < x_max)){
/*      if ((ls_kts[i] > DBL_EPSILON) && (1.0 - ls_kts[i] > DBL_EPSILON)){ */
  st_knots[j++] = ls_kts[i];
      }
    }
    if (j == 0){
      j = 1;
      st_knots[0] = my_kts[1];
    }
  }
  *k = j;

/*
  free(xxxx);
  free(yyyy);
*/

  free(ls_y);
  free(ls_ip);
  free(ls_ad);
  free(ls_coef);
  free(ls_dp);
  free(ls_logl);
  free(ls_kts);
  free(my_kts);
}



/* END Knot utilities */

/* BEGIN Generic utilities */

void swap(double **a, double **b){
  double *c;

  c=*a;
  *a=*b;
  *b=c;
}

int binary_search(double candidate, double *array, int size){
  int first=0;
  int last=size-1;
  int mid;

  while(first < last){
    mid=(first+last)/2;
    if( array[mid]<= candidate &&  candidate <= array[mid+1]){
      return mid;
    }
    if(candidate>array[mid]){
      first=mid+1;
    }
    else{
      last=mid-1;
    }
  }
  return first;
}

/* END Generic utilities */

/* BEGIN X Values utilities */

void normalize(int n, double *xin, double *xout, double *xmin, double *xmax){


  int nn = n,i;
  double *x,*xx,x0,xn;

  x0 = *xin;
  xn = *xin;
  x = xin;
  xx = xout;
  while(nn--){
    if((*x) > xn) xn = *x;
    if((*x) < x0) x0 = *x;
    x++;
  }
  nn = n;
  x = xin;
  while(nn--){
    *xx++ = ((*x++) - x0)/(xn - x0);
  }
  for(i = 1;i<n;i++){
    if((xout[i] + DBL_EPSILON) <= xout[i-1]){
      mexPrintf("X values must be strictly increasing.\n");
      mexPrintf("x[%i] = %1.20le, x[%i] = %1.20le\n",i-1,xin[i-1],i,xin[i]);
      mexErrMsgTxt("\n");
    } else
    if(xout[i] <= xout[i-1]){
      printf("warning: multiple x values\n");
      xout[i] = xout[i-1] + DBL_EPSILON;
    }
  }
  *xmin = x0;
  *xmax = xn;
}

void duplicity(int n, double *x, double *xout, int *cout,
           int *nout){
  /* memory allocation should be done prior to entry */
  int changed = 0,i,*itmp,i2,n2;
  double *dtmp;

  n2 = n;
  dtmp = tssdvec(n);
  itmp = tssivec(n);
  memcpy(xout,x,sizeof(double)*n);
  for(i=0;i<n;i++){
    cout[i] = 1;
  }
  for(i=1;i<n2;i++){
    if (xout[i] <= xout[i-1]) {
      if (i < (n2 - 1)){
    i2 = n2 - i - 1;
    memcpy(dtmp,xout + i + 1,sizeof(double)*i2);
    memcpy(itmp,cout + i + 1,sizeof(int)*i2);
    memcpy(xout + i,dtmp,sizeof(double)*i2);
    memcpy(cout + i,itmp,sizeof(int)*i2);
      }
      cout[i-1]++;
      n2--;
      i--;
    }
  }
  free(dtmp);
  free(itmp);
  *nout = n2;
}

void normalize_bin(int n, double *bwout, double *xin,
           double *xout, double *xmin, double *xmax){


  int nn = n,i;
  double *x,*xx,x0,xn,bw;

  x0 = *xin;
  xn = *xin;
  x = xin;
  xx = xout;
  while(nn--){
    if((*x) > xn) xn = *x;
    if((*x) < x0) x0 = *x;
    x++;
  }

  bw = (xn - x0)/((double)(n - 1));

  for(i=1;i<n;i++){
    if ((xin[i] - xin[(i-1)] - bw) > (10.0 * DBL_EPSILON)){
      mexErrMsgTxt("Check the input data. The bins must have equal widths.\n");
    }
  }


/*  printf("x0 = %lf, xn = %lf\n",x0,xn); */
  x0 -= 0.5 * bw;
  xn += 0.5 * bw;
  *bwout = bw;

/* printf("%lf %lf %lf %lf\n",xin[0],xin[(n-1)],x0,xn); */
  nn = n;
  x = xin;
  while(nn--){
    *xx++ = ((*x++) - x0)/(xn - x0);
  }
/*  printf("%lf %lf\n",xout[0],xout[(n-1)]);  */
  for(i = 1;i<n;i++){
    if((xout[i] + DBL_EPSILON) <= xout[i-1]){
      mexPrintf("X values must be strictly increasing.\n");
      mexPrintf("x[%i] = %1.20le, x[%i] = %1.20le\n",i-1,xin[i-1],i,xin[i]);
      mexErrMsgTxt("\n");
    }
    if(xout[i] <= xout[i-1]){
      printf("warning: multiple x values\n");
      xout[i] = xout[i-1] + DBL_EPSILON;
    }
  }
  *xmin = x0;
  *xmax = xn;
}

void lin_interp(double *x0, double *x1, int n, int *m,
         double *c11){
  double *firstx = x1;
  int j = 0;
/* j = 0 inside, so that m always has RELATIVE position */
  while(n--){
    j = 0;
    while (*x1 < *x0) {x1++; j++;}
    if (x1 > firstx) {x1--; j--;}
    *m++ = j;
    if (*x1 < *x0){
      *c11++ = (*(x1+1) - *x0)/(*(x1+1) - *x1);
    } else {
      *c11++ = 1.0;
    }
    x0++;
  }
}

void mymult(double *t_basis, double *a_basis, int n0, int napp, int p,
        int *mm1, double *c11){
int i,j,n,*m2;
double *c2,*a2;

while(p--){
  c2 = c11;
  m2 = mm1;
  a2 = a_basis;
  n = n0;
  while(n--){
    a2+=(*m2++);
    *t_basis++ = (*c2) * (*a2) + (1.0 - (*c2)) * (*(a2+1));
    c2++;
  }
  a_basis+=napp;
}
}

/* END X Values utilities */


/* BEGIN Basic Math */

void scale(int n, double a, double *x){
  while(n--){
    *x++ *= a;
  }
}


void affine(int n, double a, double *x, double *b){
  /* x -> ax + b, a is a scalar. */

  while(n--){
    *x *= a;
    *x++ += *b++;
  }
}

void transpose_matrix (int r, int c, double *mat_in, double *mat_out){
  int i,j;
  for(i=0;i < r;i++){
    for(j=0;j<c;j++){
      mat_out[c*i + j] = mat_in[i + r*j];
    }
  }
}

double squared_norm(double *z, int n){
  double r = 0.0;
  while(n--){
    r += pow((*z),2);
    z++;
  }
  return r;
}



/* END Basic Math */

/* BEGIN IO functions */

void write_matrix(FILE *f, double *b, int nr, int nc){
  int i,j;
  for(i=0;i<nr;i++){
    for(j=0;j<nc;j++){
      fprintf(f,"%.12le ", b[i + nr * j]);
    }
    fprintf(f,"\n");
  }
  fprintf(f,"\n");
}

void write_vector(FILE *f, double *b, int nr){
  int i;
  for(i=0;i<nr;i++){
    fprintf(f,"%.12le ", b[i]);
  }
  fprintf(f,"\n");
}

/* END IO functions */

/* BEGIN Poisson Regression functions */

double glmp_getsig2(int n, int p, double *y, double *mu){
  double v = 0.0;
  int m = n-p;
  while(n--){
    v += (pow((*y) - (*mu),2) / (*mu));
    y++;
    mu++;
  }
  return (v/m);
}

double glmp_getsigma(int n, int p, double *y, double *mu){
  return sqrt(glmp_getsig2(n,p,y,mu));
}

void glmp_getz(int n, double *g, double *y, double *mu){
  /*
  int n1 = n;
  double *g1,*y1,*mu1;
  g1 = g;
  y1 = y;
  mu1 = mu;
*/
  while(n--){
    *g++ = log(*mu) + ((*y++) - (*mu))/(*mu);
/*    *g++ = log(*mu) + ((*y++)/(*mu)) - 1; */
    mu++;
  }
}

void glmp_geth(int n, int p, double *x, double *mu, double *h){

  int n0;
  double *mu1;

  /* memcpy(h,x,(n*p)); */
  while(p--){
    n0 = n;
    mu1 = mu;
    while(n0--){
      *h++ = (*x++) * (*mu1++);
    }
  }
}

void glmp_getmu(int n, double *theta, double *mu){
  while(n--){
    *mu++ = exp(*theta++);
  }
}

void exb(int n, int p, double *x, double *b, double *y){
  double *v,dm;
  int n2=n;

  dm = log(0.1 * DBL_MAX);
  v = y;
  while(n2--){
    *v++ = 0.0;
  }
  while(p--){
    n2 = n;
    v = y;
    while(n2--){
      *v++ += ((*x++) * (*b));
    }
    b++;
  }
  while(n--){
    *y = ((*y) > dm) ? (0.1 * DBL_MAX) : exp(*y);
/*    *y = exp(*y);
    if ((*y) > 100000) *y = 1000000.0;
*/
    y++;
  }
}

double glmp_loglik(int n, double *theta, double *mu, double *y){

  double v = 0.0;
  while(n--){
    v += (((*y++) * (*theta++)) - (*mu++));
  }
  return v;
}

double glmp_loglik_3(int n, double *mu, double *y){
  double v = 0.0,w = 0.0;
  while(n--){
    w = log(*mu);
    v += (((*y++) * w) - (*mu++));
  }
  return v;
}


double glmp_fitloglik(int n, double *smu, double *y, int iter){
  double v = 0.0;
  while(n--){
    v += (((*y++) * (log(*smu) - log(iter))) - ((*smu)/iter));
    smu++;
  }
  return v;
}

/* END Poisson Regression functions */

/* BEGIN Parameter functions */

double get_max(int n, double *x){
  double maxx;

  maxx = *x;
  while(n--){
    if(*x > maxx){
      maxx = *x;
    }
    x++;
  }
  return maxx;
}
double get_min(int n, double *x){
  double minx;

  minx = *x;
  while(n--){
    if(*x < minx){
      minx = *x;
    }
    x++;
  }
  return minx;
}

double get_mode(int n, double *x, double *y, double x1, double xn){
  double *vy,*vx,vz;
  vy = y;
  vx = x;
  while(n--){
    if (*y > *vy){
      vx = x;
      vy = y;
    }
    x++;
    y++;
  }
  vz = x1 + ((*vx) * (xn - x1));
  return vz;
}

void get_mustats(int n, double *x, double *y, double x1,
         double xn, double *z){
  double *mxy,*mxx,*mny,*mnx;
  mxx = x;
  mnx = x;
  mxy = y;
  mny = y;
  while(n--){
    if (*y > *mxy){
      mxy = y;
      mxx = x;
    }
    if (*y < *mny){
      mny = y;
      mnx = x;
    }
    x++;
    y++;
  }
  *z++ = x1 + ((*mnx) * (xn - x1));
  *z++ = *mny;
  *z++ = x1 + ((*mxx) * (xn - x1));
  *z = *mxy;
}

void get_mode_stats(int n, double *x, double *y, double x1,
         double xn, double *z){
  double *mxy,*mxx;
  mxx = x;
  mxy = y;
  while(n--){
    if (*y > *mxy){
      mxy = y;
      mxx = x;
    }
    x++;
    y++;
  }
  *z++ = x1 + ((*mxx) * (xn - x1));
  *z = *mxy;
}
void get_mode_stats0(int n, double *x, double *y,
         double *z){
  double *mxy,*mxx;
  mxx = x;
  mxy = y;
  while(n--){
    if (*y > *mxy){
      mxy = y;
      mxx = x;
    }
    x++;
    y++;
  }
  *z++ = *mxx;
  *z = *mxy;
}

void get_mustats2(int n, double *x, double *y, double x1,
          double xn, double *z, int bps, int trials){
  double mfac = (double) bps / (double) trials;
  get_mustats(n,x,y,x1,xn,z);
  z[1] *= mfac;
  z[3] *= mfac;
}

void get_mode_stats2(int n, double *x, double *y, double x1,
             double xn, double *z, int bps, int trials){
  double mfac = (double) bps / (double) trials;
  get_mode_stats(n,x,y,x1,xn,z);
  z[1] *= mfac;
}

void cumsum(int n, double *result, double *x){
  while(n--)
    *result++ += *x++;
}

double getNorm(double *A, int n){
  double nout = 0.0,cval;
  int n0 = n,n1;
  while(n0--){
    cval = 0.0;
    n1 = n;
    while(n1--){
      cval += fabs(*A++);
    }
    if (cval > nout) nout = cval;
  }
  return nout;
}


void L2_delta_one(int n, double h, double *y, double *mu, double *del){
  *del = 0.0;
  *del -= 0.5 * pow((*y) - (*mu),2);
  while(n--){
    *del += pow((*y++) - (*mu++),2);
  }
  *del -= 0.5 * pow((*(y-1)) - (*(mu-1)),2);
  *del *= h;
  *del = sqrt(*del);
/* uses trapezoid method for integration. */
}

void L2_deltas(int n, int p, double h, double *ymat, double *mu, double *del){
  while(p--){
    L2_delta_one(n,h,ymat,mu,del++);
    ymat += n;
  }
}



void get_post_mean(int n, int p, double *ymat, double *mu){
  double *m0;
  int n0 = n,p0 = p;
  m0 = mu;
  while(n0--){
    *m0++ = 0.0;
  }
  while(p0--){
    n0 = n;
    m0 = mu;
    while(n0--){
      *m0++ += *ymat++;
    }
  }
  while(n--){
    *mu++ /= ((double)p);
  }
/*
  m0 = mu;
  while(n--){
    *m0 = (*m0)/ ((double)p);
    m0++;
  }
*/
}

void sum_stats(int n, double *X, double *results){
  int i;
  double sx = 0,sxx = 0,n2;
  n2 = (double)n;
  for(i=0;i<n;i++){
    sx += *X;
    sxx += ((*X) * (*X));
    X++;
  }
  sxx = (sxx - (sx * sx / n2))/(n2-1.0);
  sx = (sx / n2);
  *results++ = sx;
  *results = sxx;
}

void select_one(double *a, int n, int k){
  /* partial sort of a such that a[0] ... a[k-1] < a[k] < a[k+1] .. a[n-1] */
  double v,t;
  int i,j,l,r,j2;
  l = 0;
  r = n-1;
  while (r > l){
    v = *(a+r);
    i = l - 1;
    j = r;
    for(;;){
      while (*(a+(++i)) < v){
      }
      j2 = 1;
      while(j2){
    if (j){
      if (*(a+(--j)) <= v) j2 = 0;
    } else {
      j2 = 0;
    }
      }
/*
      while (*(a+(--j)) > v){
      }
*/
      if (i >= j) break;
      t = *(a+i);
      *(a+i) = *(a+j);
      *(a+j) = t;
    }
    t = *(a+i);
    *(a+i) = *(a+r);
    *(a+r) = t;
    if (i >= k) r = i - 1;
    if (i <= k) l = i+1;
  }
}

void conf_int95(double *a, int n, double *ci){
  double m,v;
  sum_stats(n,a,ci);
  m = *ci;
  v = *(ci+1);
  v = sqrt(v) * 1.96;
  *ci++ += v;
  *ci = m - v;
}

void conf_int(double *a, int n, double perc, double *ci){
  double low0,up0,perc2;
  int low1,low2,up1,up2;

  perc2 = (1.0 - perc) * 0.5;

  low0 = perc2 * (double)(n-1);
  up0 = (1.0 - perc2) * (double)(n-1);
  low1 = (int) low0;
  low2 = low1 + 1;
  up1 = (int) up0;
  up2 = up1 + 1;

  if (up2 >= n) up2 = n-1;
  select_one(a,n,low2);
  select_one(a,low2,low2 - 1);
/*  select_one(a+low2+1,n-low1,up1-low2-1); */
  select_one(a+low2+1,n-low2-1,up1-low2-1);
  select_one(a+up2,n-up2,0);
  *ci++ = a[low1] * ((double)low2 - low0) + a[low2] * (low0 - (double)low1);
  *ci = a[up1] * ((double)up2 - up0) + a[up2] * (up0 - (double)up1);
}

void conf_int_one(double *a, int n, double perc, double *ci){
  double up0;
  int up1,up2;

  up0 = perc * (double)(n-1);
  up1 = (int) up0;
  up2 = up1 + 1;
  select_one(a,n,up1);
  select_one(a+up2,n - up2,0);
  *ci = a[up1] * ((double)up2 - up0) + a[up2] * (up0 - (double)up1);
}

/* END Parameter functions */


/* BEGIN Interpolating spline functions */

double fspline(double x){
  return (pow(x,3) - x);
}

double fprime(double x){
  return (3.0*x*x - 1.0);
}

double nstep_spline(int n, double v, double *x, double *y,
            double *u, double *p){
  double t,r1,r2;
  int i;
  if (v < x[0] + DBL_EPSILON) return 0.0;
  else if (v > x[(n-1)] - DBL_EPSILON) return 0.0;
  else {
    i = binary_search(v, x, n);
    t = (v - x[i])/u[i];
    r1 = ((y[i+1] - y[i])/u[i]) +
      u[i] * (fprime(t) * p[i+1] - (fprime((1.0-t)) * p[i])) / 6.0;
    r2 = (t * p[i+1] + (1.0-t) * p[i]);
    return r1/r2;
  }
}

void makespline1(int n, double *x, double *y,
        double *u, double *d, double *w, double *p){
  int i;
  d[0] = 0.0;
  w[0] = 0.0;
  u[0] = x[1] - x[0];
  p[0] = 0.0;
  p[n-1] = 0.0;
  for(i=1;i<n-1;i++){
    d[i] = 2.0 * (x[i+1] - x[i-1]);
    u[i] = x[i+1] - x[i];
    w[i] = 6.0 * ((y[i+1] - y[i])/u[i] - (y[i] - y[i-1])/u[i-1]);
  }
  for(i=1;i<n-2;i++){
    w[i+1] -= w[i]*u[i]/d[i];
    d[i+1] -= u[i]*u[i]/d[i];
  }
  for(i=n-2;i > 0;i--){
    p[i] = (w[i] -u[i] * p[i+1])/d[i];
  }
}


double evalspline(int n, double v, double *x, double *y,double *u, double *p){
  double t;
  int i;
  if (v < x[0] + DBL_EPSILON) return y[0];
  else if (v > x[(n-1)] - DBL_EPSILON) return y[(n-1)];
  else {
    i = binary_search(v, x, n);
    t = (v - x[i])/u[i];
    return t * y[i+1] + (1.0-t)*y[i] +
      u[i]*u[i] * (fspline(t)*p[i+1] + fspline(1.0-t)*p[i])/6.0;
  }
}



void max_spline(int n, double *x, double *y, double *xy_start,
           double *u, double *d, double *w, double *p){
/* forms an interpolating spline based on x and y,
   then maximizes it based on the starting value given in xy_start.
   u,d,w,p are for workspace and store the spline information.
   new max values are put into xy_start.
*/


  double x0,h,y0;
  int i=0;
  makespline1(n,x,y,u,d,w,p);
  x0 = *xy_start;
  h = 1.0;
  while((i++ < 20) && (fabs(h) > 0.000001)){

   h = nstep_spline(n,x0,x,y,u,p);
   x0 -= h;
  }
  if (x0 < x[0] + DBL_EPSILON) {
    x0 = x[0];
    y0 = y[0];
  } else if (x0 > x[(n-1)] - DBL_EPSILON){
    x0 = x[(n-1)];
    y0 = y[(n-1)];
  } else {
    y0 = evalspline(n,x0,x,y,u,p);
  }
  if (y0 > *(xy_start + 1)){
    *xy_start = x0;
    *(xy_start + 1) = y0;
  }
}


void vevalspline(int n, double *x, double *y, double *u, double *p,
         int m, double *vx, double *vy){
  /* assumes that x[0] <= vx[0], x[n] >= vx[m], and both x and vx are
     sorted, without repeats */
  double t,x0,xn,y0,yn;
  x0 = x[0] + DBL_EPSILON;
  y0 = y[0];
  xn = x[(n-1)] - DBL_EPSILON;
  yn = y[(n-1)];


  while(m--){

    if (*vx < x0){
      vx++;
      *vy++ = y0;
    } else if (*vx > xn){
      vx++;
      *vy++ = yn;
    } else {
      while((*vx > *(x+1)) && (n>2)){
  n--;
  x++;
  y++;
  u++;
  p++;
      }
      t = (*vx++ - *x)/(*u);
      *vy++ = t * (*(y+1)) + (1.0-t)*(*y) + (*u)* (*u) *
  (fspline(t)*(*(p+1)) + fspline(1.0 - t) * (*p))/6.0;
    }
  }
}

double d1evalspline(int n, double v, double *x, double *y,
            double *u, double *p){
/* returns the first derivative at v */
  double t;
  int i;
  if (v < x[0] + DBL_EPSILON) return y[0];
  else if (v > x[(n-1)] - DBL_EPSILON) return y[(n-1)];
  else {
    i = binary_search(v, x, n);
    t = (v - x[i])/u[i];
    return y[i+1] - y[i] + u[i]*u[i] * (fprime(t) * p[i+1] - (fprime(1.0 - t) * p[i]))/6.0;
  }
}

double d2evalspline(int n, double v, double *x, double *y,
            double *u, double *p){
/* returns the second derivative at v */
  double t;
  int i;
  if (v < x[0] + DBL_EPSILON) return y[0];
  else if (v > x[(n-1)] - DBL_EPSILON) return y[(n-1)];
  else {
    i = binary_search(v, x, n);
    t = (v - x[i])/u[i];
    return u[i] * u[i] * (t*p[i+1] + (1.0 - t)*p[i]);
  }
}

/* END Interpolating Spline functions */



/* BEGIN Basis functions */

void outer_plus(int dx, int *Y,int dy, int *xy){
  int i, j;
  for (i=0; i<dx; i++){
    for(j=0; j<dy; j++){
      xy[i*dy+j]=i+Y[j];
    }
  }
}

void basis_cw_trans(double *o_basis, double *n_basis, int *index,
            int order, int n){
  int i,j,a;

  /* new basis is column_wise, old basis is row_wise.  */

  for(i=0;i<n;i++){
    for(j=0;j<order;j++){
      a= i + n * index[j*n+i];
      n_basis[a]=o_basis[i*order + j];
    }
  }
}








/* The following routines are from Bates and Venables, 1998. */

/*  Routines for manipulating B-splines.  These are intended for use with
 *  S or S-PLUS or R.
 *
 *     Copyright (C) 1998 Douglas M. Bates and William N. Venables.
 *
 *
 * This program is free software; you can distribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * These functions are distributed in the hope that they will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * The text of the GNU General Public License, version 2, is available
 * as http://www.gnu.org/copyleft or by writing to the Free Software
 * Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 * The routines are loosely based on the pseudo-code in Schumaker (Wiley,
 * 1981) and the CMLIB library DBSPLINES.
 */

static double *ldel, *rdel;
static int orderm1;


static void
diff_table(double *ti, double x, int ndiff)
{
  register double *r = rdel, *l = ldel, *dpt = ti;

  while (ndiff--) {
    *r++ = *dpt++ - x;
    *l++ = x - *--ti;
  }
}

static void
basis_funcs(double *ti, double x, double *b)
{
  int j, r;
  double saved, term;

  diff_table(ti, x, orderm1);
  b[0] = 1.;
  for (j = 1; j <= orderm1; j++) {
    saved = 0.;
    for (r = 0; r < j; r++) {
      term = b[r]/(rdel[r] + ldel[j - 1 - r]);
      b[r] = saved + rdel[r] * term;
      saved = ldel[j - 1 - r] * term;
    }
    b[j] = saved;
  }
}

static double
evaluate(double *ti, double x, double *a, int nder)
{
  register double *lpt, *rpt, *apt;
  register int inner;
  int outer = orderm1;

  while(nder--) {
    for(inner = outer, apt = a, lpt = ti - outer; inner--; apt++, lpt++)
      *apt = outer * (*(apt + 1) - *apt)/(*(lpt + outer) - *lpt);
    outer--;
  }
  diff_table(ti, x, (int) outer);
  while(outer--)
    for(apt = a, lpt = ldel + outer, rpt = rdel, inner = outer + 1;
    inner--; lpt--, rpt++, apt++)
      *apt = (*(apt + 1) * *lpt + *apt * *rpt)/(*rpt + *lpt);
  return(*a);
}

void
spline_value(double *knots, double *coeff, int *ncoeff,
         int *order, double *x, int *nx, int *derivv,
         double *y)
{
  int n = *nx;
  double *a, *last = knots + *ncoeff;

  a = (double *) calloc((long) *order, sizeof(double));
  orderm1 = *order - 1L;    /* allocate difference tables */
  rdel = (double *) calloc((long) orderm1, sizeof(double));
  ldel = (double *) calloc((long) orderm1, sizeof(double));

  knots += *order;      /* First *order knots must be <= all x's */
  while(n--) {
    while(knots <= last && *knots <= *x) {knots++; coeff++;}
    memcpy(a, coeff, *order);
    *y++ = evaluate(knots, *x++, a, (int) *derivv);
  }
  free(ldel); free(rdel); free(a);
}

void spline_basis(double *knots, int *ncoeff, int *order,
          double *xvals, int *derivs, int *nx,
          double *basis, int *offsets)
{               /* evaluate the non-zero B-spline basis */
                /* functions (or their derivatives) at */
                /* xvals.  */
  int n = *nx, i, j;
  double *dpt, *coeff, *last = knots + *ncoeff;

  orderm1 = *order - 1L;
  rdel = (double *) calloc((long) orderm1, sizeof(double));
  ldel = (double *) calloc((long) orderm1, sizeof(double));
  coeff = (double *) calloc((long) *order, sizeof(double));
  dpt = (knots += *order);  /* first *order knots must be <= all xvals */
  for( ; n--; xvals++, derivs++) {
    while(dpt < last && *dpt <= *xvals) dpt++;
    if (*derivs) {      /* slow method for derivatives */
      for(i = 0; i < *order; i++) {
    for(j = 0; j < *order; j++) coeff[j] = 0;
    coeff[i] = 1;
    *basis++ = evaluate(dpt, *xvals, coeff, (int) *derivs);
      }
    }
    else {
      basis_funcs(dpt, *xvals, basis); /* fast method for value */
      basis += *order;
    }
    *offsets++ = (int)(dpt - knots);
  }
  free(ldel); free(rdel); free(coeff);
}

/* END Basis functions */
