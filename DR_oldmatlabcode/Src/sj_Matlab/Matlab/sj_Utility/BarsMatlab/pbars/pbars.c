/* $Revision: 1.1 $*/
// pbars.c
// Automatically generated by Matlab AppWizard version 1.0
//
// This is the gateway routine for a MATLAB Math/Graphics Library-based
// C MATLAB MEX File.


/* This mex file provides a wrapper to run Bayesian Adaptive Regression
   Splines (BARS) as put forth in:

   Dimatteo, I., Genovese, C.R., and Kass, R.E. (2001)  Bayesian
   curve-fitting with free-knot splines, Biometrika, 88: 1055-1071.

   The associated source code uses Charles Kooperberg's implementation
   of LOGSPLINE, routines for manipulating B-Splines written by Bates
   and Venables (included in the release of R, 2003), and for random
   number generation obtained from Ranlib (Brown and Lovato, 1996).

   This file was adapted from source code written by Garrick Wallstrom,
   Jeffrey Liebner and Robert E. Kass in order to compile in Windows. It
   incorporates calls to BLAS / LAPACK libraries that are bundled with
   Matlab. For details see:

   Wallstrom, G., Liebner, J., and Kass, R.E. (2005) An implementation
   of Bayesian Adaptive Regression Splines (BARS) with S and R wrappers,
   under revision for Journal of Statistical Software.

   John Curtis 2005
 */

#include "mex.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "ranlib.h"
#include "barsP_utils.h"
#include "barsP_funcs.h"

#define X_DAT 0
#define Y_DAT 1
#define LGSP  2
#define IKNT  3
#define PR    4
#define PRP   5
#define BRNIN 6
#define SMS   7
#define TU    8
#define CNT   9
#define FTS  10
#define CNF  11
#define TRL  12
#define BNS  13
#define VRB  14


int VERBOSE;

/* see barslib.h for struct definitions */
/*****************************************************************************/
void filed(double *x, int *y,int *trials,int *n){
  int i,j;
  FILE *f;
  j = *n;
  f=fopen("bars_points","w");
  fprintf(f,"%i %i\n", *n,*trials);
  for (i=0;i<j;i++){
    fprintf(f,"%1.12f %i\n",x[i],y[i]);
  }
 fclose(f);
}

/*****************************************************************************/
void filed2(int *x, double *y,int *n){
  int i,j;
  FILE *f;
  j = *n;
  f=fopen("prior_file","w");
  for (i=0;i<j;i++){
    fprintf(f,"%i %1.12f\n",x[i],y[i]);
  }
 fclose(f);
}

/*****************************************************************************/
void priorsetup(int *burnin, int *sims, int *iknots,
      double *tau,double *c,double *conf,int *grid){
  FILE *f;
  f=fopen("bars_params","w");
  fprintf(f,"SET burn-in_iterations = %i\n",*burnin);
  fprintf(f,"SET sample_iterations = %i\n", *sims);
  fprintf(f,"SET initial_number_of_knots = %i\n",*iknots);
  fprintf(f,"SET beta_iterations = 3\n");
  fprintf(f,"SET beta_threshhold = -10.0\n");
  fprintf(f,"SET proposal_parameter_tau = %8.4f\n",*tau);
  fprintf(f,"SET reversible_jump_constant_c = %1.4f\n",*c);
  fprintf(f,"SET confidence_level = %8.4f\n",*conf);
  fprintf(f,"SET number_of_grid_points = %i\n",*grid);
  fprintf(f,"SET sampled_knots_file = samp_knots\n");
  fprintf(f,"SET sampled_mu_file = samp_mu\n");
  fprintf(f,"SET sampled_mu-grid_file = samp_mugrid\n");
  fprintf(f,"SET sampled_params_file = samp_params\n");
  fprintf(f,"SET summary_mu_file = summ_mu\n");
  fprintf(f,"SET summary_mu-grid_file = summ_mugrid\n");
  fprintf(f,"SET summary_params_file = summ_params\n");
  fclose(f);
}

/*****************************************************************************/
void paramp(double *priorparam){
  FILE *f;
  f=fopen("bars_params","a");
  fprintf(f,"SET prior_form = Poisson\n");
  fprintf(f,"SET Poisson_parameter_lambda = %8.4f\n",*priorparam);
  fclose(f);
}

/*****************************************************************************/
void paramu(int *upper,int *lower){
  FILE *f;
  f=fopen("bars_params","a");
  fprintf(f,"SET prior_form = Uniform\n");
  fprintf(f,"SET Uniform_parameter_U = %i\n",*upper);
  fprintf(f,"SET Uniform_parameter_L = %i\n",*lower);
  fclose(f);
}

/*****************************************************************************/
void paramuser(){
  FILE *f;
  f=fopen("bars_params","a");
  fprintf(f,"SET prior_form = User\n");
  fclose(f);
}

/*****************************************************************************/
void write_user_priorparam(double *priorparam, int nprm){
  int ndx;
  FILE *f;
  f = fopen("user_prior_params","w");
  for (ndx = 0; ndx < nprm; ndx++)
  {
    fprintf(f, "%i %.8f", (int)priorparam[ndx], priorparam[ndx + nprm]);
  }
  fclose(f);
}

/*****************************************************************************/
void write_datfile(int *bns,int *trl, double *xin, __int32 *yin, int ndat){
  int ndx;
  FILE *f;
  f = fopen("input_data","w");
  fprintf(f, "%i %i\n", *bns, *trl);
  for (ndx = 0; ndx < ndat; ndx++)
  {
    fprintf(f, "%.15f %i\n", xin[ndx], yin[ndx]);
  }
  fclose(f);
}

/*****************************************************************************/
void logfalse(){
  FILE *f;
  f=fopen("bars_params","a");
  fprintf(f,"SET Use_Logspline = false\n");
  fprintf(f,"SET verbose = false\n");
  fclose(f);
}

/*****************************************************************************/
void logtrue(){
  FILE *f;
  f=fopen("bars_params","a");
  fprintf(f,"SET Use_Logspline = true\n");
  fprintf(f,"SET verbose = false\n");
  fclose(f);
}

/*****************************************************************************/
void setBinnedInfo(Bars_BinnedData *bd,double bwidth, int tr){
  double sf = bwidth * (double)tr;
  (*bd).bin_width = bwidth;
  (*bd).trials = tr;
  (*bd).scale_factor = (sf > 0.0) ? (1.0/sf) : 0.0;
}

/*****************************************************************************/
int read_params_file(char **param_name, char **param_value, char *fname){
  int fval, pnum = 0, isset=0;
  FILE *the_file;
  char *prefix;


  the_file = fopen(fname,"r");
  prefix = tsscvec(1000);
  do{
    fval = fscanf(the_file,"%s %s = %s\n",prefix, param_name[pnum],param_value[pnum]);
    isset = 0;
    if (fval == 3){
      isset = streq(prefix,"SET");
      if (isset) pnum++;
    }
  } while((fval != EOF) && (pnum < MAX_PARAMS));
  if (VERBOSE){
    printf("%i parameters read.\n",pnum);
  }
  if (pnum == MAX_PARAMS)
    printf("WARNING: maximum number of parameters read before end of file.\n");
  fclose(the_file);
  free(prefix);
  return pnum;
}

/*****************************************************************************/
Bars_BinnedData *read_data_file(char *fname, int nf){
  FILE *data_file;
  Bars_X_Grid *xg;
  Bars_BinnedData *bd;
  int n,tr,i,fval;
  double bw,*y,*x,*x_norm,xmin,xmax,delta,d0;

/* can add a parameter to allow for non-central bin markers.
   parameter should be stored in binned_data and passed
   to normalize_bin to adjust the min and max appropriately.
*/

  data_file = fopen(fname,"r");
  fval = fscanf(data_file,"%i %i\n",&n,&tr);
  if (fval != 2) {
	  fclose(data_file);
	  mexErrMsgTxt("Data file format error: 1st line should contain two integers: number of bins, and the number of trials.\n\n");
  }
  bd = New_BinnedData(n,nf);
  y = (*bd).y;
  x = (*bd).x_raw;
  for(i=0;i<n;i++){
    fval = fscanf(data_file,"%lf %lf\n",x+i,y+i);
    if (fval != 2){
      if (fval >= 0) {
        fclose(data_file);
        mexPrintf("Data file format error: Data for bin %i contains format errors.\n",(i+1));
        mexErrMsgTxt("\n");
      }
      if (fval == EOF) {
        fclose(data_file);
        mexPrintf("Data file format error: EOF reached before data for all %i bins could be read.\n",n);
        mexErrMsgTxt("\n");
      }
    } else {
      if (y[i] < 0.0){
        fclose(data_file);
        mexPrintf("Data file error: Count for bin %i is negative.\n",(i+1));
        mexErrMsgTxt("\n");
      }
    }
  }
  fclose(data_file);
  xg = (*bd).xg;
  x_norm = (*xg).x;
  normalize_bin(n,&bw,x,x_norm,&xmin,&xmax);
  setBinnedInfo(bd,bw,tr);
  (*bd).x_rawmin = xmin;
  (*bd).x_rawmax = xmax;
  delta = x_norm[(n-1)] - x_norm[0];
  d0 = x_norm[0];
  xg = (*bd).fg;
  x_norm = (*xg).x;
  for(i=0;i<nf;i++){
    x_norm[i] = (((double)i)/((double)(nf - 1)) * delta) + d0;
  }
  return bd;
}

/*****************************************************************************/
void printParams(Bars_BarsParams *bp, Bars_PriorParams *pp){
  int i;
  if ((*bp).use_logspline){
    printf("Initial knots: from Logspline\n");
  } else {
    printf("Initial knots: %i, equally spaced\n",(*bp).k);
  }
  printf("Prior: ");
  switch((*pp).prior_id){
  case POISSON:
    printf("Poisson(%lf)\n",(*pp).dparams[0]);
    break;
  case UNIFORM:
    printf("Uniform(%i,..,%i)\n",(*pp).iparams[0],(*pp).iparams[1]);
    break;
  case USER:
    printf("User specified:\n");
    for(i=(*pp).iparams[0];i<=(*pp).iparams[1];i++){
      printf("pi( k = %i ) = %le\n",i,(*pp).dparams[i]);
    }
    printf("\n");
    break;
  default:
    printf("Unknown\n");
    break;
  }
  fflush(stdout);
}

/*****************************************************************************/
void updateParams(Bars_BarsParams *bp, Bars_PriorParams *pp,
          int n_pairs, char **param_names, char **param_values, double *dparin, int *iparin){
  int i,j,*ipar,ni = 0,nd = 0;
  enum PriorForm pid = UNSPECIFIED;
  double d,*dpar;
  char *pn, *pv;
  /* First check for user-specified prior_form parameter */
  for(i=0;i<n_pairs;i++){
    pn = param_names[i];
    pv = param_values[i];
    if ((strcmp(pn,"prior_form") == 0) && (strcmp(pv,"Poisson") == 0)){
      pid = POISSON;
      nd = 1;
      dpar = tssdvec(nd);
      dpar[0] = 6.0;
    }
    if ((streq(pn,"prior_form")) && (streq(pv,"Uniform"))){
      pid = UNIFORM;
      ni = 2;
      ipar = tssivec(ni);
      ipar[0] = 1;
      ipar[1] = MAXKNOTS;
    }
    if ((streq(pn,"prior_form")) && (streq(pv,"User"))){
      pid = USER;
      ni = 2;
      ipar = tssivec(ni);
      ipar[0] = iparin[0];
      ipar[1] = iparin[1];
      nd = MAXKNOTS + 1;
      d = 0.0;
      for(j=1;j<nd;j++){
        d += dparin[j];
      }
      if (d < (1.0 - sqrt(DBL_EPSILON))) {
        /* something is wrong - probably no prior file */
        mexErrMsgTxt("Please check that a valid prior file was supplied.\n");
      }
      dpar = tssdvec(nd);
      memcpy(dpar,dparin,sizeof(double)*nd);
    }
  }
  /* Check for other user-specified parameters */
  for(i=0;i<n_pairs;i++){
    pn = param_names[i];
    pv = param_values[i];
    if (streq(pn,"verbose")){
      VERBOSE = strbool(pv);
    } else
    if (streq(pn,"Use_Logspline")){
      (*bp).use_logspline = strbool(pv);
    } else
    if (strcmp(pn,"burn-in_iterations") == 0){
      sscanf(pv,"%i",&j);
      (*bp).burn_iter = j;
    } else
    if (strcmp(pn,"sample_iterations") == 0){
      sscanf(pv,"%i",&j);
      (*bp).samp_iter = j;
    } else
    if (strcmp(pn,"initial_number_of_knots") == 0){
      sscanf(pv,"%i",&j);
      (*bp).k = j;
    } else
    if (strcmp(pn,"beta_iterations") == 0){
      sscanf(pv,"%i",&j);
      (*bp).beta_iter = j;
    } else
    if (strcmp(pn,"beta_threshhold") == 0){
      sscanf(pv,"%lf",&d);
      (*bp).threshhold = d;
    } else
    if (strcmp(pn,"reversible_jump_constant_c") == 0){
      sscanf(pv,"%lf",&d);
      (*bp).probbd = d;
    } else
    if (strcmp(pn,"proposal_parameter_tau") == 0){
      sscanf(pv,"%lf",&d);
      (*bp).tau = d;
    } else
    if (strcmp(pn,"confidence_level") == 0){
      sscanf(pv,"%lf",&d);
      (*bp).conf_level = d;
    } else
    if (strcmp(pn,"number_of_grid_points") == 0){
      sscanf(pv,"%i",&j);
      (*bp).nf = j;
    } else
    if (strcmp(pn,"sampled_knots_file") == 0){
      strcpy((*bp).iter_knots_fname,pv);
      (*bp).use_iter_knots = isNotNone(pv);
    } else
    if (strcmp(pn,"sampled_mu_file") == 0){
      strcpy((*bp).iter_mu_fname,pv);
      (*bp).use_iter_mu = isNotNone(pv);
    } else
    if (strcmp(pn,"sampled_mu-grid_file") == 0){
      strcpy((*bp).iter_mufine_fname,pv);
      (*bp).use_iter_mufine = isNotNone(pv);
    } else
    if (strcmp(pn,"sampled_params_file") == 0){
      strcpy((*bp).iter_params_fname,pv);
      (*bp).use_iter_params = isNotNone(pv);
    } else
    if (strcmp(pn,"summary_mu_file") == 0){
      strcpy((*bp).summ_mu_fname,pv);
      (*bp).use_summ_mu = isNotNone(pv);
    } else
    if (strcmp(pn,"summary_mu-grid_file") == 0){
      strcpy((*bp).summ_mufine_fname,pv);
      (*bp).use_summ_mufine = isNotNone(pv);
    } else
    if (strcmp(pn,"summary_params_file") == 0){
      strcpy((*bp).summ_params_fname,pv);
      (*bp).use_summ_params = isNotNone(pv);
    } else
    if ((strcmp(pn,"Poisson_parameter_lambda") == 0) && (pid == POISSON)){
      sscanf(pv,"%lf",dpar);
    } else
    if ((streq(pn,"Uniform_parameter_L")) && (pid == UNIFORM)){
        sscanf(pv,"%i",ipar);
    } else
    if ((streq(pn,"Uniform_parameter_U")) && (pid == UNIFORM)){
        sscanf(pv,"%i",ipar+1);
    }
  }
  imposeParamConstraints(bp,pp);
  Set_PriorParams(pp,bp,pid,dpar,ipar);
  if (ni > 0) free(ipar);
  if (nd > 0) free(dpar);
}

/*****************************************************************************/
void setReadParamValues(Bars_BarsParams *bp, Bars_PriorParams *pp,
            char *c, double *dpar, int *ipar){
  int param_pairs;
  char **param_names, **param_values;
  param_names = tsscmat(MAX_PARAMS,MAX_VAL_LEN);
  param_values = tsscmat(MAX_PARAMS,MAX_VAL_LEN);
  param_pairs = read_params_file(param_names,param_values,c);
  updateParams(bp,pp,param_pairs,param_names,param_values,dpar,ipar);
  free_cmat(param_values,MAX_PARAMS);
  free_cmat(param_names,MAX_PARAMS);
}

/*****************************************************************************/
void setReadPriorValues(double *dpar, int *ipar,
      char *c){

  FILE *prior_file;
  int i,j,fval,firstnz,lastnz;
  double pii,epspi;

  epspi = sqrt(DBL_EPSILON);
  fval = 2;
  for(i=0;i<=MAXKNOTS;i++){
    dpar[i] = 0.0;
  }
  prior_file = fopen(c,"r");
  while(fval == 2){
    fval = fscanf(prior_file,"%i %lf\n",&i,&pii);
    if (fval == 2){
      if ((i >= 1) && (i <= MAXKNOTS)){
        dpar[i] = pii;
      }
    } else {
      if (fval != EOF){
        fclose(prior_file);
        mexErrMsgTxt("Prior file format error.\n");
      }
    }
  }
  fclose(prior_file);
  pii = 0.0;
  firstnz = -1;
  lastnz = -1;
  for(i=0;i<=MAXKNOTS;i++){
    if (dpar[i] < DBL_EPSILON) dpar[i] = 0.0;
    pii += dpar[i];
  }
  if (pii >= (MAXKNOTS * DBL_EPSILON)){
    for(i=0;i<=MAXKNOTS;i++){
      if ((dpar[i] >= DBL_EPSILON) && (firstnz < 0)) firstnz = i;
      j = MAXKNOTS - i;
      if ((dpar[j] >= DBL_EPSILON) && (lastnz < 0)) lastnz = j;
    }
    if ((firstnz < 1) || (firstnz > MAXKNOTS) ||
        (lastnz < 1) || (lastnz > MAXKNOTS) ||
        (firstnz > lastnz)){
        mexErrMsgTxt("error while processing prior file.\n");
    } else {
      pii = 0.0;
      for(i=firstnz;i<=lastnz;i++){
        if (dpar[i] < epspi) dpar[i] = epspi;
        pii += dpar[i];
      }
      for(i=firstnz;i<=lastnz;i++){
        dpar[i] /= pii;
      }
      ipar[0] = firstnz;
      ipar[1] = lastnz;
    }
  }
}

/************************* matlab entry point **************************/
void mexFunction(
  int nlhs,              // Number of left hand side (output) arguments
  mxArray *plhs[],       // Array of left hand side arguments
  int nrhs,              // Number of right hand side (input) arguments
  const mxArray *prhs[]  // Array of right hand side arguments
)
{
  int strlen, nusrprms, numdat, ik, br, sm, ft, bn, prpu, prpl, tr;
  char *inbuf;
  double *xd, *priorp, *tu, *cnt, *cnf;
  int *iknt, *brnin, *sms, *fts, *bns, *prparu, *prparl, *trl;
  __int32 *yd;
  bool *vrb, vr, *lgs, lg;

  FILE *data_file, *output_file, *fit_file, *param_file, *iter_file;
  int n,nf,np,*ipar;
  double *dpar;
  Bars_X_Grid *xg;
  Bars_BinnedData *bd;
  Bars_BarsParams *bp;
  Bars_PriorParams *pp;
  Bars_Model *m1,*m2;
  Bars_WorkSpace *ws;
  Bars_SampStat *ss;
  Bars_OutputStat *os;
  time_t t1,t2;

  time(&t1);
  setseed(&t1);
  VERBOSE = 0;

  /*******************input error check*********************/
  if (nrhs != 15) mexErrMsgTxt("Fifteen inputs required.");
  if ((mxGetN(prhs[X_DAT]) != 1) | (mxGetN(prhs[Y_DAT]) !=1 ))
  {
    mexErrMsgTxt("x and y must be column vectors");
  }
  if (mxGetM(prhs[X_DAT]) != mxGetM(prhs[Y_DAT]))
  {
    mexErrMsgTxt("x and y must have same length");
  }
  if (!mxIsInt32(prhs[Y_DAT]))
  {
    mexErrMsgTxt("y data type must be int32.");
  }
  if ((!mxIsLogical(prhs[LGSP])) || ((mxGetN(prhs[LGSP]) * mxGetM(prhs[LGSP])) != 1))
  {
    mexErrMsgTxt("initial must be a scalar logical.");
  }
  if (mxGetNumberOfElements(prhs[IKNT]) != 1)
  {
    mexErrMsgTxt("iknots must be scalar");
  }
  if (mxGetN(prhs[BRNIN]) * mxGetM(prhs[BRNIN]) != 1)
  {
    mexErrMsgTxt("burnin must be scalar");
  }
  if (mxGetN(prhs[SMS]) * mxGetM(prhs[SMS]) != 1)
  {
    mexErrMsgTxt("sims must be scalar");
  }
  if (mxGetN(prhs[TU]) * mxGetM(prhs[TU]) != 1)
  {
    mexErrMsgTxt("tau must be scalar");
  }
  if (mxGetN(prhs[CNT]) * mxGetM(prhs[CNT]) != 1)
  {
    mexErrMsgTxt("c must be scalar");
  }
  if (mxGetN(prhs[FTS]) * mxGetM(prhs[FTS]) != 1)
  {
    mexErrMsgTxt("fits must be scalar");
  }
  if (mxGetN(prhs[CNF]) * mxGetM(prhs[CNF]) != 1)
  {
    mexErrMsgTxt("conf must be scalar");
  }
  if (mxGetN(prhs[BNS]) * mxGetM(prhs[BNS]) != 1)
  {
    mexErrMsgTxt("bins must be scalar");
  }
  if (mxGetN(prhs[TRL]) * mxGetM(prhs[TRL]) != 1)
  {
    mexErrMsgTxt("trials must be scalar");
  }
  if ((!mxIsLogical(prhs[VRB])) || ((mxGetN(prhs[VRB]) * mxGetM(prhs[VRB])) != 1))
  {
    mexErrMsgTxt("verbose must a scalar logical.");
  }
  if (!mxIsChar(prhs[PR]))
  {
    mexErrMsgTxt("prior must be a string.");
  }else
  {
    strlen = (mxGetN(prhs[PR]) * mxGetM(prhs[PR])) + 1;
    inbuf = mxCalloc(strlen, sizeof(char));
    if (mxGetString(prhs[PR],inbuf,strlen) == 0)
    {
      if (strcmp(inbuf, "Poisson") == 0)
      {
        if((mxGetN(prhs[PRP]) * mxGetM(prhs[PRP])) != 1)
        {
          mxFree(inbuf);
          mexErrMsgTxt("if prior = Poisson, priorparam must be scalar");
        }
        priorp = mxGetPr(prhs[PRP]);
      }
      else if (strcmp(inbuf, "Uniform") == 0)
      {
        if((mxGetN(prhs[PRP]) * mxGetM(prhs[PRP])) != 2)
        {
          mxFree(inbuf);
          mexErrMsgTxt("if prior = Uniform, priorparam must have 2 elements.");
        }
        priorp = mxGetPr(prhs[PRP]);
    if (priorp[0] >= priorp[1])
    {
      mxFree(inbuf);
      mexErrMsgTxt("when prior = Uniform, priorparam(1) must be < priorparam(2)");
    }
    prpl = (int)priorp[0]; prparl = &prpl;
    prpu = (int)priorp[1]; prparu = &prpu;
      }
      else if (strcmp(inbuf, "User") == 0)
      {
        if (mxGetM(prhs[PRP]) != 2)
        {
          mxFree(inbuf);
          mexErrMsgTxt("priorparam must have 2 columns.");
        }
        nusrprms = mxGetN(prhs[PRP]);
        priorp = mxGetPr(prhs[PRP]);
        write_user_priorparam(priorp,nusrprms);
      }
      else
      {
        mxFree(inbuf);
        mexErrMsgTxt("prior must be Poisson, Uniform, or User.");
      }
    }
  }

  /***************get input data and parameters*****************/
  numdat = mxGetM(prhs[X_DAT]);
  xd = mxGetPr(prhs[X_DAT]);
  yd = mxGetData(prhs[Y_DAT]);
  tr = (int)*mxGetPr(prhs[TRL]); trl = &tr;
  ik = (int)*mxGetPr(prhs[IKNT]); iknt = &ik;
  br = (int)*mxGetPr(prhs[BRNIN]); brnin = &br;
  sm = (int)*mxGetPr(prhs[SMS]); sms = &sm;
  ft = (int)*mxGetPr(prhs[FTS]); fts = &ft;
  bn = (int)*mxGetPr(prhs[BNS]); bns = &bn;
  tu = mxGetPr(prhs[TU]);
  cnt = mxGetPr(prhs[CNT]);
  cnf = mxGetPr(prhs[CNF]);
  lgs = mxGetData(prhs[LGSP]); lg = *lgs;
  vrb = mxGetData(prhs[VRB]); vr = *vrb;

  /***************write parameter and data files***************/
  priorsetup(brnin,sms,iknt,tu,cnt,cnf,fts);
  if (strcmp(inbuf, "Poisson") == 0) paramp(priorp);
  if (strcmp(inbuf, "Uniform") == 0) paramu(prparu,prparl);
  if (strcmp(inbuf, "User") == 0) paramuser();
  if (lg) logtrue(); else logfalse();
  if (vr) VERBOSE = 1;
  write_datfile(bns,trl,xd,yd,numdat);


  np = 3;
  bp = New_BarsParams();
  pp = New_PriorParams();
  setDefaultParamValues(bp,pp);
  dpar = tssdvec0(MAXKNOTS + 1);
  ipar = tssivec(2);
  ipar[0] = 1;
  ipar[1] = MAXKNOTS;
  if (strcmp(inbuf, "User") == 0)
  {
  setReadPriorValues(dpar,ipar,"user_prior_params");
  }
  setReadParamValues(bp,pp,"bars_params",dpar,ipar);
  nf = (*bp).nf;
  bd = read_data_file("input_data",nf);
  n = (*bd).n;
  m1 = New_Model(n,nf);
  m2 = New_Model(n,nf);
  ws = New_WorkSpace(n,nf);
  ss = New_SampStat(n, nf,np,(*bp).samp_iter);
  os = New_OutputStat(np);
  if (VERBOSE) printParams(bp,pp);

  /*******************begin running mcmc*********************/
  Bars_MCMC(m1,m2,bd,ws,bp,pp,ss,os);
  if (VERBOSE){
    printf("\nProportion of Birth Moves Accepted: %lf\n",(*os).accept[0]);
    printf("Proportion of Death Moves Accepted: %lf\n",(*os).accept[1]);
    printf("Proportion of Relocation Moves Accepted: %lf\n",(*os).accept[2]);
    printf("Overall Proportion of Moves Accepted: %lf\n",(*os).accept[3]);
  }

  Free_OutputStat(os);
  Free_SampStat(ss);
  Free_WorkSpace(ws);
  Free_Model(m2);
  Free_Model(m1);
  Free_PriorParams(pp);
  Free_BarsParams(bp);
  Free_BinnedData(bd);
  mxFree(inbuf);

  if (VERBOSE){
    time(&t2);
    printf("\nTotal run time: %lf seconds\n",difftime(t2,t1));
  }
}
