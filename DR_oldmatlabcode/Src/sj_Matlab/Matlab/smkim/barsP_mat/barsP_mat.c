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

int VERBOSE;

/* see barslib.h for struct definitions */

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


void paramp(double *priorparam){
  FILE *f;
  f=fopen("bars_params","a");
  fprintf(f,"SET prior_form = Poisson\n");
  fprintf(f,"SET Poisson_parameter_lambda = %8.4f\n",*priorparam);
  fclose(f);
}

void paramu(int *upper,int *lower){
  FILE *f;
  f=fopen("bars_params","a");
  fprintf(f,"SET prior_form = Uniform\n");
  fprintf(f,"SET Uniform_parameter_U = %i\n",*upper);
  fprintf(f,"SET Uniform_parameter_L = %i\n",*lower);
  fclose(f);
}

void paramuser(){
  FILE *f;
  f=fopen("bars_params","a");
  fprintf(f,"SET prior_form = User\n");
  fclose(f);
}


void logfalse(){
  FILE *f;
  f=fopen("bars_params","a");
  fprintf(f,"SET Use_Logspline = false\n");
  fprintf(f,"SET verbose = false\n");
  fclose(f);
}

void logtrue(){
  FILE *f;
  f=fopen("bars_params","a");
  fprintf(f,"SET Use_Logspline = true\n");
  fprintf(f,"SET verbose = false\n");
  fclose(f);
}

void setBinnedInfo(Bars_BinnedData *bd,double bwidth, int tr){
  double sf = bwidth * (double)tr;
  (*bd).bin_width = bwidth;
  (*bd).trials = tr;
  (*bd).scale_factor = (sf > 0.0) ? (1.0/sf) : 0.0;
}

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
  mxFree(prefix);
  return pnum;
}

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
    printf("Data file format error: The first line should contain two integer entries: the number of bins, and the number of trials.\n\n");
    exit(1);
  }
  bd = New_BinnedData(n,nf);
  y = (*bd).y;
  x = (*bd).x_raw;
  for(i=0;i<n;i++){
    fval = fscanf(data_file,"%lf %lf\n",x+i,y+i);
    if (fval != 2){
      if (fval >= 0) printf("Data file format error: Data for bin %i contains format errors.\n\n",(i+1));
      if (fval == EOF) printf("Data file format error: EOF reached before data for all %i bins could be read.\n\n",n);
      fclose(data_file);
      exit(1);
    } else {
      if (y[i] < 0.0){
	printf("Data file error: Count for bin %i is negative.\n\n",(i+1));
	fclose(data_file);
	exit(1);
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

Bars_BinnedData *read_mat_struct(double * data, int tr, int n, int nf)
{
  Bars_X_Grid *xg;
  Bars_BinnedData *bd;
  int i,fval;
  double bw,*y,*x,*x_norm,xmin,xmax,delta,d0;

/* can add a parameter to allow for non-central bin markers.
   parameter should be stored in binned_data and passed
   to normalize_bin to adjust the min and max appropriately. 
*/

  bd = New_BinnedData(n,nf);
  y = (*bd).y;
  x = (*bd).x_raw;
  for(i=0;i<n;i++)
  {
    x[i] = data[i];
    y[i] = data[i+n];
      
    if (y[i] < 0.0)
    {
        printf("Data file error: Count for bin %i is negative.\n\n",i);
        exit(1);
    }
  }

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
void updateParams(Bars_BarsParams *bp, Bars_PriorParams *pp, 
          int n_pairs, char **param_names, char **param_values, double *dparin, int *iparin){
  int i,j,*ipar,ni = 0,nd = 0;
  enum PriorForm pid = UNSPECIFIED;
  double d,*dpar;
  char *pn, *pv;
  /* First check for user-specified prior_form parameter */

  fprintf(stderr, "before checking for prior form\n");
  fprintf(stderr, "ipar[0] = %i\n",ipar[0]);
  fprintf(stderr, "\n");
  fprintf(stderr, "ipar[1] = %i\n",ipar[1]);


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
	printf("Please check that a valid prior file was supplied.\n");
	exit(1);
      }
      dpar = tssdvec(nd);
      memcpy(dpar,dparin,sizeof(double)*nd);
    }
  }

  fprintf(stderr, "after checking prior form, before checking values\n");
  fprintf(stderr, "ipar[0] = %i\n",ipar[0]);
  fprintf(stderr, "ipar[1] = %i\n",ipar[1]);
  fprintf(stderr, "\n");


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

  fprintf(stderr, "after checking other params\n");
  fprintf(stderr, "ipar[0] = %i\n",ipar[0]);
  fprintf(stderr, "ipar[1] = %i\n",ipar[1]);
  fprintf(stderr, "\n");

  if (ni > 0) mxFree(ipar);
  if (nd > 0) mxFree(dpar);
}


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
	printf("Prior file format error.\n");
	fclose(prior_file);
	exit(1);
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
      if ((dpar[i] >= DBL_EPSILON) && (firstnz < 0))
	firstnz = i;
      j = MAXKNOTS - i;
      if ((dpar[j] >= DBL_EPSILON) && (lastnz < 0))
	lastnz = j;
    }
    if ((firstnz < 1) || (firstnz > MAXKNOTS) ||
	(lastnz < 1) || (lastnz > MAXKNOTS) ||
	(firstnz > lastnz)){
      printf("error while processing prior file.\n");
      exit(1);
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

void mexFunction(
	int nlhs,              /* Number of left hand side (output) arguments */
	mxArray *plhs[],       /* Array of left hand side arguments */
	int nrhs,              /* Number of right hand side (input) arguments */
	const mxArray *prhs[]  /* Array of right hand side arguments */ )
{
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
    int trials;
    
    time(&t1);
    setseed(&t1);
    VERBOSE = 1;

    np = 3;
    bp = New_BarsParams();
    pp = New_PriorParams();
    setDefaultParamValues(bp,pp);
    dpar = tssdvec0(MAXKNOTS + 1);
    ipar = tssivec(2);
    ipar[0] = 1;
    ipar[1] = MAXKNOTS;

    if (nrhs == 1)
        trials = 1;
    else
        trials = (int)mxGetScalar(prhs[1]);

	if (nrhs > 3)
    {
        setReadPriorValues(dpar,ipar,mxArrayToString(prhs[3]));        
    }
    if (nrhs > 2)
    {
        setReadParamValues(bp,pp,mxArrayToString(prhs[2]),dpar,ipar);
    }

    nf = (*bp).nf;
/*   bd = read_data_file("example.data",nf);*/
    bd = read_mat_struct(mxGetPr(prhs[0]),trials,mxGetM(prhs[0]),nf);
   
   n = (*bd).n;

    m1 = New_Model(n,nf);
    m2 = New_Model(n,nf);
    ws = New_WorkSpace(n,nf);
    ss = New_SampStat(n, nf,np,(*bp).samp_iter);

    os = New_OutputStat(np);
    if (VERBOSE) printParams(bp,pp);

    Bars_MCMC(m1,m2,bd,ws,bp,pp,ss,os);
    if (VERBOSE)
    {
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
  
    if (VERBOSE)
    {
        time(&t2);
        printf("\nTotal run time: %lf seconds\n",difftime(t2,t1));
    }
}