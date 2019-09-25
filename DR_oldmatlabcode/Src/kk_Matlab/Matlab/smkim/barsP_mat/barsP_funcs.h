/* barsP_funcs.h */
#ifndef _barsPfuncs
#define _barsPfuncs
/* iterations for glm fitting. */
#define MAXIT_GLM 20

/* constants for memory allocation */
/* N_FINE is a default value for grid size. User can override default at runtime. */
#define MAXKNOTS 80
#define N_FINE 500
#define MAX_PARAMS 500
#define MAX_VAL_LEN 200
#define EPS_FACTOR 500.0


enum PriorForm {UNSPECIFIED,POISSON,UNIFORM,USER};

typedef struct X_Grid {
  int n;
  double *x;
} Bars_X_Grid;

typedef struct BinnedData {
  int n, nf,trials;
  Bars_X_Grid *xg,*fg;
  double *x_raw, x_rawmax,x_rawmin;
  double *y, bin_width, scale_factor;
} Bars_BinnedData;

typedef struct Model {
  int n,nf,k,fit_info,fit_fine;
  double *knots_interior, *knots_all, *basis, *fitted_values;
  double loglikelihood,*J,*beta,sigma,*basis_fine;
  double accept_prob,*random_mu_fine,*random_mu,full_bic,full_loglikelihood;
} Bars_Model;

typedef struct BarsParams {
  int burn_iter, samp_iter, k,nf,use_logspline,beta_iter;
  double probbd, tau, conf_level,threshhold;
  char *iter_knots_fname,*iter_mu_fname,*iter_mufine_fname,*iter_params_fname;
  char *summ_mu_fname,*summ_mufine_fname,*summ_params_fname;
  FILE *iter_knots_file,*iter_mu_file,*iter_mufine_file,*iter_params_file;
  FILE *summ_mu_file,*summ_mufine_file,*summ_params_file;
  int use_iter_knots, use_iter_mu, use_iter_mufine, use_iter_params;
  int use_summ_mu, use_summ_mufine, use_summ_params;
} Bars_BarsParams;

typedef struct PriorParams {
  enum PriorForm prior_id;
  int num_dparams,num_iparams;
  int *iparams;
  double *dparams;
} Bars_PriorParams;

typedef struct SampStat {
  int np,sim;
  double *random_mu_mean, *random_mufine_mean;
  double *random_mu_mode, *random_mufine_mode;
  double **params;
  double *means,*modes;
  double *full_bic,max_full_bic;
} Bars_SampStat;

typedef struct OutputStat {
/*
     each row of params is a vector of length four:
     ci_low, ci_high, mean, mode
     */
  int np;
  double **params,*accept;
} Bars_OutputStat;

typedef struct WorkSpace {
  int n_max; 
  double *conwork,*con_basis,*t_con_basis,*tmp_basis;
  double *x_con,*work,*tau,*tmp_knots;
  int *coniwork,*deriv,*off,*index_m,*deriv_con;
} Bars_WorkSpace;

Bars_X_Grid *New_X_Grid (int);
void Free_X_Grid(Bars_X_Grid *);
Bars_BarsParams *New_BarsParams();
void Free_BarsParams(Bars_BarsParams *);
Bars_PriorParams *New_PriorParams();
void Set_PriorParams(Bars_PriorParams *, Bars_BarsParams *, enum PriorForm, double *, int *);
void Free_PriorParams(Bars_PriorParams *);
Bars_BinnedData *New_BinnedData(int, int);
void Free_BinnedData(Bars_BinnedData *);
Bars_Model *New_Model(int, int);
void Free_Model(Bars_Model *);
Bars_WorkSpace *New_WorkSpace(int,int);
void Free_WorkSpace(Bars_WorkSpace *);
Bars_OutputStat *New_OutputStat(int);
void Free_OutputStat(Bars_OutputStat *);
Bars_SampStat *New_SampStat (int, int, int, int);
void Free_SampStat(Bars_SampStat *);
void Bars_MCMC(Bars_Model *, Bars_Model *, Bars_BinnedData *, 
           Bars_WorkSpace *, Bars_BarsParams *, 
           Bars_PriorParams *,Bars_SampStat *,Bars_OutputStat *);    
void setDefaultParamValues(Bars_BarsParams *, Bars_PriorParams *);         
           
#endif












